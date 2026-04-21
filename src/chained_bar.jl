using Statistics

function logsumexp(x::AbstractVector{<:Real})
    m = maximum(x)
    return m + log(sum(exp.(x .- m)))
end

"""
Estimates free energy differences between adjacent temperature states using
chained backward FEP (exponential averaging).

For each adjacent pair (k-1, k), samples from state k must have log_norm
recorded at both beta_collapses[k-1] and beta_collapses[k], i.e.
beta_collapses[k-1] must lie within the measurement window of state k.
The backward FEP estimate is:

  ΔF_k = F_k - F_{k-1} = log E_k[exp(2·(log_norm(β_{k-1}) - log_norm(β_k)))]

Note: 2·(log_norm(β_{k-1}) - log_norm(β_k)) > 0 because imaginary time
cooling accumulates increasingly negative log norms with growing beta.

This is one-sided rather than bidirectional BAR because samples from the
lower state cannot evaluate the upper state's reduced potential.

Arguments
---------
all_measurements : all_measurements[k][i] = measurement Dict for METTS sample i
                   from collapse temperature k; must contain "log_norm" key
all_betas        : all_betas[k] = Vector of beta values measured for state k
beta_collapses   : K collapse temperatures, sorted ascending; beta_collapses[k-1]
                   must lie within the measurement window of state k

Returns
-------
delta_f       : Vector of K-1 free energy differences, delta_f[k] = F_{k+1} - F_k
free_energies : Vector of K cumulative free energies with F_1 = 0 as reference
"""
function chained_bar(all_measurements, all_betas, beta_collapses)
    K = length(beta_collapses)
    @assert issorted(beta_collapses) "beta_collapses must be sorted ascending"
    @assert length(all_measurements) == K && length(all_betas) == K

    delta_f = zeros(K - 1)

    for k in 2:K
        β_k   = beta_collapses[k]
        β_km1 = beta_collapses[k - 1]
        betas_k = all_betas[k]
        meas_k  = all_measurements[k]
        N_k = length(meas_k)

        idx_km1 = findfirst(b -> isapprox(b, β_km1; atol=1e-8), betas_k)
        idx_k   = findfirst(b -> isapprox(b, β_k;   atol=1e-8), betas_k)
        isnothing(idx_km1) && error("beta_collapses[$(k-1)] = $β_km1 not found in " *
                                    "betas of state $k — extend the measurement window")
        isnothing(idx_k)   && error("beta_collapses[$k] = $β_k not found in betas of state $k")

        log_ratios = [2 * (meas_k[i]["log_norm"][idx_km1] - meas_k[i]["log_norm"][idx_k])
                      for i in 1:N_k]
        delta_f[k - 1] = logsumexp(log_ratios) - log(N_k)
    end

    return delta_f, [0.0; cumsum(delta_f)]
end


"""
Refine free energy estimates by iterating the MBAR self-consistency equations.

Starting from an initial guess (e.g. from chained_bar), iterates:

  f_j ← -log Σ_{(k,i): β_j ∈ window_k} exp(2·log_norm_{k,i}(β_j) - log_denom_{k,i})

where log_denom_{k,i} = logsumexp_j'[ log(N_j') + f_j' + 2·log_norm_{k,i}(β_j') ]

For adjacent pairs with overlapping windows this converges to the BAR (Bennett
Acceptance Ratio) estimator, which uses samples from both sides and has much
lower variance than the one-sided backward FEP used by chained_bar.

Arguments
---------
all_measurements : same as chained_bar
all_betas        : same as chained_bar
beta_collapses   : same as chained_bar
f_init           : initial free energy vector (length K, f_init[1] = 0)
max_iter         : maximum number of iterations
tol              : convergence threshold on max |f_new - f|

Returns
-------
free_energies : length-K vector of self-consistent free energies (F_1 = 0)
"""
function mbar_free_energies(all_measurements, all_betas, beta_collapses, f_init;
                            max_iter=500, tol=1e-10)
    K  = length(beta_collapses)
    Nk = [length(all_measurements[k]) for k in 1:K]
    f  = copy(f_init)

    avail_states = [Int[] for _ in 1:K]
    avail_idx    = [Int[] for _ in 1:K]
    for k in 1:K, j in 1:K
        idx = findfirst(b -> isapprox(b, beta_collapses[j]; atol=1e-8), all_betas[k])
        if !isnothing(idx)
            push!(avail_states[k], j)
            push!(avail_idx[k], idx)
        end
    end

    for _ in 1:max_iter
        log_denom = [[logsumexp([log(Nk[j]) + f[j] +
                                 2*all_measurements[k][i]["log_norm"][avail_idx[k][jj]]
                                 for (jj, j) in enumerate(avail_states[k])])
                      for i in 1:Nk[k]]
                     for k in 1:K]

        f_new = zeros(K)
        for j in 1:K
            terms = Float64[]
            for k in 1:K
                jj = findfirst(==(j), avail_states[k])
                isnothing(jj) && continue
                for i in 1:Nk[k]
                    push!(terms, 2*all_measurements[k][i]["log_norm"][avail_idx[k][jj]] -
                                 log_denom[k][i])
                end
            end
            f_new[j] = isempty(terms) ? f[j] : -logsumexp(terms)
        end
        f_new .-= f_new[1]

        maximum(abs.(f_new .- f)) < tol && return f_new
        f .= f_new
    end
    return f
end


"""
MBAR reweighting of an observable using chained BAR free energies.

For each beta covered by any state's window, collects all available samples
and assigns MBAR weights:

  log w_i(β) = 2·log_norm_i(β)
             - logsumexp_j[ log(N_j) + f_j + 2·log_norm_i(β_cj) ]

where the sum runs over all states j whose collapse temperature β_cj was
recorded in sample i's measurement window. The free energies f_j from
chained_bar enter the denominator so that samples from different states
are combined consistently.

Arguments
---------
all_measurements : all_measurements[k][i] = measurement Dict for METTS sample i
all_betas        : all_betas[k] = Vector of beta values for state k
beta_collapses   : K collapse temperatures, sorted ascending
free_energies    : length-K vector returned by chained_bar (F_1 = 0)
observable       : key string present in each measurement Dict

Returns
-------
betas_out : beta values sorted ascending (each unique beta appears once)
obs_mean  : MBAR estimate of ⟨observable⟩ at each beta
obs_neff  : effective sample size 1/Σ_i w_i² (reweighting quality diagnostic)
"""
function mbar_reweight_observable(all_measurements, all_betas, beta_collapses,
                                  free_energies, observable::String)
    K  = length(beta_collapses)
    Nk = [length(all_measurements[k]) for k in 1:K]

    # For each state k, find which collapse temperatures appear in its window
    # and store their indices within all_betas[k]
    avail_states = [Int[] for _ in 1:K]
    avail_idx    = [Int[] for _ in 1:K]
    for k in 1:K
        for j in 1:K
            idx = findfirst(b -> isapprox(b, beta_collapses[j]; atol=1e-8), all_betas[k])
            if !isnothing(idx)
                push!(avail_states[k], j)
                push!(avail_idx[k],    idx)
            end
        end
    end

    # Precompute log denominator for every sample:
    # log_denom[k][i] = logsumexp_j[ log(N_j) + f_j + 2·log_norm_i(β_cj) ]
    log_denom = [[logsumexp([log(Nk[j]) + free_energies[j] +
                             2 * all_measurements[k][i]["log_norm"][avail_idx[k][jj]]
                             for (jj, j) in enumerate(avail_states[k])])
                  for i in 1:Nk[k]]
                 for k in 1:K]

    # Collect all unique beta values
    all_beta_vals = sort(unique(vcat(all_betas...)))

    betas_out = Float64[]
    obs_mean  = Float64[]
    obs_neff  = Float64[]

    for beta in all_beta_vals
        log_w = Float64[]
        obs   = Float64[]

        for k in 1:K
            beta_idx = findfirst(b -> isapprox(b, beta; atol=1e-8), all_betas[k])
            isnothing(beta_idx) && continue
            for i in 1:Nk[k]
                ln_num = 2 * all_measurements[k][i]["log_norm"][beta_idx]
                push!(log_w, ln_num - log_denom[k][i])
                push!(obs,   all_measurements[k][i][observable][beta_idx])
            end
        end

        isempty(log_w) && continue
        log_w .-= maximum(log_w)
        w = exp.(log_w)
        w ./= sum(w)

        push!(betas_out, beta)
        push!(obs_mean,  sum(w .* obs))
        push!(obs_neff,  1 / sum(w .^ 2))
    end

    return betas_out, obs_mean, obs_neff
end


"""
Same as mbar_reweight_observable but also returns bootstrap standard errors for
both the observable and the free energies.

Each bootstrap resample draws N_k samples with replacement from state k,
re-runs the chained BAR free energy estimation followed by MBAR self-consistency
iterations, and re-estimates the observable and free energies. The standard
deviation over bootstrap replicates is the error estimate.

Note: METTS samples are Markov-chain correlated, so bootstrap errors
underestimate the true statistical error. Use n_bootstrap >= 500.

Returns
-------
betas_out       : beta values at which the observable is estimated
obs_mean        : MBAR estimate of ⟨observable⟩ at each beta
obs_err         : bootstrap standard error of obs_mean
obs_neff        : effective sample size at each beta
free_energies   : self-consistent MBAR free energies (length K, F_1 = 0)
free_energy_err : bootstrap standard error of free_energies (free_energy_err[1] = 0)
"""
function bootstrap_mbar_reweight_observable(all_measurements, all_betas, beta_collapses,
                                            observable::String; n_bootstrap::Int=500)
    # Full-data estimate: start from chained BAR, refine with MBAR iteration
    _, f_bar       = chained_bar(all_measurements, all_betas, beta_collapses)
    free_energies  = mbar_free_energies(all_measurements, all_betas, beta_collapses, f_bar)
    betas_out, obs_mean, obs_neff =
        mbar_reweight_observable(all_measurements, all_betas, beta_collapses,
                                 free_energies, observable)

    n_betas  = length(betas_out)
    K        = length(beta_collapses)
    boot_obs_mat = zeros(n_bootstrap, n_betas)
    boot_f_mat   = zeros(n_bootstrap, K)

    for b in 1:n_bootstrap
        meas_boot = [all_measurements[k][rand(1:length(all_measurements[k]),
                                              length(all_measurements[k]))]
                     for k in 1:K]
        _, f_boot      = chained_bar(meas_boot, all_betas, beta_collapses)
        f_boot_mbar    = mbar_free_energies(meas_boot, all_betas, beta_collapses, f_boot)
        _, boot_obs, _ = mbar_reweight_observable(meas_boot, all_betas, beta_collapses,
                                                   f_boot_mbar, observable)
        boot_obs_mat[b, :] = boot_obs
        boot_f_mat[b, :]   = f_boot_mbar
    end

    obs_err         = [std(boot_obs_mat[:, j]) for j in 1:n_betas]
    free_energy_err = [std(boot_f_mat[:, j])   for j in 1:K]
    return betas_out, obs_mean, obs_err, obs_neff, free_energies, free_energy_err
end
