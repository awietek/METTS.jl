using Printf

"""
Time evolution using TDVP with initial basis extension, in-situ measurements,
and a final collapse.

Performs imaginary time evolution to max(t_max, t_collapse) on a grid of
integer multiples of tau. The initial phase extends the basis up to tau0 (using
nsubdiv logarithmically subdivided 1-site sweeps), then takes a single 2-site
sweep from tau0 to tau if tau0 < tau. All subsequent steps are of size tau.
Measurements are taken at each grid point in [t_min, t_max].
At t_collapse the collapse function is called on psi (unmodified).

Note: t_min, t_max, t_collapse are imaginary times, i.e. beta/2 for inverse
temperature beta.

Arguments
---------
H              : Hamiltonian MPO
psi0           : initial product-state MPS
t_min          : lower bound of measurement window (inclusive), = beta_min/2
t_max          : upper bound of measurement window (inclusive), = beta_max/2; must be > t_min
t_collapse     : imaginary time at which collapse is performed = beta_collapse/2;
                 must be an integer multiple of tau
measure        : psi::MPS -> Dict{String, Float64}
collapse_func  : psi::MPS -> Vector{Int}; must not modify psi

Returns
-------
psi            : evolved MPS at max(t_max, t_collapse), not collapsed
log_norm       : accumulated log norm from imaginary time evolution
measurements   : Dict{String, Vector{Float64}} with same keys as returned by
                 measure, each mapped to values at successive grid times in
                 [t_min, t_max]; also contains a "log_norm" key with
                 the cumulative log norm at each measurement point
collapsed_state: Vector{Int} from collapse_func applied to psi at t_collapse
times          : Vector{Float64} of imaginary times at which measurements were taken
"""
function timeevo_tdvp_extend_measurements(
    H::MPO, psi0::MPS,
    t_min::Number, t_max::Number, t_collapse::Number,
    measure::Function, collapse_func::Function;
    tau::Number=0.1, cutoff::Float64=1e-6,
    maxm::Int64=1000, tau0::Float64=0.05,
    nsubdiv::Int64=4, kkrylov::Int64=3,
    normalize::Bool=true, silent=false,
    solver_backend::AbstractString="applyexp",
    shift::Real=0.)

    @assert t_min < t_max "t_min must be strictly less than t_max"
    @assert tau0 <= tau "tau0 must be <= tau"
    n_collapse = round(Int, t_collapse / tau)
    @assert isapprox(n_collapse * tau, t_collapse; rtol=1e-8, atol=1e-10) "t_collapse must be an integer multiple of tau"

    log_norm = 0.0
    N = length(psi0)

    # Logarithmically subdivided 1-site steps summing to tau0
    times_init = [tau0 / 2^(nsubdiv - 1)]
    for e in (nsubdiv-1):-1:1
        push!(times_init, tau0 / 2^e)
    end

    psi = copy(psi0)
    current_time = 0.0
    measurements = Dict{String, Vector{Float64}}()
    times_meas = Float64[]

    function maybe_record!(psi, t)
        t_min-1e-8 <= t <= t_max+1e-8 || return
        result = measure(psi)
        push!(times_meas, t)
        push!(get!(measurements, "log_norm", Float64[]), log_norm)
        for (k, v) in result
            push!(get!(measurements, k, Float64[]), v)
        end
    end

    # --- Basis extension: 0 -> tau0 (nsubdiv logarithmically subdivided 1-site sweeps) ---
    for isub in 1:nsubdiv
        l1 = maxlinkdim(psi)
        t_ext = @elapsed begin
            if !silent
                println("    Performing basis extension ...")
                flush(stdout)
            end
            psi = basis_extend(psi, H; extension_krylovdim=kkrylov,
                extension_cutoff=1e-12)
        end
        l2 = maxlinkdim(psi)
        if !silent
            @printf("    Basis extension, maxm: %4d -> %4d, time: %.5f secs\n",
                l1, l2, t_ext)
        end
        t = @elapsed begin
            psi = tdvp(H, -times_init[isub], psi;
                updater_backend=solver_backend, nsweeps=1, nsite=1,
                cutoff=cutoff, maxdim=maxm, normalize=false)
            log_norm += log(norm(psi))
            normalize!(psi)
        end
        current_time += times_init[isub]
        if !silent
            @printf("    1TDVP sweep, tau: %.5f, SvN: %.4f, maxm: %6d, time: %.5f secs\n",
                times_init[isub], entropy_von_neumann(psi, N ÷ 2), maxlinkdim(psi), t)
            flush(stdout)
        end
    end
    GC.gc()

    # --- One 2-site sweep: tau0 -> tau  (skipped when tau0 == tau) ---
    if tau0 < tau
        step = tau - tau0
        t = @elapsed begin
            psi = tdvp(H, -step, psi;
                updater_backend=solver_backend, nsweeps=1, nsite=2,
                cutoff=cutoff, maxdim=maxm, normalize=false)
            log_norm += log(norm(psi))
            normalize!(psi)
        end
        current_time = tau  # first grid point
        maybe_record!(psi, current_time)
        if !silent
            @printf("    2TDVP sweep, tau: %.5f, SvN: %.4f, maxm: %6d, time: %.5f secs\n",
                step, entropy_von_neumann(psi, N ÷ 2), maxlinkdim(psi), t)
            flush(stdout)
        end
        GC.gc()
    else
        current_time = tau  # tau0 == tau
        maybe_record!(psi, current_time)
    end

    # --- Bulk segment 1: tau -> t_collapse  (n_collapse - 1 steps of tau) ---
    for step_idx in 2:n_collapse
        nsite = maxlinkdim(psi) < maxm ? 2 : 1
        t = @elapsed begin
            psi = tdvp(H, -tau, psi;
                updater_backend=solver_backend, nsweeps=1, nsite=nsite,
                cutoff=cutoff, maxdim=maxm, normalize=false)
            log_norm += log(norm(psi))
            normalize!(psi)
        end
        current_time = step_idx * tau
        maybe_record!(psi, current_time)
        if !silent
            @printf("    %dTDVP sweep, tau: %.5f, SvN: %.4f, maxm: %6d, time: %.5f secs\n",
                nsite, Float64(tau), entropy_von_neumann(psi, N ÷ 2), maxlinkdim(psi), t)
            flush(stdout)
        end
        GC.gc()
    end
    # current_time == t_collapse

    collapsed_state = collapse_func(psi)

    # --- Bulk segment 2: t_collapse -> floor(t_max/tau)*tau ---
    if t_max > t_collapse
        n_extra = floor(Int, (t_max - t_collapse) / tau + 1e-10)
        for extra_idx in 1:n_extra
            nsite = maxlinkdim(psi) < maxm ? 2 : 1
            t = @elapsed begin
                psi = tdvp(H, -tau, psi;
                    updater_backend=solver_backend, nsweeps=1, nsite=nsite,
                    cutoff=cutoff, maxdim=maxm, normalize=false)
                log_norm += log(norm(psi))
                normalize!(psi)
            end
            current_time = t_collapse + extra_idx * tau
            maybe_record!(psi, current_time)
            if !silent
                @printf("    %dTDVP sweep, tau: %.5f, SvN: %.4f, maxm: %6d, time: %.5f secs\n",
                    nsite, Float64(tau), entropy_von_neumann(psi, N ÷ 2), maxlinkdim(psi), t)
                flush(stdout)
            end
            GC.gc()
        end
    end

    return psi, log_norm, measurements, collapsed_state, times_meas
end
