using Printf

""" computes the number of steps, and the remainder of the time step if not commensurate """
function n_steps_remainder(time::Number, tau::Real)
    n_steps = floor(abs(time) / tau)
    time_tau = n_steps * tau * time / abs(time)
    remainder = time - time_tau
    return n_steps, remainder, time_tau
end

"""
    entropy_von_neumann(psi::MPS, b::Int) -> Float64

Von Neumann entanglement entropy `S = -Σ p log p` of the bipartition between
sites `b` and `b+1`. `psi` is not modified and need not be normalized.
"""
function entropy_von_neumann(psi::MPS, b::Int)
    N = length(psi)
    1 <= b < N || error("bond b=$b must lie in 1:$(N - 1)")
    psi = orthogonalize(psi, b)
    linds = b == 1 ? (siteind(psi, b),) : (linkind(psi, b - 1), siteind(psi, b))
    _, S, _ = svd(psi[b], linds...)
    p = [S[n, n]^2 for n in 1:dim(S, 1)]
    p ./= sum(p)
    # singular values that are zero (or numerically so) do not contribute
    return -sum((x * log(x) for x in p if x > 1e-16); init=0.0)
end

"""
Basic time evolution using TDVP

starting out with 2TDVP until a maximal bond dimension is achieved
then it switches to 1TDVP
a commensurate step size is automatically chosen
"""
function timeevo_tdvp(H::MPO, psi0::MPS, time::Number;
    tau::Number=0.1, cutoff::Float64=1e-6,
    maxm::Int64=1000,
    silent=false, solver_backend::AbstractString="applyexp")
    # nothing to evolve; n_steps_remainder would divide 0 by 0
    iszero(time) && return copy(psi0), 0.0
    N = length(psi0)
    n_steps, remainder, time_tau = n_steps_remainder(time, tau)
    psi = copy(psi0)

    log_norm = 0.0

    for step in 1:n_steps
        if maxlinkdim(psi) < maxm
            t = @elapsed begin
                psi = tdvp(H, time_tau / n_steps, psi;
                    updater_backend=solver_backend,
                    nsweeps=1,
                    nsite=2,
                    cutoff=cutoff,
                    maxdim=maxm,
                    normalize=false)

                log_norm += 2.0*log(norm(psi))
                normalize!(psi)
            end

            svn = entropy_von_neumann(psi, N ÷ 2)
            linkdim = maxlinkdim(psi)
            if !silent
                @printf("    2TDVP sweep, tau: %.5f, SvN: %.4f, maxm: %6d, time: %.5f secs\n",
                    tau, svn, linkdim, t)
                flush(stdout)
            end
        else
            t = @elapsed begin
                psi = tdvp(H, time_tau / n_steps, psi;
                    updater_backend=solver_backend,
                    nsweeps=1,
                    nsite=1,
                    cutoff=cutoff,
                    maxdim=maxm,
                    normalize=false)

                log_norm += 2.0*log(norm(psi))
                normalize!(psi)
            end
            svn = entropy_von_neumann(psi, N ÷ 2)
            linkdim = maxlinkdim(psi)
            if !silent
                @printf("    1TDVP sweep, tau: %.5f, SvN: %.4f, maxm: %6d, time: %.5f secs\n",
                    tau, svn, linkdim, t)
                flush(stdout)
            end
        end


        GC.gc()

    end
    if !isapprox(remainder, 0; rtol=1e-6, atol=1e-6)
        if maxlinkdim(psi) < maxm
            t = @elapsed begin
                psi = tdvp(H, remainder, psi;
                    updater_backend=solver_backend,
                    nsweeps=1,
                    nsite=2,
                    cutoff=cutoff,
                    maxdim=maxm,
                    normalize=false)
                
                log_norm += 2.0*log(norm(psi))
                normalize!(psi)
            end
            svn = entropy_von_neumann(psi, N ÷ 2)
            linkdim = maxlinkdim(psi)
            if !silent
                @printf("    2TDVP sweep, tau: %.5f, SvN: %.4f, maxm: %6d, time: %.5f secs\n",
                    abs(remainder), svn, linkdim, t)
                flush(stdout)
            end
        else
            t = @elapsed begin
                psi = tdvp(H, remainder, psi;
                    updater_backend=solver_backend,
                    nsweeps=1,
                    nsite=1,
                    cutoff=cutoff,
                    maxdim=maxm,
                    normalize=false)

                log_norm += 2.0*log(norm(psi))
                normalize!(psi)
            end
            svn = entropy_von_neumann(psi, N ÷ 2)
            linkdim = maxlinkdim(psi)
            if !silent
                @printf("    1TDVP sweep, tau: %.5f, SvN: %.4f, maxm: %6d, time: %.5f secs\n",
                    abs(remainder), svn, linkdim, t)
                flush(stdout)
            end
        end

        GC.gc()
    end

    return psi, log_norm
end


"""
Time evolution using TDVP with intitial basis extension

starting out with 2TDVP until a maximal bond dimension is achieved
then it switches to 1TDVP
a commensurate step size is automatically chosen

Returns the evolved state, always normalized, and `log_norm = log ‖exp(time H) psi0‖²`
for a normalized `psi0`; for METTS with `time = -β/2` and a product state
`|σ⟩` this is `log ⟨σ|exp(-βH)|σ⟩`. Weight lost to truncation is included.
"""
function timeevo_tdvp_extend(H::MPO, psi0::MPS, time::Number;
    tau::Number=0.1, cutoff::Float64=1e-6,
    maxm::Int64=1000, tau0::Float64=0.05,
    nsubdiv::Int64=4, kkrylov::Int64=3,
    silent=false,
    solver_backend::AbstractString="applyexp")

    iszero(time) && return copy(psi0), 0.0
    log_norm = 0.0

    N = length(psi0)
    tau_init = min(tau0, abs(time))
    time_init = time / abs(time) * tau_init

    # distribute intial times logarithmically
    times_init = [time_init / (2^(nsubdiv - 1))]
    for exp in (nsubdiv-1):-1:1
        append!(times_init, time_init / (2^exp))
    end

    psi = copy(psi0)

    for isub in 1:nsubdiv
        l1 = maxlinkdim(psi)
        t = @elapsed begin
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
                l1, l2, t)
        end
        t = @elapsed begin
            psi = tdvp(H, times_init[isub], psi;
                updater_backend=solver_backend,
                nsweeps=1,
                nsite=1,
                cutoff=cutoff,
                maxdim=maxm,
                normalize=false)
            
            log_norm += 2.0*log(norm(psi))
            normalize!(psi)
        end

        svn = entropy_von_neumann(psi, N ÷ 2)
        linkdim = maxlinkdim(psi)
        if !silent
            @printf("    1TDVP sweep, tau: %.5f, SvN: %.4f, maxm: %6d, time: %.5f secs\n",
                abs(times_init[isub]), svn, linkdim, t)
            flush(stdout)
        end
    end
    GC.gc()

    # Perform the bulk time evolution
    time_bulk = time - time_init
    psi, norm_tmp = timeevo_tdvp(H, psi, time_bulk; tau=tau, cutoff=cutoff, maxm=maxm,
        silent=silent, solver_backend=solver_backend)

    log_norm += norm_tmp

    return psi, log_norm
end
