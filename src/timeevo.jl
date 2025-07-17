using Printf

""" computes the number of steps, and the remainder of the time step if not commensurate """
function n_steps_remainder(time::Number, tau::Real)
    n_steps = floor(abs(time) / tau)
    time_tau = n_steps * tau * time / abs(time)
    remainder = time - time_tau
    return n_steps, remainder, time_tau
end

function entropy_von_neumann(psi::MPS, b::Int)
    nrm = norm(psi)
    s = siteinds(psi)  
    orthogonalize!(psi, b)
    _,S = svd(psi[b], (linkind(psi, b-1), s[b]))
    SvN = 0.0
    for n in 1:dim(S, 1)
        p = (S[n,n] / nrm)^2
        SvN -= p * log(p)
    end
    return SvN
end

"""
Basic time evolution using TDVP

starting out with 2TDVP until a maximal bond dimension is achieved
then it switches to 1TDVP
a commensurate step size is automatically chosen
"""
function timeevo_tdvp(H::MPO, psi0::MPS, time::Number;
                      tau::Number=0.1, cutoff::Float64=1e-6,
                      maxm::Int64=1000, normalize::Bool=true,
                      silent=false, solver_backend::AbstractString="applyexp",
                      shift::Real = 0.)
    N = length(psi0)
    n_steps, remainder, time_tau = n_steps_remainder(time, tau)
    psi = copy(psi0)

    for step in 1:n_steps        
        if maxlinkdim(psi) < maxm
            t = @elapsed begin
                # psi = tdvp(H, psi, time_tau / n_steps;
                #            nsweeps=1, nsite=2, cutoff=cutoff,
                #            maxdim=maxm, normalize=normalize,
                #            solver_backend=solver_backend, shift=shift)

                psi = tdvp(H, time_tau / n_steps, psi;
                           updater_backend=solver_backend,
                           nsweeps=1,
                           nsite=2,
                           cutoff=cutoff,
                           maxdim=maxm,
                           normalize=normalize)

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
                # psi = tdvp(H, psi, time_tau / n_steps;
                #            nsweeps=1, nsite=1, cutoff=cutoff,
                #            maxdim=maxm, normalize=normalize,
                #            solver_backend=solver_backend, shift=shift)
                psi = tdvp(H, time_tau / n_steps, psi;
                           updater_backend=solver_backend,
                           nsweeps=1,
                           nsite=1,
                           cutoff=cutoff,
                           maxdim=maxm,
                           normalize=normalize)
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
                # psi = tdvp(H, psi, remainder;
                #            nsweeps=1, nsite=2, cutoff=cutoff,
                #            maxdim=maxm, normalize=normalize,
                #            solver_backend=solver_backend, shift=shift)

                psi = tdvp(H, remainder, psi;
                           updater_backend=solver_backend,
                           nsweeps=1,
                           nsite=2,
                           cutoff=cutoff,
                           maxdim=maxm,
                           normalize=normalize)
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
                # psi = tdvp(H, psi, remainder;
                #            nsweeps=1, nsite=1, cutoff=cutoff,
                #            maxdim=maxm, normalize=normalize,
                #            solver_backend=solver_backend, shift=shift)
                
                psi = tdvp(H, remainder, psi;
                           updater_backend=solver_backend,
                           nsweeps=1,
                           nsite=1,
                           cutoff=cutoff,
                           maxdim=maxm,
                           normalize=normalize)
            end
            if !silent
                @printf("    1TDVP sweep, tau: %.5f, SvN: %.4f, maxm: %6d, time: %.5f secs\n",
                        abs(remainder), svn, linkdim, t)
                flush(stdout)
            end
        end

        GC.gc()
    end
    
    return psi
end


"""
Time evolution using TDVP with intitial basis extension

starting out with 2TDVP until a maximal bond dimension is achieved
then it switches to 1TDVP
a commensurate step size is automatically chosen
"""
function timeevo_tdvp_extend(H::MPO, psi0::MPS, time::Number;
                             tau::Number=0.1, cutoff::Float64=1e-6,
                             maxm::Int64=1000, tau0::Float64=0.05,
                             nsubdiv::Int64=4, kkrylov::Int64=3,
                             normalize::Bool=true, silent=false,
                             solver_backend::AbstractString="applyexp",
                             shift::Real = 0.)

    N = length(psi0)
    tau_init = min(tau0, abs(time))
    time_init = time / abs(time) * tau_init

    # distribute intial times logarithmically
    times_init = [time_init / (2 ^ (nsubdiv - 1))]
    for exp in (nsubdiv-1):-1:1
        append!(times_init, time_init / (2 ^ exp))
    end

    psi = copy(psi0)

    for isub in 1:nsubdiv
        l1 = maxlinkdim(psi)
        t = @elapsed begin
            println("before basis_extend")
            flush(stdout)
            psi = basis_extend(psi, H; extension_krylovdim=kkrylov,
                               extension_cutoff=1e-12)
            println("after basis_extend")
            flush(stdout)
        end
        l2 = maxlinkdim(psi)
        if !silent
        @printf("    Basis extension, maxm: %4d -> %4d, time: %.5f secs\n",
                l1, l2, t)
        end
        t = @elapsed begin
            # psi = tdvp(H, psi, times_init[isub];
            #            nsweeps=1, nsite=1, cutoff=1e-16,
            #            maxdim=maxm, normalize=normalize,
            #            solver_backend=solver_backend, shift=shift)
            psi = tdvp(H, times_init[isub], psi;
                       updater_backend=solver_backend,
                       nsweeps=1,
                       nsite=1,
                       cutoff=cutoff,
                       maxdim=maxm,
                       normalize=normalize)
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
    psi = timeevo_tdvp(H, psi, time_bulk; tau=tau, cutoff=cutoff, maxm=maxm,
                        normalize=normalize, silent=silent,
                        solver_backend=solver_backend, shift=shift)   


    return psi
end