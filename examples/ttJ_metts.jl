using ITensors
using ITensorCorrelators
using Printf
using METTS
using LinearAlgebra
using Random
using Serialization
using HDF5
using Dumper 



function parse_args(args::Vector{String})
    if length(args) < 13
        error("""
        Usage:
          julia ttJ_metts.jl L W J filling t t_prime T tau maxD cutoff seed NMETTS Nwarm
        """)
    end
    L       = parse(Int64,   args[1])
    W       = parse(Int64,   args[2])
    J       = parse(Float64, args[3])
    filling = parse(Float64, args[4])
    t       = parse(Float64, args[5])
    t_prime = parse(Float64, args[6])
    T       = parse(Float64, args[7])
    tau     = parse(Float64, args[8])
    maxD    = parse(Int64,   args[9])
    cutoff  = parse(Float64, args[10])
    seed    = parse(Int64,   args[11])
    nmetts  = parse(Int64,   args[12])  
    Nwarm   = parse(Int64,   args[13])
    return L, W, J, filling, t, t_prime, T, tau, maxD, cutoff, seed, nmetts, Nwarm
end

L, W, J, filling, t, t_prime, T, tau, maxD, cutoff, seed, nmetts, Nwarm = parse_args(ARGS)



function configure_threading()

    ITensors.Strided.set_num_threads(1)
    BLAS.set_num_threads(1)
    ITensors.enable_threaded_blocksparse()

    println("=== Thread configuration ===")
    println("Julia threads                = $(Threads.nthreads())")
    println("ITensors Strided threads     = 1")
    println("BLAS threads                 = $(BLAS.get_num_threads())")
    println("Threaded block-sparse        = ENABLED")
    println("============================")
end

###############################################################################
# Helper Functions
###############################################################################

function avg_err(v::Vector{Float64})
    N = length(v)
    avg  = sum(v) / N
    avg2 = sum(v.^2) / N
    return avg, sqrt(max(avg2 - avg^2, 0.0) / N)
end


function parse_int_list(s::AbstractString)
    out = Int[]
    for m in eachmatch(r"-?\d+", s)
        push!(out, parse(Int, m.match))
    end
    return out
end

function load_samples(filename::AbstractString)
    if !isfile(filename)
        return Dict{Int,Vector{Int}}(), 0
    end
    samples = Dict{Int,Vector{Int}}()
    last_key = 0
    open(filename, "r") do file
        for line in eachline(file)
            line = strip(line)
            isempty(line) && continue
            parts = split(line, ":"; limit=2)
            length(parts) != 2 && continue
            key = parse(Int, strip(parts[1]))
            vec = parse_int_list(strip(parts[2]))
            isempty(vec) && continue
            samples[key] = vec
            last_key = max(last_key, key)
        end
    end
    return samples, last_key
end

function save_samples(filename::AbstractString, samples::Dict{Int,Vector{Int}})
    open(filename, "w") do file
        for (k, v) in sort(collect(samples); by = x -> x[1])
            write(file, "$k: $(v)\n")
        end
    end
end

@inline occ_to_label(x::Int) = (x == 1 ? "Emp" : (x == 2 ? "Up" : "Dn"))


function make_product_mps(sites, state::Vector{String})
 
    try
        return productMPS(sites, state)
    catch
        try
            return MPS(sites, state)
        catch
            return randomMPS(sites, state)
        end
    end
end

###############################################################################
# Extended t-J Hamiltonian
###############################################################################

function tJ_cylinder_mpo(L, W, t, t_prime, J, sites)
    lattice = square_lattice(L, W; yperiodic = true)
    os = OpSum()

    # Nearest-neighbor terms
    for b in lattice
        # Hopping
        os += -t, "Cdagup", b.s1, "Cup", b.s2
        os += -t, "Cdagdn", b.s1, "Cdn", b.s2
        os += -t, "Cdagup", b.s2, "Cup", b.s1
        os += -t, "Cdagdn", b.s2, "Cdn", b.s1

        # Heisenberg J
        os += 0.5 * J, "S+",   b.s1, "S-",   b.s2
        os += 0.5 * J, "S-",   b.s1, "S+",   b.s2
        os += J,       "Sz",   b.s1, "Sz",   b.s2
        os += -0.25*J, "Ntot", b.s1, "Ntot", b.s2
    end

    # Next-nearest-neighbor hopping (t')
    for i in 1:(L-1)
        for j in 1:W
            s1 = (i-1)*W + j

            # right-diagonal neighbor
            s2_r = i*W + mod(j, W) + 1
            os += -t_prime, "Cdagup", s1,   "Cup", s2_r
            os += -t_prime, "Cdagdn", s1,   "Cdn", s2_r
            os += -t_prime, "Cdagup", s2_r, "Cup", s1
            os += -t_prime, "Cdagdn", s2_r, "Cdn", s1

            # left-diagonal neighbor
            s2_l = i*W + mod(j-2, W) + 1
            os += -t_prime, "Cdagup", s1,   "Cup", s2_l
            os += -t_prime, "Cdagdn", s1,   "Cdn", s2_l
            os += -t_prime, "Cdagup", s2_l, "Cup", s1
            os += -t_prime, "Cdagdn", s2_l, "Cdn", s1
        end
    end

    return MPO(os, sites)
end

###############################################################################
# Four-Fermionic Correlators
###############################################################################

function precompute_all_nn_indices(L, W)
    lattice = square_lattice(L, W; yperiodic = true)
    valid_combinations = NTuple{4,Int}[]
    for b1 in lattice
        for b2 in lattice
            push!(valid_combinations, (b1.s1, b1.s2, b2.s1, b2.s2))
        end
    end
    return valid_combinations
end

function compute_four_fermionic_correlators(psi::MPS, valid_indices, outfile::DumpFile, step_total::Int, Nwarm::Int)
    isempty(valid_indices) && return

    correlator_types = [("Cdagup", "Cdagdn", "Cup", "Cdn")]
    results = Dict{Tuple{NTuple{4,String}, NTuple{4,Int}}, Any}()

    for ctype in correlator_types
        for indices in valid_indices
            perms = [
                indices,
                (indices[2], indices[1], indices[3], indices[4]),
                (indices[1], indices[2], indices[4], indices[3]),
                (indices[2], indices[1], indices[4], indices[3])
            ]
            for perm in perms
                C = correlator(psi, ctype, [perm])
                if haskey(C, perm)
                    int_indices = ntuple(i -> Int(perm[i]), 4)
                    results[(ctype, int_indices)] = C[perm]
                end
            end
        end
    end

    meas = step_total - Nwarm

    # Save results to HDF5
    h5open(outfile.filename, "r+") do h5file
        for ((ctype, i4), val) in results
            corr_type_str = join(ctype, "")
            corr_folder   = joinpath("fourfermioncorr", corr_type_str, "step_$(meas)")
            if !haskey(h5file, corr_folder)
                HDF5.create_group(h5file, corr_folder)
            end
            idx_str  = "idx_$(i4[1])_$(i4[2])_$(i4[3])_$(i4[4])"
            corr_key = joinpath(corr_folder, idx_str)
            if !haskey(h5file, corr_key)
                dump!(outfile, corr_key, val)
            end
        end
    end
end


@inline function unwrap_tdvp_result(res)
    return res isa Tuple ? res[1] : res
end

###############################################################################
# Main Simulation Function
###############################################################################

function main(; L::Int, W::Int, J::Float64, filling::Float64, t::Float64, t_prime::Float64,
               T::Float64, tau::Float64, maxD::Int, cutoff::Float64, seed::Int,
               NMETTS::Int, Nwarm::Int)

    configure_threading()

    filename = @sprintf(
        "/metts.cylinder.tj.tjp/outfiles.metts/L.%d.W.%d/J.%.4f/t.%.0f/t_prime.%.4f/filling.%.5f/T.%.5f/D.%.0f/tau.%.2f.cutoff.%.1e.seed.%d/outfile.h5",
        L, W, J, t, t_prime, filling, T, maxD, tau, cutoff, seed
    )

    println("=== METTS Script ===")
    println("Geometry: Square Lattice Cylinder, L=$L, W=$W => Nsites=$(L*W)")
    println("Filling: $filling => doping δ=$(1.0 - filling)")
    println("Temperature: $T => beta=$(1/T)")
    println("Nwarm = $Nwarm, NMETTS(total incl warmup) = $NMETTS")
    println("Using random seed = $seed")
    println("Output file: $filename")

    mkpath(dirname(filename))
    outfile = DumpFile(filename)

    Random.seed!(seed)

    beta = 1.0 / T
    N    = L * W
    println("beta = $beta, N = $N")

    # Build Hamiltonian
    sites = siteinds("tJ", N; conserve_qns = true)
    H     = tJ_cylinder_mpo(L, W, t, t_prime, J, sites)

   
    samples_filename = joinpath(dirname(filename), "samples.txt")
    samples, last_meas = load_samples(samples_filename)
    @show last_meas

    last_total = (last_meas > 0) ? (Nwarm + last_meas) : 0
    first_step_total = last_total + 1
    println("Restart: last_total_step = $last_total, starting at step_total = $first_step_total")

    # Precompute indices for four-fermion correlators
    valid_indices = precompute_all_nn_indices(L, W)

    ############################################################################
    # 1) INITIAL PRODUCT STATE
    ############################################################################

    state = fill("Emp", N)

    if last_meas > 0
        last_sample = samples[last_meas]
        length(last_sample) == N || error("Saved sample length $(length(last_sample)) != N=$N")
        for i in 1:N
            state[i] = occ_to_label(last_sample[i])
        end
        println("Initialized from saved sample meas=$last_meas (total step=$last_total).")
    else
        
        Ne = round(Int, filling * N)
        isodd(Ne) && error("Ne must be even for equal # up and down electrons; got Ne=$Ne")

        perm = Random.randperm(N)
        up_indices   = perm[1:(Ne ÷ 2)]
        down_indices = perm[(Ne ÷ 2 + 1):Ne]

        for idx in up_indices
            state[idx] = "Up"
        end
        for idx in down_indices
            state[idx] = "Dn"
        end
        println("Initialized fresh random product state at target filling.")
    end

    psi = make_product_mps(sites, state)

    ############################################################################
    # 2) DMRG WARMUP (only if starting from scratch)
    ############################################################################

    if last_total == 0
        local maxdim::Vector{Int}
        local nsweeps::Int

        if T < 0.055
            println("Performing 60-sweep DMRG warmup for T < 0.055 ...")
            maxdim = [
                10, 10, 20, 20, 50, 50, 50, 50,
                100, 100, 100, 100, 100, 100,
                200, 200, 200, 200, 200,
                500, 500, 500, 500, 500, 500,
                1000, 1000, 1000, 1000, 1000,
                1500, 1500, 1500, 1500, 1500, 1500, 1500, 1500,
            ]
            nsweeps = 60
        elseif 0.055 <= T <= 0.10
            println("Performing 30-sweep DMRG warmup for 0.055 ≤ T ≤ 0.10 ...")
            maxdim = [
                10, 10, 20, 20, 50, 50, 50, 50,
                100, 100, 100, 100,
                200, 200, 200, 200,
                500, 500, 500, 500, 500,
                1000
            ]
            nsweeps = 30
        else
            println("Performing ~15-sweep DMRG warmup for T > 0.10 ...")
            maxdim  = [10, 10, 50, 50, 50, 50, 100, 100, 100, 100, 100, 500, 500, 500]
            nsweeps = 20
        end

        cutoff_vec = [cutoff]
        energy_dmrg, psi = dmrg(H, psi; nsweeps = nsweeps, maxdim = maxdim, cutoff = cutoff_vec)
        println("DMRG warmup energy = $energy_dmrg")

        println("Sampling the post-DMRG wavefunction ...")
        samp0 = collapse_with_qn!(psi, "Z")
        @show samp0

        new_state = [occ_to_label(samp0[j]) for j in 1:N]
        psi = make_product_mps(sites, new_state)
    else
        println("Skipping DMRG warmup since restart was detected.")
    end

    ############################################################################
    # 3) METTS Loop
    ############################################################################

    energies = Float64[]

    function evolve_beta_half(psi_in::MPS, beta::Float64, H)
        if beta <= 1e-14
            return psi_in
        end
        res = timeevo_tdvp_extend(
            H, psi_in, -(beta/2);
            tau            = tau,
            normalize      = true,
            solver_backend = "applyexp",
            maxm           = maxD,
            kkrylov        = 2,
            tau0           = 0.02,
            nsubdiv        = 2,
        )
        return unwrap_tdvp_result(res)
    end

    for step_total in first_step_total:NMETTS
        println("\n--- METTS iteration (total) = $step_total ---")

        if step_total <= Nwarm
            println("  step_total <= Nwarm => DRY run (no storing).")
            psi = evolve_beta_half(psi, beta, H)
            @assert psi isa MPS "psi became $(typeof(psi)) after TDVP; expected MPS."

            samp = collapse_with_qn!(psi, "X")
            @show samp

            new_state = [occ_to_label(samp[j]) for j in 1:N]
            psi = make_product_mps(sites, new_state)

        else
            meas = step_total - Nwarm
            metts_time = @elapsed begin
                println("  step_total > Nwarm => evolve, measure & store data... (meas=$meas)")
                psi = evolve_beta_half(psi, beta, H)
                @assert psi isa MPS "psi became $(typeof(psi)) after TDVP; expected MPS."

                # Energy
                println("measure energy")
                @time energy = real(ITensors.inner(psi', H, psi))
                push!(energies, energy)
                @printf("  Energy of METTS (meas=%d, total=%d) = %.10f\n", meas, step_total, energy)
                a_E, err_E = avg_err(energies)
                @printf("  Running mean (this run) = %.10f ± %.10f  [%.10f, %.10f]\n",
                        a_E, err_E, a_E - err_E, a_E + err_E)

                # Entanglement entropy
                svn = entropy_von_neumann(psi, div(N, 2))

                # Local observables
                println("measure Sz")
                @time szs = expect(psi, "Sz")

                println("measure Ntot")
                @time Ntotal = expect(psi, "Ntot")
                density = sum(Ntotal) / (L * W)
                println("Average density per site: $density")

                # Correlation matrices
                println("measure SzSz")
                @time szszs = correlation_matrix(psi, "Sz", "Sz"; ishermitian = true)

                println("measure SpSm")
                @time spsms = correlation_matrix(psi, "S+", "S-"; ishermitian = false)

                println("measure SmSp")
                @time smsps = correlation_matrix(psi, "S-", "S+"; ishermitian = false)

                println("measure NtotNtot")
                @time ntotntot = correlation_matrix(psi, "Ntot", "Ntot"; ishermitian = true)

                println("measure CdagupCup")
                @time cdagupcup = correlation_matrix(psi, "Cdagup", "Cup"; ishermitian = true)

                println("measure CdagdnCdn")
                @time cdagdncdn = correlation_matrix(psi, "Cdagdn", "Cdn"; ishermitian = true)

                println("measure Total S.S correlations")
                @time SS = szszs .+ 0.5 .* (smsps .+ spsms)

             
                dump!(outfile, "energy",    energy)
                dump!(outfile, "SZ",        szs)
                dump!(outfile, "Ntotal",    Ntotal)
                dump!(outfile, "SZSZ",      szszs)
                dump!(outfile, "SPSM",      spsms)
                dump!(outfile, "SMSP",      smsps)
                dump!(outfile, "NtotNtot",  ntotntot)
                dump!(outfile, "SS",        SS)
                dump!(outfile, "cdagupcup", cdagupcup)
                dump!(outfile, "cdagdncdn", cdagdncdn)
                dump!(outfile, "svn",       svn)
            end
            println("Time for METTS iteration with two-body observables (total=$step_total): $metts_time seconds")

            # Four-fermion correlators schedule (same logic; meas = step_total - Nwarm)
            if T <= 0.04 && mod(meas, 1) == 0
                correlator_time = @elapsed compute_four_fermionic_correlators(psi, valid_indices, outfile, step_total, Nwarm)
                println("Time for four-fermionic correlators at meas=$meas: $correlator_time seconds")
                GC.gc()
            elseif T > 0.04 && T <= 0.08 && mod(meas, 2) == 0
                correlator_time = @elapsed compute_four_fermionic_correlators(psi, valid_indices, outfile, step_total, Nwarm)
                println("Time for four-fermionic correlators at meas=$meas: $correlator_time seconds")
                GC.gc()
            elseif T > 0.08 && T <= 0.2 && mod(meas, 3) == 0
                correlator_time = @elapsed compute_four_fermionic_correlators(psi, valid_indices, outfile, step_total, Nwarm)
                println("Time for four-fermionic correlators at meas=$meas: $correlator_time seconds")
                GC.gc()
            elseif T > 0.2 && mod(meas, 5) == 0
                correlator_time = @elapsed compute_four_fermionic_correlators(psi, valid_indices, outfile, step_total, Nwarm)
                println("Time for four-fermionic correlators at meas=$meas: $correlator_time seconds")
                GC.gc()
            end

            # Collapse + store sample under measurement index meas
            samp = collapse_with_qn!(psi, "X")
            @show samp

            samples[meas] = samp
            save_samples(samples_filename, samples)

            # Re-init MPS from collapsed product state
            new_state = [occ_to_label(samp[j]) for j in 1:N]
            psi = make_product_mps(sites, new_state)
        end
    end

    println("All METTS steps done.")
    return nothing
end


main(L=L, W=W, J=J, filling=filling, t=t, t_prime=t_prime,
     T=T, tau=tau, maxD=maxD, cutoff=cutoff, seed=seed,
     NMETTS=nmetts, Nwarm=Nwarm)