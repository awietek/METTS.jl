using ITensors
using ITensorMPS
using METTS
using LinearAlgebra
using Random
using Printf
using Serialization
using HDF5

if length(ARGS) < 11
    println("Usage: julia --project examples/hubbard_cylinder_metts.jl ",
            "L W U filling t T tau maxD seed nmetts Nwarm")
    exit(1)
end

L       = parse(Int,    ARGS[1])
W       = parse(Int,    ARGS[2])
U       = parse(Float64,ARGS[3])
filling = parse(Float64,ARGS[4])
t       = parse(Float64,ARGS[5])
T       = parse(Float64,ARGS[6])
tau     = parse(Float64,ARGS[7])
maxD    = parse(Int,    ARGS[8])
seed    = parse(Int,    ARGS[9])
nmetts  = parse(Int,    ARGS[10])
Nwarm   = parse(Int,    ARGS[11])

if nmetts <= Nwarm
    @warn "nmetts <= Nwarm -> no measured samples will be stored."
end

# --------------------------- Helpers ----------------------------
# running mean +/- stderr
function avg_err(v::Vector{<:Real})
    N = length(v)
    avg  = sum(v) / N
    avg2 = sum(x->x*x, v) / N
    return avg, sqrt(max(0, avg2 - avg^2)/N)
end


function load_samples(filename::AbstractString)
    if isfile(filename)
        samples = Dict{Int,Vector{Int}}()
        last = 0
        open(filename, "r") do io
            for line in eachline(io)
                step_str, sample_str = split(line, ":")
                step = parse(Int, strip(step_str))
                sample_vec = eval(Meta.parse(strip(sample_str)))
                samples[step] = sample_vec
                last = max(last, step)
            end
        end
        return samples, last
    else
        return Dict{Int,Vector{Int}}(), 0
    end
end

function save_samples(filename::AbstractString, samples::Dict{Int,<:Any})
    open(filename, "w") do io
        for (step, sample) in sort!(collect(samples); by=first)
            write(io, "$(step): $(sample)\n")
        end
    end
end

# Hubbard MPO on an LxW cylinder
function fermi_hubbard_cylinder_mpo(L, W, t, U, sites)
    N = L*W
    lat = square_lattice(L, W; yperiodic=true)
    os = OpSum()
    for b in lat
        os += -t, "Cdagup", b.s1, "Cup", b.s2
        os += -t, "Cdagdn", b.s1, "Cdn", b.s2
        os += -t, "Cdagup", b.s2, "Cup", b.s1
        os += -t, "Cdagdn", b.s2, "Cdn", b.s1
    end
    for n in 1:N
        os += U, "Nupdn", n
    end
    return MPO(os, sites)
end

# HDF5: create file and write metadata once
function h5_init(outfile; L,W,U,filling,t,T,tau,maxD,seed,Nwarm)
    h5open(outfile, "w") do f
        g = create_group(f, "meta")
        g["L"] = L; g["W"] = W; g["U"] = U
        g["filling"] = filling; g["t"] = t; g["T"] = T
        g["tau"] = tau; g["maxD"] = maxD;
        g["seed"] = seed; g["Nwarm"] = Nwarm
    end
end

# HDF5: write one measured step into /step_<k>/*
function h5_write_step(outfile::AbstractString, step::Int; energy, Ntotal, dob, szszs, spsms, smsps, ntnt, SS)
    h5open(outfile, "r+") do f
        gname = @sprintf("step_%06d", step)
        g = haskey(f, gname) ? f[gname] : create_group(f, gname)
        g["energy"]    = energy
        g["Ntotal"]    = Ntotal
        g["DO"]        = dob
        g["SZSZ"]      = szszs
        g["SPSM"]      = spsms
        g["SMSP"]      = smsps
        g["NtotNtot"]  = ntnt
        g["SS"]        = SS
    end
end

# --------------------------- Main -------------------------------
function main()

    outroot = get(ENV, "METTS_OUTDIR", joinpath(pwd(), "outfiles.metts"))
    outdir  = @sprintf("%s/L.%d.W.%d/U.%.3f/filling.%.5f/T.%.5f/D.%d/tau.%.2f.seed.%d",
                       outroot, L, W, U, filling, T, maxD, tau, seed)
    mkpath(outdir)
    outfile = joinpath(outdir, "outfile.h5")
    smpfile = joinpath(outdir, "samples.txt")

    println("=== METTS Hubbard cylinder ===")
    println("LxW = $(L)x$(W)  (N=$(L*W))  U=$U  t=$t")
    println("filling = $filling  (doping=$(1.0 - filling))  T=$T (beta=$(1/T))")
    println("tau = $tau  maxD = $maxD")
    println("seed = $seed  nmetts = $nmetts  Nwarm = $Nwarm")
    println("outdir = $outdir")

    Random.seed!(seed)

    N     = L*W
    beta  = 1/T
    sites = siteinds("Electron", N; conserve_qns=true)
    H     = fermi_hubbard_cylinder_mpo(L, W, t, U, sites)

    # init HDF5 with metadata
    h5_init(outfile; L,W,U,filling,t,T,tau,maxD,seed,Nwarm)

    # resume if samples exist
    samples, start_step = load_samples(smpfile)
    println("start_step = $start_step")

    # initial product state
    decode_occ(x) = x == 1 ? "Emp" : x == 2 ? "Up" : x == 3 ? "Dn" : "UpDn"
    state = Vector{String}(undef, N)
    if start_step > 0
        last = samples[start_step]
        @assert length(last) == N
        @inbounds for i in 1:N
            state[i] = decode_occ(last[i])
        end
    else
        Ne = round(Int, filling*N)
        @assert iseven(Ne) "Ne must be even for equal up/down"
        perm = Random.randperm(N)
        up   = perm[1:Ne÷2]
        dn   = setdiff(perm, up)[1:Ne÷2]
        fill!(state, "Emp")
        @inbounds for i in up; state[i] = "Up"; end
        @inbounds for i in dn
            state[i] = (state[i] == "Up") ? "UpDn" : "Dn"
        end
    end
    psi = randomMPS(sites, state)

    # DMRG warmup only if no resume
    if start_step == 0
        println("DMRG warmup ...")
        nsweeps = 10
        maxdim  = [10, 20, 20, 50, 50, 100, 100, 100, 200, 200]
        cutoff = 1e-6
        energy_dmrg, psi_dmrg = dmrg(H, psi; nsweeps, maxdim, cutoff=cutoff)
        println("Sampling post-DMRG wavefunction ...")
        samp0 = collapse_with_qn!(psi_dmrg, "Z")
        new_state = [decode_occ(samp0[j]) for j in 1:N]
        psi = randomMPS(sites, new_state)
        samples[0] = samp0; save_samples(smpfile, samples)
    else
        println("Skipping warmup (resumed run).")
    end

    # TDVP (imaginary time beta/2)
    function evo!(psi_in, beta_in)
        if beta_in > 1e-14
            return timeevo_tdvp_extend(H, psi_in, -(beta_in/2);
                                       tau=tau, normalize=true,
                                       solver_backend="applyexp",
                                       maxm=maxD, kkrylov=2, tau0=0.02, nsubdiv=2)
        else
            return psi_in
        end
    end

    energies = Float64[]
    first_step = start_step + 1
    for step in first_step:nmetts
        println("\n--- METTS step $step ---")
        if step <= Nwarm
            println("Warm-up step (no storage).")
            psi = evo!(psi, beta)
            samp = collapse_with_qn!(psi, "X")
            psi = randomMPS(sites, [decode_occ(samp[j]) for j in 1:N])
        else
            t_iter = @elapsed begin
                println("Evolving to beta ...")
                psi = evo!(psi, beta)

                println("Measuring ...")
                energy = inner(psi', H, psi)
                Ntotal = expect(psi, "Ntot")
                dob    = expect(psi, "Nupdn")
                szszs  = correlation_matrix(psi, "Sz","Sz"; ishermitian=true)
                spsms  = correlation_matrix(psi, "S+","S-"; ishermitian=false)
                smsps  = correlation_matrix(psi, "S-","S+"; ishermitian=false)
                ntnt   = correlation_matrix(psi, "Ntot","Ntot"; ishermitian=true)
                SS     = szszs + 0.5*(smsps + spsms)

                push!(energies, energy)
                aE, eE = avg_err(energies)
                @printf("Energy[%d] = %.6f   running mean(E) = %.6f +/- %.6f\n", step, energy, aE, eE)
                dens = sum(Ntotal)/(L*W)
                @printf("density = %.6f\n", dens)

                # HDF5 write for this step
                h5_write_step(outfile, step;
                              energy=energy, Ntotal=Ntotal, dob=dob,
                              szszs=szszs, spsms=spsms, smsps=smsps,
                              ntnt=ntnt, SS=SS)
            end
            @printf("Step %d meas+dump time: %.3f s\n", step, t_iter)

            # Collapse, save sample, reinit
            samp = collapse_with_qn!(psi, "X")
            samples[step - Nwarm] = samp
            save_samples(smpfile, samples)
            psi = randomMPS(sites, [decode_occ(samp[j]) for j in 1:N])
            GC.gc()
        end
    end

    println("Done. Output: $outfile ; samples: $smpfile")
end

main()
