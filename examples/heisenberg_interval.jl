using ITensors
using ITensorMPS
using METTS
using Dumper
using Random

# Parameters can be overridden from the command line:
#   julia heisenberg_interval.jl <beta_collapse> <beta_min> <beta_max> <outfile>
let
    N = 12
    nmetts = 100
    random_seed = 42

    beta_collapse = length(ARGS) >= 1 ? parse(Float64, ARGS[1]) : 4.0
    beta_min      = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : beta_collapse / 2
    beta_max      = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : beta_collapse
    outfile_path  = length(ARGS) >= 4 ? ARGS[4] : "metts_bc$(beta_collapse).h5"

    # Time evolution parameters
    tau    = 0.1
    cutoff = 1e-6
    maxm   = 100
    tau0   = 0.05
    nsubdiv    = 4
    kkrylov    = 3
    normalize  = true
    silent     = false
    solver_backend = "applyexp"
    shift  = 0.

    outfile = DumpFile(outfile_path)

    # Store run metadata once
    outfile["beta_collapse"] = beta_collapse
    outfile["beta_min"]      = beta_min
    outfile["beta_max"]      = beta_max
    outfile["tau"]           = tau

    # Create Heisenberg Chain model (open b.c.)
    ops = OpSum()
    for s in 1:(N-1)
        ops += 0.5, "S+", s, "S-", s+1
        ops += 0.5, "S-", s, "S+", s+1
        ops += "Sz", s, "Sz", s+1
    end

    sites = siteinds("S=1/2", N; conserve_sz=true)
    outfile["local_states"] = local_state_strings(sites[1])

    product_state = random_product_state(sites, random_seed; nup=N÷2)
    psi = MPS(sites, product_state)
    H = MPO(ops, sites)

    measure       = psi -> Dict("energy" => inner(psi', H, psi))
    collapse_func = psi -> collapse_with_qn!(psi, "X")

    for step in 1:nmetts
        psi, log_norm, measurements, product_state, times =
            timeevo_tdvp_extend_measurements(
                H, psi, beta_min/2, beta_max/2, beta_collapse/2, measure, collapse_func;
                tau=tau, cutoff=cutoff, maxm=maxm, tau0=tau0, nsubdiv=nsubdiv,
                kkrylov=kkrylov, normalize=normalize, silent=silent,
                solver_backend=solver_backend, shift=shift)

        if step == 1
            dump!(outfile, "betas", 2 .* times)
        end
        dump!(outfile, "energy",       measurements["energy"])
        dump!(outfile, "log_norm",     measurements["log_norm"])
        dump!(outfile, "product_state", product_state)

        psi = MPS(sites, product_state)
    end
end
