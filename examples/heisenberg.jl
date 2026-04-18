using ITensors
using ITensorMPS
using METTS
using Dumper
using Random

let
    N = 12
    beta = 10
    nmetts = 100
    random_seed = 42

    # Time evolution parameters
    tau=0.2                   # TDVP time step
    cutoff=1e-6               # TDVP cutofff
    maxm=100                  # maximal bond dimension
    tau0=0.05                 # initial timestep for basis extension
    nsubdiv=4                 # subdivisions for basis extension
    kkrylov=3                 # number of krylov steps for basis ext
    normalize=true            # flag if result is normalized
    silent=false              # print some output to console
    solver_backend="applyexp" # efficient exponentiation backend
    shift=0.                  # no shift in energy

    # Create outfile
    outfile = DumpFile("outfile.h5")
        
    # Create Heisenberg Chain model (open b.c.)
    ops = OpSum()
    for s in 1:(N-1)
        ops += 0.5, "S+", s, "S-", s+1
        ops += 0.5, "S-", s, "S+", s+1
        ops += "Sz", s, "Sz", s+1
    end

    # Set up initial state and Hamiltonian MPO
    sites = siteinds("S=1/2", N; conserve_sz=true)
    dump!(outfile, "local_states", local_state_strings(sites[1]))

    product_state = random_product_state(sites, random_seed; nup=N÷2)
    psi = MPS(sites, product_state)
    H = MPO(ops, sites)

    # main METTS loop
    for step in 1:nmetts
        # time evolution
        psi, log_norm = timeevo_tdvp_extend(H, psi, beta/2;
            tau=tau, cutoff=cutoff, maxm=maxm, tau0=tau0, nsubdiv=nsubdiv,
            kkrylov=kkrylov, normalize=normalize, silent=silent,
            solver_backend=solver_backend, shift=shift)

        # measurements
        energy = inner(psi', H, psi)
        dump!(outfile, "energy", energy)

        # Collapse
        dump!(outfile, "product_state", product_state)
        product_state = collapse_with_qn!(psi, "X")
    end
end
