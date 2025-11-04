@doc raw"""
metts_single_temperature(measure_func, psi0, beta_sample, H, τ, cutoff, maxBond, tau0, nsubdiv; args...)

Performs the imaginary-time evolution of a single product state up to the desired inverse temperature.
The resulting METTS is used to perform measurments on observables defined by the user. 

Arguments
---------
measure_func::Function
    A user-defined function that measures observables on the evolved MPS state. 
    It should accept at least two arguments: a `MPS` (the evolved state) and an `MPO` (the Hamiltonian). 
    Additional positional arguments can be passed via `args...`.

psi0::MPS
    The initial matrix product state from which the imaginary time evolution starts.

beta_sample::Real
    The inverse temperature for which the METTS sampling is performed.

H::MPO
    The Hamiltonian of the system, represented as a MPO.

τ::Real
    The imaginary time step used for evolving the state.

cutoff::Real
    The truncation cutoff for singular value decompositions during MPS manipulations.

maxBond::Int
    Maximum allowed bond dimension for the MPS during evolution.

tau0::Real
    Initial time step for subdividing the imaginary time evolution.

nsubdiv::Int
    Number of subdivisions used.

args...
    Optional extra positional arguments forwarded to `measure_func`.

Returns
-------
result
    The output of `measure_func` applied to the METTS-evolved state.

phi::MPS
    The MPS after imaginary-time evolution to the target inverse temperature.

Example
-------
measure_energy(mps, H) = expect(H, mps)
result, phi = metts_single_temperature(measure_energy, psi0, 1.0, H, 0.1, 1e-6, 512, 0.1, 1)
This example measures the energy of the system at inverse temperature β = 1.0.
"""
function metts_single_temperature(measure_func::Function, psi0::MPS, beta_sample::Real, H::MPO, τ::Real, cutoff::Real, maxBond::Int, tau0::Real, nsubdiv::Int; args...)

    ### Time evolution from product state
    phi, norm_1 = timeevo_tdvp_extend(H, psi0, -beta_sample / 2.0;  # Note: current implementation may not allow certain options
        tau=τ,
        cutoff=cutoff,
        normalize=true,
        maxm=maxBond,
        tau0=tau0,
        nsubdiv=nsubdiv,
        solver_backend="applyexp",
        shift=0.0)

    return measure_func(phi, H, args...), phi
end


@doc raw"""
metts_interval_temperature(measure_func, psi0, beta_sample, beta_array, H, τ, cutoff, maxBond, tau0, nsubdiv; args...)

Performs the imaginary-time evolution of a single product state along a specified inverse temperature array 
and measures observables at each of this temperatures.

Arguments
---------
measure_func::Function
    A user-defined function that measures observables on the evolved MPS state. 
    It must accept at least two arguments: a `MPS` (the evolved state) and an `MPO` (the Hamiltonian). 
    Additional positional arguments can be passed via `args...`.
    It must return a Dictionary containing the computed observables.

psi0::MPS
    The initial matrix product state from which the imaginary time evolution starts.

beta_sample::Real
    The target inverse temperature for the METTS collapse. Must appear exactly once in `beta_array`.

beta_array::Array{<:Real}
    An array of inverse temperatures defining the time evolution steps.

H::MPO
    The Hamiltonian of the system, represented as a matrix product operator (MPO).

τ::Real
    The imaginary time step used for evolving the state.

cutoff::Real
    The truncation cutoff for singular value decompositions during MPS manipulations.

maxBond::Int
    Maximum allowed bond dimension for the MPS during evolution.

tau0::Real
    Initial time step for subdividing the imaginary time evolution (used in TDVP extension).

nsubdiv::Int
    Number of subdivisions used in the Trotter or other evolution schemes.

args...
    Optional extra positional arguments forwarded to `measure_func`.

Returns
-------
logweight::Vector{Float64}
    Logarithm of the normalization weights accumulated at each time step. Used to avoid numerical overflow.

observables::Vector{Dict}
    A vector of observables measured by `measure_func` at each time step.

psi_collapse::MPS
    The MPS at the target inverse temperature `beta_sample` which can be used for collapse in METTS.

Example
-------
measure_energy(mps, H) = expect(H, mps)
logw, observables, psi_c = metts_interval_temperature(
    measure_energy, psi0, 1.0, collect(2.0:0.1:3.0), H, 0.01, 1e-10, 100, 0.01, 10
)
This example evolves the system along inverse temperatures from 2.0 to 3.0 and measures the energy at each step.
"""
function metts_interval_temperatures(measure_func::Function, psi0::MPS, beta_sample::Real, beta_array::Array, H::MPO, τ::Real, cutoff::Real, maxBond::Int, tau0::Real, nsubdiv::Int; args...)

    #### create array for measurements
    logweight = zeros(Float64, length(beta_array))
    observables = Dict[]

    #### get index of beta_sample
    collapse_idx = findall(x -> x == beta_sample, beta_array)

    @assert length(collapse_idx) == 1  "beta_array should contain beta_sample exactly once"
    
    collapse_idx = collapse_idx[1]

    #### compute time intervals
    diffs_in_beta = diff(beta_array)
    time_intervals = vcat(beta_array[1], diffs_in_beta) ./ 2.0

    ############ First time step ####################
    psi0, norm_1 = timeevo_tdvp_extend(H, psi0, -time_intervals[1];
        tau=τ,
        cutoff=cutoff,
        normalize=true,
        maxm=maxBond,
        tau0=tau0,
        nsubdiv=nsubdiv,
        solver_backend="applyexp", shift=0.0)

    logweight[1] = norm_1
    normalize!(psi0)
    push!(observables, measure_func(psi0,H,args...))

    for i in 2:collapse_idx
        psi0, norm_1 = timeevo_tdvp(H, psi0, -time_intervals[i];
            tau=τ,
            cutoff=cutoff,
            normalize=false,
            maxm=maxBond,
            solver_backend="applyexp", shift=0.0)
        
        logweight[i] = norm_1 + logweight[i-1]
        normalize!(psi0)
        push!(observables, measure_func(psi0,H,args...))
    end

    #### store MPS at target temperature
    psi_collapse = copy(psi0)

    ############ Evolution after target temperature ####################
    for i in (collapse_idx+1):length(time_intervals)
        psi0, norm_1 = timeevo_tdvp(H, psi0, -time_intervals[i];
            tau=τ,
            cutoff=cutoff,
            normalize=false,
            maxm=maxBond,
            solver_backend="applyexp", shift=0.0)

        logweight[i] = norm_1 + logweight[i-1]
        normalize!(psi0)
        push!(observables, measure_func(psi0,H,args...))
    end

    return logweight, observables, psi_collapse
end
