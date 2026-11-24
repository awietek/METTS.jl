# ---------------------------------------------------------------------------
# Sets of product states in HDF5, e.g. initial states for METTS chains:
#
#   /local_states   String (nlocal,)            local basis state names
#   /states         UInt8  (nsites, nstates)    0-based into local_states
#
# the METTSLibrary convention for states, with local_states as a dataset.
# In memory a product state is a Vector{Int} of 1-based ITensors state
# indices, as returned by random_product_state and collapse.
# ---------------------------------------------------------------------------

# local state names of sites, which must all have the same site type
function _local_states(sites)
    isempty(sites) && error("empty sites")
    ls = local_state_strings(sites[1])
    all(local_state_strings(s) == ls for s in sites) || error("all sites must have the same site type")
    return ls
end

# local_state_strings is ordered by the ITensors index, so encoding is a shift
function _encode(state::AbstractVector{<:Integer}, nlocal::Integer)
    all(x -> 1 <= x <= nlocal, state) ||
        error("state indices must lie in 1:$nlocal, got the range $(extrema(state))")
    return UInt8.(state .- 1)
end
_decode(x::AbstractVector{<:Integer}) = Int.(x) .+ 1

function _check_local_states(stored, sites, filename)
    expected = _local_states(sites)
    stored == expected ||
        error("'$filename' stores the local states $stored, the sites have $expected")
    return nothing
end

"""
    write_product_states(filename, sites, states)

Write the product states `states` (a vector of 1-based ITensors state index
vectors, as returned by `random_product_state` or `sample`) to the HDF5 file
`filename`, as `UInt8` dataset `states` of size `(nsites, nstates)` holding
0-based indices into the string dataset `local_states`. Errors if `filename`
already exists.
"""
function write_product_states(filename::AbstractString, sites,
                              states::AbstractVector{<:AbstractVector{<:Integer}})
    isfile(filename) && error("'$filename' already exists")
    isempty(states) && error("no product states given")
    N = length(sites)
    ls = _local_states(sites)
    M = Matrix{UInt8}(undef, N, length(states))
    for (j, s) in enumerate(states)
        length(s) == N || error("product state $j has $(length(s)) sites, expected $N")
        M[:, j] = _encode(s, length(ls))
    end
    mkpath(dirname(abspath(filename)))
    h5open(filename, "w") do f
        f["local_states"] = ls
        f["states"] = M
    end
    return filename
end

"""
    read_product_states(filename, sites) -> Vector{Vector{Int}}

Read the product states written by `write_product_states` as 1-based ITensors
state indices for `sites`. Errors if the stored local states do not match the
site type of `sites`.
"""
function read_product_states(filename::AbstractString, sites)
    h5open(filename, "r") do f
        _check_local_states(read(f, "local_states"), sites, filename)
        M = read(f, "states")
        size(M, 1) == length(sites) ||
            error("'$filename' stores states on $(size(M, 1)) sites, expected $(length(sites))")
        return [_decode(M[:, j]) for j in 1:size(M, 2)]
    end
end

"""
    random_product_states(sites, nstates; rng=Random.default_rng(), nup=nothing, ndn=nothing)
        -> Vector{Vector{Int}}

`nstates` independent random product states, see `random_product_state`.
"""
function random_product_states(sites, nstates::Integer; rng::AbstractRNG=Random.default_rng(),
                               nup=nothing, ndn=nothing)
    return [random_product_state(rng, sites; nup, ndn) for _ in 1:nstates]
end

"""
    sample_product_states(psi::MPS, nstates; rng=Random.default_rng()) -> Vector{Vector{Int}}

`nstates` independent samples of product states `σ` in the z basis with
probability `|⟨σ|psi⟩|²`, as 1-based ITensors state indices. `psi` is not
modified and need not be normalized or orthogonalized. A QN conserving `psi`
gives samples in its quantum number sector.
"""
function sample_product_states(psi::MPS, nstates::Integer; rng::AbstractRNG=Random.default_rng())
    phi = orthogonalize(psi, 1)
    normalize!(phi)
    return [sample(rng, phi) for _ in 1:nstates]
end
