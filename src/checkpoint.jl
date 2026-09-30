# ---------------------------------------------------------------------------
# Checkpointed METTS chains. One HDF5 file holds one Markov chain:
#
#   /local_states    String (nlocal,)
#   /rng_seed        UInt64 (nsteps + 1,)         seed of the Xoshiro random
#                    number generator step k uses
#   /product_state   UInt8  (nsites, nsteps + 1)  0-based into local_states;
#                    column k is the product state evolved in step k
#   /<observable>    (nsteps,) for numbers, (size..., nsteps) for arrays
#
# Entry k of rng_seed and product_state is the start of step k, and
# observable[..., k] was measured on the METTS of product_state[:, k]. The
# last entries of rng_seed and product_state are the start of the next step,
# so a chain can be continued as soon as metts_start has written it, even if
# no step completed.
#
# Writing is ordered so that product_state always comes last: metts_start
# writes the first seed, then the initial state; step k writes its
# observables at index k, then the next seed and the next state at index k+1.
# With nstates columns of product_state, nstates - 1 steps are completed.
# Every entry is written at its step index rather than appended, so whatever
# a crash left of an unfinished step is simply overwritten when the step is
# redone. The file belongs to the checkpoint: every dataset other than
# local_states counts as per-step data.
#
# At the end of every step the random number generator draws the seed of a
# fresh Xoshiro for the next step. Storing that seed is all it takes to
# continue a chain exactly as it would have run uninterrupted.
# ---------------------------------------------------------------------------

const _PRODUCT_STATE = "product_state"
const _RNG_SEED = "rng_seed"
const _LOCAL_STATES = "local_states"
const _RESERVED = (_PRODUCT_STATE, _RNG_SEED, _LOCAL_STATES)

_entry(x) = x isa AbstractArray ? x : fill(x)       # numbers as 0-dim arrays
_colons(n) = ntuple(_ -> Colon(), n)
_nentries(ds) = size(ds)[end]

# write x as entry k along the last dimension of dataset `name`, which then
# holds exactly k entries; entries beyond k, left by a crash, are dropped
function _write_entry!(f, name::AbstractString, k::Integer, x)
    a = _entry(x)
    sz = size(a)
    if !haskey(f, name)
        chunk = (sz..., clamp(2^16 ÷ max(1, prod(sz)), 1, 1024))
        create_dataset(f, name, datatype(eltype(a)),
                       dataspace((sz..., 0); max_dims=(sz..., -1)); chunk)
    end
    ds = f[name]
    n = _nentries(ds)
    n >= k - 1 || error("cannot write entry $k of '$name', which has only $n entries")
    HDF5.set_extent_dims(ds, (sz..., k))
    ds[_colons(length(sz))..., k:k] = reshape(a, sz..., 1)
    return nothing
end

"""
    valid_checkpoint(filename) -> Bool

Whether `filename` holds a METTS chain that `metts_resume` can continue,
i.e. `metts_start` has completed on it. `false` if the file does not exist or
`metts_start` was interrupted; the chain is then started with `metts_start`.
Errors if `filename` exists but is not an HDF5 file, so that nothing
unreadable is ever overwritten.

# Example
```julia
if !valid_checkpoint(outfile)
    rng = Xoshiro(seed)
    initial_states = read_product_states(initfile, sites)
    metts_start(outfile, sites, initial_states[rand(rng, 1:length(initial_states))], rng)
end
state, rng, nsteps_done = metts_resume(outfile, sites)
for step in (nsteps_done + 1):nsteps
    psi, _ = timeevo_tdvp_extend(H, MPS(sites, state), -beta / 2)
    energy = real(inner(psi', H, psi))
    state = collapse_with_qn(psi, "X"; rng)
    metts_dump_step!(outfile, (; energy), state, rng)
end
```
"""
function valid_checkpoint(filename::AbstractString)
    isfile(filename) || return false
    # metts_start writes the initial state last, so it has completed iff
    # product_state exists; any inconsistency after that is left to
    # metts_resume to report, so that a damaged chain is never overwritten
    return h5open(filename, "r") do f
        haskey(f, _PRODUCT_STATE) && _nentries(f[_PRODUCT_STATE]) >= 1
    end
end

"""
    metts_start(filename, sites, initial_state, rng)

Start a new METTS chain in `filename` whose first step evolves the product
state `initial_state` (1-based ITensors state indices) and uses a `Xoshiro`
random number generator seeded from `rng`. Continue it with `metts_resume`,
which returns that state and generator. A file left by an interrupted
`metts_start` is replaced; a `valid_checkpoint` is never overwritten.
"""
function metts_start(filename::AbstractString, sites, initial_state::AbstractVector{<:Integer},
                     rng::AbstractRNG)
    valid_checkpoint(filename) &&
        error("'$filename' already holds a METTS chain; continue it with metts_resume or remove it")
    length(initial_state) == length(sites) ||
        error("initial state has $(length(initial_state)) sites, expected $(length(sites))")
    ls = _local_states(sites)
    encoded = _encode(initial_state, length(ls))
    seed = rand(rng, UInt64)
    mkpath(dirname(abspath(filename)))
    h5open(filename, "w") do f
        f[_LOCAL_STATES] = ls
        _write_entry!(f, _RNG_SEED, 1, seed)
        _write_entry!(f, _PRODUCT_STATE, 1, encoded)  # last: marks the start as completed
    end
    return nothing
end

"""
    metts_resume(filename, sites) -> (state, rng, nsteps_done)

Continue the METTS chain in `filename`, which must be a `valid_checkpoint`.
Returns the product state the next step starts from (1-based ITensors state
indices), its random number generator, and the number of completed steps (0
right after `metts_start`). With the returned `rng` the chain continues
exactly as it would have run uninterrupted. The file is not modified;
entries left by an unfinished step are overwritten when it is redone.
"""
function metts_resume(filename::AbstractString, sites)
    valid_checkpoint(filename) ||
        error("'$filename' holds no METTS chain; start it with metts_start")
    return h5open(filename, "r") do f
        _check_local_states(read(f, _LOCAL_STATES), sites, filename)
        ps = f[_PRODUCT_STATE]
        size(ps, 1) == length(sites) ||
            error("'$filename' stores states on $(size(ps, 1)) sites, expected $(length(sites))")
        nstates = _nentries(ps)
        # the write order guarantees these; anything else is a damaged file
        haskey(f, _RNG_SEED) || error("'$filename' has product states but no '$_RNG_SEED'")
        for name in keys(f)
            name in (_PRODUCT_STATE, _LOCAL_STATES) && continue
            expected = name == _RNG_SEED ? nstates : nstates - 1
            _nentries(f[name]) >= expected ||
                error("'$filename': '$name' has $(_nentries(f[name])) entries, expected $expected")
        end
        state = _decode(ps[:, nstates])
        seed = f[_RNG_SEED][nstates:nstates][1]
        return state, Xoshiro(seed), nstates - 1
    end
end

"""
    metts_dump_step!(filename, observables, next_state, rng::Xoshiro)

Write one completed METTS step to the chain in `filename` (started with
`metts_start`): the `observables` measured on the METTS of the current
product state and the collapsed `next_state` (1-based ITensors state indices)
the next step starts from.

`rng` is the generator the step used, as returned by `metts_resume`. It draws
the seed for the next step and is then re-seeded with it in place, so the
caller simply keeps using `rng`; the seed is stored so that `metts_resume`
continues with the same generator. If the step cannot be written, `rng` is
left unchanged.

`observables` maps names to numbers or arrays, e.g. a `NamedTuple` or `Dict`;
every step must record the same names with the same sizes; the element type
is fixed by the first step (later values are converted). Everything is
checked before anything is written, and `next_state` is written last, marking
the step as completed.
"""
function metts_dump_step!(filename::AbstractString, observables,
                          next_state::AbstractVector{<:Integer}, rng::Xoshiro)
    seed = rand(copy(rng), UInt64)                   # rng itself changes only on success
    h5open(filename, "r+") do f
        haskey(f, _PRODUCT_STATE) ||
            error("'$filename' is not a METTS chain; start it with metts_start")
        ps = f[_PRODUCT_STATE]
        length(next_state) == size(ps, 1) ||
            error("next_state has $(length(next_state)) sites, expected $(size(ps, 1))")
        encoded = _encode(next_state, length(read(f, _LOCAL_STATES)))
        k = _nentries(ps)                            # the step being completed

        obs = [String(name) => _entry(v) for (name, v) in pairs(observables)]
        names = first.(obs)
        any(in(_RESERVED), names) && error("the observable names $(_RESERVED) are reserved")
        existing = sort!([name for name in keys(f) if !(name in _RESERVED)])
        if k > 1
            existing == sort(names) ||
                error("'$filename': this step records $(sort(names)), previous steps $existing")
            for (name, a) in obs
                sz = size(f[name])[1:end-1]
                sz == size(a) || error("'$name' has size $(size(a)), previous steps $sz")
            end
        end

        if k == 1
            # no step completed yet: anything here was left by a crashed
            # attempt at step 1, which may even have recorded other observables
            foreach(name -> delete_object(f, name), existing)
        end
        for (name, a) in obs
            _write_entry!(f, name, k, a)
        end
        _write_entry!(f, _RNG_SEED, k + 1, seed)
        _write_entry!(f, _PRODUCT_STATE, k + 1, encoded)   # last: marks the step as completed
    end
    Random.seed!(rng, seed)                          # the generator of the next step
    return nothing
end

"""
    read_metts(filename) -> (; nsteps, local_states, product_states, observables)

Read the completed steps of the METTS chain in `filename`, without modifying
it. `product_states` is a `Matrix{Int}` of size `(nsites, nsteps)` with
1-based ITensors state indices, column `k` being the product state of step `k`
(the state the next step would start from is not included). `observables`
maps each name to its values, the last dimension running over the steps.
"""
function read_metts(filename::AbstractString)
    h5open(filename, "r") do f
        haskey(f, _PRODUCT_STATE) || error("'$filename' is not a METTS chain")
        ps = read(f, _PRODUCT_STATE)
        n = size(ps, 2) - 1
        observables = Dict{String,Array}()
        for name in keys(f)
            name in _RESERVED && continue
            a = read(f[name])
            observables[name] = collect(selectdim(a, ndims(a), 1:min(n, size(a, ndims(a)))))
        end
        return (; nsteps=n, local_states=read(f, _LOCAL_STATES),
                product_states=Int.(ps[:, 1:n]) .+ 1, observables)
    end
end

# ---------------------------------------------------------------------------
# Earlier interface, for files with a Dumper-written product state dataset.
# Note that for chains written by metts_dump_step!, count_existing_steps
# returns nsteps + 1, as the initial state is stored too.
# ---------------------------------------------------------------------------

"""
    count_existing_steps(filename::String, tag::String) -> Int

Returns the number of completed METTS steps already stored in the HDF5 file
by counting the entries in the `product_state` dataset. Returns 0 if the file
does not exist. Throws an error if it is corrupted or has no
`product_state` dataset.
"""
function count_existing_steps(filename::String, product_state_name::String)::Int
    isfile(filename) || return 0
    try
        result = h5open(filename, "r") do f
            if !haskey(f, product_state_name)
                existing_datasets = keys(f)

                if isempty(existing_datasets)
                    @warn "Checkpoint file exists but is empty" filename
                    return 0
                else
                    found_keys_str = join(existing_datasets, ", ")
                    error("Checkpoint file '$filename' contains data ($found_keys_str) but is missing the mandatory '$product_state_name' dataset.")
                end
            end

            # Get dimensions
            dataset_dims = size(f[product_state_name])

            if length(dataset_dims) ∉ (1, 2)
                error("Data shape mismatch in '$filename': expected a 1D or 2D array for '$product_state_name', but got dimensions $(dataset_dims).")
            elseif length(dataset_dims) == 2
                return dataset_dims[2]
            else
                return 1
            end
            return result
        end
    catch e
        @error "Failed to count existing steps due to an error:" filename exception = (e, catch_backtrace())
        rethrow(e)
    end
end


"""
    read_last_product_state(filename::String)
Reads the last collapsed product state from the HDF5 file. Errors if
the file cannot be read or the dataset is empty.
"""
function read_last_product_state(filename::String, product_state_name::String)
    isfile(filename) || error("File '$filename' not found.")
    try
        product_state = h5open(filename, "r") do f
            if !haskey(f, product_state_name)
                error("No '$product_state_name' dataset found in '$filename'.")
            end

            dataset = f[product_state_name]
            dataset_dims = size(dataset)

            if 0 in dataset_dims || isempty(dataset_dims)
                error("Empty '$product_state_name' dataset in '$filename'.")
            end

            if length(dataset_dims) ∉ (1, 2)
                error("Data shape mismatch in '$filename': expected a 1D or 2D array for '$product_state_name', but got dimensions $(dataset_dims).")
            elseif length(dataset_dims) == 2
                return dataset[:, end]
            else
                return dataset[:]
            end
        end
        return product_state
    catch e
        @error "Failed to read product state due to an error:" filename exception = (e, catch_backtrace())
        rethrow(e)
    end
end
