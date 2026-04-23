"""
    random_product_state(sites; nup=nothing, ndn=nothing)

Return a `Vector{Int}` of state indices for a random product state for the given sites.
The integers represent basis state indices as defined by ITensors (1-based).

Supported site types: `"S=1/2"`, `"tJ"`, `"Electron"`.

- `seed`: random seed of random number generator (optional)
- `nup`: number of up-spin electrons (optional)
- `ndn`: number of down-spin electrons (optional; ignored for `"S=1/2"`)

For `"S=1/2"`, `ndn` is inferred as `N - nup`; any explicit `ndn` is ignored with a warning.
For `"tJ"` and `"Electron"`, both `nup` and `ndn` must be specified together (or both omitted).
For `"Electron"`, double occupancy is allowed when up/down electrons share the same site.
If neither `nup` nor `ndn` is specified, each site is assigned a uniformly random state.
"""
function random_product_state(sites, seed=42; nup=nothing, ndn=nothing)
    rng = MersenneTwister(seed)
    
    N = length(sites)
    if N == 0
        error("empty sites")
    end
    if nup !== nothing && nup < 0
        error("nup must be nonnegative, got nup=$(nup)")
    end
    if ndn !== nothing && ndn < 0
        error("ndn must be nonnegative, got ndn=$(ndn)")
    end

    if hastags(sites[1], "S=1/2") || hastags(sites[1], "Spinhalf")
        if ndn !== nothing
            @warn "ndn is ignored for S=1/2 sites; it is inferred as N - nup"
        end
        up = local_state_index(sites[1], "Up")
        dn = local_state_index(sites[1], "Dn")
        if nup === nothing
            return [rand(rng, Bool) ? up : dn for _ in 1:N]
        else
            if nup > N
                error("nup=$(nup) exceeds N=$(N) for S=1/2")
            end
            states = fill(dn, N)
            states[randperm(rng, N)[1:nup]] .= up
            return states
        end

    elseif hastags(sites[1], "tJ")
        if (nup === nothing) != (ndn === nothing)
            error("For tJ sites, both nup and ndn must be specified together, or both omitted")
        end
        emp = local_state_index(sites[1], "Emp")
        up  = local_state_index(sites[1], "Up")
        dn  = local_state_index(sites[1], "Dn")
        if nup === nothing
            return [rand(rng, [emp, up, dn]) for _ in 1:N]
        else
            if nup + ndn > N
                error("nup + ndn must be ≤ N=$(N) for tJ, got nup=$(nup), ndn=$(ndn)")
            end
            states = fill(emp, N)
            perm = randperm(rng, N)
            states[perm[1:nup]] .= up
            states[perm[nup+1:nup+ndn]] .= dn
            return states
        end

    elseif hastags(sites[1], "Electron")
        if (nup === nothing) != (ndn === nothing)
            error("For Electron sites, both nup and ndn must be specified together, or both omitted")
        end
        emp  = local_state_index(sites[1], "Emp")
        up   = local_state_index(sites[1], "Up")
        dn   = local_state_index(sites[1], "Dn")
        updn = local_state_index(sites[1], "UpDn")
        if nup === nothing
            return [rand(rng, [emp, up, dn, updn]) for _ in 1:N]
        else
            if nup > N || ndn > N
                error("nup and ndn must be ≤ N=$(N) for Electron")
            end
            up_sites = Set(randperm(rng, N)[1:nup])
            dn_sites = Set(randperm(rng, N)[1:ndn])
            states = Vector{Int}(undef, N)
            for i in 1:N
                if i in up_sites && i in dn_sites
                    states[i] = updn
                elseif i in up_sites
                    states[i] = up
                elseif i in dn_sites
                    states[i] = dn
                else
                    states[i] = emp
                end
            end
            return states
        end

    else
        error("Unsupported site type: $(tags(sites[1]))")
    end
end
