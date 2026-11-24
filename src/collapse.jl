"""
    x_rotation(s::Index) -> ITensor

Single-site unitary `U` rotating the x-basis spin states onto the z basis,
`U|+x⟩ = |Up⟩` and `U|-x⟩ = |Dn⟩`, with `|±x⟩ = (|Up⟩ ± |Dn⟩)/√2`. This is the
Hadamard gate acting on the spin states `Up`, `Dn`; `Emp` and `UpDn` are left
unchanged, so the particle number is conserved while Sz is not.

Supported site types: `"S=1/2"`, `"tJ"`, `"Electron"`. `s` must not carry
quantum numbers, since `U` changes Sz.
"""
function x_rotation(s::Index)
    hasqns(s) && error("x_rotation changes Sz and cannot act on a QN conserving index; use dense(s)")
    local_state_strings(s)                   # errors for unsupported site types
    up = local_state_index(s, "Up")
    dn = local_state_index(s, "Dn")
    c = 1 / sqrt(2)
    M = Matrix{Float64}(I, dim(s), dim(s))
    M[up, up] = c
    M[up, dn] = c
    M[dn, up] = c
    M[dn, dn] = -c
    return ITensor(M, s', dag(s))
end

# Sample phi in the given basis; phi is a copy that may be modified.
function _collapse!(rng::AbstractRNG, phi::MPS, basis::AbstractString)
    if basis == "X"
        phi = apply([x_rotation(s) for s in siteinds(phi)], phi)
    elseif basis != "Z"
        error("Invalid basis \"$basis\" specified in METTS collapse, must be \"Z\" or \"X\".")
    end
    # sample requires the orthogonality center at site 1 and norm 1
    orthogonalize!(phi, 1)
    normalize!(phi)
    return sample(rng, phi)
end

"""
    collapse(psi::MPS, basis="Z"; rng=Random.default_rng()) -> Vector{Int}

Sample a product state from `|⟨σ|psi⟩|²` in the `"Z"` or `"X"` basis and
return it as 1-based ITensors state indices. `psi` is not modified and need not
be normalized or orthogonalized.

In the `"X"` basis the returned index of `Up` (`Dn`) means the spin state
`|+x⟩` (`|-x⟩`) on that site, while `Emp` and `UpDn` keep their meaning, see
`x_rotation`. The `"X"` basis requires `psi` without quantum numbers; use
`collapse_with_qn` for QN conserving states.
"""
function collapse(psi::MPS, basis::AbstractString="Z"; rng::AbstractRNG=Random.default_rng())
    basis == "X" && hasqns(psi) &&
        error("collapse in the X basis does not conserve Sz; use collapse_with_qn for QN conserving MPS")
    return _collapse!(rng, copy(psi), basis)
end

"""
    collapse_with_qn(psi::MPS, basis="Z"; rng=Random.default_rng()) -> Vector{Int}

Like `collapse`, for a `psi` that conserves quantum numbers. In the `"X"`
basis the quantum numbers are removed before rotating, and the sample is
returned with the same labels as `collapse`: index of `Up` (`Dn`) means `|+x⟩`
(`|-x⟩`). `psi` is not modified.

Continuing the METTS chain with the returned indices read as a z-basis product
state (as needed to keep conserving quantum numbers) is equivalent to the true
x-basis state only if the Hamiltonian is invariant under global spin rotations.
Sz then fluctuates between steps while the particle number stays conserved.
"""
function collapse_with_qn(psi::MPS, basis::AbstractString="Z"; rng::AbstractRNG=Random.default_rng())
    phi = basis == "X" ? dense(psi) : copy(psi)
    return _collapse!(rng, phi, basis)
end
