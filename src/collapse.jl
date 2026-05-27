function collapse(psi::MPS, basis::AbstractString=nothing)
    if basis === nothing || basis == "Z"
        return sample(psi)
    elseif basis == "X"
        s = siteinds(psi)
        N = length(psi)
        H_gates = [op("H", s[n]) for n in 1:N]
        psi2 = apply(H_gates, psi)
        orthogonalize!(psi2, 1)
        return sample(psi2)
    else
        error("Invalid basis specified in METTS collapse.")
    end
end

function collapse_with_qn(psi::MPS, basis::AbstractString=nothing)
    if basis === nothing || basis == "Z"
        return sample(psi)
    elseif basis == "X"
        # psi_noqn = removeqns(psi)
        psi_noqn = dense(psi)
        s = siteinds(psi_noqn)
        N = length(psi_noqn)
        H_gates = [op("H", s[n]) for n in 1:N]
        psi_noqn = apply(H_gates, psi_noqn)
        return sample!(psi_noqn)
    else
        error("Invalid basis specified in METTS collapse.")
    end
end
