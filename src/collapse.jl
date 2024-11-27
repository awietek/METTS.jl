function ITensors.op(::OpName"Ry", ::SiteType"Spinhalf")
    s2 = 1.0 / sqrt(2)
    return [
        s2  s2 
        -s2 s2 
    ]
end

function ITensors.op(::OpName"Ry", ::SiteType"tJ")
    s2 = 1.0 / sqrt(2)
    return [
        1.0 0.0 0.0
        0.0 s2  s2 
        0.0 -s2 s2 
    ]
end

function ITensors.op(::OpName"Ry", ::SiteType"Electron")
    s2 = 1.0 / sqrt(2)
    return [
        1.0 0.0 0.0 0.0
        0.0 s2  s2  0.0
        0.0 -s2 s2  0.0
        0.0 0.0 0.0 1.0
    ]
end

function collapse!(psi::MPS, basis::AbstractString=nothing)
    if basis === nothing || basis == "Z"
        return sample!(psi)
    elseif basis == "X"
        s = siteinds(psi)
        N = length(psi)
        Ry_gates = [op("Ry", s[n]) for n in 1:N]
        # Ry_gates = ops([("Ry", n, (θ=π / 2,)) for n in 1:N], s)
        psi = apply(Ry_gates, psi)
        return sample!(psi)
    else
        error("Invalid basis specified in METTS collapse.")
    end
end



function collapse_with_qn!(psi::MPS, basis::AbstractString=nothing)
    if basis === nothing || basis == "Z"
        return sample!(psi)
    elseif basis == "X"
        # psi_noqn = removeqns(psi)
        psi_noqn = dense(psi)
        s = siteinds(psi_noqn)
        N = length(psi_noqn)
        Ry_gates = [op("Ry", s[n]) for n in 1:N]
        psi_noqn = apply(Ry_gates, psi_noqn)
        return sample!(psi_noqn)
    else
        error("Invalid basis specified in METTS collapse.")
    end
end
