using ITensors, ITensorMPS
using LinearAlgebra

@testset "timeevo_tdvp_extend log_norm" begin
    N = 6
    s = siteinds("S=1/2", N)
    os = OpSum()
    for j in 1:(N - 1)
        os += 0.5, "S+", j, "S-", j + 1
        os += 0.5, "S-", j, "S+", j + 1
        os += "Sz", j, "Sz", j + 1
    end
    H = MPO(os, s)
    psi0 = MPS(s, [isodd(j) ? "Up" : "Dn" for j in 1:N])

    # exact log ⟨σ|exp(-βH)|σ⟩ from dense matrices
    C = combiner(s...)
    Cp = combiner(prime.(s)...)
    Hm = array(prod(H) * C * Cp, combinedind(Cp), combinedind(C))
    v = array(prod(psi0) * C, combinedind(C))
    exact(β) = log(v' * exp(-β * Symmetric(Hm)) * v)

    @testset "β = $β" for β in (0.5, 2.0)
        psi, log_norm = timeevo_tdvp_extend(H, psi0, -β / 2;
                                            tau=0.05, cutoff=1e-12, maxm=64, silent=true)
        @test isapprox(log_norm, exact(β); atol=1e-3)
        @test norm(psi) ≈ 1
    end

    @testset "zero time" begin
        psi, log_norm = timeevo_tdvp_extend(H, psi0, 0.0; silent=true)
        @test log_norm == 0.0
        @test inner(psi, psi0) ≈ 1
        psi, log_norm = timeevo_tdvp(H, psi0, 0.0; silent=true)
        @test log_norm == 0.0
    end

    @testset "time shorter than tau0 (no bulk evolution)" begin
        β = 0.04
        psi, log_norm = timeevo_tdvp_extend(H, psi0, -β / 2;
                                            tau=0.05, tau0=0.05, cutoff=1e-12, maxm=64, silent=true)
        @test isfinite(log_norm)
        @test isapprox(log_norm, exact(β); atol=1e-4)
    end
end
