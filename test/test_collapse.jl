using ITensors, ITensorMPS
using LinearAlgebra
using Random

@testset "collapse" begin

    @testset "x_rotation is unitary and conserves the particle number: $st" for st in ("S=1/2", "tJ", "Electron")
        s = siteind(st)
        U = array(x_rotation(s), s', s)
        @test U' * U ≈ I
        for name in local_state_strings(s)
            name in ("Up", "Dn") && continue
            n = local_state_index(s, name)
            @test U[:, n] == [k == n ? 1.0 : 0.0 for k in 1:dim(s)]
        end
        @test_throws ErrorException x_rotation(siteind(st; conserve_qns=true))
    end

    @testset "X basis labels: Up means |+x⟩, Dn means |-x⟩ (S=1/2)" begin
        N = 6
        s = siteinds("S=1/2", N)
        up, dn = local_state_index(s[1], "Up"), local_state_index(s[1], "Dn")
        @test collapse(MPS(s, fill("X+", N)), "X") == fill(up, N)
        @test collapse(MPS(s, fill("X-", N)), "X") == fill(dn, N)
        @test collapse(MPS(s, fill("X+", N)), "X") == collapse_with_qn(MPS(s, fill("X+", N)), "X")
    end

    @testset "X basis labels: $st" for st in ("tJ", "Electron")
        # x_rotation is its own inverse, so it maps |Up⟩ to |+x⟩ and |Dn⟩ to |-x⟩
        names = local_state_strings(siteind(st))
        N = 2 * length(names)
        s = siteinds(st, N)
        config = [names[mod1(n, length(names))] for n in 1:N]
        phi = apply([x_rotation(si) for si in s], MPS(s, config))
        @test collapse(phi, "X") == [local_state_index(s[n], config[n]) for n in 1:N]
    end

    @testset "Z basis, QN conserving, product state: $st" for st in ("S=1/2", "tJ", "Electron")
        names = local_state_strings(siteind(st))
        N = 2 * length(names)
        s = siteinds(st, N; conserve_qns=true)
        config = [names[mod1(n, length(names))] for n in 1:N]
        idx = [local_state_index(s[n], config[n]) for n in 1:N]
        @test collapse_with_qn(MPS(s, config), "Z") == idx
        @test collapse(MPS(s, config)) == idx
    end

    @testset "psi is not modified, orthogonality center and norm are arbitrary" begin
        N = 8
        s = siteinds("Electron", N; conserve_qns=true)
        psi = random_mps(s, [isodd(n) ? "Up" : "Dn" for n in 1:N]; linkdims=8)
        psi = orthogonalize(psi, N)             # center at the last site
        psi[N] *= 3.0                           # norm 3
        psi_copy = copy(psi)
        for basis in ("Z", "X")
            smp = collapse_with_qn(psi, basis; rng=Xoshiro(1))
            @test length(smp) == N
            # the particle number is conserved in both bases
            labels = [local_state_string(s[n], smp[n]) for n in 1:N]
            @test sum(l == "UpDn" ? 2 : l == "Emp" ? 0 : 1 for l in labels) == N
        end
        @test all(norm(psi[n] - psi_copy[n]) == 0 for n in 1:N)
        @test orthocenter(psi) == N
        @test_throws ErrorException collapse(psi, "X")
        @test_throws ErrorException collapse_with_qn(psi, "Y")
    end

    @testset "rng makes the collapse reproducible" begin
        N = 8
        s = siteinds("S=1/2", N; conserve_qns=true)
        psi = random_mps(s, [isodd(n) ? "Up" : "Dn" for n in 1:N]; linkdims=8)
        for basis in ("Z", "X")
            @test collapse_with_qn(psi, basis; rng=Xoshiro(42)) ==
                  collapse_with_qn(psi, basis; rng=Xoshiro(42))
        end
    end

    @testset "X basis sampling probabilities (S=1/2)" begin
        # single site |ψ⟩ = cos θ |Up⟩ + sin θ |Dn⟩: P(+x) = (cos θ + sin θ)^2 / 2
        s = siteinds("S=1/2", 2)
        θ = 0.3
        up1, dn1 = local_state_index(s[1], "Up"), local_state_index(s[1], "Dn")
        v = zeros(2)
        v[up1], v[dn1] = cos(θ), sin(θ)
        psi = MPS(ITensor(v, s[1]) * state(s[2], "Up"), s)
        up = local_state_index(s[1], "Up")
        rng = Xoshiro(7)
        nsamples = 20000
        frac = count(_ -> collapse(psi, "X"; rng)[1] == up, 1:nsamples) / nsamples
        @test isapprox(frac, (cos(θ) + sin(θ))^2 / 2; atol=0.02)
    end
end

@testset "entropy_von_neumann" begin
    N = 8
    s = siteinds("S=1/2", N; conserve_qns=true)

    # product state: zero entropy at every bond, including b = 1
    prod = MPS(s, [isodd(n) ? "Up" : "Dn" for n in 1:N])
    for b in 1:(N - 1)
        @test entropy_von_neumann(prod, b) == 0.0
    end

    # two-site singlet: S = log 2, independent of the norm
    sd = siteinds("S=1/2", 2)
    singlet = MPS(ITensor([0.0 1.0; -1.0 0.0] / sqrt(2), sd[1], sd[2]), sd)
    @test entropy_von_neumann(singlet, 1) ≈ log(2)
    @test entropy_von_neumann(5.0 * singlet, 1) ≈ log(2)

    # does not modify psi
    psi = random_mps(s, [isodd(n) ? "Up" : "Dn" for n in 1:N]; linkdims=8)
    psi = orthogonalize(psi, 1)
    S = entropy_von_neumann(psi, N ÷ 2)
    @test orthocenter(psi) == 1
    @test isfinite(S) && S >= 0

    @test_throws ErrorException entropy_von_neumann(psi, 0)
    @test_throws ErrorException entropy_von_neumann(psi, N)
end
