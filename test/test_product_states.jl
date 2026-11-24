using ITensors, ITensorMPS
using HDF5
using Random

@testset "product_states" begin

    @testset "write/read round trip: $st" for st in ("S=1/2", "tJ", "Electron")
        N = 7
        sites = siteinds(st, N; conserve_qns=true)
        states = random_product_states(sites, 5; rng=Xoshiro(1))
        mktempdir() do dir
            filename = joinpath(dir, "sub", "states.h5")      # directories are created
            write_product_states(filename, sites, states)
            @test read_product_states(filename, sites) == states
            h5open(filename, "r") do f
                @test read(f, "local_states") == local_state_strings(sites[1])
                M = read(f, "states")
                @test eltype(M) == UInt8
                @test size(M) == (N, 5)
                @test M[:, 2] == UInt8.(states[2] .- 1)
            end
            @test_throws ErrorException write_product_states(filename, sites, states)
        end
    end

    @testset "errors" begin
        mktempdir() do dir
            sites = siteinds("Electron", 4)
            filename = joinpath(dir, "states.h5")
            @test_throws ErrorException write_product_states(filename, sites, Vector{Int}[])
            @test_throws ErrorException write_product_states(filename, sites, [[1, 2, 3]])
            @test_throws ErrorException write_product_states(filename, sites, [[1, 2, 3, 5]])
            @test_throws ErrorException write_product_states(filename, [siteind("Electron"), siteind("tJ")], [[1, 2]])
            write_product_states(filename, sites, [[1, 2, 3, 4]])
            @test_throws ErrorException read_product_states(filename, siteinds("tJ", 4))
            @test_throws ErrorException read_product_states(filename, siteinds("Electron", 5))
        end
    end

    @testset "random_product_state with an rng" begin
        sites = siteinds("Electron", 10)
        # the seed method is the rng method with a fresh MersenneTwister
        @test random_product_state(sites, 5; nup=4, ndn=3) ==
              random_product_state(MersenneTwister(5), sites; nup=4, ndn=3)
        states = random_product_states(sites, 50; rng=Xoshiro(3), nup=4, ndn=3)
        @test states == random_product_states(sites, 50; rng=Xoshiro(3), nup=4, ndn=3)
        @test length(unique(states)) > 1
        for s in states
            labels = [local_state_string(sites[n], s[n]) for n in 1:10]
            @test count(in(("Up", "UpDn")), labels) == 4
            @test count(in(("Dn", "UpDn")), labels) == 3
        end
    end

    @testset "sample_product_states" begin
        N = 8
        sites = siteinds("Electron", N; conserve_qns=true)
        config = [isodd(n) ? "Up" : "Dn" for n in 1:N]
        prod = MPS(sites, config)
        @test sample_product_states(prod, 3) == fill([local_state_index(sites[n], config[n]) for n in 1:N], 3)

        psi = random_mps(sites, config; linkdims=8)
        psi = orthogonalize(psi, N)
        psi[N] *= 2.0
        psi_copy = copy(psi)
        smps = sample_product_states(psi, 20; rng=Xoshiro(2))
        @test smps == sample_product_states(psi, 20; rng=Xoshiro(2))
        for s in smps
            labels = [local_state_string(sites[n], s[n]) for n in 1:N]
            @test count(in(("Up", "UpDn")), labels) == N ÷ 2        # QN sector kept
            @test count(in(("Dn", "UpDn")), labels) == N ÷ 2
        end
        @test all(norm(psi[n] - psi_copy[n]) == 0 for n in 1:N)
    end
end
