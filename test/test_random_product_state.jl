using ITensors
using ITensorMPS

@testset "random_product_state" begin

    N = 10

    @testset "S=1/2" begin
        sites = siteinds("S=1/2", N)
        up = local_state_index(sites[1], "Up")
        dn = local_state_index(sites[1], "Dn")

        @testset "no nup/ndn" begin
            state = random_product_state(sites)
            @test length(state) == N
            @test all(s in (up, dn) for s in state)
        end

        @testset "nup specified" begin
            for nup in [0, 3, N÷2, N]
                state = random_product_state(sites; nup=nup)
                @test length(state) == N
                @test all(s in (up, dn) for s in state)
                @test count(==(up), state) == nup
                @test count(==(dn), state) == N - nup
            end
        end

        @testset "ndn warns and is ignored" begin
            @test_logs (:warn,) random_product_state(sites; nup=3, ndn=7)
            state = random_product_state(sites; nup=4, ndn=99)
            @test count(==(up), state) == 4
        end

        @testset "errors" begin
            @test_throws ErrorException random_product_state(sites; nup=N+1)
            @test_throws ErrorException random_product_state(sites; nup=-1)
            @test_throws ErrorException random_product_state(sites; ndn=-1)
        end
    end

    @testset "tJ" begin
        sites = siteinds("tJ", N)
        emp = local_state_index(sites[1], "Emp")
        up  = local_state_index(sites[1], "Up")
        dn  = local_state_index(sites[1], "Dn")
        valid = (emp, up, dn)

        @testset "no nup/ndn" begin
            state = random_product_state(sites)
            @test length(state) == N
            @test all(s in valid for s in state)
        end

        @testset "nup and ndn specified" begin
            for (nup_val, ndn_val) in [(3, 4), (0, 5), (5, 0), (0, 0), (N, 0), (0, N)]
                state = random_product_state(sites; nup=nup_val, ndn=ndn_val)
                @test length(state) == N
                @test all(s in valid for s in state)
                @test count(==(up), state) == nup_val
                @test count(==(dn), state) == ndn_val
                @test count(==(emp), state) == N - nup_val - ndn_val
            end
        end

        @testset "errors" begin
            @test_throws ErrorException random_product_state(sites; nup=6, ndn=6)
            @test_throws ErrorException random_product_state(sites; nup=3)
            @test_throws ErrorException random_product_state(sites; ndn=3)
            @test_throws ErrorException random_product_state(sites; nup=-1, ndn=3)
            @test_throws ErrorException random_product_state(sites; nup=3, ndn=-1)
        end
    end

    @testset "Electron" begin
        sites = siteinds("Electron", N)
        emp  = local_state_index(sites[1], "Emp")
        up   = local_state_index(sites[1], "Up")
        dn   = local_state_index(sites[1], "Dn")
        updn = local_state_index(sites[1], "UpDn")
        valid = (emp, up, dn, updn)

        @testset "no nup/ndn" begin
            state = random_product_state(sites)
            @test length(state) == N
            @test all(s in valid for s in state)
        end

        @testset "nup and ndn specified" begin
            for (nup_val, ndn_val) in [(3, 4), (0, 5), (5, 0), (0, 0), (N, N)]
                state = random_product_state(sites; nup=nup_val, ndn=ndn_val)
                @test length(state) == N
                @test all(s in valid for s in state)
                @test count(s -> s == up || s == updn, state) == nup_val
                @test count(s -> s == dn || s == updn, state) == ndn_val
            end
        end

        @testset "errors" begin
            @test_throws ErrorException random_product_state(sites; nup=N+1, ndn=0)
            @test_throws ErrorException random_product_state(sites; nup=0, ndn=N+1)
            @test_throws ErrorException random_product_state(sites; nup=3)
            @test_throws ErrorException random_product_state(sites; ndn=3)
            @test_throws ErrorException random_product_state(sites; nup=-1, ndn=3)
            @test_throws ErrorException random_product_state(sites; nup=3, ndn=-1)
        end
    end

    @testset "unsupported site type" begin
        sites = siteinds("Boson", N)
        @test_throws ErrorException random_product_state(sites)
    end

    @testset "empty sites error" begin
        @test_throws ErrorException random_product_state(Index[])
    end

end
