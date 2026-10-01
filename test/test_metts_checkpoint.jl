using ITensors, ITensorMPS
using HDF5
using Random

@testset "metts_checkpoint" begin
    N = 6
    sites = siteinds("Electron", N; conserve_qns=true)
    initial = random_product_state(Xoshiro(0), sites; nup=3, ndn=3)

    # a stand-in for a METTS step that uses the rng like a collapse does
    function fake_step(state, rng)
        observables = (; energy=sum(state) + rand(rng), n=float.(state), m=rand(rng, 2, 3))
        next_state = random_product_state(rng, sites; nup=3, ndn=3)
        return observables, next_state
    end

    # the pattern of a METTS script: start the chain if needed, then continue it
    function run_chain!(filename, nsteps; seed=11)
        valid_checkpoint(filename) ||
            metts_start(filename, sites, initial, Xoshiro(seed); parameters=(; seed, T=0.5, init="dmrg"))
        state, rng, ndone = metts_resume(filename, sites)
        for _ in (ndone + 1):nsteps
            obs, state = fake_step(state, rng)
            metts_dump_step!(filename, obs, state, rng)
        end
        return ndone
    end

    # extend datasets by one entry, as a step interrupted after writing them
    function fake_crash!(filename, names)
        h5open(filename, "r+") do f
            for name in names
                ds = f[name]
                sz = size(ds)
                HDF5.set_extent_dims(ds, (sz[1:end-1]..., sz[end] + 1))
            end
        end
    end

    @testset "metts_start" begin
        mktempdir() do dir
            filename = joinpath(dir, "sub", "chain.h5")      # directories are created
            @test !valid_checkpoint(filename)
            metts_start(filename, sites, initial, Xoshiro(1))
            @test valid_checkpoint(filename)                  # continuable right away
            seed = rand(Xoshiro(1), UInt64)
            h5open(filename, "r") do f
                @test read(f, "local_states") == local_state_strings(sites[1])
                @test read(f, "product_state") == reshape(UInt8.(initial .- 1), N, 1)
                @test read(f, "rng_seed") == [seed]
            end
            @test metts_resume(filename, sites) == (initial, Xoshiro(seed), 0)
            @test read_metts(filename).nsteps == 0
            @test read_metts(filename).parameters == Dict{String,Any}()
        end
    end

    @testset "parameters" begin
        mktempdir() do dir
            filename = joinpath(dir, "chain.h5")
            run_chain!(filename, 3)
            r = read_metts(filename)
            @test r.parameters == Dict{String,Any}("seed" => 11, "T" => 0.5, "init" => "dmrg")
            @test !haskey(r.observables, "parameters")
            h5open(filename, "r") do f
                @test read(f["parameters/T"]) == 0.5
            end
            @test_throws ErrorException metts_dump_step!(filename, (; parameters=1.0), initial, Xoshiro(1))
            # a chain is never overwritten
            @test_throws ErrorException metts_start(filename, sites, initial, Xoshiro(1))
        end
    end

    @testset "layout: entry k of rng_seed and product_state is the start of step k" begin
        mktempdir() do dir
            filename = joinpath(dir, "chain.h5")
            metts_start(filename, sites, initial, Xoshiro(2))
            state, rng, _ = metts_resume(filename, sites)
            states = [state]
            rngs = [copy(rng)]
            for _ in 1:4
                obs, state = fake_step(state, rng)
                @test metts_dump_step!(filename, obs, state, rng) === nothing
                push!(states, state)
                push!(rngs, copy(rng))
            end
            r = read_metts(filename)
            @test r.nsteps == 4
            @test r.local_states == local_state_strings(sites[1])
            @test r.product_states == reduce(hcat, states[1:4])
            @test size(r.observables["energy"]) == (4,)
            @test size(r.observables["n"]) == (N, 4)
            @test size(r.observables["m"]) == (2, 3, 4)
            @test r.observables["n"] == float.(r.product_states)
            h5open(filename, "r") do f
                @test size(f["product_state"]) == (N, 5)
                seeds = read(f, "rng_seed")
                @test eltype(seeds) == UInt64 && size(seeds) == (5,)
                @test all(Xoshiro(seeds[k]) == rngs[k] for k in 1:5)
            end
            @test metts_resume(filename, sites) == (states[5], rngs[5], 4)
        end
    end

    @testset "a continued chain equals an uninterrupted one" begin
        mktempdir() do dir
            straight = joinpath(dir, "straight.h5")
            pieces = joinpath(dir, "pieces.h5")
            @test run_chain!(straight, 7) == 0
            @test run_chain!(pieces, 0) == 0                 # started, no step
            @test run_chain!(pieces, 3) == 0
            # the seed is only used to start the chain
            @test run_chain!(pieces, 5; seed=999) == 3
            @test run_chain!(pieces, 7; seed=999) == 5
            @test run_chain!(pieces, 7) == 7                 # nothing left to do
            a, b = read_metts(straight), read_metts(pieces)
            @test a.product_states == b.product_states
            @test a.observables == b.observables
        end
    end

    @testset "a crash during a step is overwritten when it is redone: $what" for (what, names) in (
            "after some observables" => ["energy", "n"],
            "after all observables and the next seed" => ["energy", "n", "m", "rng_seed"])
        mktempdir() do dir
            straight = joinpath(dir, "straight.h5")
            crashed = joinpath(dir, "crashed.h5")
            run_chain!(straight, 4)
            run_chain!(crashed, 2)
            fake_crash!(crashed, names)
            bytes = read(crashed)
            @test metts_resume(crashed, sites)[3] == 2
            @test read(crashed) == bytes                     # resuming does not modify the file
            @test run_chain!(crashed, 4) == 2
            a, b = read_metts(straight), read_metts(crashed)
            @test a.product_states == b.product_states
            @test a.observables == b.observables
            h5open(crashed, "r") do f
                @test size(f["energy"]) == (4,)
                @test size(f["rng_seed"]) == (5,)
            end
        end
    end

    @testset "a crash during the first step is continued" begin
        mktempdir() do dir
            straight = joinpath(dir, "straight.h5")
            crashed = joinpath(dir, "crashed.h5")
            run_chain!(straight, 2)
            run_chain!(crashed, 0)
            h5open(crashed, "r+") do f
                METTS._write_entry!(f, "energy", 1, 1.0)     # half-written first step
            end
            @test valid_checkpoint(crashed)
            @test run_chain!(crashed, 2) == 0
            @test read_metts(crashed).observables == read_metts(straight).observables
        end
    end

    @testset "a crashed first attempt may have recorded other observables" begin
        mktempdir() do dir
            straight = joinpath(dir, "straight.h5")
            crashed = joinpath(dir, "crashed.h5")
            run_chain!(straight, 2)
            run_chain!(crashed, 0)
            h5open(crashed, "r+") do f
                METTS._write_entry!(f, "old", 1, [1.0, 2.0])  # an observable since dropped
                METTS._write_entry!(f, "n", 1, zeros(N + 3))  # another size before
            end
            @test run_chain!(crashed, 2) == 0
            @test read_metts(crashed).observables == read_metts(straight).observables
        end
    end

    @testset "an interrupted metts_start starts over" begin
        mktempdir() do dir
            straight = joinpath(dir, "straight.h5")
            crashed = joinpath(dir, "crashed.h5")
            run_chain!(straight, 2)
            h5open(crashed, "w") do f                   # died before the initial state
                f["local_states"] = local_state_strings(sites[1])
                METTS._write_entry!(f, "rng_seed", 1, UInt64(5))
            end
            @test !valid_checkpoint(crashed)
            @test run_chain!(crashed, 2) == 0
            @test read_metts(crashed).product_states == read_metts(straight).product_states
        end
    end

    @testset "read_metts ignores an unfinished step without modifying the file" begin
        mktempdir() do dir
            filename = joinpath(dir, "chain.h5")
            run_chain!(filename, 2)
            fake_crash!(filename, ["energy"])
            @test length(read_metts(filename).observables["energy"]) == 2
            h5open(filename, "r") do f
                @test size(f["energy"]) == (3,)
            end
        end
    end

    @testset "errors" begin
        mktempdir() do dir
            filename = joinpath(dir, "chain.h5")
            @test_throws Exception metts_dump_step!(filename * ".missing", (; e=1.0), initial, Xoshiro(1))
            @test_throws ErrorException metts_resume(filename * ".missing", sites)

            # a file that is not HDF5 is never taken for a missing checkpoint
            notes = joinpath(dir, "notes.h5")
            write(notes, "not an HDF5 file")
            @test_throws Exception valid_checkpoint(notes)

            metts_start(filename, sites, initial, Xoshiro(1))
            state, rng, _ = metts_resume(filename, sites)
            metts_dump_step!(filename, (; energy=1.0, n=zeros(N)), state, rng)
            rng_before = copy(rng)
            # changing observables, sizes, reserved names, bad states
            @test_throws ErrorException metts_dump_step!(filename, (; energy=1.0), state, rng)
            @test_throws ErrorException metts_dump_step!(filename, (; energy=1.0, n=zeros(N + 1)), state, rng)
            @test_throws ErrorException metts_dump_step!(filename, (; energy=1.0, n=zeros(N), rng_seed=1.0), state, rng)
            @test_throws ErrorException metts_dump_step!(filename, (; energy=1.0, n=zeros(N)), state[1:end-1], rng)
            @test_throws ErrorException metts_dump_step!(filename, (; energy=1.0, n=zeros(N)), fill(5, N), rng)
            @test read_metts(filename).nsteps == 1                  # nothing was written
            @test rng == rng_before                                 # and rng is unchanged
            # only the generator type metts_resume returns is accepted
            @test_throws MethodError metts_dump_step!(filename, (; energy=1.0, n=zeros(N)), state, MersenneTwister(1))
            # a different site type or system size
            @test_throws ErrorException metts_resume(filename, siteinds("tJ", N))
            @test_throws ErrorException metts_resume(filename, siteinds("Electron", N + 1))
            # a damaged chain is reported, never taken for a missing one
            damaged = joinpath(dir, "damaged.h5")
            run_chain!(damaged, 2)
            h5open(damaged, "r+") do f
                HDF5.set_extent_dims(f["rng_seed"], (2,))       # 3 states, 2 seeds
            end
            @test valid_checkpoint(damaged)
            @test_throws ErrorException metts_resume(damaged, sites)
            h5open(damaged, "r+") do f
                delete_object(f, "rng_seed")
            end
            @test valid_checkpoint(damaged)
            @test_throws ErrorException metts_resume(damaged, sites)
            # a wrong initial state
            @test_throws ErrorException metts_start(joinpath(dir, "new.h5"), sites, initial[1:end-1], Xoshiro(1))
            @test_throws ErrorException metts_start(joinpath(dir, "new.h5"), sites, fill(5, N), Xoshiro(1))
        end
    end
end
