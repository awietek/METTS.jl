using Test
using HDF5
using ITensors

@testset "checkpoint.jl" begin

    @testset "count_existing_steps" begin

        @testset "returns 0 if file does not exist" begin
            @test count_existing_steps("nonexistent.hdf5") == 0
        end

        @testset "returns 0 if file is empty" begin
            mktempdir() do dir
                filename = joinpath(dir, "empty.hdf5")
                h5open(filename, "w") do f end
                @test count_existing_steps(filename) == 0
            end
        end

        @testset "deletes file and returns 0 if product_state dataset is missing" begin
            mktempdir() do dir
                filename = joinpath(dir, "no_product_state.hdf5")
                h5open(filename, "w") do f
                    f["energy"] = [1.0, 2.0, 3.0]
                end
                @test_warn "product_state" count_existing_steps(filename) == 0
                @test !isfile(filename)
            end
        end

        @testset "deletes file and returns 0 if file is corrupt" begin
            mktempdir() do dir
                filename = joinpath(dir, "corrupt.hdf5")
                write(filename, "this is not a valid hdf5 file")
                @test_warn "Failed to read" count_existing_steps(filename) == 0
                @test !isfile(filename)
            end
        end
        @testset "returns correct number of steps" begin
            mktempdir() do dir
                filename = joinpath(dir, "valid.hdf5")
                h5open(filename, "w") do f
                    f["product_state"] = ["Up", "Dn", "Up", "Dn", "Up"]
                end
                @test count_existing_steps(filename) == 5
            end
        end


    end

    @testset "read_last_product_state" begin

        @testset "errors if file does not exist" begin
            @test_throws ErrorException read_last_product_state("nonexistent.hdf5")
        end

        @testset "errors if product_state dataset is missing" begin
            mktempdir() do dir
                filename = joinpath(dir, "no_product_state.hdf5")
                h5open(filename, "w") do f
                    f["energy"] = [1.0, 2.0, 3.0]
                end
                @test_throws ErrorException read_last_product_state(filename)
            end
        end

        @testset "errors if product_state dataset is empty" begin
            mktempdir() do dir
                filename = joinpath(dir, "empty_product_state.hdf5")
                h5open(filename, "w") do f
                    f["product_state"] = String[]
                end
                @test_throws ErrorException read_last_product_state(filename)
            end
        end

        @testset "returns last product state" begin
            mktempdir() do dir
                filename = joinpath(dir, "valid.hdf5")
                states = ["Up", "Dn", "Up", "Dn", "Up"]
                h5open(filename, "w") do f
                    f["product_state"] = states
                end
                @test read_last_product_state(filename) == "Up"
            end
        end

    end

end