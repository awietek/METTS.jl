using Test
using Dumper
using HDF5
using ITensors

@testset "checkpoint.jl" begin
    using Test
    using HDF5

    @testset "count_existing_steps" begin

        @testset "returns 0 if file does not exist" begin
            @test count_existing_steps("nonexistent.hdf5", "product_state") == 0
        end

        @testset "returns 0 and warns if file is empty" begin
            mktempdir() do dir
                filename = joinpath(dir, "empty.hdf5")
                h5open(filename, "w") do f
                end

                @test_logs (:warn, r"Checkpoint file exists but is empty") begin
                    @test count_existing_steps(filename, "product_state") == 0
                end
            end
        end

        @testset "returns correct number of steps with one step" begin
            mktempdir() do dir
                filename = joinpath(dir, "valid.hdf5")
                h5open(filename, "w") do f
                    f["product_state"] = [1, 1, 2, 1, 2]
                end
                @test count_existing_steps(filename, "product_state") == 1
            end
        end

        @testset "returns correct number of steps with multiple steps" begin
            mktempdir() do dir
                filename = joinpath(dir, "valid.hdf5")
                outfile = DumpFile(filename)

                dump!(outfile, "product_state", [1, 1, 2, 1, 2])
                dump!(outfile, "product_state", [1, 1, 2, 1, 2])
                dump!(outfile, "product_state", [1, 2, 2, 1, 2])

                @test count_existing_steps(filename, "product_state") == 3
            end
        end

        @testset "throws a specific error if product_state does not exist but other datasets do" begin
            mktempdir() do dir
                filename = joinpath(dir, "test.hdf5")
                outfile = DumpFile(filename)

                dump!(outfile, "energy", [0.0])
                expected_err = "Checkpoint file '$filename' contains data (energy) but is missing the mandatory 'product_state' dataset."
                @test_throws expected_err count_existing_steps(filename, "product_state")
            end
        end

        @testset "throws a specific error if the array is not 2D" begin
            mktempdir() do dir
                filename = joinpath(dir, "invalid_shape.hdf5")

                h5open(filename, "w") do f
                    f["product_state"] = zeros(Int, (3, 3, 3))
                end

                expected_err = r"Data shape mismatch.*expected a 1D or 2D array"
                @test_throws expected_err count_existing_steps(filename, "product_state")
            end
        end

        @testset "logs an error and rethrows if file is corrupted (catch block)" begin
            mktempdir() do dir
                filename = joinpath(dir, "corrupted.hdf5")

                write(filename, "Test corrupted.")

                @test_logs (:error, r"Failed to count existing steps due to an error:") match_mode = :any begin
                    @test_throws Exception count_existing_steps(filename, "product_state")
                end
            end
        end

    end

    @testset "read_last_product_state" begin

        @testset "returns the correct vector with only one array saved" begin
            mktempdir() do dir
                filename = joinpath(dir, "test.hdf5")
                h5open(filename, "w") do f
                    f["product_state"] = [1, 1, 2, 1, 2]
                end
                @test read_last_product_state(filename, "product_state") == [1, 1, 2, 1, 2]
            end
        end

        @testset "returns the correct vector with multiple arrays saved" begin
            mktempdir() do dir
                filename = joinpath(dir, "test.hdf5")
                outfile = DumpFile(filename)

                dump!(outfile, "product_state", [1, 1, 1, 1, 2])
                dump!(outfile, "product_state", [1, 2, 2, 2, 2])
                dump!(outfile, "product_state", [1, 2, 1, 2, 1])

                @test read_last_product_state(filename, "product_state") == [1, 2, 1, 2, 1]
            end
        end
    end


    @testset "read_last_product_state" begin

        @testset "throws error if file does not exist" begin
            @test_throws "File 'nonexistent.hdf5' not found." read_last_product_state("nonexistent.hdf5", "product_state")
        end

        @testset "throws error if dataset is missing" begin
            mktempdir() do dir
                filename = joinpath(dir, "missing_data.hdf5")
                outfile = DumpFile(filename)

                dump!(outfile, "energy", [0.0])

                expected_err = "No 'product_state' dataset found in '$filename'."
                @test_throws expected_err read_last_product_state(filename, "product_state")
            end
        end

        @testset "throws error if dataset is empty" begin
            mktempdir() do dir
                filename = joinpath(dir, "empty_data.hdf5")
                h5open(filename, "w") do f
                    f["product_state"] = Int[]
                end

                expected_err = "Empty 'product_state' dataset in '$filename'."
                @test_throws expected_err read_last_product_state(filename, "product_state")
            end
        end

        @testset "reads correct state for a single step" begin
            mktempdir() do dir
                filename = joinpath(dir, "single_step.hdf5")
                outfile = DumpFile(filename)

                state = [1, 1, 2, 1, 2]
                dump!(outfile, "product_state", state)

                @test read_last_product_state(filename, "product_state") == state
            end
        end

        @testset "reads correct last state for multiple steps (2D array)" begin
            mktempdir() do dir
                filename = joinpath(dir, "multiple_steps.hdf5")
                outfile = DumpFile(filename)

                step1 = [1, 1, 2, 1, 2]
                step2 = [2, 2, 3, 2, 3]
                step3 = [0, 1, 0, 1, 0]

                dump!(outfile, "product_state", step1)
                dump!(outfile, "product_state", step2)
                dump!(outfile, "product_state", step3)

                @test read_last_product_state(filename, "product_state") == step3
            end
        end

        @testset "throws error if array shape is greater than 2D" begin
            mktempdir() do dir
                filename = joinpath(dir, "invalid_shape.hdf5")
                h5open(filename, "w") do f
                    f["product_state"] = zeros(Int, 2, 2, 2)
                end

                expected_err = r"Data shape mismatch.*expected a 1D or 2D array"
                @test_throws expected_err read_last_product_state(filename, "product_state")
            end
        end

        @testset "logs error and rethrows if file is corrupted (catch block)" begin
            mktempdir() do dir
                filename = joinpath(dir, "corrupted.hdf5")
                write(filename, "Not a valid HDF5 file.")

                @test_logs (:error, r"Failed to read product state due to an error:") match_mode = :any begin
                    @test_throws Exception read_last_product_state(filename, "product_state")
                end
            end
        end

    end
end