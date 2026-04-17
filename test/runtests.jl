using METTS
using Test

@testset "METTS.jl" begin
    include("test_local_state.jl")
    include("test_random_product_state.jl")
end
