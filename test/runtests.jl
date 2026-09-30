using METTS
using Test

@testset "METTS.jl" begin
    include("test_local_state.jl")
    include("test_random_product_state.jl")
    include("test_checkpoint.jl")
    include("test_collapse.jl")
    include("test_timeevo.jl")
    include("test_product_states.jl")
    include("test_metts_checkpoint.jl")
end
