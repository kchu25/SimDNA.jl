module SimDNA

using Distributions, CUDA

include("src/ground_truth_setup.jl")
include("src/sample.jl")
include("src/helpers.jl")
include("src/save.jl")
include("src/sim_wrapper.jl")




end
