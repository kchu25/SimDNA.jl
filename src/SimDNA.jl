module SimDNA

using Distributions, CUDA

export Sim_DNA

include("ground_truth_setup.jl")
include("sample.jl")
include("helpers.jl")
include("save.jl")
include("sim_wrapper.jl")



end
