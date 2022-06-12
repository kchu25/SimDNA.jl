module SimDNA

using CUDA, Distributions

export Sim_DNA, 
       single_block_motif,
       gapped_k_block_motif,
       mixture_k_block_motif,
       mixture_gapped_k_block_motif


# println("ho")

include("ground_truth_setup.jl")
include("sample.jl")
include("helpers.jl")
include("save.jl")
include("sim_wrapper.jl")



end
