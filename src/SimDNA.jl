module SimDNA

using CUDA, Distributions, SeqShuffle, DataFrames, CSV

export Sim_DNA, 
       single_block_motif,
       gapped_k_block_motif,
       mixture_k_block_motif,
       mixture_gapped_k_block_motif,
       motif_type,
       sim_dna_str_w_motif,
       simulate_two_block_overlapping_motifs,
       simulate_two_block_gapped_motifs


# println("ho")

include("ground_truth_setup.jl")
include("sample.jl")
include("helpers.jl")
include("save.jl")
include("sim_wrapper.jl")
include("create.jl")



end
