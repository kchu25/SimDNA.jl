
function data_2_dummy(dna_sim_data_vec::Vector{sim_dna_str_w_motif}; F=Float32)
    how_many_strings = length(dna_sim_data_vec);
    @assert F <: Real "input F must be subtype of Real"
    @assert how_many_strings != 0 "There aren't DNA strings found in the input";
    len = length(dna_sim_data_vec[1].str); # length of each dna string in data    
    # note that for the simulated data here, each string always has 
    # the same length
    S = Array{F, 2}(undef, (4*len, how_many_strings));
    for i = 1:how_many_strings
        @inbounds S[:, i] = dna2dummy(dna_sim_data_vec[i].str);
    end
    return S
end