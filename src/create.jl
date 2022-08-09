
function simulate_two_block_overlapping_motifs(output_folder::String; 
                                               bcn=8:22, 
                                               N=2000, 
                                               len=100, 
                                               bern_prob=1.0,
                                               save=true
                                               )
    # bcn: block column number
    # nm: number of modes
    c1, c2 = rand(bcn), rand(bcn)
    overlapped_columns = Int(ceil(min(c1,c2)*rand()))
    c2_start = c1-overlapped_columns+1;
    # [1:c1, c2_start:c2_start+c2-1]
    block_modes = [1:c1, c2_start:c2_start+c2-1]
    _m_ = mixture_k_block_motif(block_modes);      
    data = Sim_DNA{int_t, dat_t}(_m_, N, output_folder, len, bern_prob; save=save);
    @info "block lengths: ($c1, $c2)"
    @info "olap: $overlapped_columns"
    @info "block modes: ($block_modes)"
    return data
end

function simulate_two_block_gapped_motifs(output_folder::String;
                                          bcn=8:22, 
                                          N=2000, 
                                          len=100)
    c1, c2 = rand(bcn), rand(bcn)
    gap = rand(0:10)
    _m_ = gapped_k_block_motif([c1,c2],[gap]);
    data = Sim_DNA{int_t, dat_t}(_m_, N, output_folder, len, bern_prob; save=save);
    @info "block lengths: ($c1, $c2)"
    @info "gap: $gap"
    return data
end