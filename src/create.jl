function simulate_two_block_overlapping_motifs(;bcn=8:22, N=2000, len=100)
    # bcn: block column number
    c1, c2 = rand(bcn), rand(bcn)
    overlap_ratio = round(rand(),digits=1)
    overlapped_columns = Int(ceil(c1*overlap_ratio))
    c2_start = c1-overlapped_columns+1;
    block_modes = [1:c1, c2_start:c2_start+c2-1]
    _m_ = mixture_k_block_motif(block_modes);      
    data = Sim_DNA{int_t, dat_t}(_m_; num_seqs=N, seq_len=len);
    @info "block lengths: ($c1, $c2)"
    @info "olap ratio: $overlap_ratio"
    @info "olap: $overlapped_columns"
    @info "block modes: ($block_modes)"
    return data
end

function simulate_two_block_gapped_motifs(output_folder::String;
                                          bcn=8:22, 
                                          N=2000, 
                                          train_test_split_ratio=1.0,
                                          len=100)
    c1, c2 = rand(bcn), rand(bcn)
    gap = rand(0:10)
    _m_ = gapped_k_block_motif([c1,c2],[gap]);
    data = Sim_DNA{int_t, dat_t}(_m_, N, output_folder, len, bern_prob;train_test_split_ratio=train_test_split_ratio, save=save);
    @info "block lengths: ($c1, $c2)"
    @info "gap: $gap"
    return data
end