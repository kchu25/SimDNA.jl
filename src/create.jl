function simulate_two_block_overlapping_motifs(output_folder::String; bcn=8:22, 
                N=2000, len=100, bern_prob=1.0, train_test_split_ratio=1.0)
    # bcn: block column number
    c1, c2 = rand(bcn), rand(bcn)
    overlap_ratio = rand(0.1:0.01:0.9)
    c1_times_or = Int(ceil(c1*overlap_ratio))
    overlapped_columns = min(c1_times_or, c2)
    c2_start = c1-c1_times_or+1;
    block_modes = [1:c1, c2_start:c2_start+c2-1]
    _m_ = mixture_k_block_motif(block_modes);      
    data = Sim_DNA{int_t, dat_t}(_m_, N, output_folder, len, bern_prob; save=true, train_test_split_ratio=train_test_split_ratio);
    @info "block lengths: ($c1, $c2)"
    @info "olap ratio: $overlap_ratio"
    @info "olap columns: $overlapped_columns"
    @info "block modes: ($block_modes)"
    return (c1,c2), overlap_ratio, overlapped_columns, block_modes, data
end

function simulate_two_block_gapped_motifs(output_folder::String;
                                          bcn=8:22, 
                                          N=2000, 
                                          train_test_split_ratio=1.0,
                                          bern_prob=1.0,
                                          gap_range=1:10,
                                          save=true,
                                          len=100)
    c1, c2 = rand(bcn), rand(bcn)
    gap = rand(gap_range)
    _m_ = gapped_k_block_motif([c1,c2],[gap]);
    data = Sim_DNA{int_t, dat_t}(_m_, N, output_folder, len, bern_prob; train_test_split_ratio=train_test_split_ratio, save=save);
    @info "block lengths: ($c1, $c2)"
    @info "gap: $gap"
    return (c1,c2), gap, data
end