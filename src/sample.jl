len(m::single_part_motif) = m.len;
# view(m::single_part_motif) = @pt :header = DNA_char_cap view(m.P); 
sample_m(m::single_part_motif)::String = join([DNA_char_cap[i] for i in rand(m.P)]);
sample_whole(m::single_part_motif) = (motif=sample_m(m), mode=1);

function sample_contiguous_subset(m::single_part_motif, r::UnitRange{Int})::String
    if maximum(r) > m.len || minimum(r) < 1
        error("unit range specified is not covered by 1:$(m.len)");
    end
    sample_m(m)[r] # make it Vector{String} for concatenation purposes later (sample fake data)
end

view(m::k_parts_motif) = for v in m.P view(v) end;
sample_m(m::k_parts_motif)::String = join([sample_m(v) for v in m.P]);
sample_vec(m::k_parts_motif) = [sample_m(v) for v in m.P];
sample_whole(m::k_parts_motif) = (motif=sample_m(m), mode=1);

# distribution on the number of gap alphbet for now is assumed to be discrete uniform (for now)
function sample_m(m::gapped_k_parts_motif)::String
    s = [];
    parts = sample_vec(m.P_motifs);
    for i = 1:(m.K-1)
        push!(s, parts[i]);        
        num_alphabet = rand(DiscreteUniform(0, m.gap_len[i]));
        push!(s, join([DNA_char_nocap[rand(m.P_gap)] for j = 1:num_alphabet]));
    end
    push!(s,parts[m.K]);
    return join(s);
end

sample_whole(m::gapped_k_parts_motif) = (motif=sample_m(m), mode=1);

#=
sample a (contiguous) subset of gapped_p_fam_k_motif_PC 
according to a unitrange r
=#

function sample_contiguous_subset(m::gapped_k_parts_motif, r::UnitRange{Int})::String
    r_max = maximum(r); r_min = minimum(r);
    if r_max > m.K || r_min < 1
        error("unit range specified is not covered by 1:$(m.len)");
    end
    s = [];
    parts = sample_vec(m.P_motifs);
    for i = r_min:(r_max-1)
        push!(s, parts[i]);
        num_alphabet = rand(DiscreteUniform(0,m.gap_len[i]));
        push!(s, join([DNA_char_nocap[rand(m.P_gap)] for j = 1:num_alphabet]));
    end
    push!(s, parts[r_max]);
    return join(s);
end

view(m::gapped_k_parts_motif) = view(m.P_motifs);

# example ----------------------------------
# a = gapped_k_parts_motif([3,3,3],[2,2]);
# sample_m(a)
# view(a)
# ------------------------------------------

# Sample the ith mode of the motif, a gapped_p_fam_k_motif

function sample_m(m::mixture_gapped_k_parts_motifs)::String
    i = rand(m.mixture_weights);
    sample_contiguous_subset(m.motif, m.modes[i]); 
end

function sample_whole(m::mixture_gapped_k_parts_motifs)
    i = rand(m.mixture_weights);
    str = sample_contiguous_subset(m.motif, m.modes[i]);
    (motif=str, mode=i);
end

# example ----------------------------------
# a=mixture_gapped_k_parts_motifs([7,7,7],[4,4],[1:2,2:3]);
# a.modes
# a.num_modes
# a.mixture_weights.p
# rand(a.mixture_weights)
# sample_m(a)

# b = mixture_gapped_k_parts_motifs([1:9, 6:12, 8:15], "uniform");
# b.modes
# b.num_modes
# b.mixture_weights.p
# rand(b.mixture_weights)
# sample_m(b)
# ------------------------------------------

#=
Sample the ith mode of the motif, a mixture_p_fam_motifs
=#
function sample_m(m::mixture_k_parts_motifs)::String
    i = rand(m.mixture_weights);
    sample_contiguous_subset(m.motif, m.modes[i]); 
end

function sample_whole(m::mixture_k_parts_motifs)
    i = rand(m.mixture_weights);
    str = sample_contiguous_subset(m.motif, m.modes[i]);
    # println("mode $i");
    (motif=str, mode=i);
end

# example ----------------------------------
# a = mixture_k_parts_motifs([1:9, 3:8, 6:15]);
# a.modes
# a.num_modes
# a.mixture_weights.p
# rand(a.mixture_weights)
# sample_m(a)
# sample_whole(a)

# b = mixture_gapped_p_fam_k_motifs([1:7, 6:8, 6:10], "uniform");
# b.modes
# b.num_modes
# b.mixture_weights.p
# rand(b.mixture_weights)
# sample_m(b)
# ------------------------------------------


# Sample a background-string of length len (default: len = data_pt_len)
sample_background_single(len::Integer = data_pt_len)::String = 
        join([DNA_char_nocap[rand(background)] for _ in 1:len]);


#=
sample a background string of length len (default: len = data_pt_len)
with a motif embedded in it; 

returns a named-tuple:
(str=dna_string, motif_where=motif_start:motif_end, mode=binding_mode), 
which has a motif-embedded background string , an UnitRange that 
indicates where motif is embedded, and the mode of binding
=#

function sample_backgound_with_motif_single(motif::motif_type,
                        len::Integer,
                        bern_appear
                        )::sim_dna_str_w_motif
    if rand(bern_appear)
        motif_realization = sample_whole(motif);

        ########### either make it complement or not ###########    
        motif_str = motif_realization.motif; motif_mode = motif_realization.mode;
        is_complement = rand() > percent_complement;
        motif_realization = is_complement ? (motif=reverse_complement(motif_str),mode=motif_mode) : motif_realization;
        ########################################################

        motif_len = length(motif_realization.motif);
        motif_len > len && error("Please use a motif that has width shorter than $len");

        avail_positions::Int = (len-motif_len+1)-sample_near_center;
        at::Int = rand(DiscreteUniform(1+sample_near_center, avail_positions));
        before_at::Int = (at-1);
        dna_str = join([sample_background_single(before_at), 
            motif_realization.motif::String,
            sample_background_single((len-(before_at+motif_len)))]);
        
        return  (str=dna_str, 
            motif_where=at:(at+motif_len-1), 
            mode=motif_realization.mode, 
            complement=is_complement);
    else
        return (str=sample_background_single(len), 
                motif_where=0:0,
                mode=0,
                complement=false);
    end
end

# example ----------------------------------
# a = gapped_k_parts_motif([3,3,3],[2,2]);
# sample_a = sample_backgound_with_motif_single(a);
# sample_a
# sample_a.mode
# sample_a[1][sample_a.motif_where]

# b = mixture_gapped_k_parts_motifs([3,3,3],[2,2],[1:2,2:3]);
# sample_b = sample_backgound_with_motif_single(b);
# sample_b
# sample_b.mode
# sample_b[1][sample_b.motif_where]
# ------------------------------------------

#=
sample multiple background string of length len (default: len = data_pt_len)
with a motif embedded in it; see sample_backgound_with_motif_single. 
Instead of returning an instance of sim_dna_str_w_motif, it returns an 
instance of Vector{sim_dna_str_w_motif}.
=#

function sample_backgound_with_motif_multiple(
                            motif::motif_type, 
                            how_many::Integer,
                            len::Integer=100, # default to 100 bp for each sequence
                            bern_prob::Real=1.0 # default to every sequence to contain motif(s)                    
                            )::Vector{sim_dna_str_w_motif}
    @assert  0 ≤ bern_prob ≤ 1 "probability for Bernoulli should be in between 0 and 1"
    bern_appear = Bernoulli(bern_prob); 
    return [sample_backgound_with_motif_single(motif, len, bern_appear)
        for _ in 1:how_many];
end

# example ----------------------------------
# b = mixture_gapped_k_parts_motifs([3,3,3],[2,2],[1:2,2:3]);
# sample_b = sample_backgound_with_motif_multiple(b, 5);
# sample_b
# ------------------------------------------
