###### create data-set for as specified in DNAdataset###
function dna2dummy(dna_string::String; F=Float32)
    v = Array{F,1}(undef, 4*length(dna_string));
    for (index, alphabet) in enumerate(dna_string)
        start = (index-1)*4+1;
        v[start:start+3] = dummy[uppercase(alphabet)];
    end
    return v
end

#=
get the set of dummy-vectors from a set of dna-strings
the dummy-vectors are all of same length (for now)
=#
function data_2_dummy(dna_strings::Vector{String}; F=Float32)
    how_many_strings = length(dna_strings);
    @assert how_many_strings != 0 "There aren't DNA strings found in the input";
    _len_ = length(dna_strings[1]); # length of each dna string in data    
    _S_ = Array{F, 2}(undef, (4*_len_, how_many_strings));
    for i = 1:how_many_strings
        length(dna_strings[i]) == _len_ && (@inbounds _S_[:, i] = dna2dummy(dna_strings[i]);)
    end
    return _S_
end

function reverse_complement(s::String)    
    join(islowercase(s[si]) ? s[si] : DNA_complement[s[si]] for si = length(s):-1:1)
end


function get_train_test_inds(dna_read, train_test_split_ratio, shuffle)
    len_dna_read = length(dna_read)
    how_may_entries_in_test = Int(floor((1-train_test_split_ratio)*len_dna_read));
    test_set_inds = nothing;
    if shuffle 
        test_set_inds = sample(1:len_dna_read, 
                      how_may_entries_in_test, 
                      replace=false)
    else
        test_set_inds = collect(len_dna_read-how_may_entries_in_test+1:len_dna_read)
    end
    train_set_inds = setdiff(1:len_dna_read, test_set_inds)
    return train_set_inds, test_set_inds
end


function get_data_matrices(dna_read; 
                           k=2, 
                           train_test_split_ratio=0.85,
                           FloatType=dat_t,
                           shuffle=true)
    # set train_test_split_ratio = 1.0 if no test set is needed
    train_set_inds, test_set_inds = get_train_test_inds(dna_read, train_test_split_ratio, shuffle)

    dna_read_train = @view dna_read[train_set_inds]
    dna_read_test = @view dna_read[test_set_inds]    

    shuffled_dna_read_train = seq_shuffle.(dna_read_train; k=k);
    data_matrix_train = data_2_dummy(dna_read_train; F=FloatType);
    data_matrix_bg_train = data_2_dummy(shuffled_dna_read_train; F=FloatType);

    shuffled_dna_read_test = seq_shuffle.(dna_read_test; k=k);
    data_matrix_test = data_2_dummy(dna_read_test; F=FloatType);
    data_matrix_bg_test = data_2_dummy(shuffled_dna_read_test; F=FloatType);

    # shuffled_dna_read = seq_shuffle.(dna_read; k=k);
    # data_matrix = data_2_dummy(dna_read; F=FloatType);
    # data_matrix_bg = data_2_dummy(shuffled_dna_read; F=FloatType);

    return data_matrix_train, 
           data_matrix_bg_train, 
           data_matrix_test, 
           data_matrix_bg_test,
           length(dna_read_train),
           length(dna_read_test)

end

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
