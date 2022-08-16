# mutable so that shape of data can be changed 
mutable struct Sim_DNA{T <: Integer, S <: Real} 
    raw_data::Vector{sim_dna_str_w_motif}
    raw_data_test::Vector{sim_dna_str_w_motif}
    motif::motif_type
    N::T
    L::T
    data_matrix::Union{Array{S,3}, Array{S,2}}
    data_matrix_bg::Union{Nothing, Array{S,3}, Array{S,2}}
    target_folder::Union{Nothing, String}
    prob_per_seq::S
    data_matrix_test::Union{Nothing, Array{S,3}, Array{S,2}}
    data_matrix_bg_test::Union{Nothing, Array{S,3}, Array{S,2}}
    N_test::T

    function Sim_DNA{T, S}(motif::motif_type, 
                        num_data_pts::Integer,
                        output_folder::String,
                        data_pt_len::Integer=100,
                        bern_prob::Real=1.0; 
                        k=2, # kmer frequency in the test set
                        train_test_split_ratio=0.85,
                        save=false,
                        shuffle=true
        ) where {T <: Integer, S <: Real}        
        # simulate the data
        raw_data = sample_backgound_with_motif_multiple(motif,num_data_pts,data_pt_len,bern_prob);
        raw_data_DNA_str = [uppercase.(i.str) for i in  raw_data];

        data_matrix, data_matrix_bg, data_matrix_test, data_matrix_bg_test, 
            N_train, N_test, train_set_inds, test_set_inds = 
            get_data_matrices(raw_data_DNA_str; k=k, 
                              train_test_split_ratio=train_test_split_ratio, 
                              FloatType=S, 
                              shuffle=shuffle
                              );
        # save the simulated data as a fasta file in the target folder
        # just save the training set
        if save 
            save_sim_data_as_fasta(output_folder, 
                                       raw_data[train_set_inds], 
                                       N_train, 
                                       motif);
            save_sim_data_as_fasta_test(output_folder, 
                                        raw_data[test_set_inds], 
                                        N_test)
        end

        L, _ = size(data_matrix);
        data_matrix = reshape(data_matrix, (L,1,N_train));
        train_test_split_ratio < 1.0 && (data_matrix_test = reshape(data_matrix_test, (L,1, N_test));)
        new(
            raw_data[train_set_inds],
            raw_data[test_set_inds],
            motif,
            T(N_train),
            T(L/4),
            data_matrix,
            data_matrix_bg,
            output_folder,
            S(bern_prob),
            data_matrix_test,
            data_matrix_bg_test,
            T(N_test)
        )
    end # update this later
    
    function Sim_DNA{T,S}(motif::motif_type; 
                          num_seqs::Integer=1000,
                          seq_len::Integer=100,
                          bern_prob::Real=1.0,
                          k=2, # kmer frequency in the test set
                          train_test_split_ratio=0.85,
                          shuffle=true
                          ) where {T <: Integer, S <: Real}     
        raw_data = sample_backgound_with_motif_multiple(motif,num_seqs,seq_len,bern_prob);
        raw_data_DNA_str = [uppercase.(i.str) for i in  raw_data];
        data_matrix, data_matrix_bg = get_data_matrices(raw_data_DNA_str; FloatType=S);

        data_matrix, data_matrix_bg, data_matrix_test, data_matrix_bg_test, 
            N_train, N_test, train_set_inds, test_set_inds = 
                get_data_matrices(raw_data_DNA_str; k=k, 
                                train_test_split_ratio=train_test_split_ratio, 
                                FloatType=S, 
                                shuffle=shuffle
                                );

        L, _ = size(data_matrix);
        data_matrix = reshape(data_matrix, (L, 1, N_train));
        train_test_split_ratio < 1.0 && (data_matrix_test = reshape(data_matrix_test, (L, 1, N_test));)
        new(
            raw_data[train_set_inds], 
            raw_data[test_set_inds],
            motif, 
            T(N_train), 
            T(L/4),
            data_matrix, 
            data_matrix_bg, 
            nothing, 
            S(bern_prob),
            data_matrix_test,
            data_matrix_bg_test,
            T(N_test)
        )
    end
end


function reshape_for_search!(s::Sim_DNA)
    s.data_matrix = reshape(s.data_matrix, (s.L,s.N_train));
    s.data_matrix_gpu = reshape(s.data_matrix_gpu, (s.L,S.N_train));
end