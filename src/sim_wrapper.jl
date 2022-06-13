struct Sim_DNA{T <: Integer, S <: Real} 
    raw_data::Vector{sim_dna_str_w_motif}
    motif::motif_type
    N::T
    L::T
    data_matrix::Array{S,3}
    data_matrix_gpu::Union{Nothing, CuArray{S,3},  CuArray{S,2}}
    data_matrix_bg::Union{Nothing, Array{S,3}, Array{S,2}}
    target_folder::Union{Nothing, String}
    prob_per_seq::S

    # function Sim_DNA{T, S}(motif::motif_type, 
    #                     num_data_pts::Integer,
    #                     output_folder::String,
    #                     data_pt_len::Integer=100,
    #                     bern_prob::Real=1.0
    #     ) where {T <: Integer, S <: Real}        
    #     # simulate the data
    #     raw_data = sample_backgound_with_motif_multiple(motif,num_data_pts,data_pt_len,bern_prob);
    #     raw_data_DNA_str = [uppercase.(i.str) for i in  raw_data];
    #     # save the simulated data as a fasta file in the target folder
    #     save_sim_data_as_fasta(output_folder, raw_data, num_data_pts, motif);
    #     data_matrix, data_matrix_bg = get_data_matrices(raw_data_DNA_str; FloatType=S);
    #     new(
    #         raw_data,
    #         motif,
    #         T(num_data_pts),
    #         T(size(data_matrix,1)/4),
    #         data_matrix,
    #         cu(data_matrix),
    #         data_matrix_bg,
    #         output_folder,
    #         S(bern_prob)            
    #     )
    # end # update this later
    
    function Sim_DNA{T,S}(motif::motif_type; 
                          num_seqs::Integer=1000,
                          seq_len::Integer=100,
                          bern_prob::Real=1.0) where {T <: Integer, S <: Real}     
        raw_data = sample_backgound_with_motif_multiple(motif,num_seqs,seq_len,bern_prob);
        raw_data_DNA_str = [uppercase.(i.str) for i in  raw_data];
        data_matrix, data_matrix_bg = get_data_matrices(raw_data_DNA_str; FloatType=S);
        L, N = size(data_matrix);
        data_matrix = reshape(data_matrix, (L,1,N));
        new(
            raw_data, motif, T(N), T(L/4),
            data_matrix, cu(data_matrix), data_matrix_bg, nothing, S(bern_prob)
        )                          
    end
end

function reshape_for_search!(s::Sim_DNA)
    s.data_matrix = reshape(s.data_matrix, (s.L,s.N));
    s.data_matrix_gpu = reshape(s.data_matrix_gpu, (s.L,S.N));
end