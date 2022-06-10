function save_sim_data_as_fasta(target_folder::String, 
                                raw_data, 
                                N::Integer,
                                motif)
    !isdir(target_folder) && mkpath(target_folder);
    save_gt_motifs(target_folder, motif);

    strings = [raw_data[i].str for i = 1:N];
    strings_all_upcase = [uppercase.(raw_data[i].str) for i = 1:N];
    
    open(target_folder*"/data_w_answer.fa","w") do file
        for (ind,s) in enumerate(strings)
            write(file, string(">sequence_", string(ind),"\n"));
            write(file, string(s,"\n"))
        end
    end
    open(target_folder*"/data.fa","w") do file
        for (ind,s) in enumerate(strings_all_upcase)
            write(file, string(">sequence_", string(ind),"\n"));
            write(file, string(s,"\n"))
        end
    end
end