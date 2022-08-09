function save_gt_motifs(target_folder::String, motif::mixture_gapped_k_block_motif)    
    for j = 1:length(motif.motif.P_motifs.P)
        gt = reduce(hcat, [motif.motif.P_motifs.P[j].P.pc[i].p for i = 1:length(motif.motif.P_motifs.P[j].P.pc)]);
        CSV.write(target_folder*"/pfm_gt$j.csv",  Tables.table(gt), writeheader=false);
    end     
end

function save_gt_motifs(target_folder::String, motif::gapped_k_block_motif)    
    for j = 1:motif.K
        gt = reduce(hcat, [motif.P_motifs.P[j].P.pc[i].p for i = 1:length(motif.P_motifs.P[j].P.pc)]);
        CSV.write(target_folder*"/pfm_gt$j.csv",  Tables.table(gt), writeheader=false);
    end
end

function save_gt_motifs(target_folder::String, motif::mixture_k_block_motif)    
    for j = 1:motif.num_modes
        gt = reduce(hcat, [motif.motif.P.pc[i].p for i = motif.modes[j]]);
        CSV.write(target_folder*"/pfm_gt$j.csv",  Tables.table(gt), writeheader=false);
    end
end

function save_gt_motifs(target_folder::String, motif::single_block_motif)    
    gt = reduce(hcat, [motif.P.pc[i].p for i = 1:length(motif.P.pc)]);
    CSV.write(target_folder*"/pfm_gt1.csv",  Tables.table(gt), writeheader=false);
end

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

function save_sim_data_as_fasta_test(target_folder::String, 
                                raw_data, 
                                N::Integer)
    !isdir(target_folder) && mkpath(target_folder);

    strings = [raw_data[i].str for i = 1:N];
    strings_all_upcase = [uppercase.(raw_data[i].str) for i = 1:N];
    
    open(target_folder*"/data_w_answer_test.fa","w") do file
        for (ind,s) in enumerate(strings)
            write(file, string(">sequence_", string(ind),"\n"));
            write(file, string(s,"\n"))
        end
    end
    open(target_folder*"/data_test.fa","w") do file
        for (ind,s) in enumerate(strings_all_upcase)
            write(file, string(">sequence_", string(ind),"\n"));
            write(file, string(s,"\n"))
        end
    end
end