%%%%%%%
traits = {'abeta_42_beta_pval','alpha_syn_beta_pval','benton_beta_pval',...
		  'ess_beta_pval', 'gco_beta_pval', 'gds_beta_pval', 'hvlt_delayed_recall_beta_pval',...
		  'hvlt_recog_disc_index_beta_pval','hvlt_retention_beta_pval','hvlt_total_recall_beta_pval', ...
		  'lns_beta_pval','moca_beta_pval','p_tau181p_beta_pval','pigd_scores_beta_pval',...
		  'quip_beta_pval','rem_beta_pval','schwab_beta_pval','scopa_beta_pval','semantic_fluency_beta_pval',...
		  'stai_beta_pval','symbol_digit_beta_pval','total_tau_beta_pval','tremor_scores_beta_pval',...
		  'updrs1_beta_pval','updrs2_beta_pval','updrs3_beta_pval','updrs4_beta_pval'}; 
cutoffs = [0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5];%0.5% top, 1% top, ...
resultsMatrix = zeros(length(cutoffs), length(traits));

if ~exist('../GPSnet_result/', 'dir')
            mkdir('../GPSnet_result/')
end

for i = 1:length(traits)
    trait = traits{i};
    for j = 1:length(cutoffs)
        cutoff = cutoffs(j);
        alpha = 0.5;
            
        load(['../Raw_module/Raw_module_', trait, '_50']);
        ppmi_gene = Cancer_Module_Calculation(trait, alpha, cutoff);
        save_f = ['../GPSnet_result/', trait, '_cutoff_', num2str(cutoff), '_alpha_', num2str(alpha), '.txt'];
        f = fopen(save_f, 'w');
        fprintf(f, '%d\n', ppmi_gene);
        fclose(f);

        resultsMatrix(j, i) = length(ppmi_gene);
    end
end

resultsTable = array2table(resultsMatrix, 'RowNames', cellstr(num2str(cutoffs')), 'VariableNames', traits);
writetable(resultsTable, '../GPSnet_result/summary_results.csv', 'WriteRowNames', true);
