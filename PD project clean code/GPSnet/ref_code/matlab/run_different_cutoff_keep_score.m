%%%%%%%
traits = {'updrs3_beta_pval'}; 
cutoffs = [0.005, 0.01];%0.5% top, 1% top, ...
resultsMatrix = zeros(length(cutoffs), length(traits));

if ~exist('../GPSnet_result_keep_score/', 'dir')
            mkdir('../GPSnet_result_keep_score/')
end

for i = 1:length(traits)
    trait = traits{i};
    for j = 1:length(cutoffs)
        cutoff = cutoffs(j);
        alpha = 0.5;
            
        load(['../Raw_module/Raw_module_', trait, '_50.mat']);
        ppmi_gene = Cancer_Module_Calculation_keep_score(trait, alpha, cutoff);
        save_f = ['../GPSnet_result_keep_score/', trait, '_cutoff_', num2str(cutoff), '_alpha_', num2str(alpha), '.txt'];
        csvwrite(save_f,ppmi_gene)
        % f = fopen(save_f, 'w');
        % fprintf(f, '%d\t%d\n', ppmi_gene);
        % fclose(f);

        resultsMatrix(j, i) = length(ppmi_gene);
    end
end

resultsTable = array2table(resultsMatrix, 'RowNames', cellstr(num2str(cutoffs')), 'VariableNames', traits);
writetable(resultsTable, '../GPSnet_result_keep_score/summary_results.csv', 'WriteRowNames', true);
