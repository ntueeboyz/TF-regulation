norm = {'without', 'ZSCORE', 'LOG'};
corr = {'Pearson', 'Spearman', 'Distance', 'MI1_sign', 'MI2_sign', 'coagreement', 'coagreement_linear', 'coagreement_sigmoid', 'coagreement_linear_new', 'coagreement_sigmoid_new'};


%% Create result table
statistic_result(1,:) = {'Data type', 'Sensitivity', 'Specificity', 'F1 score for upregulation', 'F1 score for downregulation','AUC', 'Precision-positive', 'Precision-negtive', 'Data type', 'Sensitivity', 'Specificity', 'F1 score for upregulation', 'F1 score for downregulation','AUC', 'Precision-positive', 'Precision-negtive'};
idx = 2; %index of row number


for j = 1:length(norm)
    for k = 1:length(corr)
        name = [norm{j},'_',corr{k},'.mat'];
        load (name)
                sensitivity_D = sum([resultsD{2:end,5}])/length([resultsD{2:end,5}]); %fraction of catched downregulation
                sensitivity_U = sum([resultsU{2:end,5}])/length([resultsU{2:end,5}]); %fraction of catched upregulation
                precision_D = sum([resultsD{2:end,5}])/(sum([resultsD{2:end,5}])+sum([resultsU{2:end,5}] == 0)); %Precision of downregulation
                precision_U = sum([resultsU{2:end,5}])/(sum([resultsU{2:end,5}])+sum([resultsD{2:end,5}] == 0)); %Precision of downregulation
                F1_score_D = 2*(precision_D*sensitivity_D/(sensitivity_D+precision_D));
                F1_score_U = 2*(precision_U*sensitivity_U/(sensitivity_U+precision_U));
                auc = AUC(sensitivity_D, sensitivity_U);
                data_name = [norm{j},'_',corr{k}];
                
                statistic_result{idx, 1} = data_name;
                statistic_result{idx, 2} = sensitivity_U;
                statistic_result{idx, 3} = sensitivity_D;
                statistic_result{idx, 4} = F1_score_U;
                statistic_result{idx, 5} = F1_score_D;
                statistic_result{idx, 6} = auc;
                statistic_result{idx, 7} = precision_U;
                statistic_result{idx, 8} = precision_D;
                
                
                % Threshold
                sensitivity_threshold = sum([resultsU_threshold{2:end,5}])/length([resultsU_threshold{2:end,5}]);
                specificity_threshold = sum([resultsD_threshold{2:end,5}])/length([resultsD_threshold{2:end,5}]);
                precision_U = sum([resultsU_threshold{2:end,5}])/(sum([resultsU_threshold{2:end,5}])+sum([resultsD_threshold{2:end,5}] == 0));
                precision_D = sum([resultsD_threshold{2:end,5}])/(sum([resultsD_threshold{2:end,5}])+sum([resultsU_threshold{2:end,5}] == 0));
                F1_score_U = 2*(precision_U*sensitivity_threshold/(precision_U+sensitivity_threshold));
                F1_score_D = 2*(precision_D*specificity_threshold/(precision_D+specificity_threshold));
                auc_threshold = AUC(sensitivity_threshold, specificity_threshold);
                data_name_threshold = ['0.05_',norm{j},'_',corr{k}];
                
                statistic_result{idx, 9} = data_name_threshold;
                statistic_result{idx, 10} = sensitivity_threshold;
                statistic_result{idx, 11} = specificity_threshold;
                statistic_result{idx, 12} = F1_score_U;
                statistic_result{idx, 13} = F1_score_D;
                statistic_result{idx, 14} = auc_threshold;
                statistic_result{idx, 15} = precision_U;
                statistic_result{idx, 16} = precision_D;
                statistic_result{idx, 17} = length([resultsU_threshold{2:end, 5}]);
                statistic_result{idx, 18} = length([resultsD_threshold{2:end, 5}]);
                
                idx = idx + 1;
                
    end
end

save('result.mat', 'statistic_result');

function auc = AUC(sensitivity, specificity)
if isnan(sensitivity) || isnan(specificity)
    auc = NaN;
else
    y = [0;sensitivity;1];
    x = [0;1-specificity;1];
    auc = trapz(x,y);
end
end
