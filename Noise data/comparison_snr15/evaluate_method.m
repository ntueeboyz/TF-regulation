function evaluate_method (gold_list, method_result)
[num,TF_regulation,~] = xlsread(gold_list);
regulation_list = TF_regulation(2:end,1:3);

resultsU(1,:)={'TF','Gene','Type of regulation (exp)', 'value','Catched?'};
resultsD(1,:)={'TF','Gene','Type of regulation (exp)', 'value','Catched?'};
idx_U = 2;
idx_D = 2;
for i = 2:size(method_result,1)
    
    if method_result{i,3} > 0
        trend = 'upregulation';
    elseif method_result{i,3} < 0
        trend = 'downregulation';
    elseif method_result{i,3} == 0
        trend = 'NAN';
    end
    
    TF = method_result{i,1};
    gene = method_result{i,2};
    idx_TF = find(strcmp(TF, regulation_list(:,1)));
    
    for j = 1:length(idx_TF)
        if strcmp(regulation_list{idx_TF(j), 2},gene)
            if strcmp(regulation_list{idx_TF(j), 3}, 'upregulation')
                resultsU{idx_U, 1} = TF;
                resultsU{idx_U, 2} = gene;
                resultsU{idx_U, 3} = regulation_list{idx_TF(j), 3};
                resultsU{idx_U, 4} = method_result{i, 3};
                resultsU{idx_U, 5} = strcmp(regulation_list{idx_TF(j), 3},trend);
                idx_U = idx_U + 1;
            elseif strcmp(regulation_list{idx_TF(j), 3}, 'downregulation')
                resultsD{idx_D, 1} = TF;
                resultsD{idx_D, 2} = gene;
                resultsD{idx_D, 3} = regulation_list{idx_TF(j), 3};
                resultsD{idx_D, 4} = method_result{i, 3};
                resultsD{idx_D, 5} = strcmp(regulation_list{idx_TF(j), 3},trend);
                idx_D = idx_D + 1;
            end
        end
    end
end

statistic_result(1,:) = {'Sensitivity', 'Specificity', 'F1 score for upregulation', 'F1 score for downregulation','AUC', 'Precision-positive', 'Precision-negtive'};

sensitivity_D = sum([resultsD{2:end,5}])/length([resultsD{2:end,5}]); %fraction of catched downregulation
sensitivity_U = sum([resultsU{2:end,5}])/length([resultsU{2:end,5}]); %fraction of catched upregulation
precision_D = sum([resultsD{2:end,5}])/(sum([resultsD{2:end,5}])+sum([resultsU{2:end,5}] == 0)); %Precision of downregulation
precision_U = sum([resultsU{2:end,5}])/(sum([resultsU{2:end,5}])+sum([resultsD{2:end,5}] == 0)); %Precision of downregulation
F1_score_D = 2*(precision_D*sensitivity_D/(sensitivity_D+precision_D));
F1_score_U = 2*(precision_U*sensitivity_U/(sensitivity_U+precision_U));
auc = AUC(sensitivity_D, sensitivity_U);

statistic_result{2, 1} = sensitivity_U;
statistic_result{2, 2} = sensitivity_D;
statistic_result{2, 3} = F1_score_U;
statistic_result{2, 4} = F1_score_D;
statistic_result{2, 5} = auc;
statistic_result{2, 6} = precision_U;
statistic_result{2, 7} = precision_D;

save('insilico_size100_1_snr15_GRNVBEM.mat', 'resultsU', 'resultsD', 'statistic_result','-append');

end

%% Calculate AUC value
function auc = AUC(sensitivity, specificity)
if isnan(sensitivity) || isnan(specificity)
    auc = NaN;
else
    y = [0;sensitivity;1];
    x = [0;1-specificity;1];
    auc = trapz(x,y);
end
end