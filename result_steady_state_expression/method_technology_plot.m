function method_technology_plot(result)

result = 'result.mat';
load (result)

data = {'microarray', 'RNAseq', 'CAGE'};
data_label = {'Microarray', 'RNA-seq', 'CAGE'};
norm = {'without', 'LOG', 'ZSCORE'};
% corr = {'Pearson', 'Spearman', 'SDCS','SDCP'};
corr = {'SDCP','SDCS'};
corr_label = {'SDC-P','SDC-S'};
statistics_properties = {'Sensitivity', 'Specificity', 'F1-Score positive', 'F1-Score negative', 'AUC', 'Precision-positive', 'Precision-negative'};
x_fontsize = 15;

y_microarray = [];
e_microarray = [];
y_rnaseq = [];
e_rnaseq = [];
y_cage = [];
e_cage = [];

figure

%% Extract data from technologies
for i = 1:length(corr)
    plot_data = {};
    plot_data(1,:) = {'norm', 'Sensitivity', 'Specificity', 'F1-Score for upregulation', 'F1-Score for downregulation', 'AUC', 'Precision-positive', 'Precision-negative'};
    idx_data = 2;
    for j = 1:length(data)
        exp_data = regexp(statistic_result(:,1) ,strcat('(?<=',data{j}, '_', ').*(?=',corr{i},')'), 'match');
        
        sens_v = [];
        spe_v =[];
        f1_p = [];
        f1_n = [];
        AUC_v = [];
        pre_p = [];
        pre_n = [];
        for k = 1:length(exp_data)
            if ~isempty(exp_data{k})
                sens_v = [sens_v, statistic_result{k,2}];
                spe_v =[spe_v, statistic_result{k,3}];
                f1_p = [f1_p, statistic_result{k,4}];
                f1_n = [f1_n, statistic_result{k,5}];
                AUC_v = [AUC_v, statistic_result{k,6}];
                pre_p = [pre_p, statistic_result{k,7}];
                pre_n = [pre_n, statistic_result{k,8}];
            end
        end
        sens_mean = mean(sens_v);
        sens_error = std(sens_v)/sqrt(length(sens_v));
        
        spe_mean = mean(spe_v);
        spe_error = std(spe_v)/sqrt(length(spe_v));
        
        f1_p_mean = mean(f1_p);
        f1_p_error = std(f1_p)/sqrt(length(f1_p));
        
        f1_n_mean = mean(f1_n);
        f1_n_error = std(f1_n)/sqrt(length(f1_n));
        
        AUC_mean = mean(AUC_v);
        AUC_error = std(AUC_v)/sqrt(length(AUC_v));
        
        pre_p_mean = mean(pre_p);
        pre_p_error = std(pre_p)/sqrt(length(pre_p));
        
        pre_n_mean = mean(pre_n);
        pre_n_error = std(pre_n)/sqrt(length(pre_n));
        
        plot_data{idx_data, 1} = data{j};
        plot_data{idx_data, 2} = sens_mean;
        plot_data{idx_data, 3} = spe_mean;
        plot_data{idx_data, 4} = f1_p_mean;
        plot_data{idx_data, 5} = f1_n_mean;
        plot_data{idx_data, 6} = AUC_mean;
        plot_data{idx_data, 7} = pre_p_mean;
        plot_data{idx_data, 8} = pre_n_mean;
        
        idx_data = idx_data + 1;
        
        plot_data{idx_data, 1} = [data{j}, '_std_error'];
        plot_data{idx_data, 2} = sens_error;
        plot_data{idx_data, 3} = spe_error;
        plot_data{idx_data, 4} = f1_p_error;
        plot_data{idx_data, 5} = f1_n_error;
        plot_data{idx_data, 6} = AUC_error;
        plot_data{idx_data, 7} = pre_p_error;
        plot_data{idx_data, 8} = pre_n_error;
        
        idx_data = idx_data + 1;
        
    end
    
    y = [];
    e = [];
    for l = 2:size(plot_data,2)
        stat_v = [plot_data{2, l} plot_data{4, l} plot_data{6, l}];
        error_v = [plot_data{3, l} plot_data{5, l} plot_data{7, l}];
        y = [y; stat_v];
        e = [e; error_v];
    end
    
    max_v = [];
    max_idx = [];
    for q = 1:size(y,1)
        m = max(y(q,:));
        max_v = [max_v, m];
        m_idx = find(y(q,:) == m);
        max_idx = [max_idx, m_idx];
    end
    subplot(2,1,i);
    
    x = 1:7;
    h=bar(x,y);
    set(gca, 'XTickLabel', statistics_properties, 'Fontsize', x_fontsize);
    ylim([0 1])
    title(corr_label{i},'Fontsize', 18)
    if i == 2
        legend(data_label,'Location', 'southoutside','Orientation','horizontal');
    end
    for t=1:7
        
        text(t+h(max_idx(t)).XOffset,max_v(t)+0.02,num2str(round(max_v(t),2)), ...
            'VerticalAlignment','bottom','horizontalalign','center','Fontsize',15)
    end
    hold on
    
    ngroups = size(y, 1);
    nbars = size(y, 2);
    groupwidth = min(0.8, nbars/(nbars + 1.5));
    for m = 1:nbars
        x = (1:ngroups) - groupwidth/2 + (2*m-1) * groupwidth / (2*nbars);
        errorbar(x, y(:,m), e(:,m), '.', 'Color', [0 0 0]);
    end
    hold off
    
end
t = suptitle('Performance of technologies across different correlations');
set(t,'FontSize',30);

end