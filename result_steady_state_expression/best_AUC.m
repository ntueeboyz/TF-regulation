function best_AUC (result)
result = 'result.mat';
load (result)

data = {'microarray', 'RNAseq', 'CAGE'};
data_title = {'Microarray', 'RNA-seq', 'CAGE'};
corr = {'Pearson', 'Spearman', 'SDCP','SDCS'};
corr_label = {'Pearson', 'Spearman', 'SDC-P', 'SDC-S'}
x_fontsize = 15;
for i = 1:length(data)
    plot_data = {};
    idx_data = 1;
    for j = 1:length(corr)
        exp_data = regexp(statistic_result(:,1) ,['(?<=',data{i}, '_', ').*(?=','_',corr{j},')'], 'match');
        
        AUC_v = [];
        norm = [];
        idx_norm = 1;
        for k = 1:length(exp_data)
            if ~isempty(exp_data{k})
                AUC_v = [AUC_v, statistic_result{k,6}];
                if strcmp(exp_data{k},'without')
                    norm{idx_norm} = 'Unnormalized';
                    idx_norm = idx_norm + 1;
                else
                    norm{idx_norm} = exp_data{k}{1};
                    idx_norm = idx_norm + 1;
                end
            end
        end
        
        max_auc = max(AUC_v);
        idx_max = find(max_auc == AUC_v);
        if strcmp(corr{j}, 'Pearson')
            plot_data{idx_data, 1} = [norm{idx_max(1)}];
            plot_data{idx_data, 2} = max_auc;
        elseif strcmp(corr{j}, 'Spearman')
            plot_data{idx_data, 1} = [norm{idx_max(1)}];
            plot_data{idx_data, 2} = max_auc;
        elseif strcmp(corr{j}, 'SDCP')
            plot_data{idx_data, 1} = [norm{idx_max(1)}];
            plot_data{idx_data, 2} = max_auc;
        elseif strcmp(corr{j}, 'SDCS')
            plot_data{idx_data, 1} = [norm{idx_max(1)}];
            plot_data{idx_data, 2} = max_auc;
        end
        
        
        idx_data = idx_data + 1;
        
        
    end
    
    y = [];
    leg = [];
    for l = 1:size(plot_data,1)
        y = [y plot_data{l,2}];
        leg{l} = plot_data{l,1};
    end
    
    subplot(1,3,i);
    
    x = 1:4;
    bar(x,y, 0.3);
    bar(x, [y(1) nan nan nan], 0.3, 'FaceColor', [0 0 0.6]);
    hold on
    bar(x, [nan y(2) nan nan], 0.3, 'FaceColor', [0 0.8 1]);
    hold on
    bar(x, [nan nan y(3) nan], 0.3, 'FaceColor', [0.4 0.6 0]);
    hold on
    bar(x, [nan nan nan y(4)], 0.3, 'FaceColor', 'y');
    hold on
    
    set(gca, 'XTickLabel', leg, 'Fontsize', 15);
    xtickangle(-15)
    for m = 1:length(y)
        text(m,y(m)+0.01,num2str(round(y(m),2)),...
            'VerticalAlignment','bottom','horizontalalign','center','Fontsize',15)
    end
    ylabel('AUC', 'FontSize', 17)
    ylim([0 0.8])
    
    title(data_title{i},'Fontsize', 18)
    if i == 2
       legend(corr_label,'Location','southoutside', 'Orientation', 'horizontal','FontSize', 15)

    end
    
    
end
  t = suptitle('The best performance of correlations across different technologies');
  set(t,'FontSize',30);

end
