function best_AUC (result)
result = 'result.mat';
load (result)

data = {'THP1_single_cell', 'insilico_size100_1', 'insilico_size100_2'};
data_title = {'THP1', 'Sim 1', 'Sim 2'};
corr = {'Pearson', 'Spearman', 'SDCP', 'SDCS','TSNI', 'GRNVBEM', 'SINCERITIES'};
corr_label = {'Pearson', 'Spearman', 'SDC-P', 'SDC-S','TSNI', 'GRNVBEM', 'SINCERITIES'};
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
        elseif strcmp(corr{j}, 'TSNI')
            if strcmp(norm{idx_max(1)}, 'log')
                plot_data{idx_data, 1} = [upper(norm{idx_max(1)})];
                plot_data{idx_data, 2} = max_auc;
            else
                plot_data{idx_data, 1} = [norm{idx_max(1)}];
                plot_data{idx_data, 2} = max_auc;
            end
        elseif strcmp(corr{j}, 'GRNVBEM')
            if strcmp(norm{idx_max(1)}, 'log')
                plot_data{idx_data, 1} = [upper(norm{idx_max(1)})];
                plot_data{idx_data, 2} = max_auc;
            elseif strcmp(norm{idx_max(1)}, 'zscore')
                plot_data{idx_data, 1} = [upper(norm{idx_max(1)})];
                plot_data{idx_data, 2} = max_auc;
            else
                plot_data{idx_data, 1} = [norm{idx_max(1)}];
                plot_data{idx_data, 2} = max_auc;
            end
        elseif strcmp(corr{j}, 'SINCERITIES')
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
    
    subplot(3,1,i);
    
    x = 1:7;
    bar(x,y, 0.3);
    bar(x, [y(1) nan nan nan nan nan nan], 0.3, 'FaceColor', [0 0.4470 0.7410]);
    hold on
    bar(x, [nan y(2) nan nan nan nan nan], 0.3, 'FaceColor', [0.8500 0.3250 0.0980]);
    hold on
    bar(x, [nan nan y(3) nan nan nan nan], 0.3, 'FaceColor', [0.9290 0.6940 0.1250]);
    hold on
    bar(x, [nan nan nan y(4) nan nan nan], 0.3, 'FaceColor', [0.4940 0.1840 0.5560]);
    hold on
    bar(x, [nan nan nan nan y(5) nan nan], 0.3, 'FaceColor', [0.4660 0.6740 0.1880]);
    hold on
    bar(x, [nan nan nan nan nan y(6) nan], 0.3, 'FaceColor', [0.6350 0.0780 0.1840]);
    hold on
    bar(x, [nan nan nan nan nan nan y(7)], 0.3, 'FaceColor', 'y');
    hold on
    
    set(gca, 'XTickLabel', leg, 'Fontsize', x_fontsize);
    for m = 1:length(y)
        text(m,y(m)+0.01,num2str(round(y(m),2)),...
            'VerticalAlignment','bottom','horizontalalign','center','Fontsize',15)
    end
    ylabel('AUC', 'FontSize', 17)
    ylim([0 1])
    title(data_title{i},'Fontsize', 18)
    if i == 3
        legend(corr_label,'Location','southoutside', 'Orientation', 'horizontal')
    end
    
    
end
 t = suptitle('The best performance of methods across different datasets');
 set(t,'FontSize',30);

end
