function noise_plot2(result)

result = 'result.mat';
load (result)

data = {'THP1_single_cell', 'insilico_size100_1', 'insilico_size100_2'};
data_label = {'THP1', 'Simulated data 1', 'Simulated data 2'};
noise = {'','snr15','snr10', 'snr5'};
noise_label = {'Without noise', 'SNR 15', 'SNR 10', 'SNR 5'};
% Best auc in THP1
corr_t = {'LOG_Pearson', 'without_Spearman', 'without_SDCP', 'without_SDCS', 'without_TSNI', 'zscore_GRNVBEM', 'without_SINCERITIES'};
% Best auc in Sim1
corr_s1 = {'LOG_Pearson', 'without_Spearman', 'without_SDCP', 'without_SDCS','log_TSNI', 'log_GRNVBEM', 'without_SINCERITIES'};
% Best auc in Sim2
corr_s2 = {'LOG_Pearson', 'without_Spearman', 'LOG_SDCP', 'LOG_SDCS','without_TSNI', 'without_GRNVBEM', 'without_SINCERITIES'};

corr = {'Pearson', 'Spearman', 'SDC-P', 'SDC-P', 'TSNI', 'GRNVBEM', 'SINCERITIES'};
x_fontsize = 18;
subplot_idx = 1;
for i = 1:length(data)
    plot_data = {};
    plot_data(1,:) = {'norm', 'AUC', 'Pre-P', 'Pre-N'};
    idx_data = 2;
    for j = 1:length(noise)
        if i == 1
            for o = 1:length(corr_t)
                if j == 1
                    exp_data = find(strcmp(statistic_result(:,1), [data{i},'_',corr_t{o}]));
                else
                    exp_data = find(strcmp(statistic_result(:,1), [data{i},'_',noise{j},'_',corr_t{o}]));
                end                
                
                AUC_mean = statistic_result{exp_data,6};
                
                pre_p_mean = statistic_result{exp_data,7};
                
                pre_n_mean = statistic_result{exp_data,8};
                
                plot_data{idx_data, 1} = [data{i},'_',noise{j}, '_',corr_t{o}];
                plot_data{idx_data, 2} = AUC_mean;
                plot_data{idx_data, 3} = pre_p_mean;
                plot_data{idx_data, 4} = pre_n_mean;
                
                idx_data = idx_data + 1;
                
            end
        elseif i == 2
            for o = 1:length(corr_s1)
                if j == 1
                    exp_data = find(strcmp(statistic_result(:,1), [data{i},'_',corr_s1{o}]));
                else
                    exp_data = find(strcmp(statistic_result(:,1), [data{i},'_',noise{j},'_',corr_s1{o}]));
                end
                
                AUC_mean = statistic_result{exp_data,6};
                
                pre_p_mean = statistic_result{exp_data,7};
                
                pre_n_mean = statistic_result{exp_data,8};
                
                plot_data{idx_data, 1} = [data{i},'_',noise{j}, '_',corr_s1{o}];
                plot_data{idx_data, 2} = AUC_mean;
                plot_data{idx_data, 3} = pre_p_mean;
                plot_data{idx_data, 4} = pre_n_mean;
                
                idx_data = idx_data + 1;
                
            end
        elseif i == 3
            for o = 1:length(corr_s2)
                if j == 1
                    exp_data = find(strcmp(statistic_result(:,1), [data{i},'_',corr_s2{o}]));
                else
                    exp_data = find(strcmp(statistic_result(:,1), [data{i},'_',noise{j},'_',corr_s2{o}]));
                end
                
                AUC_mean = statistic_result{exp_data,6};
                
                pre_p_mean = statistic_result{exp_data,7};
                
                pre_n_mean = statistic_result{exp_data,8};
                
                plot_data{idx_data, 1} = [data{i},'_',noise{j}, '_',corr_s2{o}];
                plot_data{idx_data, 2} = AUC_mean;
                plot_data{idx_data, 3} = pre_p_mean;
                plot_data{idx_data, 4} = pre_n_mean;
                
                idx_data = idx_data + 1;
                
            end
        end     
    end
    
    y_a_p = [];
    y_pp_p = [];
    y_pn_p = [];
    
    y_a_s = [];
    y_pp_s = [];
    y_pn_s = [];
    
    y_a_sdcp = [];
    y_pp_sdcp = [];
    y_pn_sdcp = [];
    
    y_a_sdcs = [];
    y_pp_sdcs = [];
    y_pn_sdcs = [];
    
    y_a_t = [];
    y_pp_t = [];
    y_pn_t = [];
    
    y_a_g = [];
    y_pp_g = [];
    y_pn_g = [];
    
    y_a_si = [];
    y_pp_si = [];
    y_pn_si = [];
    
    idx_pd = 2;
    for r = 1:4 % start from without, snr15, snr10, snr5.
        
        y_a_p = [y_a_p plot_data{idx_pd, 2}];
        y_pp_p = [y_pp_p plot_data{idx_pd, 3}];
        y_pn_p = [y_pn_p plot_data{idx_pd, 4}];
        
        y_a_s = [y_a_s plot_data{idx_pd+1, 2}];
        y_pp_s = [y_pp_s plot_data{idx_pd+1, 3}];
        y_pn_s = [y_pn_s plot_data{idx_pd+1, 4}];
        
        y_a_sdcp = [y_a_sdcp plot_data{idx_pd+2, 2}];
        y_pp_sdcp = [y_pp_sdcp plot_data{idx_pd+2, 3}];
        y_pn_sdcp = [y_pn_sdcp plot_data{idx_pd+2, 4}];
        
        y_a_sdcs = [y_a_sdcs plot_data{idx_pd+3, 2}];
        y_pp_sdcs = [y_pp_sdcs plot_data{idx_pd+3, 3}];
        y_pn_sdcs = [y_pn_sdcs plot_data{idx_pd+3, 4}];
        
        y_a_t = [y_a_t plot_data{idx_pd+4, 2}];
        y_pp_t = [y_pp_t plot_data{idx_pd+4, 3}];
        y_pn_t = [y_pn_t plot_data{idx_pd+4, 4}];
        
        y_a_g = [y_a_g plot_data{idx_pd+5, 2}];
        y_pp_g = [y_pp_g plot_data{idx_pd+5, 3}];
        y_pn_g = [y_pn_g plot_data{idx_pd+5, 4}];
        
        y_a_si = [y_a_si plot_data{idx_pd+6, 2}];
        y_pp_si = [y_pp_si plot_data{idx_pd+6, 3}];
        y_pn_si = [y_pn_si plot_data{idx_pd+6, 4}];
        
        idx_pd = idx_pd + 7;
    end
    
    subplot(3,3,subplot_idx) % Precision-positive row in subplot
    x = 1:4;
    h = plot(x,y_a_p,x,y_a_s,x,y_a_sdcp,x,y_a_sdcs,x,y_a_t,x,y_a_g,x,y_a_si,'LineWidth',3);
    set(h,{'Color'},{[0 0.4470 0.7410];[0.8500 0.3250 0.0980];[0.9290 0.6940 0.1250];[0.4940 0.1840 0.5560];[0.4660 0.6740 0.1880];[0.6350 0.0780 0.1840];[1 1 0]},{'LineStyle'},{':';'--';'-.';'-.';'-';'-';'-'})
    xticks([1 2 3 4])
    xticklabels(noise_label);
    a = get(gca,'xTickLabel');
    set(gca,'xTickLabel',a,'Fontsize',15);
    ylabel('AUC','Fontsize', x_fontsize)
    if subplot_idx == 1 || subplot_idx == 2 || subplot_idx == 3
        ylim([0 1])
    elseif subplot_idx == 7 || subplot_idx == 8 || subplot_idx == 9
        ylim([0.4 0.8])
    else
        ylim([0.4 1])
    end
    title(data_label{i},'Fontsize', 18)
    
    subplot_idx = subplot_idx + 1;
    
    subplot(3,3,subplot_idx); % Precision-positive row in subplot
    h1 = plot(x,y_pp_p,x,y_pp_s,x,y_pp_sdcp,x,y_pp_sdcs,x,y_pp_t,x,y_pp_g,x,y_pp_si,'LineWidth',3);
    set(h1,{'Color'},{[0 0.4470 0.7410];[0.8500 0.3250 0.0980];[0.9290 0.6940 0.1250];[0.4940 0.1840 0.5560];[0.4660 0.6740 0.1880];[0.6350 0.0780 0.1840];[1 1 0]},{'LineStyle'},{':';'--';'-.';'-.';'-';'-';'-'})
    xticks([1 2 3 4])
    xticklabels(noise_label)
    a = get(gca,'xTickLabel');
    set(gca,'xTickLabel',a,'Fontsize',15);
    ylabel('Precision-positive','Fontsize', x_fontsize)
    if subplot_idx == 1 || subplot_idx == 2 || subplot_idx == 3
        ylim([0 1])
    elseif subplot_idx == 7 || subplot_idx == 8 || subplot_idx == 9
        ylim([0.4 0.8])
    else
        ylim([0.4 1])
    end
    title(data_label{i},'Fontsize', 18)
    if i == 3
        legend(corr,'Location','southoutside','Orientation','horizontal','Fontsize', 16)
    end
    
    subplot_idx = subplot_idx + 1;
    
    
    subplot(3,3,subplot_idx); % Precision-negative row in subplot
    h2 = plot(x,y_pn_p,x,y_pn_s,x,y_pn_sdcp,x,y_pn_sdcs,x,y_pn_t,x,y_pn_g,x,y_pn_si,'LineWidth',3);
    set(h2,{'Color'},{[0 0.4470 0.7410];[0.8500 0.3250 0.0980];[0.9290 0.6940 0.1250];[0.4940 0.1840 0.5560];[0.4660 0.6740 0.1880];[0.6350 0.0780 0.1840];[1 1 0]},{'LineStyle'},{':';'--';'-.';'-.';'-';'-';'-'})
    xticks([1 2 3 4])
    xticklabels(noise_label)
    a = get(gca,'xTickLabel');
    set(gca,'xTickLabel',a,'Fontsize',15);
    ylabel('Precision-negative','Fontsize', x_fontsize)
    if subplot_idx == 1 || subplot_idx == 2 || subplot_idx == 3
        ylim([0 1])
    elseif subplot_idx == 7 || subplot_idx == 8 || subplot_idx == 9
        ylim([0.4 0.8])
    else
        ylim([0.4 1])
    end
    title(data_label{i},'Fontsize', 18)
    
    
    subplot_idx = subplot_idx + 1;
end
t = suptitle('Noise data');
set(t,'FontSize',30);
end