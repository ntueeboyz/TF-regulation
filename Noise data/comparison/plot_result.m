%% THP1

load THP1.mat

t_p_auc = mean([statistic_result{2,6} statistic_result{3,6} statistic_result{4,6}]);

t_s_auc = mean([statistic_result{5,6} statistic_result{6,6} statistic_result{7,6}]);

t_d_auc = mean([statistic_result{8,6} statistic_result{9,6} statistic_result{10,6}]);

t_p_e = std([statistic_result{2,6} statistic_result{3,6} statistic_result{4,6}])/sqrt(3);

t_s_e = std([statistic_result{5,6} statistic_result{6,6} statistic_result{7,6}])/sqrt(3);

t_d_e = std([statistic_result{8,6} statistic_result{9,6} statistic_result{10,6}])/sqrt(3);


load THP1_TSNI.mat

t_tsni_auc = statistic_result{2,5};

load THP1_SINCERITIES.mat

t_sincerities_auc = statistic_result{2,5};

load THP1_GRNVBEM.mat

t_grnvbem_auc = statistic_result{2,5};

%% simulated data 1

load insilico_size100_1.mat

s1_p_auc = mean([statistic_result{2,6} statistic_result{3,6} statistic_result{4,6}]);

s1_s_auc = mean([statistic_result{5,6} statistic_result{6,6} statistic_result{7,6}]);

s1_d_auc = mean([statistic_result{8,6} statistic_result{9,6} statistic_result{10,6}]);

s1_p_e = std([statistic_result{2,6} statistic_result{3,6} statistic_result{4,6}])/sqrt(3);

s1_s_e = std([statistic_result{5,6} statistic_result{6,6} statistic_result{7,6}])/sqrt(3);

s1_d_e = std([statistic_result{8,6} statistic_result{9,6} statistic_result{10,6}])/sqrt(3);

load insilico_size100_1_TSNI.mat

s1_tsni_auc = statistic_result{2,5};

load insilico_size100_1_SINCERITIES.mat

s1_sincerities_auc = statistic_result{2,5};

load insilico_size100_1_GRNVBEM.mat

s1_grnvbem_auc = statistic_result{2,5};

%%

load insilico_size100_2.mat

s2_p_auc = mean([statistic_result{2,6} statistic_result{3,6} statistic_result{4,6}]);

s2_s_auc = mean([statistic_result{5,6} statistic_result{6,6} statistic_result{7,6}]);

s2_d_auc = mean([statistic_result{8,6} statistic_result{9,6} statistic_result{10,6}]);

s2_p_e = std([statistic_result{2,6} statistic_result{3,6} statistic_result{4,6}])/sqrt(3);

s2_s_e = std([statistic_result{5,6} statistic_result{6,6} statistic_result{7,6}])/sqrt(3);

s2_d_e = std([statistic_result{8,6} statistic_result{9,6} statistic_result{10,6}])/sqrt(3);

load insilico_size100_2_TSNI.mat

s2_tsni_auc = statistic_result{2,5};

load insilico_size100_2_SINCERITIES.mat

s2_sincerities_auc = statistic_result{2,5};

load insilico_size100_2_GRNVBEM.mat

s2_grnvbem_auc = statistic_result{2,5};

%% Plot
method = {'Pearson', 'Spearman', 'Distance', 'TSNI', 'SINCERITIES', 'GRNVBEM'};
width = 0.3;
xfontsize = 15;
yfontsize = 15;
tfontsize = 18;

figure

subplot(3,1,1);
x = 1:6;
y_t = [t_p_auc t_s_auc t_d_auc t_tsni_auc t_sincerities_auc t_grnvbem_auc];
e_t = [t_p_e t_s_e t_d_e 0 0 0];
bar(x,y_t,width);
set(gca, 'XTickLabel', method, 'Fontsize', xfontsize);
ylabel('AUC', 'Fontsize', yfontsize)
title('THP1','Fontsize', tfontsize)
ylim([0 1])
yticks([0.2 0.4 0.6 0.8 1])
for m = 1:length(y_t)
    text(m-0.06,y_t(m)+0.1,num2str(round(y_t(m),2)),'Fontsize', 12)
end
hold on
errorbar(x,y_t,e_t,'.', 'Color', [0 0 0])
hold off

subplot(3,1,2);
x = 1:6;
y_s1 = [s1_p_auc s1_s_auc s1_d_auc s1_tsni_auc s1_sincerities_auc s1_grnvbem_auc];
e_s1 = [s1_p_e s1_s_e s1_d_e 0 0 0];
bar(x,y_s1,width);
set(gca, 'XTickLabel', method, 'Fontsize', xfontsize);
ylabel('AUC', 'Fontsize', yfontsize)
title('Simulated data 1','Fontsize', tfontsize)
ylim([0 1])
yticks([0.2 0.4 0.6 0.8 1])
for m = 1:length(y_s1)
    text(m-0.06,y_s1(m)+0.1,num2str(round(y_s1(m),2)),'Fontsize', 12)
end
hold on
errorbar(x,y_s1,e_s1,'.', 'Color', [0 0 0])
hold off

subplot(3,1,3);
x = 1:6;
y_s2 = [s2_p_auc s2_s_auc s2_d_auc s2_tsni_auc s2_sincerities_auc s2_grnvbem_auc];
e_s2 = [s2_p_e s2_s_e s2_d_e 0 0 0];
bar(x,y_s2,width);
set(gca, 'XTickLabel', method, 'Fontsize', xfontsize);
ylabel('AUC', 'Fontsize', yfontsize)
title('Simulated data 2','Fontsize', tfontsize)
ylim([0 1])
yticks([0.2 0.4 0.6 0.8 1])
for m = 1:length(y_s2)
    text(m-0.06,y_s2(m)+0.1,num2str(round(y_s2(m),2)),'Fontsize', 12)
end
hold on
errorbar(x,y_s2,e_s2,'.', 'Color', [0 0 0])
hold off

t = suptitle('Method comparison');
set(t,'FontSize',30);