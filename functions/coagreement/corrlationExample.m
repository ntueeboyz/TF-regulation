%% Data gerneration
c_1 = copularnd('Gaussian',0.99999,500); %Pear = 1
c_08 = copularnd('Gaussian',0.8,500); %Pear = 0.8
c_04 = copularnd('Gaussian',0.4,500); %Pear = 0.4
c_00 = copularnd('Gaussian',0,500); %Pear = 0
c_m04 = copularnd('Gaussian',-0.4,500); %Pear = -0.4
c_m08 = copularnd('Gaussian',-0.8,500); %Pear = -0.8
c_m1 = copularnd('Gaussian',-0.99999,500); %Pear = -1

%% Plot
subplot(1,7,1)
plot_corrrelations(c_1(:,1),c_1(:,2));

subplot(1,7,2)
plot_corrrelations(c_08(:,1),c_08(:,2));

subplot(1,7,3)
plot_corrrelations(c_04(:,1),c_04(:,2));

subplot(1,7,4)
plot_corrrelations(c_00(:,1),c_00(:,2));

subplot(1,7,5)
plot_corrrelations(c_m04(:,1),c_m04(:,2));

subplot(1,7,6)
plot_corrrelations(c_m08(:,1),c_m08(:,2));

subplot(1,7,7)
plot_corrrelations(c_m1(:,1),c_m1(:,2));


%% Plot function
function plot_corrrelations(x,y)
pear = corr(x,y); %Pearson
spear = corr(x,y,'type','sp'); %Spearman
kndl = corr(x,y,'type','kendall'); %Kendall
dist = signed_dist(x,y); %Distance
[mi, nmi1, nmi2] = mi_gg(x,y); %NMI1 suggested by Carlo, NMI2 suggested by Claudio
coag = coagreement(x,y); %Normal coagreement
coag_linear = coagreement_nn_18July(x,y,'linear'); %Final version of coagreement
coag_sigmoid = coagreement_nn_18July(x,y,'sigmoide'); %Final version of coagreement

hold on
scatter(x,y,'.');
axis off
title_p = ['Pearson: ',num2str(round(pear,2))];
title_sp = ['Spearman: ',num2str(round(spear,2))];
title_k = ['kendall: ',num2str(round(kndl,2))];
title_d = ['Distance: ',num2str(round(dist,2))];
title_nmi1 = ['NMI1: ',num2str(round(nmi1,2))];
title_nmi2 = ['NMI2: ',num2str(round(nmi2,2))];
title_coag = ['Coagreement: ',num2str(round(coag,2))];
title_coag_linear = ['Coagreement linear: ',num2str(round(coag_linear,2))];
title_coag_sigmoid = ['Coagreement sigmoid: ',num2str(round(coag_sigmoid,2))];
title({title_p,title_sp,title_k,title_d,title_nmi1,title_nmi2,title_coag,title_coag_linear,title_coag_sigmoid});
end