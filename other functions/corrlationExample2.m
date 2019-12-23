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
[mi, nmi2] = mi_gg(x,y); %MI2 normalized by entropy
nmi1 = mi/(mi +1); %MI1 normalized by I/(I+1)
coag = coagreement(x,y); %Normal coagreement
coag_linear = coagreement_nn_17July(x,y,'linear'); %Final version of coagreement
coag_sigmoid = coagreement_nn_17July(x,y,'sigmoide'); %Final version of coagreement

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


%% Signed distance correlation
function c = signed_dist(x,y)
dist_corr = distcor(x, y);
if isnan(dist_corr)
    c = NaN;
else
    distance=[];
    for i = 1:length(x)-1
        for j = i+1:length(x)
            deltaX = x(i) - x(j);
            deltaY = y(i) - y(j);
            twoPoint = [x(i),y(i); x(j),y(j)];
            dist = pdist(twoPoint, 'euclidean');
            value = deltaX*deltaY*dist;
            distance = [distance, value];
        end
    end
    c = sign(mean(distance))*dist_corr;
end
end

function dcor = distcor(x,y)
% This function calculates the distance correlation between x and y.
% Reference: http://en.wikipedia.org/wiki/Distance_correlation
% Date: 18 Jan, 2013
% Author: Shen Liu (shen.liu@hotmail.com.au)
% Check if the sizes of the inputs match

if size(x,1) ~= size(y,1)
    error('Inputs must have the same number of rows')
end
% Delete rows containing unobserved values
N = any([isnan(x) isnan(y)],2);
x(N,:) = [];
y(N,:) = [];
% Calculate doubly centered distance matrices for x and y
a = pdist2(x,x);
mcol = mean(a);
mrow = mean(a,2);
ajbar = ones(size(mrow))*mcol;
akbar = mrow*ones(size(mcol));
abar = mean(mean(a))*ones(size(a));
A = a - ajbar - akbar + abar;
b = pdist2(y,y);
mcol = mean(b);
mrow = mean(b,2);
bjbar = ones(size(mrow))*mcol;
bkbar = mrow*ones(size(mcol));
bbar = mean(mean(b))*ones(size(b));
B = b - bjbar - bkbar + bbar;
% Calculate squared sample distance covariance and variances
dcov = sum(sum(A.*B))/(size(mrow,1)^2);
dvarx = sum(sum(A.*A))/(size(mrow,1)^2);
dvary = sum(sum(B.*B))/(size(mrow,1)^2);
% Calculate the distance correlation
dcor = sqrt(dcov/sqrt(dvarx*dvary));
end


%% Mutual information
% I is unnormalized
% NI is normalized by joint entropy
function [I , NI] = mi_gg(x, y, biascorrect, demeaned)
% MI_GG Mutual information (MI) between two Gaussian variables in bits
%   I = mi_gg(x,y) returns the MI between two (possibly multidimensional)
%   Gassian variables, x and y, with bias correction.
%   If x and/or y are multivariate rows must correspond to samples, columns
%   to dimensions/variables. (Samples first axis) 
%
%   biascorrect : true / false option (default true) which specifies whether
%   bias correction should be applied to the esimtated MI.
%   demeaned : false / true option (default false) which specifies whether the
%   input data already has zero mean (true if it has been copula-normalized)

% ensure samples first axis for vectors
if isvector(x)
    x = x(:);
end
if isvector(y)
    y = y(:);
end
if ndims(x)~=2 || ndims(y)~=2
    error('mi_gg: input arrays should be 2d')
end
Ntrl = size(x,1);
Nvarx = size(x,2);
Nvary = size(y,2);

if size(y,1) ~= Ntrl
    error('mi_gg: number of trials do not match')
end

% default option values
if nargin<3
    biascorrect = true;
end
if nargin<4
    demeaned = false;
end

% demean data if required
if ~demeaned
    x = bsxfun(@minus,x,sum(x,1)/Ntrl);
    y = bsxfun(@minus,y,sum(y,1)/Ntrl);
end

% joint variable
xy = [x y];
Cxy = (xy'*xy) / (Ntrl - 1);
% submatrices of joint covariance
Cx = Cxy(1:Nvarx,1:Nvarx);
ystart = Nvarx + 1;
Nvarxy = Nvarx + Nvary;
Cy = Cxy(ystart:Nvarxy,ystart:Nvarxy);

chCx = chol(Cx);
chCy = chol(Cy);
chCxy = chol(Cxy);

% entropies in nats
% normalisations cancel for information
HX = sum(log(diag(chCx))); % + 0.5*Nvarx*log(2*pi*exp(1));
HY = sum(log(diag(chCy))); % + 0.5*Nvary*log(2*pi*exp(1));
HXY = sum(log(diag(chCxy))); % + 0.5*(Nvarx+Nvary)*log(2*pi*exp(1));

ln2 = log(2);
if biascorrect
    psiterms = psi((Ntrl - (1:Nvarxy))/2) / 2;
    dterm = (ln2 - log(Ntrl-1)) / 2;
    HX = (HX - Nvarx*dterm - sum(psiterms(1:Nvarx)));
    HY = (HY - Nvary*dterm - sum(psiterms(1:Nvary)));
    HXY = (HXY - Nvarxy*dterm - sum(psiterms));
end

% convert to bits
I = (HX + HY - HXY) / ln2;

% Normalization
NI = -I/HXY;

end


%% Normal coagreement
function [c,p,r] = coagreement(x,y,it) 
    % x: first variable
    % y: second variable
    % it: number of iteration for null model (suggestion:1000)
    % r: type of rule, 1 is mean/mean, 2 is median/median, 3 is mean/median and 4 is median/mean

    % co-agreement
    [c,n,r]=cannistraci(x,y);


    %null model
    if nargin == 3
        m = zeros(1,it+1);
        m(1) = c;
        parfor i=2:it+1
            x_perm = x(randperm(n));
            y_perm = y(randperm(n));
            m(i) = coagreement(x_perm, y_perm);
        end
        positionC = find(sort(m, 'ascend') == c);

        if c >= 0
            p = 1 - (positionC(1)-1)/(it+1);
        else
            p = positionC(end)/(it+1);
        end
    else
        p = nan;
    end

end

function  [z,n,r] = cannistraci(x,y)

    n=length(x);

    z1 = 1/n*sum(abs(sign(x-mean(x))+sign(y-mean(y))))-1;
    z2 = 1/n*sum(abs(sign(x-median(x))+sign(y-median(y))))-1;
    z3 = 1/n*sum(abs(sign(x-mean(x))+sign(y-median(y))))-1;
    z4 = 1/n*sum(abs(sign(x-median(x))+sign(y-mean(y))))-1;
    z=max([z1 z2 z3 z4]);
    r=find([z1 z2 z3 z4]==z,1);

end


%% Coagreement with linear and sigmoid type
 function [c,p,r] = coagreement_nn_17July(x,y,type,it)
 
% x: first variable
% y: second variable
% type: 'linear' or 'sigmoide'(default is sigmoide)
% it: number of iteration for null model (suggestion:1000)
% r: type of rule, 1 is mean/mean, 2 is median/median, 3 is mean/median and 4 is median/mean

% co-agreement
[c,n,r]=coagrenet(x,y,type);


%null model
if nargin == 4
                m = zeros(1,it+1);
                m(1) = c;
                parfor i=2:it+1
                    x_perm = x(randperm(n));
                    y_perm = y(randperm(n));
                    m(i) = coagreement(x_perm, y_perm);
                end
                positionC = find(sort(m, 'ascend') == c);
                
                if c >= 0
                    p = 1 - (positionC(1)-1)/(it+1); 
                else
                    p = positionC(end)/(it+1); 
                end
            else
                p = nan;
end

end


function  [z,n,r] = coagrenet(x,y,type)

           
%%% centering rule: mean/mean
           n=length(x);
           % binaly classifier layer
           xs=sign(x-mean(x)); ys=sign(y-mean(y));
           % coagreement layer
           cal=abs(xs+ys)-1;
           % first step output
           o = 1/n*sum(cal);
           disp(o);
           % learning weights update 
           if o==0
           z1=o; %final output case o==0
           else
           
           if o>0
           xs(xs==0)=1; ys(ys==0)=1;
           elseif o<0
           xs(xs==0)=-1; ys(ys==0)=1;
           end
           cal=abs(xs+ys)-1;   
           weight=ones(1,n);
           if strcmp(type,'linear')
           weight(cal==sign(-o))=1-abs(o);
           else
           %weight(cal==sign(-o))=1-(1./(1+exp(-(14*abs(o)-7)))+abs(o))/2; 
           weight(cal==sign(-o))=1-1./(1+exp(-(14*abs(o)-7)));
           end
           z1 = 1/n*sum(cal.*weight);  %final output case o<>0
           end
           
%%% centering rule: median/median
           n=length(x);
           % binaly classifier layer
           xs=sign(x-median(x)); ys=sign(y-median(y));
           % coagreement layer
           cal=abs(xs+ys)-1;  
           % first step output
           o = 1/n*sum(cal); 
           % learning weights update 
           if o==0
           z2=o; %final output case o==0
           else
           
           if o>0
           xs(xs==0)=1; ys(ys==0)=1;
           elseif o<0
           xs(xs==0)=-1; ys(ys==0)=1;
           end
           cal=abs(xs+ys)-1;   
           
           weight=ones(1,n);
           if strcmp(type,'linear')
           weight(cal==sign(-o))=1-abs(o);   
           else
           %weight(cal==sign(-o))=1-(1./(1+exp(-(14*abs(o)-7)))+abs(o))/2; 
           weight(cal==sign(-o))=1-1./(1+exp(-(14*abs(o)-7)));
           end
           z2 = 1/n*sum(cal.*weight);  %final output case o<>0
           end
           
%%% centering rule: mean/median
           n=length(x);
           % binaly classifier layer
           xs=sign(x-mean(x)); ys=sign(y-median(y));
           % coagreement layer
           cal=abs(xs+ys)-1;  
           % first step output
           o = 1/n*sum(cal); 
           % learning weights update 
           if o==0
           z3=o; %final output case o==0
           else
           
           if o>0
           xs(xs==0)=1; ys(ys==0)=1;
           elseif o<0
           xs(xs==0)=-1; ys(ys==0)=1;
           end
           cal=abs(xs+ys)-1;   
           
           weight=ones(1,n);
           if strcmp(type,'linear')
           weight(cal==sign(-o))=1-abs(o);   
           else
           %weight(cal==sign(-o))=1-(1./(1+exp(-(14*abs(o)-7)))+abs(o))/2; 
           weight(cal==sign(-o))=1-1./(1+exp(-(14*abs(o)-7)));
           end
           z3 = 1/n*sum(cal.*weight);  %final output case o<>0
           end
           
%%% centering rule: median/mean
           n=length(x);
           % binaly classifier layer
           xs=sign(x-median(x)); ys=sign(y-mean(y));
           % coagreement layer
           cal=abs(xs+ys)-1;  
           % first step output
           o = 1/n*sum(cal); 
           % learning weights update 
           if o==0
           z4=o; %final output case o==0
           else
           
           if o>0
           xs(xs==0)=1; ys(ys==0)=1;
           elseif o<0
           xs(xs==0)=-1; ys(ys==0)=1;
           end
           cal=abs(xs+ys)-1;   
           
           weight=ones(1,n);
           if strcmp(type,'linear')
           weight(cal==sign(-o))=1-abs(o);   
           else
           %weight(cal==sign(-o))=1-(1./(1+exp(-(14*abs(o)-7)))+abs(o))/2; 
           weight(cal==sign(-o))=1-1./(1+exp(-(14*abs(o)-7)));
           end
           z4 = 1/n*sum(cal.*weight);  %final output case o<>0
           end
          
z=max([z1 z2 z3 z4]);
r=find([z1 z2 z3 z4]==z,1);
n=length(x);

end