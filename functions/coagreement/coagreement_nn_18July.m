 function [c,v,p,r,w,a] = coagreement_nn_18July(x,y,type,it)

% INPUT
% x: first variable
% y: second variable
% type: 'linear' or 'sigmoide'(default is sigmoide)
% it: number of iteration for null model (suggestion:1000)
% r: type of rule, 1 is mean/mean, 2 is median/median, 3 is mean/median and 4 is median/mean

% OUTPUT
% c: coagreement_nn (after learning) coefficient
% v: coagreement coefficient
% p-value of c
% r: winner centering rule for c
% w: learned weightes for synapses of coagreement_nn
% a: winner centering rule values for the two naive binary classifiers of the first layer



if isrow(x), x=x'; end %we use column vectors
if isrow(y), y=y'; end %we use column vectors


% co-agreement
[c,v,n,r,w,a]=coagre_net(x,y,type);


%null model
if nargin == 4
                m = zeros(1,it+1);
                m(1) = c;
                parfor i=2:it+1
                    x_perm = x(randperm(n));
                    y_perm = y(randperm(n));
                    m(i) = coagre_net(x_perm, y_perm,type);
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


function  [z,v,n,r,w,a] = coagre_net(x,y,type)
 
n=length(x);
           
%%% centering rule: mean/mean

           % binaly classifier layer
           xs=sign(x-mean(x)); ys=sign(y-mean(y));
           % coagreement layer
           cal=abs(xs+ys)-1;  
           % first step output
           o = 1/n*sum(cal); 
           % learning weights update 
           weight=ones(size(cal));
           if o==0
           z1=o; %final output case o==0
           v1=o;
           else
           
           if o>0
           xs(xs==0)=1; ys(ys==0)=1;
           elseif o<0
           xs(xs==0)=-1; ys(ys==0)=1;
           end
           cal=abs(xs+ys)-1;   
           v1 = 1/n*sum(cal); 
                 
           if strcmp(type,'linear')
           weight(cal==sign(-o))=1-abs(o);   
           else
           %weight(cal==sign(-o))=1-(1./(1+exp(-(14*abs(o)-7)))+abs(o))/2; 
           weight(cal==sign(-o))=1-1./(1+exp(-(14*abs(o)-7)));
           end
           z1 = 1/n*sum(cal.*weight);  %final output case o<>0
           end
           w=weight;
           a=[mean(x) mean(y)];
           
%%% centering rule: median/median
           
           % binaly classifier layer
           xs=sign(x-median(x)); ys=sign(y-median(y));
           % coagreement layer
           cal=abs(xs+ys)-1;  
           % first step output
           o = 1/n*sum(cal); 
           % learning weights update 
           weight=ones(size(cal));
           if o==0
           z2=o; %final output case o==0
           v2=o;
           else
           
           if o>0
           xs(xs==0)=1; ys(ys==0)=1;
           elseif o<0
           xs(xs==0)=-1; ys(ys==0)=1;
           end
           cal=abs(xs+ys)-1;   
           v2 = 1/n*sum(cal); 
           
           if strcmp(type,'linear')
           weight(cal==sign(-o))=1-abs(o);   
           else
           %weight(cal==sign(-o))=1-(1./(1+exp(-(14*abs(o)-7)))+abs(o))/2; 
           weight(cal==sign(-o))=1-1./(1+exp(-(14*abs(o)-7)));
           end
           z2 = 1/n*sum(cal.*weight);  %final output case o<>0
           end
           w=[w weight];
           a=[a; median(x) median(y)];
           
%%% centering rule: mean/median
           
           % binaly classifier layer
           xs=sign(x-mean(x)); ys=sign(y-median(y));
           % coagreement layer
           cal=abs(xs+ys)-1;  
           % first step output
           o = 1/n*sum(cal); 
           % learning weights update 
           weight=ones(size(cal));
           if o==0
           z3=o; %final output case o==0
           v3=o;
           else
           
           if o>0
           xs(xs==0)=1; ys(ys==0)=1;
           elseif o<0
           xs(xs==0)=-1; ys(ys==0)=1;
           end
           cal=abs(xs+ys)-1;   
           v3 = 1/n*sum(cal); 
           
           if strcmp(type,'linear')
           weight(cal==sign(-o))=1-abs(o);   
           else
           %weight(cal==sign(-o))=1-(1./(1+exp(-(14*abs(o)-7)))+abs(o))/2; 
           weight(cal==sign(-o))=1-1./(1+exp(-(14*abs(o)-7)));
           end
           z3 = 1/n*sum(cal.*weight);  %final output case o<>0
           end
           w=[w weight];
           a=[a; mean(x) median(y)];
           
%%% centering rule: median/mean
           
           % binaly classifier layer
           xs=sign(x-median(x)); ys=sign(y-mean(y));
           % coagreement layer
           cal=abs(xs+ys)-1;  
           % first step output
           o = 1/n*sum(cal); 
           % learning weights update 
           weight=ones(size(cal));
           if o==0
           z4=o; %final output case o==0
           v4=o;
           else
           
           if o>0
           xs(xs==0)=1; ys(ys==0)=1;
           elseif o<0
           xs(xs==0)=-1; ys(ys==0)=1;
           end
           cal=abs(xs+ys)-1;   
           v4 = 1/n*sum(cal); 
           
           if strcmp(type,'linear')
           weight(cal==sign(-o))=1-abs(o);   
           else
           %weight(cal==sign(-o))=1-(1./(1+exp(-(14*abs(o)-7)))+abs(o))/2; 
           weight(cal==sign(-o))=1-1./(1+exp(-(14*abs(o)-7)));
           end
           z4 = 1/n*sum(cal.*weight);  %final output case o<>0
           end
           w=[w weight];
           a=[a; median(x) mean(y)];

v=[v1 v2 v3 v4];  
[~,r]=max(abs(v));
v=v(r);

z=[z1 z2 z3 z4];
[~,r]=max(abs(z));
z=z(r);
w=w(:,r);
a=a(r,:);

end



