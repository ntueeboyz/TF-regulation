function [c,p,r] = coagreement_nn(x,y,type,it)
 
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
           cal=abs(xs+ys); 
           % first learning step output
           o = 1/n*sum(cal)-1; 
           
           % weights update
           if strcmp(type,'linear')
               u=1-abs(o);
           else
               u=abs(o);
               u=1-(1./(1+exp(-(14*u-7)))+u)/2;
           end
           weight=ones(1,n); weight(cal==0)=u;
           n=n-sum(cal==0)*(1-u);
           
           %second learning step output              
           z1 = 1/n*sum(cal.*weight)-1;
           
%%% centering rule: median/median
           n=length(x);
           % binaly classifier layer
           xs=sign(x-median(x)); ys=sign(y-median(y));
           % coagreement layer
           cal=abs(xs+ys); 
           % first learning step output
           o = 1/n*sum(cal)-1; 
           
           % weights update
           if strcmp(type,'linear')
           u=1-abs(o);
           else
           u=abs(o);
           u=1-(1./(1+exp(-(14*u-7)))+u)/2;
           end
           weight=ones(1,n); weight(cal==0)=u;
           n=n-sum(cal==0)*(1-u); 
           
           %second learning step output              
           z2 = 1/n*sum(cal.*weight)-1;
           
%%% centering rule: mean/median
           n=length(x);
           % binaly classifier layer
           xs=sign(x-mean(x)); ys=sign(y-median(y));
           % coagreement layer
           cal=abs(xs+ys); 
           % first learning step output
           o = 1/n*sum(cal)-1; 
           
           % weights update
           if strcmp(type,'linear')
           u=1-abs(o);
           else
           u=abs(o);
           u=1-(1./(1+exp(-(14*u-7)))+u)/2;
           end
           weight=ones(1,n); weight(cal==0)=u;
           n=n-sum(cal==0)*(1-u); 
           
           %second learning step output              
           z3 = 1/n*sum(cal.*weight)-1;
           
%%% centering rule: median/mean
           n=length(x);
           % binaly classifier layer
           xs=sign(x-median(x)); ys=sign(y-mean(y));
           % coagreement layer
           cal=abs(xs+ys); 
           % first learning step output
           o = 1/n*sum(cal)-1; 
           
           % weights update
           if strcmp(type,'linear')
           u=1-abs(o);
           else
           u=abs(o);
           u=1-(1./(1+exp(-(14*u-7)))+u)/2;
           end
           weight=ones(1,n); weight(cal==0)=u;
           n=n-sum(cal==0)*(1-u); 
           
           %second learning step output              
           z4 = 1/n*sum(cal.*weight)-1;
          
z=max([z1 z2 z3 z4]);
r=find([z1 z2 z3 z4]==z,1);
n=length(x);

end



