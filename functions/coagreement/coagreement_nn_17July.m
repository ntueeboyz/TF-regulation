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



