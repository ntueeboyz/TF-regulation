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



