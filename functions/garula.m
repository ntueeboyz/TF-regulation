function z = garula(x,y,op)
    % op = 1; median centering
    % op = 2; mean centering
    % op = 3; distribution mode centering
    switch op
        case 1
            z = 1/(length(x))*sum(abs(sign(x-median(x))+sign(y-median(y))))-1;
        case 2
            z = 1/(length(x))*sum(abs(sign(x-mean(x))+sign(y-mean(y))))-1;
        case 3
            z = 1/(length(x))*sum(abs(sign(x-dmode(x))+sign(y-dmode(y))))-1;
    end
end

% mode based on the random variable distribution
function mode = dmode(x)
    mode = [];
    for i=1:size(x,2)
        [f,xi] = ksdensity(x(:,i));
        [~, ind] = max(f);
        mode = horzcat(mode, xi(ind));
    end
end


% Formulas that did not work
% z = 1/(length(x)-1)*sum(abs(sign(zscore(x))+sign(zscore(y))))-1;
% z = 2/(sum(sign(abs(zscore(x)+zscore(y)))))*sum(abs(zscore(x)+zscore(y)))-1;