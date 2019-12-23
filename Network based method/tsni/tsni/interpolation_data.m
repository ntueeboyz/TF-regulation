function [scss,polynomial_coefs_spline,confidence_lower,confidence_upper,smoothed_values] = interpolation_data(x,y,smoothing_parameter,sliding_window_size)
% This function is called by TSNI.m
% The code is to do make the smoothing spline with the given
% smoothinbg_parameter and sliding window size to calculate the weight on
% each time. Once the weight are calculated, it also calculated the
% confidence interval for the spline and do the interpolation at the given
% points. The smoothed curve is forced to pass from the the origin by
% making a natural cubic spline between 0 and first smoothed point. 

weight(1:length(x)) = 1;
[scss,p] = csaps(x,y,smoothing_parameter);
first_der = fnder(scss); %%% First derivative
second_der = fnder(fnder(scss)); %%% Second derivative
%%% To make the smoothed curve pass from zero. It is done by making the
%%% natural cubic spline between 0 and first smoothed point. 
v_secondpoint_zeroder = fnval(scss,x(2));  %%% Value at second point
v_second_point_firstder = fnval(first_der,x(2)); %%% First derivatve at second point
v_second_point_secondder = fnval(second_der,x(2)); %%% Second derivative at second point
mat = [x(2) x(2)^2 x(2)^3; 1 2*x(2) 3*x(2)^2; 0 2 6*x(2)]; 
points = [v_secondpoint_zeroder; v_second_point_firstder;v_second_point_secondder];
pp = inv(mat'*mat)*mat' * points;
polynomial_coefs_spline = [pp(3) pp(2) pp(1) 0];

a = min(x) - (max(x) - min(x))/(2*length(x));
b = max(x) - (max(x) - min(x))/(2*length(x));
co = (pi^4/length(x))*(b-a)^(-4);
trA = 0;
for i = 3:length(x)
    trA = trA + (1+co*smoothing_parameter*(i-1.5)^4)^(-1);
end
trA = trA + 2;
smoothed_values = fnval(scss,x);
smoothed_values(1) = 0; %%% Take care of this value I am making this pass from 0
for loop = 1:1   %%% Loop to reiterate to compute weight, but we stop at first iteration. 
    variance = 0;
    for i=1:length(x)
        variance = variance + weight(i)*(y(i) - smoothed_values(i))^2;
    end
    variance = variance/(length(x)-trA);
    for i = 1:length(x)
        residual(i) = sqrt(weight(i))*(y(i) - smoothed_values(i))/(sqrt(variance)*sqrt(1-length(x)^(-1)*trA));
    end
    
    for i = 1:length(x)
        ni = min(length(x),i+sliding_window_size);
        mi = max(1,i-sliding_window_size);
        invw(i) = 0;
        for j = mi:ni
            invw(i) = invw(i)+(ni-mi+1)^(-1)*weight(i)^(-1)*residual(j)^2;
        end
        weight(i) = 1/invw(i);
    end
end

%%% To calculate confidence interval
for i = 1:length(x)
    confidence_lower(i) = smoothed_values(i) - sqrt(variance)*smoothing_parameter^(-1/8)*weight(i)^(-3/8)*2^(-3/4);
    calculated_variance(i) = sqrt(variance)*smoothing_parameter^(-1/8)*weight(i)^(-3/8)*2^(-3/4);
end
%%%%%
confidence_upper = 2*smoothed_values - confidence_lower;
