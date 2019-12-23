function [A,B,total_states] = tsni(data,time_points,perturbation,principle_components_to_use)

% [A,B,total_states] = tsni(data,time_points,perturbation,principle_components)
% data: no_genes X time_points format, in which first time points is
% scaled to zero. 
% time_points: Time points at which measurement are done. It is a row
% vector.
% perturbation: It is a row vector and it is an input to the system. Size
% is p X no_time_points where p is number of different perturbation.
% Usually p = 1
% principle_components: Number of principle componets to
% use. If you dont know how many principle components to use, put it to 0.
% A: Recovered network. It is a fully connected network.
% B: Targets of perturbations. It is the score of each gene to be the
% target of perturbation.
% total_states: Number of principle components used.
% There is a parameter in code which is smoothing parameter. Default
% value is 0.8. It can be changed to see how the profiles of the gene
% changes. It shouldn't be too low. 

s = size(data);
no_data_points = s(2);
no_genes = s(1);
u_input = perturbation;

sliding_window_size = 5;   %%% Sliding window size to calculate the weight on each time point.
smoothing_parameter = .8;  %%% Smoothing parameter  

x = time_points;
to_interpolate(1:2:2*length(x)-1) = x;
to_interpolate(2:2:2*length(x)-1) = x(1:length(x)-1) + diff(x)/2; 

for i = 1:size(perturbation,1)
    pert = perturbation(i,:);
    % this part is just to handle constant perturbation. Will change it in
    % future to handle it properly. 
    if pert(2:end)==pert(2)
        u_model(i,1) = pert(1);
        u_model(i,2:2*length(pert)-1) = pert(2);
    end

    u_model(i,1:2:2*length(pert)-1) = pert;
    u_model(i,2:2:2*length(pert)-1) = pert(1:length(perturbation)-1) + diff(pert)/2;
end
count_polynomial_coefs_spline = 0;
for j = 1:length(to_interpolate)
    if to_interpolate(j) < (x(2) - x(1))
        count_polynomial_coefs_spline = count_polynomial_coefs_spline + 1;
    end
end
clear dd sss sd
for i = 1:no_genes
    y = data(i,:);
    [scss,polynomial_coefs_spline,confidence_lower,confidence_upper,smoothed_values] = interpolation_data(x,y,smoothing_parameter,sliding_window_size); %%% Smoothing of data and interpolating it.
    dd(i,1:count_polynomial_coefs_spline) = polyval(polynomial_coefs_spline,to_interpolate(1:count_polynomial_coefs_spline));
    dd(i,count_polynomial_coefs_spline+1:length(to_interpolate)) = fnval(scss,to_interpolate(count_polynomial_coefs_spline+1:length(to_interpolate)));
    sss(i,:) = confidence_upper - smoothed_values;
    sss(i,1) = 0;
    sd(i,:) = spline(x,sss(i,:),to_interpolate);
end

datas = dd;

[A,B,total_states] = discretization(x,datas,u_model,to_interpolate,principle_components_to_use); %%% Calls discretization function to infer network.



