function [A_svd,B_svd,total_states] = discretization(x,datas,u_model,to_interpolate,principle_components)

% This function is called by tsni.m 
% A_svd:          Network inferred 
% B_svd:          Targets inferred
% total_states:   Number of states for principle components

no_inputs = size(u_model,1);    %%% Checks how many input are there. 
s = size(datas);              
no_time_points = s(2);          %%% Find how many time points are there in data.
no_of_genes = s(1);             %%% Find how many genes are in network.  
clear s

Y1 = datas(1:no_of_genes,2:no_time_points);                             %%%  Builds the data matrix for Right hand side.
ph = [datas(:,1:no_time_points-1)' u_model(:,1:no_time_points-1)']';    %%% Builds the data matrix for left hand. 

[s1,s2,s3] = svd(ph);
diagonal = diag(s2);

states = 0;
for j =  1:length(diagonal)				%%% Counts the number of Principle Components (PC) based on singular value ratio.
    condition = s2(1,1)/s2(j,j);
    if condition <= 100
        states = states +1;
    end
end

clear diagonal condition
Z = s1'*ph;
if principle_components == 0				%%% If number of PC is not specified by user then selects then use the calculated PC's.  
    total_states = states;
else
    total_states = principle_components;
end

Z = Z(1:total_states,:);

for i=1:no_of_genes
    result1 = Y1(i,:)*pinv(Z);
    A1 = result1.^2;
    result1(1,total_states+1:no_of_genes+no_inputs) = 0;
    result(i,:) = result1*s1';
end

%%% Results
result_svd = result;
Ad_svd = result_svd(:,1:no_of_genes);   %%% Network A in discrete space
Bd_svd = result_svd(:,no_of_genes+1:no_of_genes+no_inputs);   %%% Targets inferred in discrete space
%%%%

A_svd = (2/(to_interpolate(2) - to_interpolate(1)))*(Ad_svd-1*eye(no_of_genes))*pinv(Ad_svd+1*eye(no_of_genes));  %%% Bilinear transformation 
B_svd = inv(Ad_svd - 1*eye(no_of_genes))*A_svd*Bd_svd; %%% Bilinear transformation

