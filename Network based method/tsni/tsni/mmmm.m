x_m = [];
idx = 0;
for i = 1:length(unique(time))
    temp = mean(x(:,idx+1:idx+10),2);
    x_m(:,i)=temp;
    idx = idx + 10;
end

time = unique(time);
p = ones(size(x_m,2),1);
save('simulated_5_snr15_mean.mat', 'x_m', 'gene', 'time', 'p');