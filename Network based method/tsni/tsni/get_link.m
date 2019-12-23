function get_link(nw_matrix, gene)

result(1,:) = {'TF', 'Gene', 'value'};
k = 2;
for i = 1: size(nw_matrix,1)
    for j = 1:size(nw_matrix,2)
        result{k, 1} = gene{i};
        result{k, 2} = gene{j};
        result{k, 3} = nw_matrix(i, j);
        k = k + 1;
    end
end

save('insilico_size100_5_snr15_TSNI.mat', 'result');
end




