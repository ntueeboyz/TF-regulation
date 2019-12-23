function get_TF (gene_data, TF_data, dataset)
TF = [];
idx = 1;
for i = 1:length(gene_data)
    if ~isempty(sum(strcmp(gene_data{i},TF_data)))
        TF{idx} = gene_data{i};
        idx = idx + 1;
    end
end
save(dataset, 'TF', '-append')
end