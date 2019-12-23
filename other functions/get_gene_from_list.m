function get_gene_from_list(upregulation_table, downregulation_table)
gene_list = [];
idx = 1;
for i = 1:length(downregulation_table(2:end,1))
    gene_list{idx} = downregulation_table{i+1, 1};
    idx = idx + 1;
end

for j = 1:length(downregulation_table(2:end,2))
    gene_list{idx} = downregulation_table{j+1, 2};
    idx = idx + 1;
end

for k = 1:length(upregulation_table(2:end,1))
    gene_list{idx} = upregulation_table{k+1, 1};
    idx = idx + 1;
end

for l = 1:length(upregulation_table(2:end,2))
    gene_list{idx} = upregulation_table{l+1, 2};
    idx = idx + 1;
end

gene_list_RNAseq = unique(gene_list);

disp(length(gene_list_RNAseq));
save('gene_from_list.mat', 'gene_list_RNAseq', '-append');
end