function meanValueOfCAGE(x, cage_peaks)
%This function save the matrix that calculate the mean across peaks of same
%gene and normalization befor calculating mean.
%% Mean only without normalization
% Get the gene name from data first
symbol = [];
for i = 1 : length(cage_peaks)
    temp = regexp(cage_peaks(i,3),'(?<=at_)\w*(?=_5end)', 'match'); % Catch the gene name in 3rd column of cage_peaks
    if isempty(temp{1})
        symbol{i} = char([]);
    else
        symbol{i} = temp{1}{1};
    end
end

gene = unique(symbol);
gene = gene(2:length(gene)); %To remove empty char element
x_mean = zeros(length(gene), size(x,2));

% Calculate the mean value across peaks of same gene
for i = 1:length(gene)
    logical_gene = strcmp(gene{i}, symbol(:));
    gene_peaks = find(logical_gene == 1);
    x_all_peaks = x(gene_peaks,:);
    x_peaks = mean(x_all_peaks,1);
    x_mean(i,:) = x_peaks;
end

%% ZSCORE -> mean
x_zscore = Normalization(x, 'ZSCORE');
x_zscore_mean = zeros(length(gene), size(x,2));

for i = 1:length(gene)
    logical_gene = strcmp(gene{i}, symbol(:));
    gene_peaks = find(logical_gene == 1);
    x_all_peaks = x_zscore(gene_peaks,:);
    x_peaks = mean(x_all_peaks,1);
    x_zscore_mean(i,:) = x_peaks;
end

%% LOG -> mean
x_log = Normalization(x, 'LOG');
x_log_mean = zeros(length(gene), size(x,2));

for i = 1:length(gene)
    logical_gene = strcmp(gene{i}, symbol(:));
    gene_peaks = find(logical_gene == 1);
    x_all_peaks = x_log(gene_peaks,:);
    x_peaks = mean(x_all_peaks,1);
    x_log_mean(i,:) = x_peaks;
end

%% LOG-ZSCORE -> mean
x_log_zscore = Normalization(x, 'double');
x_log_zscore_mean = zeros(length(gene), size(x,2));

for i = 1:length(gene)
    logical_gene = strcmp(gene{i}, symbol(:));
    gene_peaks = find(logical_gene == 1);
    x_all_peaks = x_log_zscore(gene_peaks,:);
    x_peaks = mean(x_all_peaks,1);
    x_log_zscore_mean(i,:) = x_peaks;
end

end

function x = Normalization(x_data, norm)
switch norm
    case 'without'
        x = x_data;
    case 'ZSCORE'
        x = zscore(x_data);
    case 'LOG'
        x = log(1+x_data);
    case 'double'
        x_temp = log(1+x_data);
        x = zscore(x_temp);
end
end