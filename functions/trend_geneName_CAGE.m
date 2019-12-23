function results =trend_geneName_CAGE(x, cage_peaks, table_regulation, genes, type_regulation)
results=cell(size(table_regulation,1)+1,6); %Create the result table.
results(1,:)={'TF','Gene','Pearson Correlation','Type of regulation (exp)','Catched?','Deviation'}; %First row of the result table.
if strcmp(type_regulation,'downregulation')
    trend = -1;
elseif strcmp(type_regulation,'upregulation')
    trend = 1;
else
    trend = 0;
end
    
cellfind = @(string)(@(cell_contents)(~isempty(strfind(cell_contents,string))));

for i=1:length(table_regulation)
    results{i+1,1}=table_regulation{i,1}; %TF name
    results{i+1,2}=table_regulation{i,2}; %Gene name
    results{i+1,4}=type_regulation;
    
   
    %Find CAGE peaks associated to transcription factor
    logical_TF = cellfun(cellfind(genes.(table_regulation{i,1}).name),cage_peaks(:,2));
    %Find CAGE peaks associated to regulated gene
    logical_gene = cellfun(cellfind(genes.(table_regulation{i,2}).name),cage_peaks(:,2));

    
    if (sum(logical_TF)< 1) || (sum(logical_gene)< 1) || strcmp(table_regulation{i,1},table_regulation{i,2})
        results{i+1,5} = [];
        results{i+1,6} = [];
        results{i+1,3} = [];
    else
        TF_peaks=find(logical_TF==1);
        gene_peaks=find(logical_gene==1);
        %Compute the expression average across the peaks
        x_all_TF = x(TF_peaks,:);
        x_TF = mean(x_all_TF,1);
        x_all_gene = x(gene_peaks,:);
        x_gene = mean(x_all_gene,1);
        
        %Compute the correlation between the expression 
        results{i+1,3} = corr(x_TF(1,:)',x_gene(1,:)');
        
        if strcmp(type_regulation,'upregulation')||  strcmp(type_regulation,'downregulation')
            results{i+1,5} = sign(results{i+1,3}) == trend;
            results{i+1,6} = abs(results{i+1,3}-trend);
        elseif strcmp(type_regulation,'not known')
            results{i+1,5}='-';
            results{i+1,6}='-';
        end
    end
end

end