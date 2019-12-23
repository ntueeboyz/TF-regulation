

data_xlsx = 'gold_standard.xlsx';
%% INPUT
%% data_xlsx => excel table with the list of TF and their regulated genes (experimentally validated)

% %% -------- TO RUN only ONCE -------- %%
% % It may take some seconds/minutes to run
% tic
% % Get automatically the latest NCBI human gene-info file
% ftpobj = ftp('ftp.ncbi.nlm.nih.gov');
% cd(ftpobj,'gene/DATA/GENE_INFO/Mammalia');
% mget(ftpobj,'Homo_sapiens.gene_info.gz');
% file = gunzip('*.gz');
% movefile('Homo_sapiens.gene_info', 'Homo_sapiens.gene_info.txt')
% toc
% % ------------------------------------- %

%% Main function
[num,TF_regulation,~] = xlsread(data_xlsx,'TF');
table = TF_regulation(2:end,1:7);

% First convert the gene name to entrez gene id
load gene_ID.mat
geneNames = unique(table(:,1:2));
[genes,geneNames_problems] = fromsymboltoentrez(geneNames);

table_mod= table;
for i=1:length(geneNames_problems)
    [r,c]= find(strcmp(table_mod(:,1:2),geneNames_problems{i}));
    table_mod(r,:)=[];
end

% Regulations experimentally validated by means, among others, of RNAseq
idx_RNAseqexp = cellfun(@(s) ~isempty(strfind(s, 'RNAseq')), table_mod(:,6)); %find the RNAseq method in the list
table_RNAseqexp = table_mod(idx_RNAseqexp,:);

% Regulations not experimentally validated by means, among others, of RNAseq
table_noRNAseqexp = table_mod;
table_noRNAseqexp(idx_RNAseqexp,:)=[]; %create the table without RNAseq
downregulation_table = table_noRNAseqexp(strcmp(table_noRNAseqexp(:,3),'downregulation'),:); %downregulation
upregulation_table = table_noRNAseqexp(strcmp(table_noRNAseqexp(:,3),'upregulation'),:); %upregulation
noInfo_regulation_table = table_noRNAseqexp(strcmp(table_noRNAseqexp(:,3),'-'),:); %not known

% Normalization and correlation
% Do the ZSCORE, LOG and LOG+ZSCORE
norm = {'without', 'ZSCORE', 'LOG'};
corr = {'Pearson','Spearman','SDCP','SDCS'};
% Create result table
statistic_result(1,:) = {'Data type', 'Sensitivity', 'Specificity', 'F1 score for upregulation', 'F1 score for downregulation','AUC', 'Precision-positive', 'Precision-negtive'};
idx_statistics_result = 2; %index of row number

%% Microarray dataset
load CCLE_microarray.mat
for i=1:length(corr)
    for j=1:length(norm)
        xmicroarray = Normalization(x_microarray, norm{j}); %Normalization
        %The function trend_entrezGeneID is using the entrez gene id
        resultsD_CCLE_microarray = trend_entrezGeneID(xmicroarray, entez_gene_id_data, downregulation_table, genes, 'downregulation', corr{i}); %Downregulation
        resultsU_CCLE_microarray = trend_entrezGeneID(xmicroarray, entez_gene_id_data, upregulation_table, genes, 'upregulation', corr{i}); %Upregulation
        
        result = ['./result_all/microarray_', norm{j},'_', corr{i}, '.mat'];
        save(result, 'resultsU_CCLE_microarray', 'resultsD_CCLE_microarray');
        
        % Sensitivity, F1 score, AUC
        sensitivity_microarray_D = sum([resultsD_CCLE_microarray{2:end,5}])/length([resultsD_CCLE_microarray{2:end,5}]); %fraction of catched downregulation
        sensitivity_microarray_U = sum([resultsU_CCLE_microarray{2:end,5}])/length([resultsU_CCLE_microarray{2:end,5}]); %fraction of catched upregulation
        precision_microarray_D = sum([resultsD_CCLE_microarray{2:end,5}])/(sum([resultsD_CCLE_microarray{2:end,5}])+sum([resultsU_CCLE_microarray{2:end,5}] == 0)); %Precision of downregulation
        precision_microarray_U = sum([resultsU_CCLE_microarray{2:end,5}])/(sum([resultsU_CCLE_microarray{2:end,5}])+sum([resultsD_CCLE_microarray{2:end,5}] == 0)); %Precision of downregulation
        F1_score_microarray_D = 2*(precision_microarray_D*sensitivity_microarray_D/(sensitivity_microarray_D+precision_microarray_D));
        F1_score_microarray_U = 2*(precision_microarray_U*sensitivity_microarray_U/(sensitivity_microarray_U+precision_microarray_U));
        auc_microarray = AUC(sensitivity_microarray_D, sensitivity_microarray_U);
        data_name = ['microarray_',norm{j},'_',corr{i}];
        
        % Upregulation is positive class 
        statistic_result{idx_statistics_result, 1} = data_name;
        statistic_result{idx_statistics_result, 2} = sensitivity_microarray_U;
        statistic_result{idx_statistics_result, 3} = sensitivity_microarray_D;
        statistic_result{idx_statistics_result, 4} = F1_score_microarray_U;
        statistic_result{idx_statistics_result, 5} = F1_score_microarray_D;
        statistic_result{idx_statistics_result, 6} = auc_microarray;
        statistic_result{idx_statistics_result, 7} = precision_microarray_U;
        statistic_result{idx_statistics_result, 8} = precision_microarray_D;
        
        idx_statistics_result = idx_statistics_result + 1;
        
    end
end
%% Atlas RNA-seq dataset
load data_Atlas_RNAseq.mat
for i=1:length(corr)
    for j=1:length(norm)
        
        xrnaseq = Normalization(x_mean, norm{j}); %Normalization
        %The function trend_geneName is a bit different for this dataset, because
        %it is using the gene names
        resultsD_Atlas_RNAseq = trend_geneName(xrnaseq,gene,downregulation_table,genes, 'downregulation',corr{i}); %Downregulation
        resultsU_Atlas_RNAseq = trend_geneName(xrnaseq,gene,upregulation_table,genes, 'upregulation', corr{i}); %Upregulation
        
        result = ['./result_all/RNAseq_', norm{j},'_', corr{i}, '.mat'];
        save(result, 'resultsU_Atlas_RNAseq', 'resultsD_Atlas_RNAseq');
        
        % Sensitivity, F1 score, AUC
        sensitivity_Atlas_D = sum([resultsD_Atlas_RNAseq{2:end,5}])/length([resultsD_Atlas_RNAseq{2:end,5}]); %fraction of catched downregulation
        sensitivity_Atlas_U = sum([resultsU_Atlas_RNAseq{2:end,5}])/length([resultsU_Atlas_RNAseq{2:end,5}]); %fraction of catched upregulation
        precision_Atlas_D = sum([resultsD_Atlas_RNAseq{2:end,5}])/(sum([resultsD_Atlas_RNAseq{2:end,5}])+sum([resultsU_Atlas_RNAseq{2:end,5}] == 0)); %Precision of downregulation
        precision_Atlas_U = sum([resultsU_Atlas_RNAseq{2:end,5}])/(sum([resultsU_Atlas_RNAseq{2:end,5}])+sum([resultsD_Atlas_RNAseq{2:end,5}] == 0)); %Precision of downregulation
        F1_score_Atlas_D = 2*(precision_Atlas_D*sensitivity_Atlas_D/(sensitivity_Atlas_D+precision_Atlas_D));
        F1_score_Atlas_U = 2*(precision_Atlas_U*sensitivity_Atlas_U/(sensitivity_Atlas_U+precision_Atlas_U));
        auc_Atlas = AUC(sensitivity_Atlas_D, sensitivity_Atlas_U);
        data_name = ['RNAseq_',norm{j},'_',corr{i}];
        
        % Upregulation is positive class
        statistic_result{idx_statistics_result, 1} = data_name;
        statistic_result{idx_statistics_result, 2} = sensitivity_Atlas_U;
        statistic_result{idx_statistics_result, 3} = sensitivity_Atlas_D;
        statistic_result{idx_statistics_result, 4} = F1_score_Atlas_U;
        statistic_result{idx_statistics_result, 5} = F1_score_Atlas_D;
        statistic_result{idx_statistics_result, 6} = auc_Atlas;
        statistic_result{idx_statistics_result, 7} = precision_Atlas_U;
        statistic_result{idx_statistics_result, 8} = precision_Atlas_D;
        
        idx_statistics_result = idx_statistics_result + 1;
        
    end
end
%% FANTOM 5 CAGE dataset
% Please create yourself the function trend_xx, as done for the previous datasets (trend_geneName/trend_entrezGeneID)
% Inside the function, of course, you will have to compute the correlation between the TF and the gene it regulates.
% Here below find an example on how to compute the correlation between GATA6 and DKK1
load FANTOM5_CAGE.mat

% This part is that compute mean first then normalize x
for i=1:length(corr)
    for j=1:length(norm)
        
        xcage = Normalization(x_mean, norm{j}); %Normalization
        resultsD_FANTOM5_CAGE = trend_geneName(xcage,gene,downregulation_table,genes, 'downregulation', corr{i}); %Downregulation
        resultsU_FANTOM5_CAGE = trend_geneName(xcage,gene,upregulation_table,genes, 'upregulation', corr{i}); %Upregulation
        
        result = ['./result_all/CAGE_', norm{j},'_', corr{i}, '.mat'];
        save(result, 'resultsU_FANTOM5_CAGE', 'resultsD_FANTOM5_CAGE');
        
        % Sensitivity, F1 score, AUC
        sensitivity_FANTOM5_CAGE_D = sum([resultsD_FANTOM5_CAGE{2:end,5}])/length([resultsD_FANTOM5_CAGE{2:end,5}]); %fraction of catched downregulation
        sensitivity_FANTOM5_CAGE_U = sum([resultsU_FANTOM5_CAGE{2:end,5}])/length([resultsU_FANTOM5_CAGE{2:end,5}]); %fraction of catched upregulation
        precision_FANTOM5_CAGE_D = sum([resultsD_FANTOM5_CAGE{2:end,5}])/(sum([resultsD_FANTOM5_CAGE{2:end,5}])+sum([resultsU_FANTOM5_CAGE{2:end,5}] == 0)); %Precision of downregulation
        precision_FANTOM5_CAGE_U = sum([resultsU_FANTOM5_CAGE{2:end,5}])/(sum([resultsU_FANTOM5_CAGE{2:end,5}])+sum([resultsD_FANTOM5_CAGE{2:end,5}] == 0)); %Precision of downregulation
        F1_score_CAGE_D = 2*(precision_FANTOM5_CAGE_D*sensitivity_FANTOM5_CAGE_D/(sensitivity_FANTOM5_CAGE_D+precision_FANTOM5_CAGE_D));
        F1_score_CAGE_U = 2*(precision_FANTOM5_CAGE_U*sensitivity_FANTOM5_CAGE_U/(sensitivity_FANTOM5_CAGE_U+precision_FANTOM5_CAGE_U));
        auc_FANTOM5_CAGE = AUC(sensitivity_FANTOM5_CAGE_D, sensitivity_FANTOM5_CAGE_U);
        
        data_name = ['CAGE_',norm{j},'_',corr{i}];
        
        % Up regulation is positive class
        statistic_result{idx_statistics_result, 1} = data_name;
        statistic_result{idx_statistics_result, 2} = sensitivity_FANTOM5_CAGE_U;
        statistic_result{idx_statistics_result, 3} = sensitivity_FANTOM5_CAGE_D;
        statistic_result{idx_statistics_result, 4} = F1_score_CAGE_U;
        statistic_result{idx_statistics_result, 5} = F1_score_CAGE_D;
        statistic_result{idx_statistics_result, 6} = auc_FANTOM5_CAGE;
        statistic_result{idx_statistics_result, 7} = precision_FANTOM5_CAGE_U;
        statistic_result{idx_statistics_result, 8} = precision_FANTOM5_CAGE_D;
        
        idx_statistics_result = idx_statistics_result + 1;
        
    end
end
save('./result_all/result.mat', 'statistic_result', '-append');

%% Convert gene symbol to gene ID by using NCBI API
function [genes,geneNames_problems] = fromsymboltoentrez(searchList)
fprintf('\n');
textprogressbar('Converting the gene symbols to entrez gene ID with NCBI API: ');
tic
baseURL = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
eutil = 'esearch.fcgi?';
dbParam = 'db=gene';
usehistoryParam = '&usehistory=y';

load gene_ID.mat
idx = [];
idx_matches =[];
idx_all = 1:length(searchList);
[ia1,ia2]=unique(round(idx_all/length(idx_all),2)*100);
for mi=1:length(idx_all)
    currentGene = searchList{mi};
    try
        genes.(searchList{mi}).entrez;
    catch
        %fprintf('\nsearchList{%d} -> %s\n', mi, searchList{mi});
        if ~isempty(strfind(currentGene, '/'))|| ~isempty(strfind(currentGene, '(')) || ~isempty(strfind(currentGene, '-'))
            idx = [idx,mi];
            continue
        else
            termParam = strcat('&term=',currentGene,'[Gene Name]) AND Homo sapiens[Organism]');
            esearchURL = [baseURL, eutil, dbParam, termParam, usehistoryParam];
            
            numberTrials = 1;
            totalTrials = 3;
            while numberTrials > 0 && numberTrials < totalTrials
                try
                    searchReport = webread(esearchURL ,weboptions('Timeout',30));
                    idMatches = regexp(searchReport ,'<Id>\d+</Id>', 'match');
                    if length(idMatches)>1
                        idx_matches=[idx_matches,mi];
                    end
                    temp = [];
                    for ix=1:length(idMatches)
                        temp{ix} = regexp(idMatches{ix},'((?<=<Id>).*(?=<\/Id>))', 'match'); %Preserve only digitals
                    end
                    if isempty(temp)|| length(temp)>1
                        idx = [idx,mi];
                    else
                        genes.(currentGene).name = currentGene;
                        genes.(currentGene).entrez = [temp{:}];
                    end
                    numberTrials = 0;
                catch ME
                    fprintf('\nOops! There is an error with ID: %s (TRIAL %d/%d)\n', currentGene, numberTrials, totalTrials);
                    fprintf('WEBREAD without success: %s\n', ME.message);
                    numberTrials = numberTrials + 1;
                end
            end
            if (numberTrials >= 3)
                genes.(currentGene).name = currentGene;
                genes.(currentGene).entrez = NaN;
            end
        end
    end
    save('gene_ID.mat', 'genes', '-append');
    if ismember(idx_all(mi),ia2)
        textprogressbar(round(idx_all(mi)/length(idx_all),2)*100);
    end
end
textprogressbar('done');
geneNames_problems=searchList(idx);
end

%% Progress bar
function textprogressbar(c)
% This function creates a text progress bar. It should be called with a
% STRING argument to initialize and terminate. Otherwise the number correspoding
% to progress in % should be supplied.
% INPUTS:   C   Either: Text string to initialize or terminate
%                       Percentage number to show progress
% OUTPUTS:  N/A
% Example:  Please refer to demo_textprogressbar.m

% Author: Paul Proteus (e-mail: proteus.paul (at) yahoo (dot) com)
% Version: 1.0
% Changes tracker:  29.06.2010  - First version

% Inspired by: http://blogs.mathworks.com/loren/2007/08/01/monitoring-progress-of-a-calculation/

%% Initialization
persistent strCR;           %   Carriage return pesistent variable

% Vizualization parameters
strPercentageLength = 10;   %   Length of percentage string (must be >5)
strDotsMaximum      = 10;   %   The total number of dots in a progress bar

%% Main

if isempty(strCR) && ~ischar(c),
    % Progress bar must be initialized with a string
    error('The text progress must be initialized with a string');
elseif isempty(strCR) && ischar(c),
    % Progress bar - initialization
    fprintf('%s',c);
    strCR = -1;
elseif ~isempty(strCR) && ischar(c),
    % Progress bar  - termination
    strCR = [];
    fprintf([c '\n']);
elseif isnumeric(c)
    % Progress bar - normal progress
    c = floor(c);
    percentageOut = [num2str(c) '%%'];
    percentageOut = [percentageOut repmat(' ',1,strPercentageLength-length(percentageOut)-1)];
    nDots = floor(c/100*strDotsMaximum);
    dotOut = ['[' repmat('.',1,nDots) repmat(' ',1,strDotsMaximum-nDots) ']'];
    strOut = [percentageOut dotOut];
    
    % Print it on the screen
    if strCR == -1,
        % Don't do carriage return during first run
        fprintf(strOut);
    else
        % Do it during all the other runs
        fprintf([strCR strOut]);
    end
    
    % Update carriage return
    strCR = repmat('\b',1,length(strOut)-1);
    
else
    % Any other unexpected input
    error('Unsupported argument type');
end
end

%% Evaluation correlation-based methods (for gene ID data)
function results = trend_entrezGeneID(x,entrez_geneID,table_regulation,genes, type_regulation,type_corr)
% results=cell(size(table_regulation,1)+1,6);
results(1,:)={'TF','Gene','Correlation','Type of regulation (exp)','Catched?', 'Deviation', 'p value of corr', 'Index TF', 'Index gene'};
if strcmp(type_regulation,'downregulation')
    trend = -1;
elseif strcmp(type_regulation,'upregulation')
    trend = 1;
else
    trend = 0;
end

j = 1; %To manage the row of the table in the loop
for i=1:length(table_regulation(:,1))
    logical_TF = strcmp(entrez_geneID, genes.(table_regulation{i,1}).entrez);
    logical_gene = strcmp(entrez_geneID,genes.(table_regulation{i,2}).entrez);
    if (sum(logical_TF)== 0) || (sum(logical_gene)== 0) || strcmp(table_regulation{i,1},table_regulation{i,2})
        continue
    else
        results{j+1,1}=table_regulation{i,1}; %TF name
        results{j+1,2}=table_regulation{i,2}; %Gene name
        results{j+1,4}=type_regulation;
        % Compute correlation
        idx_TF = find(logical_TF);
        idx_gene = find(logical_gene);
        c = corr_v2(x(idx_TF,:)', x(idx_gene,:)', type_corr); %Correlation
        results{j+1,3}=c;

        if strcmp(type_regulation,'upregulation')||  strcmp(type_regulation,'downregulation')
            results{j+1,5} = sign(results{j+1,3}) == trend;
            results{j+1,6} = abs(results{j+1,3}-trend);
        elseif strcmp(type_regulation,'not known')
            results{j+1,5}='-';
            results{j+1,6}='-';
        end
        results{j+1, 8} = idx_TF;
        results{j+1, 9} = idx_gene;
        j = j + 1;
    end
end
end

%% Evaluation correlation-based methods (for gene symbol data)
function results = trend_geneName(x,gene_name,table_regulation,genes, type_regulation, type_corr)
% results=cell(size(table_regulation,1)+1,6);
results(1,:)={'TF','Gene','Correlation','Type of regulation (exp)','Catched?','Deviation', 'p value of corr', 'Index TF', 'Index gene'};
if strcmp(type_regulation,'downregulation')
    trend = -1;
elseif strcmp(type_regulation,'upregulation')
    trend = 1;
else
    trend = 0;
end

j = 1; %To manage the row of the table in the loop
for i=1:length(table_regulation(:,1))
    logical_TF = strcmp(gene_name,genes.(table_regulation{i,1}).name);
    logical_gene = strcmp(gene_name,genes.(table_regulation{i,2}).name);

    if (sum(logical_TF)== 0) || (sum(logical_gene)== 0) || strcmp(table_regulation{i,1},table_regulation{i,2})
        continue
    else
        results{j+1,1}=table_regulation{i,1}; %TF name
        results{j+1,2}=table_regulation{i,2}; %Gene name
        results{j+1,4}=type_regulation;
        % Compute correlation
        idx_TF = find(logical_TF);
        idx_gene = find(logical_gene);

        c = corr_v2(x(idx_TF,:)', x(idx_gene,:)', type_corr); %Correlation
        results{j+1,3}=c;

        if strcmp(type_regulation,'upregulation')||  strcmp(type_regulation,'downregulation')
            results{j+1,5} = sign(results{j+1,3}) == trend;
            results{j+1,6} = abs(results{j+1,3}-trend);
        elseif strcmp(type_regulation,'not known')
            results{j+1,5}='-';
            results{j+1,6}='-';
        end
        results{j+1, 8} = idx_TF;
        results{j+1, 9} = idx_gene;
        j = j + 1;
    end
end
end

%% Apply cut-off
function results = trend_threshold(result_table, threshold)
results(1,:) = {'TF', 'Gene', 'Correlation', 'Type of regulation', 'Catched', 'Deviation', 'p value'};
j = 2;
for i=2:length(result_table(:,1))
    if result_table{i, 7} < threshold
        results{j, 1} = result_table{i, 1};
        results{j, 2} = result_table{i, 2};
        results{j, 3} = result_table{i, 3};
        results{j, 4} = result_table{i, 4};
        results{j, 5} = result_table{i, 5};
        results{j, 6} = result_table{i, 6};
        results{j, 7} = result_table{i, 7};
        j = j + 1;
    else
        continue
    end
end
end

%% Calculate AUC value
function auc = AUC(sensitivity, specificity)
if isnan(sensitivity) || isnan(specificity)
    auc = NaN;
else
    y = [0;sensitivity;1];
    x = [0;1-specificity;1];
    auc = trapz(x,y);
end
end

%% All normalization method
function x = Normalization(x_data, norm)
switch norm
    case 'without'
        x = x_data;
    case 'ZSCORE'
        x = zscore(x_data);
    case 'LOG'
        x = log(1+x_data);
end
end

%% All correlation-based method
function c = corr_v2(x,y,type)

switch type
    case 'Pearson'
        c = corr(x, y);
        
    case 'Spearman'
        c = corr(x, y , 'Type', 'Spearman');
        
    case 'SDCP'
        [d_x, d_y] = increment(x,y);
        c = corr(d_x,d_y);
    case 'SDCS'
        [d_x, d_y] = increment(x,y);
        c = corr(d_x,d_y,'Type','Spearman');
end
end

function [d_x, d_y] = increment(x,y)
d_x = [];
d_y = [];
for i = 1:length(x)-1
    for j = i+1:length(x)
        deltaX = x(i) - x(j);
        deltaY = y(i) - y(j);
        d_x = [d_x; deltaX];
        d_y = [d_y; deltaY];
    end
end
end

%% Signed distance correlation
function c = signed_dist(x,y)
dist_corr = distcor(x, y);
if isnan(dist_corr)
    c = NaN;
else
    distance=[];
    for i = 1:length(x)-1
        for j = i+1:length(x)
            deltaX = x(i) - x(j);
            deltaY = y(i) - y(j);
            twoPoint = [x(i),y(i); x(j),y(j)];
            dist = pdist(twoPoint, 'euclidean');
            value = deltaX*deltaY*dist;
            distance = [distance, value];
        end
    end
    c = sign(mean(distance))*dist_corr;
end
end

function dcor = distcor(x,y)
% This function calculates the distance correlation between x and y.
% Reference: http://en.wikipedia.org/wiki/Distance_correlation
% Date: 18 Jan, 2013
% Author: Shen Liu (shen.liu@hotmail.com.au)
% Check if the sizes of the inputs match

if size(x,1) ~= size(y,1)
    error('Inputs must have the same number of rows')
end
% Delete rows containing unobserved values
N = any([isnan(x) isnan(y)],2);
x(N,:) = [];
y(N,:) = [];
% Calculate doubly centered distance matrices for x and y
a = pdist2(x,x);
mcol = mean(a);
mrow = mean(a,2);
ajbar = ones(size(mrow))*mcol;
akbar = mrow*ones(size(mcol));
abar = mean(mean(a))*ones(size(a));
A = a - ajbar - akbar + abar;
b = pdist2(y,y);
mcol = mean(b);
mrow = mean(b,2);
bjbar = ones(size(mrow))*mcol;
bkbar = mrow*ones(size(mcol));
bbar = mean(mean(b))*ones(size(b));
B = b - bjbar - bkbar + bbar;
% Calculate squared sample distance covariance and variances
dcov = sum(sum(A.*B))/(size(mrow,1)^2);
dvarx = sum(sum(A.*A))/(size(mrow,1)^2);
dvary = sum(sum(B.*B))/(size(mrow,1)^2);
% Calculate the distance correlation
dcor = sqrt(dcov/sqrt(dvarx*dvary));
end

%% Mutual information 
function [MI , NMI1, NMI2] = mi_gg(x, y, biascorrect, demeaned)
% MI_GG Mutual information (MI) between two Gaussian variables in bits
%   I = mi_gg(x,y) returns the MI between two (possibly multidimensional)
%   Gassian variables, x and y, with bias correction.
%   If x and/or y are multivariate rows must correspond to samples, columns
%   to dimensions/variables. (Samples first axis) 
%
%   biascorrect : true / false option (default true) which specifies whether
%   bias correction should be applied to the esimtated MI.
%   demeaned : false / true option (default false) which specifies whether the
%   input data already has zero mean (true if it has been copula-normalized)

%   The code is from Robin A.A. Ince, "A statistical framework for neuroimaging 
%   data analysis based on mutual information estimated via a gaussian copula", Human
%   Brain Mapping (2016).

% ensure samples first axis for vectors
if isvector(x)
    x = x(:);
end
if isvector(y)
    y = y(:);
end
if ndims(x)~=2 || ndims(y)~=2
    error('mi_gg: input arrays should be 2d')
end
Ntrl = size(x,1);
Nvarx = size(x,2);
Nvary = size(y,2);

if size(y,1) ~= Ntrl
    error('mi_gg: number of trials do not match')
end

% default option values
if nargin<3
    biascorrect = true;
end
if nargin<4
    demeaned = false;
end

% demean data if required
if ~demeaned
    x = bsxfun(@minus,x,sum(x,1)/Ntrl);
    y = bsxfun(@minus,y,sum(y,1)/Ntrl);
end

% joint variable
xy = [x y];
Cxy = (xy'*xy) / (Ntrl - 1);
% submatrices of joint covariance
Cx = Cxy(1:Nvarx,1:Nvarx);
ystart = Nvarx + 1;
Nvarxy = Nvarx + Nvary;
Cy = Cxy(ystart:Nvarxy,ystart:Nvarxy);

chCx = chol(Cx);
chCy = chol(Cy);
chCxy = chol(Cxy);

% entropies in nats
% normalisations cancel for information
HX = sum(log(diag(chCx))); % + 0.5*Nvarx*log(2*pi*exp(1));
HY = sum(log(diag(chCy))); % + 0.5*Nvary*log(2*pi*exp(1));
HXY = sum(log(diag(chCxy))); % + 0.5*(Nvarx+Nvary)*log(2*pi*exp(1));

ln2 = log(2);
if biascorrect
    psiterms = psi((Ntrl - (1:Nvarxy))/2) / 2;
    dterm = (ln2 - log(Ntrl-1)) / 2;
    HX = (HX - Nvarx*dterm - sum(psiterms(1:Nvarx)));
    HY = (HY - Nvary*dterm - sum(psiterms(1:Nvary)));
    HXY = (HXY - Nvarxy*dterm - sum(psiterms));
end

% convert to bits
MI = (HX + HY - HXY) / ln2;

% Normalization

% The normalization suggested by Carlo Cannistraci (July 2019)
NMI1 = MI/(MI+1);

% The normalization suggested by Claudio Duran (July 2019)
NMI2 = -MI/HXY;

end

%% Coagreement
function [c,p,r] = coagreement(x,y,it) 
    % x: first variable
    % y: second variable
    % it: number of iteration for null model (suggestion:1000)
    % r: type of rule, 1 is mean/mean, 2 is median/median, 3 is mean/median and 4 is median/mean

    % co-agreement
    [c,n,r]=cannistraci(x,y);


    %null model
    if nargin == 3
        m = zeros(1,it+1);
        m(1) = c;
        parfor i=2:it+1
            x_perm = x(randperm(n));
            y_perm = y(randperm(n));
            m(i) = coagreement(x_perm, y_perm);
        end
        positionC = find(sort(m, 'ascend') == c);

        if c >= 0
            p = 1 - (positionC(1)-1)/(it+1);
        else
            p = positionC(end)/(it+1);
        end
    else
        p = nan;
    end

end


function  [z,n,r] = cannistraci(x,y)

    n=length(x);

    z1 = 1/n*sum(abs(sign(x-mean(x))+sign(y-mean(y))))-1;
    z2 = 1/n*sum(abs(sign(x-median(x))+sign(y-median(y))))-1;
    z3 = 1/n*sum(abs(sign(x-mean(x))+sign(y-median(y))))-1;
    z4 = 1/n*sum(abs(sign(x-median(x))+sign(y-mean(y))))-1;
    z=max([z1 z2 z3 z4]);
    r=find([z1 z2 z3 z4]==z,1);

end

%% Coagreement neuron network
 function [c,v,p,r,w,a] = coagreement_nn_18July(x,y,type,it)

% INPUT
% x: first variable
% y: second variable
% type: 'linear' or 'sigmoide'(default is sigmoide)
% it: number of iteration for null model (suggestion:1000)
% r: type of rule, 1 is mean/mean, 2 is median/median, 3 is mean/median and 4 is median/mean

% OUTPUT
% c: coagreement_nn (after learning) coefficient
% v: coagreement coefficient
% p-value of c
% r: winner centering rule for c
% w: learned weightes for synapses of coagreement_nn
% a: winner centering rule values for the two naive binary classifiers of the first layer



if isrow(x), x=x'; end %we use column vectors
if isrow(y), y=y'; end %we use column vectors


% co-agreement
[c,v,n,r,w,a]=coagre_net(x,y,type);


%null model
if nargin == 4
                m = zeros(1,it+1);
                m(1) = c;
                parfor i=2:it+1
                    x_perm = x(randperm(n));
                    y_perm = y(randperm(n));
                    m(i) = coagre_net(x_perm, y_perm,type);
                end
                positionC = find(sort(m, 'ascend') == c);
                
                if c >= 0
                    p = 1 - (positionC(1)-1)/(it+1); 
                else
                    p = positionC(end)/(it+1); 
                end
            else
                p = nan;
end

 end

function  [z,v,n,r,w,a] = coagre_net(x,y,type)
 
n=length(x);
           
%%% centering rule: mean/mean

           % binaly classifier layer
           xs=sign(x-mean(x)); ys=sign(y-mean(y));
           % coagreement layer
           cal=abs(xs+ys)-1;  
           % first step output
           o = 1/n*sum(cal); 
           % learning weights update 
           weight=ones(size(cal));
           if o==0
           z1=o; %final output case o==0
           v1=o;
           else
           
           if o>0
           xs(xs==0)=1; ys(ys==0)=1;
           elseif o<0
           xs(xs==0)=-1; ys(ys==0)=1;
           end
           cal=abs(xs+ys)-1;   
           v1 = 1/n*sum(cal); 
                 
           if strcmp(type,'linear')
           weight(cal==sign(-o))=1-abs(o);   
           else
           %weight(cal==sign(-o))=1-(1./(1+exp(-(14*abs(o)-7)))+abs(o))/2; 
           weight(cal==sign(-o))=1-1./(1+exp(-(14*abs(o)-7)));
           end
           z1 = 1/n*sum(cal.*weight);  %final output case o<>0
           end
           w=weight;
           a=[mean(x) mean(y)];
           
%%% centering rule: median/median
           
           % binaly classifier layer
           xs=sign(x-median(x)); ys=sign(y-median(y));
           % coagreement layer
           cal=abs(xs+ys)-1;  
           % first step output
           o = 1/n*sum(cal); 
           % learning weights update 
           weight=ones(size(cal));
           if o==0
           z2=o; %final output case o==0
           v2=o;
           else
           
           if o>0
           xs(xs==0)=1; ys(ys==0)=1;
           elseif o<0
           xs(xs==0)=-1; ys(ys==0)=1;
           end
           cal=abs(xs+ys)-1;   
           v2 = 1/n*sum(cal); 
           
           if strcmp(type,'linear')
           weight(cal==sign(-o))=1-abs(o);   
           else
           %weight(cal==sign(-o))=1-(1./(1+exp(-(14*abs(o)-7)))+abs(o))/2; 
           weight(cal==sign(-o))=1-1./(1+exp(-(14*abs(o)-7)));
           end
           z2 = 1/n*sum(cal.*weight);  %final output case o<>0
           end
           w=[w weight];
           a=[a; median(x) median(y)];
           
%%% centering rule: mean/median
           
           % binaly classifier layer
           xs=sign(x-mean(x)); ys=sign(y-median(y));
           % coagreement layer
           cal=abs(xs+ys)-1;  
           % first step output
           o = 1/n*sum(cal); 
           % learning weights update 
           weight=ones(size(cal));
           if o==0
           z3=o; %final output case o==0
           v3=o;
           else
           
           if o>0
           xs(xs==0)=1; ys(ys==0)=1;
           elseif o<0
           xs(xs==0)=-1; ys(ys==0)=1;
           end
           cal=abs(xs+ys)-1;   
           v3 = 1/n*sum(cal); 
           
           if strcmp(type,'linear')
           weight(cal==sign(-o))=1-abs(o);   
           else
           %weight(cal==sign(-o))=1-(1./(1+exp(-(14*abs(o)-7)))+abs(o))/2; 
           weight(cal==sign(-o))=1-1./(1+exp(-(14*abs(o)-7)));
           end
           z3 = 1/n*sum(cal.*weight);  %final output case o<>0
           end
           w=[w weight];
           a=[a; mean(x) median(y)];
           
%%% centering rule: median/mean
           
           % binaly classifier layer
           xs=sign(x-median(x)); ys=sign(y-mean(y));
           % coagreement layer
           cal=abs(xs+ys)-1;  
           % first step output
           o = 1/n*sum(cal); 
           % learning weights update 
           weight=ones(size(cal));
           if o==0
           z4=o; %final output case o==0
           v4=o;
           else
           
           if o>0
           xs(xs==0)=1; ys(ys==0)=1;
           elseif o<0
           xs(xs==0)=-1; ys(ys==0)=1;
           end
           cal=abs(xs+ys)-1;   
           v4 = 1/n*sum(cal); 
           
           if strcmp(type,'linear')
           weight(cal==sign(-o))=1-abs(o);   
           else
           %weight(cal==sign(-o))=1-(1./(1+exp(-(14*abs(o)-7)))+abs(o))/2; 
           weight(cal==sign(-o))=1-1./(1+exp(-(14*abs(o)-7)));
           end
           z4 = 1/n*sum(cal.*weight);  %final output case o<>0
           end
           w=[w weight];
           a=[a; median(x) mean(y)];

v=[v1 v2 v3 v4];  
[~,r]=max(abs(v));
v=v(r);

z=[z1 z2 z3 z4];
[~,r]=max(abs(z));
z=z(r);
w=w(:,r);
a=a(r,:);

end