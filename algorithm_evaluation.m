function algorithm_evaluation(data_xlsx, exp_data, gene_name_type)
%parpool(8);
% data_xlsx = 'gold_standard.xlsx';
% exp_data = 'CCLE_microarray.mat';
% gene_name_type = 'gene ID';
%% INPUT
% data_xlsx => gold standard excel table with TF column, regulated gene
% column, and regulation type column (upregulation or downregulation).
% exp_data => .mat file with gene symbol/ID and expression vaue and should
% contain exactly variable name: x (expression value with no_gene X
% time-points/samples) and gene (gene symbol/ID)
% gene_name_type => if gene in exp_data is gene ID then enter 'gene ID'
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
[num,TF_regulation,~] = xlsread(data_xlsx);
table = TF_regulation(2:end,1:3);

% First convert the gene name to entrez gene id
if nargin > 2
    if strcmp(gene_name_type, 'gene ID')
        load gene_ID.mat
        geneNames = unique(table(:,1:2));
        [genes,geneNames_problems] = fromsymboltoentrez(geneNames);
    end
    
    
    table_mod= table;
    for i=1:length(geneNames_problems)
        [r,c]= find(strcmp(table_mod(:,1:2),geneNames_problems{i}));
        table_mod(r,:)=[];
    end
    
else
    table_mod = table;
end

downregulation_table = table_mod(strcmp(table_mod(:,3),'downregulation'),:); %downregulation
upregulation_table = table_mod(strcmp(table_mod(:,3),'upregulation'),:); %upregulation
noInfo_regulation_table = table_mod(strcmp(table_mod(:,3),'-'),:); %not known

% Normalization and correlation
norm = {'without','LOG','ZSCORE'};
corr = {'Pearson','Spearman','SDCP','SDCS'};
% Create result table
statistic_result(1,:) = {'Data type', 'Sensitivity', 'Specificity', 'F1 score for upregulation', 'F1 score for downregulation','AUC', 'Precision-positive', 'Precision-negtive'};
idx = 2; %index of row number

%% Dataset
load (exp_data)
load gene_ID.mat

for i=1:length(corr)
    for j=1:length(norm)
        x_norm = Normalization(x, norm{j}); %Normalization
        %The function trend_entrezGeneID is using the entrez gene id
        if nargin == 2
            resultsD = trend_geneName(x_norm, gene, downregulation_table, genes, 'downregulation', corr{i}, 'symbol'); %Downregulation
            resultsU = trend_geneName(x_norm, gene, upregulation_table, genes, 'upregulation', corr{i}, 'symbol'); %Upregulation
        elseif nargin > 2
            resultsD = trend_geneName(x_norm, gene, downregulation_table, genes, 'downregulation', corr{i}, 'id'); %Downregulation
            resultsU = trend_geneName(x_norm, gene, upregulation_table, genes, 'upregulation', corr{i}, 'id'); %Upregulation
        end
%         resultsD = trend_geneName(x_norm, gene, downregulation_table, genes, 'downregulation', corr{i}, 'id'); %Downregulation
%         resultsU = trend_geneName(x_norm, gene, upregulation_table, genes, 'upregulation', corr{i}, 'id'); %Upregulation


        file = regexp(exp_data,'(\w*(?=.mat))', 'match');
        result = [['./result_', file{1}, '/'],file{1},'_',norm{j},'_',corr{i}, '.mat'];
        save(result, 'resultsU', 'resultsD');
        
        % Sensitivity, F1 score, AUC
        % Upregulation is positive class
        sensitivity_D = sum([resultsD{2:end,5}])/length([resultsD{2:end,5}]); %fraction of catched downregulation
        sensitivity_U = sum([resultsU{2:end,5}])/length([resultsU{2:end,5}]); %fraction of catched upregulation
        precision_D = sum([resultsD{2:end,5}])/(sum([resultsD{2:end,5}])+sum([resultsU{2:end,5}] == 0)); %Precision of downregulation
        precision_U = sum([resultsU{2:end,5}])/(sum([resultsU{2:end,5}])+sum([resultsD{2:end,5}] == 0)); %Precision of downregulation
        F1_score_D = 2*(precision_D*sensitivity_D/(sensitivity_D+precision_D));
        F1_score_U = 2*(precision_U*sensitivity_U/(sensitivity_U+precision_U));
        auc = AUC(sensitivity_D, sensitivity_U);
        data_name = [file{1},'_',norm{j},'_',corr{i}];
        
        statistic_result{idx, 1} = data_name;
        statistic_result{idx, 2} = sensitivity_U;
        statistic_result{idx, 3} = sensitivity_D;
        statistic_result{idx, 4} = F1_score_U;
        statistic_result{idx, 5} = F1_score_D;
        statistic_result{idx, 6} = auc;
        statistic_result{idx, 7} = precision_U;
        statistic_result{idx, 8} = precision_D;
        
        idx = idx + 1;
        
    end
end

save(['./result_', file{1}, '/', 'result.mat'], 'statistic_result');
end

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

%% Evaluation correlation-based methods (for gene symbol data)
function results = trend_geneName(x,gene_name,table_regulation,genes, type_regulation, type_corr, type_gene)
% results=cell(size(table_regulation,1)+1,6);
results(1,:)={'TF','Gene','Correlation','Type of regulation (exp)','Catched?','Deviation'};
if strcmp(type_regulation,'downregulation')
    trend = -1;
elseif strcmp(type_regulation,'upregulation')
    trend = 1;
else
    trend = 0;
end


if strcmp(type_gene,'symbol')
    j = 1; %To manage the row of the table in the loop
    for i=1:length(table_regulation(:,1))
        logical_TF = strcmp(gene_name,table_regulation{i,1});
        logical_gene = strcmp(gene_name,table_regulation{i,2});
        
        if (sum(logical_TF)== 0) || (sum(logical_gene)== 0) || strcmp(table_regulation{i,1},table_regulation{i,2})
            continue
        else
            results{j+1,1}=table_regulation{i,1}; %TF name
            results{j+1,2}=table_regulation{i,2}; %Gene name
            results{j+1,4}=type_regulation;
            % Compute correlation
            idx_TF = find(logical_TF);
            idx_gene = find(logical_gene);
            try
                c = corr_v2(x(idx_TF,:)', x(idx_gene,:)', type_corr); %Correlation
                results{j+1,3} = c;
                
            catch
                results{j+1,3}=NaN;
                results{j+1,7}=NaN;
            end

            if strcmp(type_regulation,'upregulation')||  strcmp(type_regulation,'downregulation')
                results{j+1,5} = sign(results{j+1,3}) == trend;
                results{j+1,6} = abs(results{j+1,3}-trend);
            elseif strcmp(type_regulation,'not known')
                results{j+1,5}='-';
                results{j+1,6}='-';
            end

            j = j + 1;
        end
        fprintf(type_corr);
        disp(i/length(table_regulation(:,1)));
    end
    
elseif strcmp(type_gene,'id')
    j = 1; %To manage the row of the table in the loop
    for i=1:length(table_regulation(:,1))
        logical_TF = strcmp(gene_name, genes.(table_regulation{i,1}).entrez);
        logical_gene = strcmp(gene_name,genes.(table_regulation{i,2}).entrez);
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
            results{j+1,3} = c;

            if strcmp(type_regulation,'upregulation')||  strcmp(type_regulation,'downregulation')
                results{j+1,5} = sign(results{j+1,3}) == trend;
                results{j+1,6} = abs(results{j+1,3}-trend);
            elseif strcmp(type_regulation,'not known')
                results{j+1,5}='-';
                results{j+1,6}='-';
            end

            j = j + 1;
        end
        fprintf(type_corr);
        disp(i/length(table_regulation(:,1)));
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

%% SDC
function [sign, dependency] = SDC(x,y,type)
[d_x, d_y] = increment(x,y);


switch type
    case 'Pearson'
        sign = corr(d_x,d_y);
        dependency = distcor(x,y,'Pearson');
    case 'Spearman'
        sign = corr(d_x,d_y,'Type','Spearman');
        dependency = distcor(x,y,'Spearman');
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

function dcor = distcor(x,y,type)
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
if strcmp(type,'Spearman')
    a = rank_matrix(a); %Rank matrix
end

mcol = mean(a);
mrow = mean(a,2);
ajbar = ones(size(mrow))*mcol;
akbar = mrow*ones(size(mcol));
abar = mean(mean(a))*ones(size(a));
A = a - ajbar - akbar + abar;


b = pdist2(y,y);
if strcmp(type,'Spearman')
    b = rank_matrix(b); %Rank matrix
end

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

function rank_x = rank_matrix(x)
r = sort(x(:));
rank_x = [];
for i = 1:size(x,1)
    for j = 1:size(x,2)
        idx = find(r == x(i,j));
        if length(idx) == 1
            rank_x(i,j) = idx;
        else
            rank_x(i,j) = sum(idx)/length(idx);
        end
    end
end
end