function algorithm_evaluation(data_xlsx, exp_data, gene_name_type)
%parpool(16);

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
if nargin == 3
    if gene_name_type == 'gene ID'
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
% Do the ZSCORE, LOG and LOG+ZSCORE
norm = {'without', 'ZSCORE', 'LOG'};
% corr = {'Pearson', 'Spearman', 'Distance','MI1_noSign', 'MI2_noSign', 'MI1_sign', 'MI2_sign', 'coagreement', 'coagreement_linear', 'coagreement_sigmoid'};
corr = {'coagreement_linear_new', 'coagreement_sigmoid_new'};
% Create result table
statistic_result(1,:) = {'Data type', 'Sensitivity', 'Specificity', 'F1 score for upregulation', 'F1 score for downregulation','AUC', 'Precision-positive', 'Precision-negtive'};
idx = 2; %index of row number

%% Microarray dataset
load (exp_data)
load gene_ID.mat

for i=1:length(corr)
    for j=1:length(norm)
        x_norm = Normalization(x, norm{j}); %Normalization
        %The function trend_entrezGeneID is using the entrez gene id
        resultsD = trend_geneName(x_norm, gene, downregulation_table, genes, 'downregulation', corr{i}); %Downregulation
        resultsU = trend_geneName(x_norm, gene, upregulation_table, genes, 'upregulation', corr{i}); %Upregulation
        resultsD_threshold = trend_threshold(resultsD,0.05); %The correlation value should be higher than 0.5 in this result
        resultsU_threshold = trend_threshold(resultsU,0.05); %The correlation value should be higher than 0.5 in this result
        file = regexp(exp_data,'(\w*(?=.mat))', 'match');
        result = [['./result_', file{1}, '/'], norm{j},'_', corr{i}, '.mat'];
        save(result, 'resultsU', 'resultsD', 'resultsD_threshold', 'resultsU_threshold');
        
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
        
        % Threshold
        sensitivity_threshold = sum([resultsU_threshold{2:end,5}])/length([resultsU_threshold{2:end,5}]);
        specificity_threshold = sum([resultsD_threshold{2:end,5}])/length([resultsD_threshold{2:end,5}]);
        precision_U = sum([resultsU_threshold{2:end,5}])/(sum([resultsU_threshold{2:end,5}])+sum([resultsD_threshold{2:end,5}] == 0));
        precision_D = sum([resultsD_threshold{2:end,5}])/(sum([resultsD_threshold{2:end,5}])+sum([resultsU_threshold{2:end,5}] == 0));
        F1_score_U = 2*(precision_U*sensitivity_threshold/(precision_U+sensitivity_threshold));
        F1_score_D = 2*(precision_D*specificity_threshold/(precision_D+specificity_threshold));
        auc_threshold = AUC(sensitivity_threshold, specificity_threshold);
        data_name_threshold = [file{1},'_0.05',norm{j},'_',corr{i}];
        
        statistic_result{idx, 1} = data_name_threshold;
        statistic_result{idx, 2} = sensitivity_threshold;
        statistic_result{idx, 3} = specificity_threshold;
        statistic_result{idx, 4} = F1_score_U;
        statistic_result{idx, 5} = F1_score_D;
        statistic_result{idx, 6} = auc_threshold; 
        statistic_result{idx, 7} = precision_U;
        statistic_result{idx, 8} = precision_D;
        
        idx = idx + 1;
    end
end

% save(['./result_', file{1}, '/', 'result.mat'], 'statistic_result');
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


if nargin ~= 3
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
            
            [c, p] = corr_v2(x(idx_TF,:)', x(idx_gene,:)', 1000, type_corr); %Correlation
            results{j+1,3}=c;
            results{j+1,7}=p;
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
        fprintf(type_corr);
        disp(i/length(table_regulation(:,1)));
    end
    
elseif nargin == 3
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
            [c, p] = corr_v2(x(idx_TF,:)', x(idx_gene,:)', 1000, type_corr); %Correlation
            results{j+1,3}=c;
            results{j+1,7}=p;
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
function [c,p] = corr_v2(x,y,n,type)

switch type
    case 'Pearson'
        c = corr(x, y);
        if isnan(c)
            p = NaN;
            c = NaN;
        else
            if nargin == 4
                m = zeros(1, n+1);
                m(1) = c;
                parfor i=2:n+1
                    m(i) = corr(x(randperm(length(x))), y(randperm(length(y))));
                end
                positionC = find(sort(m, 'ascend') == c);
                if c >= 0
                    p = 1 - ((positionC(1)-1)/(n+1));
                else
                    p = positionC(length(positionC))/(n+1);
                end
            else
                p = nan;
            end
        end
    case 'Spearman'
        c = corr(x, y , 'Type', 'Spearman');
        if isnan(c)
            p = NaN;
            c = NaN;
        else
            if nargin == 4
                m = zeros(1,n+1);
                m(1) = c;
                parfor i=2:n+1
                    m(i) = corr(x(randperm(length(x))), y(randperm(length(y))), 'Type', 'Spearman');
                end
                positionC = find(sort(m, 'ascend') == c);
                if c >= 0
                    p = 1 - ((positionC(1)-1)/(n+1));
                else
                    p = positionC(length(positionC))/(n+1);
                end
            else
                p = nan;
            end
        end
    case 'Distance'
        dist_corr = distcor(x, y);
        if isnan(dist_corr)
            p = NaN;
            c = NaN;
        else
            c = signed_dist(x, y);
            
            if nargin == 4
                m = zeros(1,n+1);
                m(1) = c;
                parfor i=2:n+1
                    x_perm = x(randperm(length(x)));
                    y_perm = y(randperm(length(y)));
                    m(i) = signed_dist(x_perm, y_perm);
                end
                positionC = find(sort(m, 'ascend') == c);
                if c >= 0
                    p = 1 - ((positionC(1)-1)/(n+1)); % If the value c more than one, choose the first position for worse case.
                else
                    p = positionC(length(positionC))/(n+1); % If the value c more than one, choose the last position for worse case.
                end
            else
                p = nan;
            end
        end
    case 'MI1_noSign' % MI1 = I/(1+I) for nalmalization
        if sum(x) == 0 || sum(y) == 0
            c = NaN;
        else
            [~, c, ~] = mi_gg(x, y);
        end
        if isnan(c)
            c = NaN;
            p = NaN;
        else
            if nargin == 4
                m = zeros(1,n+1);
                m(1) = c;
                parfor i=2:n+1
                    [~, I, ~] = mi_gg(x(randperm(length(x))), y(randperm(length(y))));
                    if I > 0
                        m(i) = I
                    else
                        m(i) = 0;
                    end
                end
                positionC = find(sort(m, 'ascend') == c);
                
                p = 1 - ((positionC(1)-1)/(n+1));
            else
                p = nan;
            end
        end
    case 'MI2_noSign' % MI2 = -I/H(x, y) for normalization
        if sum(x) == 0 || sum(y) == 0
            c = NaN;
        else
            [I, ~, NI] = mi_gg(x, y);
            if I > 0
                c = abs(NI);
            else
                c = 0;
            end
        end
        if isnan(c)
            c = NaN;
            p = NaN;
        else
            if nargin == 4
                m = zeros(1,n+1);
                m(1) = c;
                parfor i=2:n+1
                    [~, ~, I]= mi_gg(x(randperm(length(x))), y(randperm(length(y))));
                    m(i) = I;
                end
                positionC = find(sort(m, 'ascend') == c);
                
                p = 1 - ((positionC(1)-1)/(n+1));
            else
                p = nan;
            end
        end
    case 'MI1_sign' % MI1 = I/(1+I) for nalmalization
        if sum(x) == 0 || sum(y) == 0
            c = NaN;
        else
            [~, I, ~] = mi_gg(x, y);
            if I > 0
                c = sign(coagreement_nn_18July(x, y, 'linear'))*I;
            else
                c = 0;
            end
        end
        if isnan(c)
            c = NaN;
            p = NaN;
        else
            if nargin == 4
                m = zeros(1,n+1);
                m(1) = c;
                parfor i=2:n+1
                    x_perm = x(randperm(length(x)));
                    y_perm = y(randperm(length(y)));
                    [~, I, ~] = mi_gg(x_perm, y_perm);
                    if I > 0
                        m(i) = sign(coagreement_nn_18July(x_perm, y_perm, 'linear'))*I;
                    else
                        m(i) = 0;
                    end
                end
                positionC = find(sort(m, 'ascend') == c);
                
                if c >= 0
                    p = 1 - ((positionC(1)-1)/(n+1)); % If the value c more than one, choose the first position for worse case.
                else
                    p = positionC(length(positionC))/(n+1); % If the value c more than one, choose the last position for worse case.
                end
            else
                p = nan;
            end
        end
    case 'MI2_sign'
        if sum(x) == 0 || sum(y) == 0
            c = NaN;
        else
            [I, ~, NI]= mi_gg(x, y);
            if I > 0
                c = sign(coagreement_nn_18July(x, y, 'linear'))*NI;
            else
                c = 0;
            end
        end
        if isnan(c)
            c = NaN;
            p = NaN;
        else
            if nargin == 4
                m = zeros(1,n+1);
                m(1) = c;
                parfor i=2:n+1
                    x_perm = x(randperm(length(x)));
                    y_perm = y(randperm(length(y)));
                    [I, ~,NI]= mi_gg(x_perm, y_perm);
                    if I > 0
                        m(i) = sign(coagreement_nn_18July(x_perm, y_perm, 'linear'))*NI;
                    else
                        m(i) = 0;
                    end
                end
                positionC = find(sort(m, 'ascend') == c);
                
                if c >= 0
                    p = 1 - ((positionC(1)-1)/(n+1)); % If the value c more than one, choose the first position for worse case.
                else
                    p = positionC(length(positionC))/(n+1); % If the value c more than one, choose the last position for worse case.
                end
            else
                p = nan;
            end
        end
    case 'coagreement'
        if sum(x) == 0 || sum(y) == 0
            c = NaN;
            p = NaN;
        else
            if nargin == 4
                [c, p, ~] = coagreement(x, y, n);
            else
                [c, p, ~] = coagreement(x, y);
            end
        end
    case 'coagreement_linear'
        if sum(x) == 0 || sum(y) == 0
            c = NaN;
            p = NaN;
        else
            if nargin == 4
                [c, ~, p, ~, ~, ~] = coagreement_nn_18July(x, y, 'linear', n);
            else
                [c, ~, ~, ~, ~, ~] = coagreement_nn_18July(x, y, 'linear');
            end
        end
    case 'coagreement_sigmoid'
        if sum(x) == 0 || sum(y) == 0
            c = NaN;
            p = NaN;
        else
            if nargin == 4
                [c, ~, p, ~, ~, ~] = coagreement_nn_18July(x, y, 'sigmoid', n);
            else
                [c, ~, ~, ~, ~, ~] = coagreement_nn_18July(x, y, 'sigmoid');
            end
        end
    case 'coagreement_linear_new'
        if sum(x) == 0 || sum(y) == 0
            c = NaN;
            p = NaN;
        else
            if nargin == 4
                [c, ~, p, ~, ~, ~] = coagreement_nn_20Aug(x, y, 'linear', n);
            else
                [c, ~, ~, ~, ~, ~] = coagreement_nn_20Aug(x, y, 'linear');
            end
        end
    case 'coagreement_sigmoid_new'
        if sum(x) == 0 || sum(y) == 0
            c = NaN;
            p = NaN;
        else
            if nargin == 4
                [c, ~, p, ~, ~, ~] = coagreement_nn_20Aug(x, y, 'sigmoid', n);
            else
                [c, ~, ~, ~, ~, ~] = coagreement_nn_20Aug(x, y, 'sigmoid');
            end
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
