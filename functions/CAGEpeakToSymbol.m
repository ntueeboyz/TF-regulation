function CAGEpeakToSymbol (cagePeaks)
idx_all = 1:length(cagePeaks);
[ia1, ia2] = unique(round(idx_all/length(idx_all), 2)*100);
geneSymbol_cage=[];
UCSCstart = 'https://genome-asia.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=';
UCSCend = '&hgsid=472748652_ZGDKaSiJ9NPyp6DznGchCkiVtQei';
    
parfor i=1:length(cagePeaks)
    peaksPosition = regexp(cagePeaks{i}, '\w*', 'match');
    chrLocus = strcat(peaksPosition{1}, peaksPosition{2}, '-', peaksPosition{3});
    position = strcat(peaksPosition{1}, '%3A', peaksPosition{2}, '-', peaksPosition{3});
    esearchURL = [UCSCstart, position, UCSCend];
        
    numberTrials = 1;
    totalTrials = 3;
    while numberTrials > 0 && numberTrials < totalTrials
        try
            searchReport = webread(esearchURL, weboptions('Timeout', 30));
            symbol = regexp(searchReport, '((?<=gtexGene&i=)\w*(?='' class))', 'match');
            if isempty(symbol)
                geneSymbol_cage{i}=NaN;
            else
                symbol = unique(symbol);
                geneSymbol_cage{i} = symbol;
            end
            numberTrials = 0;
        catch 
            fprintf('\nOops! There is an error with peak: %s (TRIAL %d/%d)\n', chrLocus, numberTrials, totalTrials);
            numberTrials = numberTrials + 1;
        end
    end
    disp(i);
end
save('geneSymbol_cage.mat','geneSymbol_cage');
end
