# TF-regulation
## Master thesis: Correlation-based detection of activation/repression transcription factor regulation from gene expression

## The folders quick guide
1. data: The original data used in this project, including Microarray, RNA-seq, CAGE, THP1 and simulated data.
2. Figures: All the figures from the thesis were included.
3. functions: The functions used in the simulation. However, the funcions are all included in the main script (MAIN.m and algorithm_evaluation.m) and no need to add.
4. Network based method: The published reverse engineering algorithm. SINCERITIES, GRNVBEM and TSNI.
5. Noise data: The THP1 and simulated data with SNR5, SNR10 and SNR15 gaussian noise. Here includes the result of the statistical measures as well.
6. result: The results of different data.

## Evaluation script
1. MAIN.m: The main script to evaluate the steady-state expression data (Microarray, RNA-seq and CAGE).
2. algorithm_evaluation.m: The function to evaluate the different correlations (Pearson, Spearman and signed distance correlation) with gold standard and expression data input.
3. evaluation_method: The function to evaluate the output result of network based method with gold standard input.

## Plot script
#### \result_steady_state_expression\ (Steady state expression part)
1. best_AUC.m: To plot "The best performance of correlations across different technologies"
2. method_technology_plot.m: To plot "Performance of technologies across different correlations"
3. normalization_technology_plot.m: To plot "Performance of technologies across different normalizations"
4. technology_method_plot.m: To plot "Performance of correlations across different technologies"
5. technology_normalization_plot.m: To plot "Performance of normalizations across different technologies"
#### \result_time_series_expression\ (Time series expression part)
6. best_AUC.m: To plot "The best performance of methods across different datasets"
7. method_technology_plot.m: To plot "Performance of datasets across different methods"
8. normalization_technology_plot.m: To plot "Performance of datasets across different normalizations"
9. technology_method_plot.m: To plot "Performance of methods across different technologies"
10. technology_normalization_plot.m: To plot "Performance of normalizations across different datasets"
#### \Noise data\SNR data\
11. noise_plot2.m: To plot "Noise data"
