data = {'THP1_single_cell', 'insilico_size100_1', 'insilico_size100_2'};
norm = {'','without','log','zscore'};
noise = {'','snr5','snr10','snr15'};
method = {'', 'TSNI','SINCERITIES', 'GRNVBEM'};

idx =2;
result_data(1, :) = {'data','sens','spe','f1_p','f1_n','auc', 'pre_p', 'pre_n'};

for i = 1:length(data)
    for j = 1:length(norm)
        for m = 1:length(noise)
            for k = 1:length(method)
                try
                    if k == 1 && j == 1 && m == 1
                        mat_name = [data{i},'.mat'];
                        load(mat_name)
                        
                        for l = 2:size(statistic_result,1)
                            result_data{idx , 1} = statistic_result{l,1};
                            result_data{idx , 2} = statistic_result{l,2};
                            result_data{idx , 3} = statistic_result{l,3};
                            result_data{idx , 4} = statistic_result{l,4};
                            result_data{idx , 5} = statistic_result{l,5};
                            result_data{idx , 6} = statistic_result{l,6};
                            result_data{idx , 7} = statistic_result{l,7};
                            result_data{idx , 8} = statistic_result{l,8};
                            idx = idx + 1;
                        end
                    elseif k == 1 && j == 1 && m ~= 1
                        mat_name = [data{i},'_',noise{m},'.mat'];
                        load(mat_name)
                        
                        for l = 2:size(statistic_result,1)
                            result_data{idx , 1} = statistic_result{l,1};
                            result_data{idx , 2} = statistic_result{l,2};
                            result_data{idx , 3} = statistic_result{l,3};
                            result_data{idx , 4} = statistic_result{l,4};
                            result_data{idx , 5} = statistic_result{l,5};
                            result_data{idx , 6} = statistic_result{l,6};
                            result_data{idx , 7} = statistic_result{l,7};
                            result_data{idx , 8} = statistic_result{l,8};
                            idx = idx + 1;
                        end
                        
                    elseif k == 1
                        continue
                    elseif j == 1 && k ~= 1
                        continue
                    elseif m == 1 && k ~= 1
                        mat_name = [data{i}, '_', norm{j},'_', method{k},'.mat'];
                        data_name = [data{i},'_', norm{j}, '_', method{k}];
                        load(mat_name)
                        
                        result_data{idx , 1} = data_name;
                        result_data{idx , 2} = statistic_result{2,1};
                        result_data{idx , 3} = statistic_result{2,2};
                        result_data{idx , 4} = statistic_result{2,3};
                        result_data{idx , 5} = statistic_result{2,4};
                        result_data{idx , 6} = statistic_result{2,5};
                        result_data{idx , 7} = statistic_result{2,6};
                        result_data{idx , 8} = statistic_result{2,7};
                
                        idx = idx + 1;
                        
                        
                    elseif k ~= 1
                        mat_name = [data{i},'_', norm{j}, '_',noise{m},'_', method{k},'.mat'];
                        data_name = [data{i},'_',noise{m},'_',norm{j},'_', method{k}];
                        load(mat_name)
                        
                        result_data{idx , 1} = data_name;
                        result_data{idx , 2} = statistic_result{2,1};
                        result_data{idx , 3} = statistic_result{2,2};
                        result_data{idx , 4} = statistic_result{2,3};
                        result_data{idx , 5} = statistic_result{2,4};
                        result_data{idx , 6} = statistic_result{2,5};
                        result_data{idx , 7} = statistic_result{2,6};
                        result_data{idx , 8} = statistic_result{2,7};
                        idx = idx + 1;
                        
                    end
                catch
                    continue
                end
            end
        end
    end
end

save('result.mat', 'result_data')
