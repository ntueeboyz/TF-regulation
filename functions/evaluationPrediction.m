function [acc,prec,sens,spec,auc_,Fscore,auc,aupr,mcc] = evaluationPrediction(modelPrediction,scores,validationLabels,classes)

for i = 1:length(classes)
    
    auxPredictedLabels = modelPrediction;
    auxValidationLabels = validationLabels;
    
    possClass = classes{i};
    negClass = 'negClass';
    
    auxPredictedLabels(~ismember(modelPrediction,possClass)) = {negClass};
    auxValidationLabels(~ismember(validationLabels,possClass)) = {negClass};
    
    tp = sum(ismember(auxPredictedLabels,possClass) & ismember(auxValidationLabels,possClass));
    tn = sum(ismember(auxPredictedLabels,negClass) & ismember(auxValidationLabels,negClass));
    fp = sum(ismember(auxPredictedLabels,possClass) & ismember(auxValidationLabels,negClass));
    fn = sum(ismember(auxPredictedLabels,negClass) & ismember(auxValidationLabels,possClass));
    
    acc(i) = (tp + tn) / (tp + tn + fp + fn);
    sens(i) = tp/(tp+fn); %recall or TPR
    spec(i) = tn/(tn+fp); %TNR
    prec(i) = tp/(tp+fp); %PPV
    
    fpr = fp/(fp+tn);
    
    x = [0;sens(i);1];
    y = [0;fpr;1];
    
    auc_(i) = trapz(y,x);
    if auc_(i) < 0.5
        auc_(i) = 1 - auc_(i);
    end
    
    Fscore(i) = (2*prec(i)*sens(i))/(prec(i)+sens(i));  
    
    if size(scores,2) > 1
        [auc(i),aupr(i)] = aupr_evaluation(auxValidationLabels,scores(:,i),possClass);
    else
        [auc(i),aupr(i)] = aupr_evaluation(auxValidationLabels,scores,possClass);
    end
    
    if (tp == 0 && fp == 0) || (tn == 0 && fn == 0)
        mcc(i) = 0;
    else
        mcc(i) = tp*tn - fp*fn / sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)) ;
    end
end