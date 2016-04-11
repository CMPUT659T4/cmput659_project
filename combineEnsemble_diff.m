function [accu accu1 accu2 accu3]=combineEnsemble_diff(limit)

load('FBIRN/finaldata_AO/features/fBIRN_AudOdd_allsites_0003_log_degrees.mat');
%limit=2;
accu=zeros(limit,1);
accu1=zeros(limit,1);
accu2=zeros(limit,1);
accu3=zeros(limit,1);

for j=1:limit
    IndicesToTest=TestInd(data,19);
    [pred1 label]=BestFeaturesPrediction_Mario_diff(IndicesToTest);
    pred2=sgmrf_glasso_mario_diff(IndicesToTest);
    pred3=region_degrees_mario_diff(IndicesToTest);
    pred=pred1+pred2+pred3;
    %=pred1+pred3;
    final_res=zeros(size(pred));
     accu1(j,1)=sum(pred1== label) / length(label);
     accu2(j,1)=sum(pred2 == label) / length(label);
     accu3(j,1)=sum(pred3 == label) / length(label);
    for i=1:length(pred)
        if pred(i)>0
            final_res(i)=1;
        else
            final_res(i)=-1;
        end
    end
    accu(j,1)=sum(final_res == label) / length(label);
end
% final_res=pred3;
% accu=sum(pred3 == label) / length(label);
fprintf('############################################\n');
fprintf('Program is repeated for %d iterations\n',limit);
fprintf('Accuracy for Best Feature Selection is %.3f with STD %.3f\n',mean(accu1),std(accu1,1));
fprintf('Accuracy for SGMRF GLASSO is %.3f with STD %.3f\n',mean(accu2),std(accu2,1));
fprintf('Accuracy for Region Degrees SVM is %.3f with STD %.3f\n',mean(accu3),std(accu3,1));
fprintf('############################################\n');
fprintf('Total Accuracy is %.3f with STD %.3f\n',mean(accu),std(accu,1));
fprintf('############################################\n');
end