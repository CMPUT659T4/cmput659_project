function [Accuracy,PErrVec,FErrVec,AnalyseMatrix]=Region_DegreesCrossValidation(limit,itterat)
AnalyseMatrix=[];
PErrVec=zeros(1,itterat);
FErrVec=zeros(1,itterat);
addpath('lasso','-end');
%load('FBIRN/finaldata_AO/features/fBIRN_AudOdd_allsites_0003_degrees.mat')
%load('FBIRN/finaldata_AO/features/fBIRN_AudOdd_allsites_0003_log_degrees.mat')
load('FBIRN/finaldata_AO/features/fBIRN_AudOdd_allsites_0003_log_degrees_Tlms_MFG_SFG_.mat')
load('concat_roi_avg.mat')

% for s = 1:size(labels, 1)
%     for r = 1:size(data, 1)
%         series = squeeze(data(r, :, s));
%         data(r, :, s) = (series - mean(series)) / std(series);
%     end
% end

data([81, 82, 128, 184, 247, 248, 249, 250], :, :) = [];

new_data = zeros(size(data, 1), size(data, 3));
for i = 1:size(data, 3)
    c = corr(squeeze(data(:, :, i))');
    c = c - diag(diag(c));
    new_data(:, i) = sum(c > 0.7);
end
data = log(0.1+new_data);
%save(['region_degrees.mat'], 'data')

s_inds = labels(:) == 1;
h_inds = labels(:) == -1;
s_data = data(:, s_inds);
h_data = data(:, h_inds);

lambdas=0.9;
method='varsel_mrf';
%a = zeros(10, 1);
for iteration = 1:itterat
    s_kfold = crossvalind('Kfold', size(s_data, 2) / 4, 10);
    s_kfold = [s_kfold, s_kfold, s_kfold, s_kfold]';
    s_kfold = s_kfold(:);

    h_kfold = crossvalind('Kfold', size(h_data, 2) / 4, 10);
    h_kfold = [h_kfold, h_kfold, h_kfold, h_kfold]';
    h_kfold = h_kfold(:);

    s_test_inds = logical((s_kfold == 1) + (s_kfold == 2));% + (s_kfold == 3));
    h_test_inds = logical((h_kfold == 1) + (h_kfold == 2));% + (h_kfold == 3));
    s_train_inds = ~s_test_inds;
    h_train_inds = ~h_test_inds;

    s_train = s_data(:, s_train_inds)';
    h_train = h_data(:, h_train_inds)';

    s_test = s_data(:, s_test_inds)';
    h_test = h_data(:, h_test_inds)';
    
    model = MRFC_learn([s_train; h_train], [ones(size(s_train, 1), 1); -1 * ones(size(h_train, 1), 1)], method, lambdas);
    predictions = MRFC_predict([s_test; h_test], model);


     t_labels = [ones(size(s_test, 1), 1); -1 * ones(size(h_test, 1), 1)];
    FPerr = size(find(predictions > 0 & t_labels < 0 ),1);
        FNerr = size(find(predictions < 0 & t_labels > 0 ),1);
    a(iteration) = sum(predictions == t_labels) / length(t_labels);
     HTot=length(find(t_labels<0));
        STot=length(find(t_labels>0));
    PErrVec(iteration)=FPerr/HTot;
    FErrVec(iteration)=FNerr/STot;
end
Accuracy=a;
min(a)
max(a)

fprintf('Running using degrees for each region\n',limit);
fprintf('Accuracy = %.2f with std=%.2f Minimum=%.2f Max=%.2f \n False positives=%.2f \n Missed Patients=%.2f\n',mean(a),std(a,1),min(a),max(a),mean(mean(PErrVec)),mean(mean(FErrVec)));


end