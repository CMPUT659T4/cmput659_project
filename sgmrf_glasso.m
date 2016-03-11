addpath('lasso','-end');
%load('FBIRN/finaldata_AO/features/fBIRN_AudOdd_allsites_0003_degrees.mat')
%load('FBIRN/finaldata_AO/features/fBIRN_AudOdd_allsites_0003_log_degrees.mat')
load('FBIRN/finaldata_AO/features/fBIRN_AudOdd_allsites_0003_log_degrees_Tlms_MFG_SFG_.mat')
s_inds = data(:, end) == 1;
h_inds = data(:, end) == -1;
s_data = data(s_inds, :);
h_data = data(h_inds, :);

s_kfold = crossvalind('Kfold', size(s_data, 1) / 4, 10);
s_kfold = [s_kfold, s_kfold, s_kfold, s_kfold]';
s_kfold = s_kfold(:);

h_kfold = crossvalind('Kfold', size(h_data, 1) / 4, 10);
h_kfold = [h_kfold, h_kfold, h_kfold, h_kfold]';
h_kfold = h_kfold(:);

s_test_inds = logical((s_kfold == 1) + (s_kfold == 2) + (s_kfold == 3));
h_test_inds = logical((h_kfold == 1) + (h_kfold == 2) + (h_kfold == 3));
s_train_inds = ~s_test_inds;
h_train_inds = ~h_test_inds;

s_train = s_data(s_train_inds, :);
h_train = h_data(h_train_inds, :);

s_test = s_data(s_test_inds, :);
h_test = h_data(h_test_inds, :);

[s_covar, s_precision, ~, ~, ~] = GraphicalLasso(s_train(:, 1:end-1), 0.1);
[h_covar, h_precision, ~, ~, ~] = GraphicalLasso(h_train(:, 1:end-1), 0.1);