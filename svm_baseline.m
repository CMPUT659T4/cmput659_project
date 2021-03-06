addpath('lasso','-end');
%load('FBIRN/finaldata_AO/features/fBIRN_AudOdd_allsites_0003_degrees.mat')
%load('FBIRN/finaldata_AO/features/fBIRN_AudOdd_allsites_0003_log_degrees.mat')
% load('FBIRN/finaldata_AO/features/fBIRN_AudOdd_allsites_0003_log_degrees_Tlms_MFG_SFG_.mat')
load('concat_roi_avg.mat')
load('holdout.mat')

% for s = 1:size(labels, 1)
%     for r = 1:size(data, 1)
%         series = squeeze(data(r, :, s));
%         data(r, :, s) = (series - mean(series)) / std(series);
%     end
% end
% 
% data([81, 82, 128, 184, 247, 248, 249, 250], :, :) = [];
% 
% new_data = zeros(size(data, 1), size(data, 1), size(data, 3));
% for i = 1:size(data, 3)
%     d = squeeze(data(:, :, i));
%     [~, precision, ~, ~, ~] = GraphicalLasso(d', 0.01);
%     new_data(:, :, i) = precision;
%     i
% end
% data = new_data;
% save(['precision_data.mat'], 'data')
load('precision_data.mat')
holdout_set = TestInd(labels, 19);
holdout_labels = labels(holdout_set);
train_set = [1:380]';
train_set(holdout_set) = [];
train_labels = labels(train_set);

holdout = data(:, :, holdout_set);
data = data(:, :, train_set);
s_inds = train_labels(:) == 1;
h_inds = train_labels(:) == -1;
s_data = data(:, :, s_inds);
h_data = data(:, :, h_inds);

% a = zeros(5, 1);
% for iteration = 1:5
%     s_kfold = crossvalind('Kfold', size(s_data, 3) / 4, 10);
%     s_kfold = [s_kfold, s_kfold, s_kfold, s_kfold]';
%     s_kfold = s_kfold(:);
% 
%     h_kfold = crossvalind('Kfold', size(h_data, 3) / 4, 10);
%     h_kfold = [h_kfold, h_kfold, h_kfold, h_kfold]';
%     h_kfold = h_kfold(:);
% 
%     s_test_inds = logical((s_kfold == 1) + (s_kfold == 2));% + (s_kfold == 3));
%     h_test_inds = logical((h_kfold == 1) + (h_kfold == 2));% + (h_kfold == 3));
%     s_train_inds = ~s_test_inds;
%     h_train_inds = ~h_test_inds;
% 
%     s_train = s_data(:, :, s_train_inds);
%     s_train = reshape(s_train, [size(s_train, 1) * size(s_train, 2), size(s_train, 3)])';
%     h_train = h_data(:, :, h_train_inds);
%     h_train = reshape(h_train, [size(h_train, 1) * size(h_train, 2), size(h_train, 3)])';
% 
%     s_test = s_data(:, :, s_test_inds);
%     h_test = h_data(:, :, h_test_inds);
% 
%     s_test = reshape(s_test, [size(s_test, 1) * size(s_test, 2), size(s_test, 3)])';
%     h_test = reshape(h_test, [size(h_test, 1) * size(h_test, 2), size(h_test, 3)])';
%     SVMStruct = svmtrain([s_train; h_train], [ones(size(s_train, 1), 1); -1 * ones(size(h_train, 1), 1)]);
%     predictions = svmclassify(SVMStruct, [s_test; h_test]);
%     t_labels = [ones(size(s_test, 1), 1); -1 * ones(size(h_test, 1), 1)];
%     %[s_covar, s_precision, ~, ~, ~] = GraphicalLasso(s_train, 0.01);
%     %[h_covar, h_precision, ~, ~, ~] = GraphicalLasso(h_train, 0.01);
%     
%     %prior = [size(s_train, 1); size(h_train, 1)] / (size(s_train, 1) + size(h_train, 1));
%     %l = unique(labels);   n_l = size(l,1);
%     %conditional(1) = struct('C', s_precision, 'mu', mean(s_train));
%     %conditional(2) = struct('C', h_precision, 'mu', mean(h_train));
%     %model = struct('nvars', size(s_train, 2), 'nlabels', n_l, 'labels', l, ...
%     %               'class_cond', conditional, 'class_prior', prior, ...
%     %               'lambdas', 0);
% 
%     %test = cat(3, s_test, h_test);
%     %t_labels = [ones(size(s_test, 3), 1); -1 * ones(size(h_test, 3), 1)];
%     %predictions = zeros(size(test, 3), 1);
%     %for i = 1:size(test, 3)
%     %    s_ll = gaussianFit(s_precision, squeeze(test(:, :, i))', 0);
%     %    h_ll = gaussianFit(h_precision, squeeze(test(:, :, i))', 0);
%     %    %[y,pyx] = MRFC_predict(squeeze(test(:, :, i))', model);
%     %    %s_ll = sum(y == 1);
%     %    %h_ll = sum(y == -1);
%     %    if s_ll > h_ll
%     %        predictions(i) = 1;
%     %    else
%     %        predictions(i) = -1;
%     %    end
%     %end
% 
%     a(iteration) = sum(predictions == t_labels) / length(t_labels);
% end
% min(a)
% max(a)
% mean(a)

data = reshape(data, [size(data, 1) * size(data, 2), size(data, 3)])';
holdout = reshape(holdout, [size(holdout, 1) * size(holdout, 2), size(holdout, 3)])';
SVMStruct = svmtrain(data, train_labels);
predictions = svmclassify(SVMStruct, holdout);
sum(predictions == holdout_labels) / length(holdout_labels)
