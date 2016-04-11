function [predict]=region_degrees_mario()
addpath('lasso','-end');
%load('FBIRN/finaldata_AO/features/fBIRN_AudOdd_allsites_0003_degrees.mat')
%load('FBIRN/finaldata_AO/features/fBIRN_AudOdd_allsites_0003_log_degrees.mat')
load('FBIRN/finaldata_AO/features/fBIRN_AudOdd_allsites_0003_log_degrees_Tlms_MFG_SFG_.mat')
load('concat_roi_avg.mat')
load('holdout.mat')

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
data = new_data;
%save(['region_degrees.mat'], 'data', 'labels')

holdout = data(:, holdout_set);
data = data(:, train_set);
s_inds = train_labels(:) == 1;
h_inds = train_labels(:) == -1;
s_data = data(:, s_inds);
h_data = data(:, h_inds);

% method='varsel_mrf';
% params = [0.7; 0.1; 0.07; 0.01; 0.007; 0.001];
% param_acc = zeros(length(params), 1);
% for p = 1:length(params)
%     p
%     a = zeros(5, 1);
%     for iteration = 1:5
%         s_kfold = crossvalind('Kfold', size(s_data, 2) / 4, 10);
%         s_kfold = [s_kfold, s_kfold, s_kfold, s_kfold]';
%         s_kfold = s_kfold(:);
% 
%         h_kfold = crossvalind('Kfold', size(h_data, 2) / 4, 10);
%         h_kfold = [h_kfold, h_kfold, h_kfold, h_kfold]';
%         h_kfold = h_kfold(:);
% 
%         s_test_inds = logical((s_kfold == 1) + (s_kfold == 2));
%         h_test_inds = logical((h_kfold == 1) + (h_kfold == 2));
%         s_train_inds = ~s_test_inds;
%         h_train_inds = ~h_test_inds;
% 
%         s_train = s_data(:, s_train_inds)';
%         h_train = h_data(:, h_train_inds)';
% 
%         s_test = s_data(:, s_test_inds)';
%         h_test = h_data(:, h_test_inds)';
% 
%         %[s_covar, s_precision, ~, ~, ~] = GraphicalLasso(s_train, params(p));
%         %[h_covar, h_precision, ~, ~, ~] = GraphicalLasso(h_train, params(p));
%         % save(['post_glasso.mat'], 's_inds', 'h_inds', 's_data', 'h_data', ...
%         %      's_kfold', 'h_kfold', 's_test_inds', 'h_test_inds', ...
%         %      's_train_inds', 'h_train_inds', 's_train', 'h_train', ...
%         %      's_test', 'h_test', 's_covar', 'h_covar', 's_precision', 'h_precision')
% 
%         %prior = [size(s_train, 1); size(h_train, 1)] / (size(s_train, 1) + size(h_train, 1));
%         %l = unique(labels);   n_l = size(l,1);
%         %conditional(1) = struct('C', s_precision, 'mu', mean(s_train));
%         %conditional(2) = struct('C', h_precision, 'mu', mean(h_train));
%         %model = struct('nvars', size(s_train, 2), 'nlabels', n_l, 'labels', l, ...
%         %               'class_cond', conditional, 'class_prior', prior, ...
%         %               'lambdas', 0);
%                    
%         model = MRFC_learn([s_train; h_train], [ones(size(s_train, 1), 1); -1 * ones(size(h_train, 1), 1)], method, params(p));
% 
%         test = [s_test; h_test];
%         t_labels = [ones(size(s_test, 1), 1); -1 * ones(size(h_test, 1), 1)];
%         predictions = zeros(size(test, 1), 1);
%         for i = 1:size(test, 1)
%             %nll(1) = gaussianFit(s_precision, test(i, :), 0);
%             %nll(2) = gaussianFit(h_precision, test(i, :), 0);
%             [y,pyx] = MRFC_predict(test(i, :), model);
%             nll = pyx;
%             if nll(1) > nll(2)
%                 predictions(i) = 1;
%             else
%                 predictions(i) = -1;
%             end
%         end
% 
%         a(iteration) = sum(predictions == t_labels) / length(t_labels);
%     end
%     param_acc(p) = mean(a);
%     [~, index] = max(a);
% end
% 
% [~, p] = max(param_acc);
% param_acc'
% params'

s_data = s_data';
h_data = h_data';

holdout = holdout';

%[~, s_precision, ~, ~, ~] = GraphicalLasso(s_data, params(p));
%[~, h_precision, ~, ~, ~] = GraphicalLasso(h_data, params(p));

%prior = [size(s_data, 1); size(h_data, 1)] / (size(s_data, 1) + size(h_data, 1));
%l = unique(holdout_labels);   n_l = size(l,1);
%conditional(1) = struct('C', s_precision, 'mu', mean(s_data));
%conditional(2) = struct('C', h_precision, 'mu', mean(h_data));
%model = struct('nvars', size(s_train, 2), 'nlabels', n_l, 'labels', l, ...
%               'class_cond', conditional, 'class_prior', prior, ...
%               'lambdas', 0);
           
% model = MRFC_learn([s_data; h_data], [ones(size(s_data, 1), 1); -1 * ones(size(h_data, 1), 1)], method, params(p));
% 
% predictions = zeros(size(holdout, 1), 1);
% for i = 1:size(holdout, 1)
%     %nll(1) = gaussianFit(s_precision, holdout(i, :), 0);
%     %nll(2) = gaussianFit(h_precision, holdout(i, :), 0);
%     [y,pyx] = MRFC_predict(holdout(i, :), model);
%     nll = pyx;
%     if nll(1) > nll(2)
%         predictions(i) = 1;
%     else
%         predictions(i) = -1;
%     end
% end

SVMStruct = svmtrain([s_data; h_data], [ones(size(s_data, 1), 1); -1 * ones(size(h_data, 1), 1)]);
predictions = svmclassify(SVMStruct, holdout);
sum(predictions == holdout_labels) / length(holdout_labels)
predict=predictions;
end