function [predict]=sgmrf_glasso_mario_diff(holdout_set);
addpath('lasso','-end');
%load('FBIRN/finaldata_AO/features/fBIRN_AudOdd_allsites_0003_degrees.mat')
%load('FBIRN/finaldata_AO/features/fBIRN_AudOdd_allsites_0003_log_degrees.mat')
% load('FBIRN/finaldata_AO/features/fBIRN_AudOdd_allsites_0003_log_degrees_Tlms_MFG_SFG_.mat')

load('concat_roi_avg.mat')
%load('holdout.mat')
holdout_labels = labels(holdout_set);

inds = 1:length(labels);
train_set = setdiff(inds, holdout_set)';
train_labels = labels(train_set);

for s = 1:size(labels, 1)
    for r = 1:size(data, 1)
        series = squeeze(data(r, :, s));
        data(r, :, s) = (series - mean(series)) / std(series);
    end
end

data([81, 82, 128, 184, 247, 248, 249, 250], :, :) = [];

holdout = data(:, :, holdout_set);
data = data(:, :, train_set);
s_inds = train_labels(:) == 1;
h_inds = train_labels(:) == -1;
s_data = data(:, :, s_inds);
h_data = data(:, :, h_inds);

method='varsel_mrf';
params = [0.007];
p=1;
% param_acc = zeros(length(params), 1);
% for p = 1:length(params)
%     p
%     a = zeros(5, 1);
%     b = zeros(5, 1);
%     c = zeros(5, 1);
%     for iteration = 1:5
%         s_kfold = crossvalind('Kfold', size(s_data, 3) / 4, 10);
%         s_kfold = [s_kfold, s_kfold, s_kfold, s_kfold]';
%         s_kfold = s_kfold(:);
% 
%         h_kfold = crossvalind('Kfold', size(h_data, 3) / 4, 10);
%         h_kfold = [h_kfold, h_kfold, h_kfold, h_kfold]';
%         h_kfold = h_kfold(:);
% 
%         s_test_inds = logical((s_kfold == 1) + (s_kfold == 2));
%         h_test_inds = logical((h_kfold == 1) + (h_kfold == 2));
%         s_train_inds = ~s_test_inds;
%         h_train_inds = ~h_test_inds;
% 
%         s_train = s_data(:, :, s_train_inds);
%         s_train = reshape(s_train, [size(s_train, 1), size(s_train, 2) * size(s_train, 3)])';
%         h_train = h_data(:, :, h_train_inds);
%         h_train = reshape(h_train, [size(h_train, 1), size(h_train, 2) * size(h_train, 3)])';
% 
%         s_test = s_data(:, :, s_test_inds);
%         h_test = h_data(:, :, h_test_inds);
% 
%         [s_covar, s_precision, ~, ~, ~] = GraphicalLasso(s_train, params(p));
%         [h_covar, h_precision, ~, ~, ~] = GraphicalLasso(h_train, params(p));
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
%         %model = MRFC_learn([s_train; h_train], [ones(size(s_train, 1), 1); -1 * ones(size(h_train, 1), 1)], method, params(p));
% 
%         test = cat(3, s_test, h_test);
%         t_labels = [ones(size(s_test, 3), 1); -1 * ones(size(h_test, 3), 1)];
%         predictions = zeros(size(test, 3), 1);
%         for i = 1:size(test, 3)
%             nll(1) = gaussianFit(s_precision, squeeze(test(:, :, i))', 0);
%             nll(2) = gaussianFit(h_precision, squeeze(test(:, :, i))', 0);
%             %[y,pyx] = MRFC_predict(squeeze(test(:, :, i))', model);
%             %nll = sum(pyx);
%             if nll(1) > nll(2)
%                 predictions(i) = 1;
%             else
%                 predictions(i) = -1;
%             end
%         end
%         b(iteration)=sum(predictions>t_labels)/sum(t_labels==-1);
%          c(iteration)=sum(predictions<t_labels)/sum(t_labels==1);
%         t_labels;
%         a(iteration) = sum(predictions == t_labels) / length(t_labels);
%     end
%     param_acc(p) = mean(a);
%      param_false_neg(p) = mean(c);
%       param_false_pos(p) = mean(b);
% end
% 
% [~, p] = max(param_acc);
% % param_acc'
% % params'
% fprintf('Fasle Positives = %.2f, Missed Patients = %.2f,    Accuracy=%.2f \n', param_false_pos(p), param_false_neg(p),param_acc(p));

s_data = reshape(s_data, [size(s_data, 1), size(s_data, 2) * size(s_data, 3)])';
h_data = reshape(h_data, [size(h_data, 1), size(h_data, 2) * size(h_data, 3)])';
[~, s_precision, ~, ~, ~] = GraphicalLasso(s_data, params(p));
[~, h_precision, ~, ~, ~] = GraphicalLasso(h_data, params(p));

%prior = [size(s_data, 1); size(h_data, 1)] / (size(s_data, 1) + size(h_data, 1));
%l = unique(holdout_labels);   n_l = size(l,1);
%conditional(1) = struct('C', s_precision, 'mu', mean(s_data));
%conditional(2) = struct('C', h_precision, 'mu', mean(h_data));
%model = struct('nvars', size(s_train, 2), 'nlabels', n_l, 'labels', l, ...
%               'class_cond', conditional, 'class_prior', prior, ...
%               'lambdas', 0);
           
%model = MRFC_learn([s_data; h_data], [ones(size(s_data, 1), 1); -1 * ones(size(h_data, 1), 1)], method, params(p));

predictions = zeros(size(holdout, 3), 1);
for i = 1:size(holdout, 3)
    nll(1) = gaussianFit(s_precision, squeeze(holdout(:, :, i))', 0);
    nll(2) = gaussianFit(h_precision, squeeze(holdout(:, :, i))', 0);
    %[y,pyx] = MRFC_predict(squeeze(holdout(:, :, i))', model);
    %nll = sum(pyx);
    if nll(1) > nll(2)
        predictions(i) = 1;
    else
        predictions(i) = -1;
    end
end

d=sum(predictions == holdout_labels) / length(holdout_labels);
e=sum(predictions > holdout_labels) / sum(holdout_labels==-1);
f=sum(predictions < holdout_labels) / sum(holdout_labels==1);
fprintf('Fasle Positives = %.2f, Missed Patients = %.2f,    Accuracy=%.2f \n', e, f,d);
predict=predictions;
end