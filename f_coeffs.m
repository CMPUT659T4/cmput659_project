addpath('lasso','-end');
%load('FBIRN/finaldata_AO/features/fBIRN_AudOdd_allsites_0003_degrees.mat')
%load('FBIRN/finaldata_AO/features/fBIRN_AudOdd_allsites_0003_log_degrees.mat')
% load('FBIRN/finaldata_AO/features/fBIRN_AudOdd_allsites_0003_log_degrees_Tlms_MFG_SFG_.mat')
load('concat_roi_avg.mat')

% for s = 1:size(labels, 1)
%     for r = 1:size(data, 1)
%         series = squeeze(data(r, :, s));
%         data(r, :, s) = (series - mean(series)) / std(series);
%     end
% end

data([81, 82, 128, 184, 247, 248, 249, 250], :, :) = [];

new_data = zeros(size(data, 1), floor(size(data, 2) / 2), size(data, 3));
for i = 1:size(data, 3)
    d = squeeze(data(:, :, i));
    d = real(fft(d, [], 2));
    d = d(:, 2:ceil(size(d, 2) / 2));
    new_data(:, :, i) = d;
end
data = new_data;
%save(['coeffs_data.mat'], 'data')

s_inds = labels(:) == 1;
h_inds = labels(:) == -1;
s_data = data(:, :, s_inds);
h_data = data(:, :, h_inds);

s = reshape(s_data, [size(s_data, 1) * size(s_data, 2), size(s_data, 3)])';
h = reshape(h_data, [size(h_data, 1) * size(h_data, 2), size(h_data, 3)])';

% p = zeros(size(s, 2), 1);
% for i = 1:size(s, 2)
%     [~, p(i)] = ttest2(s(:, i)', h(:, i)');
% end
% [B, I] = sort(p);
% s = s(:, I(1:200));
% h = h(:, I(1:200));

a = zeros(10, 1);
for iteration = 1:10
    s_kfold = crossvalind('Kfold', size(s, 1) / 4, 10);
    s_kfold = [s_kfold, s_kfold, s_kfold, s_kfold]';
    s_kfold = s_kfold(:);

    h_kfold = crossvalind('Kfold', size(h, 1) / 4, 10);
    h_kfold = [h_kfold, h_kfold, h_kfold, h_kfold]';
    h_kfold = h_kfold(:);

    s_test_inds = logical((s_kfold == 1) + (s_kfold == 2));
    h_test_inds = logical((h_kfold == 1) + (h_kfold == 2));
    s_val_inds = logical((s_kfold == 3) + (s_kfold == 4));
    h_val_inds = logical((h_kfold == 3) + (h_kfold == 4));
    s_train_inds = ~(s_test_inds + s_val_inds);
    h_train_inds = ~(h_test_inds + h_val_inds);

    s_train = s(s_train_inds, :);
    h_train = h(h_train_inds, :);

    s_val = s(s_val_inds, :);
    h_val = h(h_val_inds, :);
    
    s_test = s(s_test_inds, :);
    h_test = h(h_test_inds, :);
    
%     p = zeros(size(s_val, 2), 1);
%     for i = 1:size(s_val, 2)
%         [~, p(i)] = ttest2(s_val(:, i)', h_val(:, i)');
%     end
%     [B, I] = sort(p);
%     s_train = s_train(:, I(1:200));
%     h_train = h_train(:, I(1:200));
%     s_test = s_test(:, I(1:200));
%     h_test = h_test(:, I(1:200));

%     [s_covar, s_precision, ~, ~, ~] = GraphicalLasso(s_train, 0.01);
%     [h_covar, h_precision, ~, ~, ~] = GraphicalLasso(h_train, 0.01);
%     
%     prior = [size(s_train, 1); size(h_train, 1)] / (size(s_train, 1) + size(h_train, 1));
%     l = unique(labels);   n_l = size(l,1);
%     conditional(1) = struct('C', s_precision, 'mu', mean(s_train));
%     conditional(2) = struct('C', h_precision, 'mu', mean(h_train));
%     model = struct('nvars', size(s_train, 2), 'nlabels', n_l, 'labels', l, ...
%                    'class_cond', conditional, 'class_prior', prior, ...
%                    'lambdas', 0);
% 
%     test = [s_test; h_test];
%     t_labels = [ones(size(s_test, 1), 1); -1 * ones(size(h_test, 1), 1)];
%     predictions = zeros(size(test, 1), 1);
%     for i = 1:size(test, 1)
%         s_ll = gaussianFit(s_precision, squeeze(test(i, :))', 0);
%         h_ll = gaussianFit(h_precision, squeeze(test(i, :))', 0);
%         %[y,pyx] = MRFC_predict(squeeze(test(:, :, i))', model);
%         %s_ll = sum(y == 1);
%         %h_ll = sum(y == -1);
%         if s_ll > h_ll
%             predictions(i) = 1;
%         else
%             predictions(i) = -1;
%         end
%     end

    t_labels = [ones(size(s_test, 1), 1); -1 * ones(size(h_test, 1), 1)];
    SVMStruct = svmtrain([s_train; h_train], [ones(size(s_train, 1), 1); -1 * ones(size(h_train, 1), 1)]);
    predictions = svmclassify(SVMStruct, [s_test; h_test]);
    a(iteration) = sum(predictions == t_labels) / length(t_labels);
end
min(a)
max(a)
mean(a)
