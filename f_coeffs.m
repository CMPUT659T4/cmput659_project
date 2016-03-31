addpath('lasso','-end');
%load('FBIRN/finaldata_AO/features/fBIRN_AudOdd_allsites_0003_degrees.mat')
%load('FBIRN/finaldata_AO/features/fBIRN_AudOdd_allsites_0003_log_degrees.mat')
% load('FBIRN/finaldata_AO/features/fBIRN_AudOdd_allsites_0003_log_degrees_Tlms_MFG_SFG_.mat')
load('concat_roi_avg.mat')
load('holdout.mat')

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

holdout = data(:, :, holdout_set);
data = data(:, :, train_set);
s_inds = train_labels(:) == 1;
h_inds = train_labels(:) == -1;
s_data = data(:, :, s_inds);
h_data = data(:, :, h_inds);

method='varsel_mrf';
params = [0.7; 0.1; 0.07; 0.01; 0.007; 0.001];
param_acc = zeros(length(params), 1);
top_100_params = zeros(100, length(params));
for p = 1:length(params)
    p
    a = zeros(5, 1);
    top_100_folds = zeros(100, 5);
    for iteration = 1:5
        s_kfold = crossvalind('Kfold', size(s_data, 3) / 4, 10);
        s_kfold = [s_kfold, s_kfold, s_kfold, s_kfold]';
        s_kfold = s_kfold(:);

        h_kfold = crossvalind('Kfold', size(h_data, 3) / 4, 10);
        h_kfold = [h_kfold, h_kfold, h_kfold, h_kfold]';
        h_kfold = h_kfold(:);

        s_test_inds = logical((s_kfold == 1) + (s_kfold == 2));
        h_test_inds = logical((h_kfold == 1) + (h_kfold == 2));
        s_val_inds = logical((s_kfold == 3) + (s_kfold == 4));
        h_val_inds = logical((h_kfold == 3) + (h_kfold == 4));
        s_train_inds = ~(s_test_inds + s_val_inds);
        h_train_inds = ~(h_test_inds + h_val_inds);
        
        s_val = s_data(:, :, s_val_inds);
        h_val = h_data(:, :, h_val_inds);
        s_val = reshape(s_val, [size(s_val, 1) * size(s_val, 2), size(s_val, 3)])';
        h_val = reshape(h_val, [size(h_val, 1) * size(h_val, 2), size(h_val, 3)])';

        s_train = s_data(:, :, s_train_inds);
        s_train = reshape(s_train, [size(s_train, 1) * size(s_train, 2), size(s_train, 3)])';
        h_train = h_data(:, :, h_train_inds);
        h_train = reshape(h_train, [size(h_train, 1) * size(h_train, 2), size(h_train, 3)])';

        s_test = s_data(:, :, s_test_inds);
        h_test = h_data(:, :, h_test_inds);
        s_test = reshape(s_test, [size(s_test, 1) * size(s_test, 2), size(s_test, 3)])';
        h_test = reshape(h_test, [size(h_test, 1) * size(h_test, 2), size(h_test, 3)])';
        
        x = zeros(size(s_val, 2), 1);
        for i = 1:size(s_val, 2)
            [~, x(i)] = ttest2(s_val(:, i)', h_val(:, i)');
        end
        [~, I] = sort(x);
        top_100_folds(:, iteration) = I(1:100);
        s_train = s_train(:, I(1:100));
        h_train = h_train(:, I(1:100));
        s_test = s_test(:, I(1:100));
        h_test = h_test(:, I(1:100));

        %[s_covar, s_precision, ~, ~, ~] = GraphicalLasso(s_train, params(p));
        %[h_covar, h_precision, ~, ~, ~] = GraphicalLasso(h_train, params(p));
        % save(['post_glasso.mat'], 's_inds', 'h_inds', 's_data', 'h_data', ...
        %      's_kfold', 'h_kfold', 's_test_inds', 'h_test_inds', ...
        %      's_train_inds', 'h_train_inds', 's_train', 'h_train', ...
        %      's_test', 'h_test', 's_covar', 'h_covar', 's_precision', 'h_precision')

        %prior = [size(s_train, 1); size(h_train, 1)] / (size(s_train, 1) + size(h_train, 1));
        %l = unique(labels);   n_l = size(l,1);
        %conditional(1) = struct('C', s_precision, 'mu', mean(s_train));
        %conditional(2) = struct('C', h_precision, 'mu', mean(h_train));
        %model = struct('nvars', size(s_train, 2), 'nlabels', n_l, 'labels', l, ...
        %               'class_cond', conditional, 'class_prior', prior, ...
        %               'lambdas', 0);
                   
        model = MRFC_learn([s_train; h_train], [ones(size(s_train, 1), 1); -1 * ones(size(h_train, 1), 1)], method, params(p));

        test = [s_test; h_test];
        t_labels = [ones(size(s_test, 1), 1); -1 * ones(size(h_test, 1), 1)];
        predictions = zeros(size(test, 1), 1);
        for i = 1:size(test, 1)
            %nll(1) = gaussianFit(s_precision, test(i, :), 0);
            %nll(2) = gaussianFit(h_precision, test(i, :), 0);
            [y,pyx] = MRFC_predict(test(i, :), model);
            nll = pyx;
            if nll(1) > nll(2)
                predictions(i) = 1;
            else
                predictions(i) = -1;
            end
        end

        a(iteration) = sum(predictions == t_labels) / length(t_labels);
    end
    param_acc(p) = mean(a);
    [~, index] = max(a);
    top_100_params(:, p) = top_100_folds(:, index);
end

[~, p] = max(param_acc);
top_100 = top_100_params(:, p);
param_acc'
params'

s_data = reshape(s_data, [size(s_data, 1) * size(s_data, 2), size(s_data, 3)])';
h_data = reshape(h_data, [size(h_data, 1) * size(h_data, 2), size(h_data, 3)])';

holdout = reshape(holdout, [size(holdout, 1) * size(holdout, 2), size(holdout, 3)])';

% x = zeros(size(s_data, 2), 1);
% for i = 1:size(s_data, 2)
%     [~, x(i)] = ttest2(s_data(:, i)', h_data(:, i)');
% end
% [~, I] = sort(x);
% s_data = s_data(:, I(1:100));
% h_data = h_data(:, I(1:100));
% holdout = holdout(:, I(1:100));
s_data = s_data(:, top_100);
h_data = h_data(:, top_100);
holdout = holdout(:, top_100);

%[~, s_precision, ~, ~, ~] = GraphicalLasso(s_data, params(p));
%[~, h_precision, ~, ~, ~] = GraphicalLasso(h_data, params(p));

%prior = [size(s_data, 1); size(h_data, 1)] / (size(s_data, 1) + size(h_data, 1));
%l = unique(holdout_labels);   n_l = size(l,1);
%conditional(1) = struct('C', s_precision, 'mu', mean(s_data));
%conditional(2) = struct('C', h_precision, 'mu', mean(h_data));
%model = struct('nvars', size(s_train, 2), 'nlabels', n_l, 'labels', l, ...
%               'class_cond', conditional, 'class_prior', prior, ...
%               'lambdas', 0);
           
model = MRFC_learn([s_data; h_data], [ones(size(s_data, 1), 1); -1 * ones(size(h_data, 1), 1)], method, params(p));

predictions = zeros(size(holdout, 1), 1);
for i = 1:size(holdout, 1)
    %nll(1) = gaussianFit(s_precision, holdout(i, :), 0);
    %nll(2) = gaussianFit(h_precision, holdout(i, :), 0);
    [y,pyx] = MRFC_predict(holdout(i, :), model);
    nll = pyx;
    if nll(1) > nll(2)
        predictions(i) = 1;
    else
        predictions(i) = -1;
    end
end

%SVMStruct = svmtrain([s_data; h_data], [ones(size(s_data, 1), 1); -1 * ones(size(h_data, 1), 1)]);
%predictions = svmclassify(SVMStruct, holdout);
sum(predictions == holdout_labels) / length(holdout_labels)
