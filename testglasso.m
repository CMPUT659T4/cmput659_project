function Accuracy = testglasso(rho,numFeatures)
addpath(genpath('FBIRN/PGMTools/'))
addpath(genpath('E:\Google Drive\University\Proabilistic Graphical Models\Project'))


addpath('FBIRN/PGMTools/SparseMRF','-end');
addpath('FBIRN/PGMTools/MRFC','-end');
addpath('FBIRN/PGMTools/SparseMRF/','-end');
addpath(genpath('FBIRN/PGMTools/SparseMRF/'))
addpath(genpath('FBIRN/PGMTools/MRFC/'))

load('test_ind2.mat');


load('FBIRN/finaldata_AO/features/fBIRN_AudOdd_allsites_0003_degrees.mat')
%load('fBIRN_AudOdd_allsites_0003_all_corr.mat')
%load('fBIRN_AudOdd_allsites_0003_log_degrees.mat')
%load('fBIRN_AudOdd_allsites_0003_log_degrees_Tlms_MFG_SFG_.mat')
%load('fBIRN_AudOdd_allsites_0003_eigenvalues.mat')
%load('fBIRN_AudOdd_allsites_0003_log_disconnection.mat')
%load('fBIRN_AudOdd_allsites_0003_log_disconnection_Tlms_MFG_SFG.mat')
%load('fBIRN_AudOdd_allsites_0003_mbi_stats.mat')
%load('fBIRN_AudOdd_allsites_0003_norm_avgrgs_degrees.mat')
%numFeatures= 1000;
%y=data(:,end);
%a=find(y==1);
%sch=a(1:38);
%a=find(y==-1);
%hel=a(1:38);
%test_ind=[hel;sch];
train = setdiff([1:380],test_ind);
[~,p] = ttest(data(train,:));
[~,I] = sort(p);
ind = I(1:numFeatures);
X = data(train,ind);
%[~,p] = ttest(data(test,:));
%[~,I] = sort(p);
X_ = data(test_ind,ind);
y = data(train,end);
y_ = data(test_ind,end);



addpath('lasso','-end');

s_inds = y(:) == 1;
h_inds = y(:) == -1;
s_train = X(s_inds,:);
h_train = X(h_inds, :);
[~, s_precision, ~, ~, ~] = GraphicalLasso(s_train, rho);
[~, h_precision, ~, ~, ~] = GraphicalLasso(h_train, rho);
predictions = zeros(size(test_ind, 1), 1);
for i = 1:length(test_ind)
    nll(1) = gaussianFit(s_precision, X_(i,:), 0);
    nll(2) = gaussianFit(h_precision, X_(i,:), 0);
    %[y,pyx] = MRFC_predict(squeeze(holdout(:, :, i))', model);
    %nll = sum(pyx);
    if nll(1) > nll(2)
        predictions(i) = 1;
    else
        predictions(i) = -1;
    end
end



% method ='varsel_mrf'; %'kataya';%'projected_grad';%'covsel';%'varsel_mrf';
% %model = MRFC_learn(X, y, method, rho);
% [y,pyx] = MRFC_predict(X_, model);

Accuracy = (sum(predictions == y_)/length(y))*100
if Accuracy>=70,
    fprintf('Holy Moly Cow!!!!');
end;
end