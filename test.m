addpath(genpath('FBIRN/PGMTools/'))
addpath(genpath('E:\Google Drive\University\Proabilistic Graphical Models\Project'))


addpath('FBIRN/PGMTools/SparseMRF','-end');
addpath('FBIRN/PGMTools/MRFC','-end');
addpath('FBIRN/PGMTools/SparseMRF/','-end');
addpath(genpath('FBIRN/PGMTools/SparseMRF/'))
addpath(genpath('FBIRN/PGMTools/MRFC/'))

load('testInd.mat');


load('FBIRN/finaldata_AO/features/fBIRN_AudOdd_allsites_0003_degrees.mat')
%load('fBIRN_AudOdd_allsites_0003_all_corr.mat')
%load('fBIRN_AudOdd_allsites_0003_log_degrees.mat')
%load('fBIRN_AudOdd_allsites_0003_log_degrees_Tlms_MFG_SFG_.mat')
%load('fBIRN_AudOdd_allsites_0003_eigenvalues.mat')
%load('fBIRN_AudOdd_allsites_0003_log_disconnection.mat')
%load('fBIRN_AudOdd_allsites_0003_log_disconnection_Tlms_MFG_SFG.mat')
%load('fBIRN_AudOdd_allsites_0003_mbi_stats.mat')
%load('fBIRN_AudOdd_allsites_0003_norm_avgrgs_degrees.mat')
numFeatures= 200;
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
method = 'varsel_mrf';
model = MRFC_learn(X, y, method, 0.7);
[y,pyx] = MRFC_predict(X_, model);

Accuracy=(sum(y==y_)/length(y))*100
if Accuracy>=70,
    fprintf('Holy Moly Cow!!!!');
end;