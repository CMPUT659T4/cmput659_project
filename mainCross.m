function result = crossTest()

addpath(genpath('FBIRN/PGMTools/'))
addpath(genpath('E:\Google Drive\University\Proabilistic Graphical Models\Project'))


addpath('FBIRN/PGMTools/SparseMRF','-end');
addpath('FBIRN/PGMTools/MRFC','-end');
addpath('FBIRN/PGMTools/SparseMRF/','-end');
addpath(genpath('FBIRN/PGMTools/SparseMRF/'))
addpath(genpath('FBIRN/PGMTools/MRFC/'))

%load('testInd.mat');


load('FBIRN/finaldata_AO/features/fBIRN_AudOdd_allsites_0003_degrees.mat')
%load('fBIRN_AudOdd_allsites_0003_all_corr.mat')
%load('fBIRN_AudOdd_allsites_0003_log_degrees.mat')
%load('fBIRN_AudOdd_allsites_0003_log_degrees_Tlms_MFG_SFG_.mat')
%load('fBIRN_AudOdd_allsites_0003_eigenvalues.mat')
%load('fBIRN_AudOdd_allsites_0003_log_disconnection.mat')
%load('fBIRN_AudOdd_allsites_0003_log_disconnection_Tlms_MFG_SFG.mat')
%load('fBIRN_AudOdd_allsites_0003_mbi_stats.mat')
%load('fBIRN_AudOdd_allsites_0003_norm_avgrgs_degrees.mat')
if size(data,1)~=380,
        error('Data format is incompatible with current test format!!');
end

rho = [0.01,0.1,0.7];
numFeatures = [4,40,400,1000];
numFolds = 5;
percentage = 20;
result = zeros(length(rho),length(numFeatures));
result_t = zeros(length(rho),length(numFeatures),numFolds);
test_ind = TestInd(data,percentage);
train_ind = setdiff([1:size(data,1)],test_ind);
data_t = data(train_ind,:);
save('test_ind1.mat','test_ind');
for i=1:length(rho),
    for j=1:length(numFeatures),
        for k=1:numFolds
            ind = crossind(data_t,numFolds);
            cross_test = find(ind==k);
            result_t(i,j,k) = testCross(data_t,rho(i),numFeatures(j),cross_test);
        end
    end
end
for i=1:length(rho),
    for j=1:length(numFeatures),
        result(i,j) = mean(result_t(i,j,:));
    end
end

end