function result = mainCrossAve()

addpath(genpath('FBIRN/PGMTools/'))
addpath(genpath('E:\Google Drive\University\Proabilistic Graphical Models\Project'))


addpath('FBIRN/PGMTools/SparseMRF','-end');
addpath('FBIRN/PGMTools/MRFC','-end');
addpath('FBIRN/PGMTools/SparseMRF/','-end');
addpath(genpath('FBIRN/PGMTools/SparseMRF/'))
addpath(genpath('FBIRN/PGMTools/MRFC/'))

%load('testInd.mat');


%load('FBIRN/finaldata_AO/features/fBIRN_AudOdd_allsites_0003_degrees.mat')
%load('fBIRN_AudOdd_allsites_0003_all_corr.mat')
load('fBIRN_AudOdd_allsites_0003_normMax_log_degrees_2.mat');
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

rho = 0.1;
numFeatures = 1000;
numLoops = 5;
percentage = 20;

result_t = zeros(numLoops,1);

test_ind = struct('ind',(percentage*380)/100);
for i=1:numLoops,
    test_ind(i).ind = TestInd(data,percentage);
    result_t(i,1) = testCross(data,rho,numFeatures,test_ind(i).ind);
end
%save('test_ind_vector.mat','test_ind');
result = mean(result_t);
save('resultAVe_Log_Degree_400.mat','result');
save('resultAVe_Log_Degree_400_cros.mat','result_t');
end