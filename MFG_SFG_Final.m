function [Accuracy,PErrVec,FErrVec,AnalyseMatrix]=MFG_SFG_Final(limit)
%limit is the number of variables we are interested in. ie, if limit=10, we
%only look at the first 10 columns of the data vector
AnalyseMatrix=[];

addpath('FBIRN/PGMTools/SparseMRF','-end');
addpath('FBIRN/PGMTools/MRFC','-end');
addpath('FBIRN/PGMTools/SparseMRF/','-end');
addpath(genpath('FBIRN/PGMTools/SparseMRF/'))
addpath(genpath('FBIRN/PGMTools/MRFC/'))


%%% generating the test and the training set
%
%
load('FBIRN/finaldata_AO/features/OurTestMFG_SFGVoxelsLogDegrees.mat');
data=OurTestMFG_SFGVoxelsLogDegrees;

%     load('FBIRN/finaldata_AO/features/OurTestMBI_stat.mat');
%     data=(OurTestMBI_stat);
%
X_train=data(:,1:limit);
Y_train=data(:,end);


lambdas=0.7;
% method='projected_gradient';
method='varsel_mrf';
%method='sinco';
%method='our_covsel';
%          model = MRFC_CVlambda_learn(X_train, Y_train, method, lambda);  %learn MRF classifier (MRFC)
%          predicted_Y = MRFC_predict(X_test,model); %test MRFC
min_err = Inf;


%for lambda2 = lambdas
% model = MRFC_learn(X_train, Y_train, method, lambda1,lambda2);  %learn MRF classifier (MRFC)
model = MRFC_learn(X_train, Y_train, method, lambdas);  %learn MRF classifier (MRFC)

load('FBIRN/finaldata_AO/features/FinalTestSFG_MFG.mat');
test=OurTestMFG_SFGVoxelsLogDegrees;
X_test=test(:,1:limit);
Y_test=test(:,end);

[predicted_Y Probs] = MRFC_predict(X_test,model); %test MRFC

FPerr = size(find(predicted_Y > 0 & Y_test < 0 ),1);
FNerr = size(find(predicted_Y < 0 & Y_test > 0 ),1);
err = FPerr +FNerr;
HTot=length(find(Y_test<0));
STot=length(find(Y_test>0));

%Just checking underfit_Overfit
%         [predicted_Y Probs] = MRFC_predict(X_train,model); %test MRFC
%
%         FPerr = size(find(predicted_Y > 0 & Y_train < 0 ),1);
%         FNerr = size(find(predicted_Y < 0 & Y_train > 0 ),1);
%         err = FPerr +FNerr;
%         HTot=length(find(Y_train<0));
%         STot=length(find(Y_train>0));
fprintf('Fasle Positives = %.2f, Missed Patients = %.2f, Error Rate for lambda=%.2f is %.2f\n',FPerr/HTot,FNerr/STot,lambdas,err/(length(Y_test)));
if err < min_err
    min_err=err;
    best_model=model;
    best_pred = predicted_Y;
end
% end


%Convert the log likelihoods to normal space
%
AnalyseMatrix=[AnalyseMatrix;Probs Y_test predicted_Y];

model=best_model;
Y_test;
predicted_Y = best_pred;
Accuracy=(1-err/length(Y_test));
PErrVec=FPerr/HTot;
FErrVec=FNerr/STot;

fprintf('****************************************************************\n \n');
fprintf('Running Voxel degree connectivity for %d features\n',limit);
%fprintf('Accuracy = %.2f with std=%.2f Minimum=%.2f Max=%.2f \n False positives=%.2f \n Missed Patients=%.2f\n',mean(Accuracy),std(Accuracy,1),min(Accuracy),max(Accuracy),mean(PErrVec),mean(FErrVec));
fprintf('Accuracy = %.2f with std=%.2f Minimum=%.2f Max=%.2f \n False positives=%.2f \n Missed Patients=%.2f\n',mean(mean(Accuracy)),std(mean(Accuracy),1),min(min(Accuracy)),max(max(Accuracy)),mean(mean(PErrVec)),mean(mean(FErrVec)));

