function [predict labels]=BestFeaturesPrediction_Mario()
%limit is the number of variables we are interested in. ie, if limit=10, we
%only look at the first 10 columns of the data vector
%This time doing the 5-fold cros validation over the data set

AnalyseMatrix=[];

limit=15;
addpath('FBIRN/PGMTools/SparseMRF','-end');
addpath('FBIRN/PGMTools/MRFC','-end');
addpath('FBIRN/PGMTools/SparseMRF/','-end');
addpath(genpath('FBIRN/PGMTools/SparseMRF/'))
addpath(genpath('FBIRN/PGMTools/MRFC/'))
%addpath('FBIRN/PGMTools/SparseMRF','-end')



%%% generating the test and the training set
%[H P H_t P_t X_train,Y_train,X_test,Y_test]=getTwoFeatureSets(5,2);
%load('FBIRN/finaldata_AO/features/OurTestAllVoxelsLogDegrees.mat');
load('FBIRN/finaldata_AO/features/fBIRN_AudOdd_allsites_0003_log_degrees.mat');
load('testInd.mat');
%IndicesToTest=TestInd(data,19);
IndicesToTest=test_ind;
size(IndicesToTest);
%data=OurTestAllVoxelsLogDegrees;
data1=data;
data1(IndicesToTest,:)=[];
size(data1);
%Select from all the data (We dont worry about balancing the training
%and the testing sets with equal number of patients and controls)
Set1=find(data1(:,end)==-1)';%healthy
Set2=find(data1(:,end)==1)'; % patients

Healthy=data1(Set1,1:end-1); %First k-1 sets are for the training set
Patients=data1(Set2,1:end-1);

X_train=data1(:,1:(end-1));
Y_train=data1(:,end);


%fprintf('size of H is %d and P is %d, total training set %d, Test set %d\n',size(H,1),size(P,1),size(X_train,1),size(X_test,1));
%diff_mu=-(mean(Healthy)-mean(Patients));%get the difference between the means of each voxel for healthy and Schizophrenic subjects

diff_mu=-(mean(Healthy)-mean(Patients));

[K index]=sort(diff_mu,'descend');%Sort the differences and get the widely changed pixel values
X_train=X_train(:,index(1:limit));

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

%load('FBIRN/finaldata_AO/features/FinalTestVoxelDegreesLog.mat');
%test=FinalTestLog;
test=data(IndicesToTest,:);
size(test)
X_test=test(:,index(1:limit));
Y_test=test(:,end);


[predicted_Y Probs] = MRFC_predict(X_test,model); %test MRFC

FPerr = size(find(predicted_Y > 0 & Y_test < 0 ),1);
FNerr = size(find(predicted_Y < 0 & Y_test > 0 ),1);
err = FPerr +FNerr;
HTot=length(find(Y_test<0));
STot=length(find(Y_test>0));

% end
AnalyseMatrix=[AnalyseMatrix;Y_test Probs Probs(:,1)-Probs(:,2)  predicted_Y];
Accuracy=(1-err/length(Y_test));
PErrVec=FPerr/HTot;
FErrVec=FNerr/STot;
fprintf('Accuracy=%.2f Fasle Positives = %.2f, Missed Patients = %.2f, Error Rate for lambda=%.2f\n',Accuracy,FPerr/HTot,FNerr/STot,lambdas)
predict=predicted_Y ;
labels=Y_test;
end