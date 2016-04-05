function [FinAccuracy,FinPErrVec,FinFErrVec]=FindLambdaRishMFG_SFG_CrossValidation(limit,itterat)
%limit is the number of variables we are interested in. ie, if limit=10, we
%only look at the first 10 columns of the data vector
addpath('FBIRN/PGMTools/SparseMRF','-end');
addpath('FBIRN/PGMTools/MRFC','-end');
addpath('FBIRN/PGMTools/SparseMRF/','-end');
addpath(genpath('FBIRN/PGMTools/SparseMRF/'))
addpath(genpath('FBIRN/PGMTools/MRFC/'))

lambdas=[0.4 0.5 0.6 0.7 0.8 0.9];
FinPErrVec=zeros(length(lambdas),1);
FinFErrVec=zeros(length(lambdas),1);
FinAccuracy=zeros(length(lambdas),1);

PErrVec=zeros(5,itterat);
FErrVec=zeros(5,itterat);
Accuracy=zeros(5,itterat);

for ghj=1:length(lambdas)
    lambdaVal=lambdas(ghj);
    for i=1:itterat
        %%% generating the test and the training set
        %
        %
        load('FBIRN/finaldata_AO/features/OurTestMFG_SFGVoxelsLogDegrees.mat');
        data=OurTestMFG_SFGVoxelsLogDegrees;
        
        %     load('FBIRN/finaldata_AO/features/OurTestMBI_stat.mat');
        %     data=(OurTestMBI_stat);
        %
        
        Set1=find(data(:,end)==-1)';%healthy
        Set2=find(data(:,end)==1)'; % patients
        
        Indices1 = crossvalind('Kfold', length(Set1)/4, 5);%Indices of healthy
        Indices2 = crossvalind('Kfold', length(Set2)/4, 5); % Indices of patientss
        % all the indices with number 5 is taken away as the testing set
        % Extend this indices to the 380 numbers
        Indices1=kron(Indices1,[1 1 1 1]');  % extend the healthy indices
        Indices2=kron(Indices2,[1 1 1 1]'); % extend the patient indices (this is similar to repmat)
        for k=1:5
            Healthy=data(Set1(find(Indices1~=k))',1:end-1); %First k-1 sets are for the training set
            Patients=data(Set2(find(Indices2~=k))',1:end-1);
            Healthy_Test=data(Set1(find(Indices1==k))',1:end-1); %Last k is the test set
            Patient_Test=data(Set2(find(Indices2==k))',1:end-1);
            
            trainInd=[Set1(find(Indices1~=k))';Set2(find(Indices2~=k))']; %Combine the indices of training indices of healthy and patients
            testInd=[Set1(find(Indices1==k))';Set2(find(Indices2==k))']; %Combine the test indices of healthy and patients
            Train=data(trainInd,:); % Extract the training data
            Test=data(testInd,:); %Extract the tes tdata
            
            X_train=Train(:,1:limit);
            Y_train=Train(:,end);
            X_test=Test(:,1:limit);
            Y_test=Test(:,end);
           
            % method='projected_gradient';
            method='varsel_mrf';
            %method='sinco';
            %method='our_covsel';     
           
            %for lambda2 = lambdas
            % model = MRFC_learn(X_train, Y_train, method, lambda1,lambda2);  %learn MRF classifier (MRFC)
            model = MRFC_learn(X_train, Y_train, method, lambdaVal);  %learn MRF classifier (MRFC)
            [predicted_Y Probs] = MRFC_predict(X_test,model); %test MRFC
            FPerr = size(find(predicted_Y > 0 & Y_test < 0 ),1);
            FNerr = size(find(predicted_Y < 0 & Y_test > 0 ),1);
            err = FPerr +FNerr;
            HTot=length(find(Y_test<0));
            STot=length(find(Y_test>0));
         
            Accuracy(k,i)=(1-err/length(Y_test));
            PErrVec(k,i)=FPerr/HTot;
            FErrVec(k,i)=FNerr/STot;
        end
        
    end
    fprintf('****************************************************************\n \n');
    fprintf('Lambda=%d',lambdaVal);
    fprintf('Accuracy = %.2f with std=%.2f Minimum=%.2f Max=%.2f \n False positives=%.2f \n Missed Patients=%.2f\n',mean(mean(Accuracy)),std(mean(Accuracy),1),min(min(Accuracy)),max(max(Accuracy)),mean(mean(PErrVec)),mean(mean(FErrVec)));
    
    FinAccuracy(ghj)=mean(mean(Accuracy));
    FinPErrVec(ghj)=mean(mean(PErrVec));
    FinFErrVec(ghj)=mean(mean(FErrVec));
end
[JHk pl]=max(FinAccuracy);
fprintf('The best value for lambda is %.6f and accuracy is %.2f\n',lambdas(pl),JHk);
end
