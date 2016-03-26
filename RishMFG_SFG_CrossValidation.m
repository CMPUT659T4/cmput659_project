function [Accuracy,PErrVec,FErrVec,Met,WrongMAt,GoodMat]=RishMFG_SFG_CrossValidation(limit,itterat)
%limit is the number of variables we are interested in. ie, if limit=10, we
%only look at the first 10 columns of the data vector
Met=[];
WrongMAt=[];
GoodMat=[];
addpath('FBIRN/PGMTools/SparseMRF','-end');
addpath('FBIRN/PGMTools/MRFC','-end');
addpath('FBIRN/PGMTools/SparseMRF/','-end');
addpath(genpath('FBIRN/PGMTools/SparseMRF/'))
addpath(genpath('FBIRN/PGMTools/MRFC/'))




PErrVec=zeros(5,itterat);
FErrVec=zeros(5,itterat);
Accuracy=zeros(5,itterat);
for i=1:itterat
    %%% generating the test and the training set
    
    
    load('FBIRN/finaldata_AO/features/OurTestMFG_SFGVoxelsLogDegrees.mat');
    data=OurTestMFG_SFGVoxelsLogDegrees;
    
    
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
        [predicted_Y Probs] = MRFC_predict(X_test,model); %test MRFC
        Probs;
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
        Met=[Probs (Probs(:,1)-Probs(:,2)) Y_test predicted_Y];
        Prob=Met;
        %     Prob(:,1)=exp(Met(:,1));
        %     Prob(:,2)=exp(Met(:,2));
        
        Prob(:,1)=exp(Met(:,1))./(exp(Met(:,1))+exp(Met(:,2)));;
        Prob(:,2)=exp(Met(:,2))./(exp(Met(:,1))+exp(Met(:,2)));;
        Met=Prob;
        WrongMAt=Met(find(Y_test~=predicted_Y),:);
        GoodMat=Met(find(Y_test==predicted_Y),:);
        model=best_model;
        Y_test;
        predicted_Y = best_pred;
        Accuracy(k,i)=(1-err/length(Y_test));
        PErrVec(k,i)=FPerr/HTot;
        FErrVec(k,i)=FNerr/STot;
    end
end
fprintf('****************************************************************\n \n');
fprintf('Running Voxel degree connectivity for %d features\n',limit);
%fprintf('Accuracy = %.2f with std=%.2f Minimum=%.2f Max=%.2f \n False positives=%.2f \n Missed Patients=%.2f\n',mean(Accuracy),std(Accuracy,1),min(Accuracy),max(Accuracy),mean(PErrVec),mean(FErrVec));
fprintf('Accuracy = %.2f with std=%.2f Minimum=%.2f Max=%.2f \n False positives=%.2f \n Missed Patients=%.2f\n',mean(mean(Accuracy)),std(mean(Accuracy),1),min(min(Accuracy)),max(max(Accuracy)),mean(mean(PErrVec)),mean(mean(FErrVec)));

