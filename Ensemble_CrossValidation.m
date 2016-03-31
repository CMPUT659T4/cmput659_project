function [Accuracy,PErrVec,FErrVec,AnalyseMatrix]=Ensemble_CrossValidation(itterat)
%limit is the number of variables we are interested in. ie, if limit=10, we
%only look at the first 10 columns of the data vector
AnalyseMatrix=[];

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
    %
    
    
    
    
    %&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& TEST 1 &&&&&&&&&&&&&&&&&&&&&&&
    load('FBIRN/finaldata_AO/features/OurTestMFG_SFGVoxelsLogDegrees.mat');
    data1=OurTestMFG_SFGVoxelsLogDegrees;
    
    %&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& TEST 2 & 3 &&&&&&&&&&&&&&&&&&&&&&&
    load('FBIRN/finaldata_AO/features/OurTestAllVoxelsLogDegrees.mat');
    data2=OurTestAllVoxelsLogDegrees;
    
    %&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& TEST 4 &&&&&&&&&&&&&&&&&&&&&&&
    load('FBIRN/finaldata_AO/features/train_concat_roi_avg.mat')
    data4=train_concat_roi_avg;
    data4([81, 82, 128, 184, 247, 248, 249, 250], :, :) = [];
    
    new_data = zeros(size(data4, 1), size(data4, 3));
    for j = 1:size(data4, 3)
        c = corr(squeeze(data4(:, :, j))');
        c = c - diag(diag(c));
        new_data(:, j) = sum(c > 0.7);
    end
    
    %data4 = log(0.1+new_data);
     data4 = new_data;
     
    limit=20;
    limit2=10;
    limit3=20;
    Set1=find(data1(:,end)==-1)';%healthy
    Set2=find(data1(:,end)==1)'; % patients
    
    Indices1 = crossvalind('Kfold', length(Set1)/4, 5);%Indices of healthy
    Indices2 = crossvalind('Kfold', length(Set2)/4, 5); % Indices of patientss
    % all the indices with number 5 is taken away as the testing set
    % Extend this indices to the 380 numbers
    Indices1=kron(Indices1,[1 1 1 1]');  % extend the healthy indices
    Indices2=kron(Indices2,[1 1 1 1]'); % extend the patient indices (this is similar to repmat)
    
    for k=1:5
        
        trainInd=[Set1(find(Indices1~=k))';Set2(find(Indices2~=k))']; %Combine the indices of training indices of healthy and patients
        testInd=[Set1(find(Indices1==k))';Set2(find(Indices2==k))']; %Combine the test indices of healthy and patients
        
        %&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& TEST 1 &&&&&&&&&&&&&&&&&&&&&&&
        Train1=data1(trainInd,:); % Extract the training data
        Test1=data1(testInd,:); %Extract the test tdata
        
        X_train1=Train1(:,1:limit);
        Y_train1=Train1(:,end);
        X_test1=Test1(:,1:limit);
        Y_test=Test1(:,end);
        
        lambdas=0.7;
        method='varsel_mrf';
        min_err = Inf;
        
        model1 = MRFC_learn(X_train1, Y_train1, method, lambdas);  %learn MRF classifier (MRFC)
        [Test1 Probs] = MRFC_predict(X_test1,model1); %test MRFC
        
        
        %&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& TEST 2 &&&&&&&&&&&&&&&&&&&&&&&
        Healthy=data2(Set1(find(Indices1~=k))',1:end-1); %First k-1 sets are for the training set
        Patients=data2(Set2(find(Indices2~=k))',1:end-1);
        
        Train2=data2(trainInd,:); % Extract the training data
        Test2=data2(testInd,:); %Extract the tes tdata
        
        X_train2=Train2(:,1:(end-1));
        Y_train2=Train2(:,end);
        X_test2=Test2(:,1:(end-1));
        
        diff_mu=-(mean(Healthy)-mean(Patients));
        
        [K index]=sort(diff_mu,'descend');%Sort the differences and get the widely changed pixel values
        
        X_train2=X_train2(:,index(1:limit2));
        X_test2=X_test2(:,index(1:limit2));
        
        lambdas=0.7;
        method='varsel_mrf';
        model2 = MRFC_learn(X_train2, Y_train2, method, lambdas);  %learn MRF classifier (MRFC)
        [Test2 Probs] = MRFC_predict(X_test2,model2); %test MRFC
        
        %&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& TEST 3 &&&&&&&&&&&&&&&&&&&&&&&
        Healthy=data2(Set1(find(Indices1~=k))',1:end-1); %First k-1 sets are for the training set
        Patients=data2(Set2(find(Indices2~=k))',1:end-1);
        
        Train3=data2(trainInd,:); % Extract the training data
        Test3=data2(testInd,:); %Extract the tes tdata
        
        X_train3=Train3(:,1:(end-1));
        Y_train3=Train3(:,end);
        X_test3=Test3(:,1:(end-1));
        
        diff_mu=(mean(Healthy)-mean(Patients));
        
        [K index]=sort(diff_mu,'descend');%Sort the differences and get the widely changed pixel values
        
        X_train3=X_train3(:,index(1:limit3));
        X_test3=X_test3(:,index(1:limit3));
        
        lambdas=0.7;
        method='varsel_mrf';
        model3 = MRFC_learn(X_train3, Y_train3, method, lambdas);  %learn MRF classifier (MRFC)
        [Test3 Probs] = MRFC_predict(X_test3,model3); %test MRFC
        
        
        
        %&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& TEST 4 &&&&&&&&&&&&&&&&&&&&&&&
            
        s_test_inds = Set2(find(Indices2==k))';% + (s_kfold == 3));
        h_test_inds = Set1(find(Indices1==k))';% + (h_kfold == 3));
        s_train_inds = Set2(find(Indices2~=k))';
        h_train_inds = Set1(find(Indices1~=k))';
        
        s_train = data4(:,s_train_inds)';
        h_train = data4(:,h_train_inds)';
              
        
        s_test = data4(:,s_test_inds)';
        h_test =  data4(:,h_test_inds)';
        indt=data4(:,testInd)';
        
        SVMStruct = svmtrain([s_train; h_train], [ones(size(s_train, 1), 1); -1 * ones(size(h_train, 1), 1)]);
        %Test4 = svmclassify(SVMStruct, [s_test; h_test]);
        Test4 = svmclassify(SVMStruct, [indt]);
        
        %&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& Ensembling &&&&&&&&&&&&&&&&&&&&&&&
        
        predicted_Y=Test4;
        
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
fprintf('Running Voxel degree connectivity for %d features\n',limit);
%fprintf('Accuracy = %.2f with std=%.2f Minimum=%.2f Max=%.2f \n False positives=%.2f \n Missed Patients=%.2f\n',mean(Accuracy),std(Accuracy,1),min(Accuracy),max(Accuracy),mean(PErrVec),mean(FErrVec));
fprintf('Accuracy = %.2f with std=%.2f Minimum=%.2f Max=%.2f \n False positives=%.2f \n Missed Patients=%.2f\n',mean(mean(Accuracy)),std(mean(Accuracy),1),min(min(Accuracy)),max(max(Accuracy)),mean(mean(PErrVec)),mean(mean(FErrVec)));

