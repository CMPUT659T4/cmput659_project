function [X_Train,Y_Train,X_Test,Y_Test]=generateBalancedDataSet(site,k,feature)
%This one balances the number of healthy vs Schizophrenia data among the
%training and the test set
% this returns the training set and test from the feature files.
% 1/k th portion of the data is selected for the test and the rest is
% selected for he training set

switch feature
    case 1
        load('FBIRN/finaldata_AO/features/fBIRN_AudOdd_allsites_0003_degrees.mat');
    case 2
        load('FBIRN/finaldata_AO/features/fBIRN_AudOdd_allsites_0003_log_degrees.mat');

    case 3
        load('FBIRN/finaldata_AO/features/fBIRN_AudOdd_allsites_0003_log_degrees_Tlms_MFG_SFG_.mat');
end
%data=data(1:20,:);
if site==6 
    %Select from all the data (We dont worry about balancing the training
    %and the testing sets with equal number of patients and controls)
    Set1=find(data(:,end)==-1)';
    Set2=find(data(:,end)==1)';
    
    Indices1 = crossvalind('Kfold', length(Set1)/4, k);
    Indices2 = crossvalind('Kfold', length(Set2)/4, k);
    % all the indices with number 5 is taken away as the testing set
    % Extend this indices to the 380 numbers
    Indices1=kron(Indices1,[1 1 1 1]');
    Indices2=kron(Indices2,[1 1 1 1]');
    %Indices=[Indices1;Indices2]
    find(Indices1<=k-1)';
    find(Indices2<=k-1)';
    Set1(find(Indices1<=k-1))';
    Set2(find(Indices2<=k-1))';
    trainInd=[Set1(find(Indices1<=k-1))';Set2(find(Indices2<=k-1))'];
    testInd=[Set1(find(Indices1>k-1))';Set2(find(Indices2>k-1))'];
   Train=data(trainInd,:);
    Test=data(testInd,:);
%     
%     Train=[data(Set1(find(Indices1<=k-1)),:);data(Set2(find(Indices2<=k-1)),:)];
%     Test=[data(Set2(find(Indices1>k-1)),:);data(Set2(find(Indices2>k-1)),:)];
%     
%     
%     
%     Train=[Set1(find(Indices1<=k-1),:);Set2(find(Indices2<=k-1),:)];
%     Test=[Set2(find(Indices1>k-1),:);Set2(find(Indices2>k-1),:)];
    X_Train=Train(:,1:(end-1));
    Y_Train=Train(:,end);
    X_Test=Test(:,1:(end-1));
    Y_Test=Test(:,end);
else
    
end
