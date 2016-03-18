function [Healthy,Patients,Healthy_Test,Patient_Test X_Train,Y_Train,X_Test,Y_Test]=getTwoFeatureSets(k,feature)
%This one returns four data arrays by dividing the patients and healthy people
% 1/k th portion of the data is kept for further processing
%for the test and the rest is
% selected for he training set

switch feature
    case 1
        load('FBIRN/finaldata_AO/features/fBIRN_AudOdd_allsites_0003_degrees.mat');
    case 2
        load('FBIRN/finaldata_AO/features/fBIRN_AudOdd_allsites_0003_log_degrees.mat');
        
    case 3
        load('FBIRN/finaldata_AO/features/fBIRN_AudOdd_allsites_0003_log_degrees_Tlms_MFG_SFG_.mat');
end

%Select from all the data (We dont worry about balancing the training
%and the testing sets with equal number of patients and controls)
data;
Set1=find(data(:,end)==-1)';%healthy
Set2=find(data(:,end)==1)';

Indices1 = crossvalind('Kfold', length(Set1)/4, k);
Indices2 = crossvalind('Kfold', length(Set2)/4, k);
% all the indices with number 5 is taken away as the testing set
% Extend this indices to the 380 numbers
Indices1=kron(Indices1,[1 1 1 1]');
Indices2=kron(Indices2,[1 1 1 1]');


Healthy=data(Set1(find(Indices1<=k-1))',:);
Patients=data(Set2(find(Indices2<=k-1))',:);
Healthy_Test=data(Set1(find(Indices1>k-1))',:);
Patient_Test=data(Set2(find(Indices2>k-1))',:);
Healthy= Healthy(:,1:end-1);
Patients=Patients(:,1:end-1);
Healthy_Test=Healthy_Test(:,1:end-1);
Patient_Test=Patient_Test(:,1:end-1);

 trainInd=[Set1(find(Indices1<=k-1))';Set2(find(Indices2<=k-1))'];
    testInd=[Set1(find(Indices1>k-1))';Set2(find(Indices2>k-1))'];
   Train=data(trainInd,:);
    Test=data(testInd,:);

    X_Train=Train(:,1:(end-1));
    Y_Train=Train(:,end);
    X_Test=Test(:,1:(end-1));
    Y_Test=Test(:,end);
