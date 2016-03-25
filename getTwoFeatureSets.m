function [Healthy,Patients,Healthy_Test,Patient_Test,X_Train,Y_Train,X_Test,Y_Test]=getTwoFeatureSets(k,feature)
%This one returns four data arrays by dividing the patients and healthy people
% 1/k th portion of the data is kept for further processing
%for the test and the rest is
% selected for the training set
% switch feature
%     case 1
%         load('FBIRN/finaldata_AO/features/fBIRN_AudOdd_allsites_0003_degrees.mat');
%     case 2
%         load('FBIRN/finaldata_AO/features/fBIRN_AudOdd_allsites_0003_log_degrees.mat');
%
%     case 3
%         load('FBIRN/finaldata_AO/features/fBIRN_AudOdd_allsites_0003_log_degrees_Tlms_MFG_SFG_.mat');
% end
switch feature
    case 1
        load('FBIRN/finaldata_AO/features/OurTestAllVoxels.mat');
        data=OurTestAllVoxels;
    case 2
        load('FBIRN/finaldata_AO/features/OurTestAllVoxelsLogDegrees.mat');
        data=OurTestAllVoxelsLogDegrees;
    case 3
        load('FBIRN/finaldata_AO/features/OurTestMFG_SFGVoxelsLogDegrees.mat');
        data=OurTestMFG_SFGVoxelsLogDegrees;
end
%Select from all the data (We dont worry about balancing the training
%and the testing sets with equal number of patients and controls)
Set1=find(data(:,end)==-1)';%healthy
Set2=find(data(:,end)==1)'; % patients

Indices1 = crossvalind('Kfold', length(Set1)/4, k);%Indices of healthy
Indices2 = crossvalind('Kfold', length(Set2)/4, k); % Indices of patientss
% all the indices with number 5 is taken away as the testing set
% Extend this indices to the 380 numbers
Indices1=kron(Indices1,[1 1 1 1]');  % extend the healthy indices
Indices2=kron(Indices2,[1 1 1 1]'); % extend the patient indices (this is similar to repmat)


Healthy=data(Set1(find(Indices1<=k-1))',:); %First k-1 sets are for the training set
Patients=data(Set2(find(Indices2<=k-1))',:);
Healthy_Test=data(Set1(find(Indices1>k-1))',:); %Last k is the test set
Patient_Test=data(Set2(find(Indices2>k-1))',:);
Healthy= Healthy(:,1:end-1); % Remove the last column
Patients=Patients(:,1:end-1);
Healthy_Test=Healthy_Test(:,1:end-1); %remove the last column
Patient_Test=Patient_Test(:,1:end-1);

trainInd=[Set1(find(Indices1<=k-1))';Set2(find(Indices2<=k-1))']; %Combine the indices of training indices of healthy and patients
testInd=[Set1(find(Indices1>k-1))';Set2(find(Indices2>k-1))']; %Combine the test indices of healthy and patients
Train=data(trainInd,:); % Extract the training data
Test=data(testInd,:); %Extract the tes tdata

X_Train=Train(:,1:(end-1));
Y_Train=Train(:,end);
X_Test=Test(:,1:(end-1));
Y_Test=Test(:,end);
