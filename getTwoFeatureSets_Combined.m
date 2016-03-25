function [Healthy,Patients,Healthy_Test,Patient_Test X_Train1 X_Train2 Y_Train X_Test1 X_Test2 Y_Test]=getTwoFeatureSets_Combined(k,feature1,feature2)
%This one returns four data arrays by dividing the patients and healthy people
% 1/k th portion of the data is kept for further processing
%for the test and the rest is
% selected for he training set

% switch feature
%     case 1
%         load('FBIRN/finaldata_AO/features/fBIRN_AudOdd_allsites_0003_degrees.mat');
%     case 2
%         load('FBIRN/finaldata_AO/features/fBIRN_AudOdd_allsites_0003_log_degrees.mat');
%
%     case 3
%         load('FBIRN/finaldata_AO/features/fBIRN_AudOdd_allsites_0003_log_degrees_Tlms_MFG_SFG_.mat');
% end
switch feature1
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

X_Train1=Train(:,1:(end-1));
Y_Train=Train(:,end);
X_Test1=Test(:,1:(end-1));
Y_Test=Test(:,end);



switch feature2
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
data;
Set1=find(data(:,end)==-1)';%healthy
Set2=find(data(:,end)==1)';
% 
% Indices1 = crossvalind('Kfold', length(Set1)/4, k);
% Indices2 = crossvalind('Kfold', length(Set2)/4, k);
% all the indices with number 5 is taken away as the testing set
% Extend this indices to the 380 numbers
% Indices1=kron(Indices1,[1 1 1 1]');
% Indices2=kron(Indices2,[1 1 1 1]');


Healthy2=data(Set1(find(Indices1<=k-1))',:);
Patients2=data(Set2(find(Indices2<=k-1))',:);
Healthy_Test2=data(Set1(find(Indices1>k-1))',:);
Patient_Test2=data(Set2(find(Indices2>k-1))',:);
Healthy2= Healthy2(:,1:end-1);
Patients2=Patients2(:,1:end-1);
Healthy_Test2=Healthy_Test2(:,1:end-1);
Patient_Test2=Patient_Test2(:,1:end-1);

trainInd2=[Set1(find(Indices1<=k-1))';Set2(find(Indices2<=k-1))'];
testInd2=[Set1(find(Indices1>k-1))';Set2(find(Indices2>k-1))'];
Train2=data(trainInd2,:);
Test2=data(testInd2,:);

X_Train2=Train2(:,1:(end-1));
Y_Train2=Train2(:,end);
X_Test2=Test2(:,1:(end-1));
Y_Test=Test(:,end);
