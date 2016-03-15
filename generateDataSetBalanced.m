
function [train_X,train_Y,test_X,test_Y]=generateBalancedDataSet(feature,k)
% feature selects which feature to use
% this returns the training set and test from the feature files.
% Balance 1/k th portion of the data is selected for the test and the rest is
% selected for the training set in way that 50% in healthy and 50% is
% schizophernic
switch feature
    case 1
        load('FBIRN/finaldata_AO/features/fBIRN_AudOdd_allsites_0003_degrees.mat');
    case 2
        load('FBIRN/finaldata_AO/features/fBIRN_AudOdd_allsites_0003_log_degrees.mat');

    case 3
        load('FBIRN/finaldata_AO/features/fBIRN_AudOdd_allsites_0003_log_degrees_Tlms_MFG_SFG_.mat');
end

    ind_sch = find(data(:,end)==1);
    ind_hel = find(data(:,end)==0); 
    L_sch = size(ind_sch,1);
    L_hel = size(ind_hel,1);
    ind_test_sch = (1:floor(L_sch/k))*k;
    ind_train_sch = setdiff((1:L_sch),ind_test_sch);
    ind_test_hel = (1:floor(L_hel/k))*k;
    ind_train_hel = setdiff((1:L_hel),ind_test_hel);
    num_test = min(length(ind_test_sch),length(ind_test_hel));
    num_train = min(length(ind_train_hel),length(ind_train_sch));
    train_X = [data(ind_train_hel(1:num_train),:) data(ind_train_sch(1:num_train),:)];
    test_X = [data(ind_test_hel(1:num_test),:) data(ind_test_sch(1:num_test),:)];
    train_Y = [data(ind_train_hel(1:num_train),end) data(ind_train_sch(1:num_train),end)];
    test_Y = [data(ind_test_hel(1:num_test),end) data(ind_test_sch(1:num_test),end)];
  
    


end