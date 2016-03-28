
function fold = generateDataSetBalanced(feature,k)
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
    data = data(1:304,:);
    for i=1:k
		fold(i) = split(data,k,i);
    end
end
function result=split(data,k,n_fold)
    ind_sch = find(data(:,end)==1);
    ind_hel= find(data(:,end)==-1);
	L_hel = length(ind_hel);
	L_sch = length(ind_sch);
    num = min(L_hel,L_sch);
    SizeofFold = floor(num/k);
    if n_fold~=k
        ind_i = (SizeofFold*(n_fold-1)+1:SizeofFold*(n_fold));
        ind_test_sch = ind_i ;
        ind_test_hel = ind_i ;
    else
        ind_test_sch = ((SizeofFold*k+1):L_sch);
        ind_test_hel = ((SizeofFold*k+1):L_hel);
    end
	result = struct('fold_hel', [],'fold_sch', [],'h_y' ,[],'s_y' ,[] ,'NumofFold', n_fold);
    a=size(data(ind_sch(ind_test_sch),1:end-1),1);
    b=size(data(ind_hel(ind_test_hel),1:end-1),1);
    temp_sch = [data(ind_sch(ind_test_sch),1:end-1)'];
    temp_sch_y = [data(ind_sch(ind_test_sch),end)'];
    temp_hel = [data(ind_hel(ind_test_hel),1:end-1)'];
    temp_hel_y = [data(ind_hel(ind_test_hel),end)'];
    result.fold_hel = temp_hel';
    result.h_y = temp_hel_y';
    result.fold_sch = temp_sch';
    result.s_y = temp_sch_y';
    if n_fold==k,
		mess=sprintf('Last fold has : %u%% of Healthy to Schizophernic or vice versa I dont know yet :)' ,floor((a/b)*100) );
		disp(mess);
    end
end
