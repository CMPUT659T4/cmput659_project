
function fold=generateDataSetBalanced(feature,k)
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

    for i=1:k
		fold(i)= split(data,k,i);
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

        %ind_train_sch = setdiff((1:L_sch),ind_test_sch);
        ind_test_sch = ind_i ;
        ind_test_hel = ind_i ;
    else
        ind_test_sch = ((SizeofFold*k+1):L_sch);
        %ind_train_sch = setdiff((1:L_sch),ind_test_sch);
        ind_test_hel = ((SizeofFold*k+1):L_hel);
    end
	%dim =size(data);
	%dim(1)= length([ind_test_hel , ind_test_sch]);
	result = struct('fold', [], 'NumofFold', n_fold);
	%DATA.folds = struct('data',size(dim));
	
	%ind_hel(ind_test_hel)
    %ind_train_hel = setdiff((1:L_hel),ind_test_hel);
    %num_test = min(length(ind_test_sch),length(ind_test_hel));
    %num_train = min(length(ind_train_hel),length(ind_train_sch));
    %train_X = [data(ind_train_hel(1:num_train),:) data(ind_train_sch(1:num_train),:)];
    %SizeofFold
	%ind_i 
	%ind_test_sch
    %ind_sch
    a=size(data(ind_sch(ind_test_sch),1:end-1),1);
    b=size(data(ind_hel(ind_test_hel),1:end-1),1);
    
    temp=[data(ind_sch(ind_test_sch),1:end-1)', data(ind_hel(ind_test_hel),1:end-1)'];
    result.fold = temp';
    if n_fold==k,
		mess=sprintf('Last fold has : %u%% of Healthy to Schizophernic or vice versa I dont know yet :)' ,floor((a/b)*100) );
		disp(mess);
	end
        
    
    %train_Y = [data(ind_train_hel(1:num_train),end) data(ind_train_sch(1:num_train),end)];
    %test_Y = [data(ind_test_hel(1:num_test),end) data(ind_test_sch(1:num_test),end)];
end
