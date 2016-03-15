%  case 'mrf'
%             model = MRFC_learn(X_train, Y_train, method, lambdas);  %learn MRF classifier (MRFC)
%             [predicted_Y,pyx] = MRFC_predict(X_test,model); %test MRFC   features_used = ranked_features(1:top_features);
%
%             outf = sprintf('results/pr_%s_top%d_model_fold%dout_of%d_%s',lower(classifier),top_features,fold,folds,in_fnm);
%             save(outf, 'pyx');
%
%             %predicted_Y = MRFC_predict_varsel(X_test,model);
%             % loglik(fold,Y_test(1)) = MRFC_loglik(X_test,Y_test,model);
%
%
%
%
% %        case 'mrf_cvlambda'

%loading the files - Copied from Neil's Code

%This is just for future use
% lowres_dir = dir('FBIRN/finaldata_AO/lowresRef');
% dirs = cell(1, length(lowres_dir) - 3);
% for folder = 3:length(lowres_dir)-1
%     f_name = lowres_dir(folder).name;
%     dirs{folder - 2} = strcat('FBIRN/finaldata_AO/lowresRef/', f_name);
% end
%
% for folder = 1:length(dirs)
%     f_name = dirs{folder};
%     fprintf([f_name,'\n']);
% end
addpath('FBIRN/PGMTools/SparseMRF','-end');
addpath('FBIRN/PGMTools/MRFC','-end');
addpath('FBIRN/PGMTools/SparseMRF/','-end');
addpath(genpath('FBIRN/PGMTools/SparseMRF/'))
%addpath('FBIRN/PGMTools/SparseMRF','-end')
load('FBIRN/finaldata_AO/features\fBIRN_AudOdd_allsites_0003_log_degrees_Tlms_MFG_SFG_.mat');
%load('FBIRN/finaldata_AO\features\fBIRN_AudOdd_allsites_0003_degrees.mat');


trainLimit=40;


%%% generating the test and the training set
Indices = crossvalind('Kfold', 95, 5);
% all the indices with number 5 is taken away as the testing set
% Extend this indices to the 380 numbers
Indices=kron(Indices,[1 1 1 1]');
Train=data(find(Indices<=4),:);
Test=data(find(Indices>4),:);
X_train=Train(:,1:(end-1)/50);
Y_train=Train(:,end);
X_test=Test(:,1:(end-1)/50);
Y_test=Test(:,end);
lambdas=[0 0.01 0.1 0.5 1 5];
%lambdas=[0.7];% 0.5 1 5];
%lambdas=[0.1:0.2:0.9];;% 0.5 1 5];
lambdas=0.7;
%method='projected_gradient';
method='varsel_mrf';
%method='sinco';



%          model = MRFC_CVlambda_learn(X_train, Y_train, method, lambda);  %learn MRF classifier (MRFC)
%          predicted_Y = MRFC_predict(X_test,model); %test MRFC
min_err = Inf;

for lambda1=lambdas
    %for lambda2 = lambdas
    % model = MRFC_learn(X_train, Y_train, method, lambda1,lambda2);  %learn MRF classifier (MRFC)
    model = MRFC_learn(X_train, Y_train, method, lambda1);  %learn MRF classifier (MRFC)
    predicted_Y = MRFC_predict(X_test,model); %test MRFC
    
    FPerr = size(find(predicted_Y > 0 & Y_test < 0 ),1);
    FNerr = size(find(predicted_Y < 0 & Y_test > 0 ),1);
    err = FPerr +FNerr;
    fprintf('Fasle positives = %d, false negatives = %d, Error Rate for lambda=%.2f is %.2f\n',FPerr,FNerr,lambda1,err/(376*.2));
    if err < min_err
        min_err=err;
        best_model=model;
        best_pred = predicted_Y;
    end
    % end
end

model=best_model;
Y_test;
predicted_Y = best_pred;




