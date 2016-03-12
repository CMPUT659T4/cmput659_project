function [Accuracy,PErrVec,FErrVec,site2]=calcRishProjGrad(limit,site,itterat)
%limit is the number of variables we are interested in. ie, if limit=10, we
%only look at the first 10 columns of the data vector
addpath('FBIRN/PGMTools/SparseMRF','-end');
addpath('FBIRN/PGMTools/MRFC','-end');
addpath('FBIRN/PGMTools/SparseMRF/','-end');
addpath(genpath('FBIRN/PGMTools/SparseMRF/'))
%addpath('FBIRN/PGMTools/SparseMRF','-end')

site2=site;
ErrVec=zeros(1,itterat);
FErrVec=zeros(1,itterat);
Accuracy=zeros(1,itterat);
for i=1:itterat
    %%% generating the test and the training set   
   
   [X_train,Y_train,X_test,Y_test]=generateDataSet(site,5);
    X_train=X_train(:,1:limit);
     X_test=X_test(:,1:limit);  
    lambdas=[0 0.01 0.1 0.5 1 5];
    %lambdas=[0.7];% 0.5 1 5];
    %lambdas=[0.1:0.2:0.9];;% 0.5 1 5];
    lambdas=0.7;
    method='projected_gradient';
    %method='varsel_mrf';
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
    Accuracy(i)=(1-err/length(Y_test));
    PErrVec(i)=FPerr;
    FErrVec(i)=FNerr;
end

