function [Accuracy,PErrVec,FErrVec,site2]=combinedPredict(limit,site,itterat)
%limit is the number of variables we are interested in. ie, if limit=10, we
%only look at the first 10 columns of the data vector

addpath('FBIRN/PGMTools/SparseMRF','-end');
addpath('FBIRN/PGMTools/MRFC','-end');
addpath('FBIRN/PGMTools/SparseMRF/','-end');
addpath(genpath('FBIRN/PGMTools/SparseMRF/'))
addpath(genpath('FBIRN/PGMTools/MRFC/'))
%addpath('FBIRN/PGMTools/SparseMRF','-end')

site2=site;
ErrVec=zeros(1,itterat);
FErrVec=zeros(1,itterat);
Accuracy=zeros(1,itterat);
for i=1:itterat
    %%% generating the test and the training set
    
    [H P H_t P_t X_train1 X_train2 Y_train X_test1 X_test2 Y_test]=getTwoFeatureSets_Combined(5,2,3);
    diff_mu=-(mean(H)-mean(P));%get the difference between the means of each voxel for healthy and Schizophrenic subjects
    [K index]=sort(diff_mu,'descend');%Sort the differences and get the widely changed pixel values
    
    
    size(X_test1);
    X_train1=X_train1(:,index(1:limit));
    X_test1=X_test1(:,index(1:limit));
    lambdas=0.7;
    % method='projected_gradient';
    method='varsel_mrf';
    %method='sinco';
    %method='our_covsel';
    X_train2=X_train2(:,1:4504);
    X_test2=X_test2(:,1:4504);
    min_err = Inf;
    
    for lambda1=lambdas
        %for lambda2 = lambdas
        % model = MRFC_learn(X_train, Y_train, method, lambda1,lambda2);  %learn MRF classifier (MRFC)
        model1 = MRFC_learn(X_train1, Y_train, method, lambda1);  %learn MRF classifier (MRFC)
        
        model2 = MRFC_learn(X_train2, Y_train, method, lambda1);
        
        [predicted_Y1 Probs1] = MRFC_predict(X_test1,model1); %test MRFC
        [predicted_Y2 Probs2] = MRFC_predict(X_test2,model2);
        
        
        ProbsTot=0.8*Probs1+0.2*Probs2;
        size(ProbsTot);
        %Just checking underfit_Overfit
        %         [predicted_Y Probs] = MRFC_predict(X_train,model); %test MRFC
        %
        %         FPerr = size(find(predicted_Y > 0 & Y_train < 0 ),1)
        %         FNerr = size(find(predicted_Y < 0 & Y_train > 0 ),1)
        %         err = FPerr +FNerr
        %         HTot=length(find(Y_train<0))
        %         STot=length(find(Y_train>0))
      
        [a,ind]=max(ProbsTot');
        a;
        predicted_Y=zeros(size(predicted_Y1));
        for kl=1:length(Y_test)
            predicted_Y(kl)=model1.labels(ind(kl));
        end
         predicted_Y- predicted_Y1;
        FPerr = size(find(predicted_Y > 0 & Y_test < 0 ),1);
        FNerr = size(find(predicted_Y < 0 & Y_test > 0 ),1);
        err = FPerr +FNerr;
        HTot=length(find(Y_test<0));
        STot=length(find(Y_test>0));
        fprintf('Fasle Positives = %.2f, Missed Patients = %.2f, Error Rate for lambda=%.2f is %.2f\n',FPerr/HTot,FNerr/STot,lambda1,err/(length(Y_test)));
        if err < min_err
            min_err=err;
            best_model1=model1;
            best_pred = predicted_Y;
        end
        % end
    end
    
    %     Met=[Probs1 (Probs1(:,1)-Probs1(:,2)) Y_test predicted_Y1];
    %    Prob=Met;
    %     Prob(:,1)=exp(Met(:,1));
    %     Prob(:,2)=exp(Met(:,2));
    %
    %     Prob(:,1)=exp(Met(:,1))./(exp(Met(:,1))+exp(Met(:,2)));;
    %     Prob(:,2)=exp(Met(:,2))./(exp(Met(:,1))+exp(Met(:,2)));;
    %     Met=Prob;
    %     WrongMAt=Met(find(Y_test~=predicted_Y1),:);
    %     GoodMat=Met(find(Y_test==predicted_Y1),:);
    
    
    model=best_model1;
    Y_test;
    predicted_Y = best_pred;
    Accuracy(i)=(1-err/length(Y_test));
    PErrVec(i)=FPerr/HTot;
    FErrVec(i)=FNerr/STot;
end
fprintf('****************************************************************\n \n');
fprintf('Running using max values of differences between -Healthy+patients for %d features\n',limit);
fprintf('Accuracy = %.2f with std=%.2f Minimum=%.2f Max=%.2f \n False positives=%.2f \n Missed Patients=%.2f\n',mean(Accuracy),std(Accuracy,1),min(Accuracy),max(Accuracy),mean(PErrVec),mean(FErrVec));

