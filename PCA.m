function [train_X,train_Y,test_X,test_Y]=PCA(internal)
[train_X,train_Y,test_X,test_Y] = generateBalancedDataSet(feature,k);

if internal==0,
    prompt = 'Please enter number of Principal component that you want to have(less than 304): ';
    principal = input(prompt);
    X_Train_pca= ppca(train_X,principal);
    X_Test_pca = ppca(X_Test,principal);
else 
    
    X_Train_pca = getPCA(train_X,k);
    X_Test_pca = getPCA(test_X,k);
end
end
function Z=getPCA(X,k)
% this is good for using on the 
% X is the feature matrix
% m is number of instances
% k is max features

Sigma = cov(X);
[U, S, V] = svd(Sigma);
U=U(:,1:k);
%for i=1:n
Z=U'*X';
%fprintf('Using k=%d we have %f of variance retained \n',k,sum(S(1:k,:))/sum(S(:)));
    %end
%plot(Z)
%Z=tZ; %(:,1:k);
end