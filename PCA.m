function [train_X_pca,test_X_pca]=PCA(internal,feature,k)
fold=generateDataSetBalanced(feature,k);
train_X = fold(2:k).fold;
train_X =train_X';
test_X = fold(1).fold;
test_X = test_X';
if internal==0,
    mini = min(size(train_X));
    prompt = sprintf('Please enter number of Principal component that you want to have(less than %d): ',mini);
    principal = input(prompt);
    train_X_pca= ppca(train_X,principal);
    test_X_pca = ppca(test_X,principal);
else 
    
    train_X_pca = getPCA(train_X,k);
    test_X_pca = getPCA(test_X,k);
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