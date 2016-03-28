function PC=PCA(internal,feature,k)
comp = struct('coeff', [],'score', [],'pcvar' ,[]);
h = comp;
s = comp;
PC = struct('h', [],'s', []);
fold=generateDataSetBalanced(feature,k);

train_hel = cat(1, fold(1).fold_hel ,fold(2:end).fold_hel);
%train_hel =train_hel;
%test_X_hel = fold(1).fold_hel;
%test_X_hel = test_X_hel';

train_sch = cat(1,fold(1).fold_sch,fold(2:end).fold_sch);
%train_sch =train_sch';
%test_X_sch = fold(1).fold_sch;
%test_X_sch = test_X_sch';
if internal~=1000,
    size(train_hel)
    mini = min(size(train_hel));
    prompt = sprintf('Please enter number of Principal component that you want to have(less than %d): ',mini);
    principal = input(prompt);
    [coeff,score,pcvar]= ppca(train_hel,principal);
    PC.h.coeff= coeff;
    PC.h.score = score;
    PC.h.pcvar = pcvar;
    %PC.test_X_hel = ppca(test_X_hel,principal);
    [coeff,score,pcvar]= ppca(train_sch,principal);
    PC.s.coeff= coeff;
    PC.s.score = score;
    PC.s.pcvar = pcvar;
    %PC.test_X_sch = ppca(test_X_sch,principal);
else 
    
    [h_coeff,score,h_pcvar] = getPCA(train_hel,k);
    PC.train_hel= [h_coeff,score,h_pcvar];
    %PC.test_X_hel = getPCA(test_X_hel,k);
    [s_coeff,s_score,s_pcvar] = getPCA(train_sch,k);
    PC.train_sch = [s_coeff,s_score,s_pcvar];
    %PC.test_X_sch = getPCA(test_X_sch,k);
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