function [y,pyx] = MRFC_predict(X, model)
% Use Gaussian Markov Random Field classification model (class-conditional
% models P(X|C) and class prior P(C)) to predict class labels of a given
% set of training instances X
% model = struct('nvars',..,'nlabels',.., 'class_cond',..,'class_prior',..);  

[n,p]=size(X);
p=model.nvars; % number of variables/dimensions
k=model.nlabels; % number of class labels
py=model.class_prior;
cond = model.class_cond; % k-dimensional array of structs with fields C and mu

% select max-likelihood label:
% max_i P(Y=i|X=x) = max_i P(X=x,Y=i) label since P(Y=i|X) = P(X=x,Y=i)/P(X=x)
CC = inv(cond(1).C);
    CC = (CC + CC')/2;
    cond(1).mu;
   mvnpdf(X,cond(1).mu, CC );
   X;
pyx = zeros(n,k);
py;
k;
% compute P(X,Y=i)=P(Y=i)*P(X=x|Y=i) for each class lable i
for class=1:k
    CC = inv(cond(class).C);
    CC = (CC + CC')/2;
    pyx(:,class) = log(py(class))+logmvnpdf(X,cond(class).mu, CC ); % multivar gaussian pdf
    %pyx(:,class) = logmvnpdf(X,cond(class).mu, CC ); % multivar gaussian pdf
end

size(pyx);
[a,ind]=max(pyx');

y=zeros(n,1);
for i=1:n
    y(i)=model.labels(ind(i));
end
