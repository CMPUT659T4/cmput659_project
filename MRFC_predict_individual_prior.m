function [y,pyx] = MRFC_predict_individual_prior(X,pre, model)
size(pre);

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

pyx = zeros(n,k);
size(pyx);
% compute P(X,Y=i)=P(Y=i)*P(X=x|Y=i) for each class lable i
for class=1:k
    CC = inv(cond(class).C);
    CC = (CC + CC')/2;
    pyx(:,class) = log(pre(:,class)')+logmvnpdf(X,cond(class).mu, CC ); % multivar gaussian pdf
    %pyx(:,class) = logmvnpdf(X,cond(class).mu, CC ); % multivar gaussian pdf
end

pyx;
[a,ind]=max(pyx');
y=zeros(n,1);
for i=1:n
    y(i)=model.labels(ind(i));
end
