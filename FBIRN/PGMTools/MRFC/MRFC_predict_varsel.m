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

pyx = zeros(n,k);
% compute P(X,Y=i)=P(Y=i)*P(X=x|Y=i) for each class lable i

% variable selection - delete independent nodes
 
C = cond(1).C;
C = C - diag(diag(C));
i1 = find(sum(C))'; % find connected nodes
clear C;

C = cond(2).C;
C = C - diag(diag(C));
i2 = find(sum(C))'; % find connected nodes
clear C;

inds = union(i1,i2);
 
for class=1:k
    C=cond(class).C(inds,inds); % use only connected variables
    mu=cond(class).mu(inds);
    XX = X(:,inds);
    
    CC = inv(C);
    CC = (CC + CC')/2;
    pyx(:,class) = py(class)*mvnpdf(XX,mu, CC ); % multivar gaussian pdf
end


[a,ind]=max(pyx');
y=zeros(n,1);
for i=1:n
    y(i)=model.labels(ind(i));
end
