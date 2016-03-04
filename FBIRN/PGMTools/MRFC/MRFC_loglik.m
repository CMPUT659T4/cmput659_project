function loglik = MRFC_loglik(X, Y, model)
% Use Gaussian Markov Random Field model to compute the likelihood
% of test samples
% we assume the test samples are of the same class

class = Y(1);

[n,p]=size(X);
p=model.nvars; % number of variables/dimensions
k=model.nlabels; % number of class labels
py=model.class_prior;
cond = model.class_cond; % k-dimensional array of structs with fields C and mu

% select max-likelihood label:
% max_i P(Y=i|X=x) = max_i P(X=x,Y=i) label since P(Y=i|X) = P(X=x,Y=i)/P(X=x)

 CC = inv(cond(class).C;
 CC = (CC + CC')/2;
 loglik = log(mvnpdf(X,cond(class).mu, CC )); % multivar gaussian pdf
 

 