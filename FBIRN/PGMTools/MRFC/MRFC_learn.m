function model = MRFC_learn(X, y,  method, lambdas)
% Learn a Markov Random Field classifier for Gaussian data X
% learn a separate multivariate Gaussain P(X|Y=i) distribution for each 
% value i of the discrete class variable Y
% % method - algorithm used to find MRF (e.g. 'covsel', 'projected_grad',
% 'our_covsel'

% lambdas vector contains particular sparsity parameter for each class

[n,p] = size(X);
labels = unique(y);   n_labels = size(labels,1);

if nargin < 4
    lambdas=ones(n_labels,1);
end
   
if nargin < 3
    method='our_covsel';
end


for j=1:n_labels
    prior(j)=size(find(y==labels(j)),1); 
end
prior = prior/n;
       size(X) ;
       size(lambdas);
       lambdas=lambdas*ones(1,n_labels);
for j=1:n_labels
    ind = find(y==labels(j));
    ind;
    conditional(j) = learn_mrf(X(ind,:),method, lambdas(j));
end

 
model = struct('nvars',p,'nlabels',n_labels, 'labels', labels, 'class_cond',conditional,'class_prior',prior,'lambdas',lambdas);  
 