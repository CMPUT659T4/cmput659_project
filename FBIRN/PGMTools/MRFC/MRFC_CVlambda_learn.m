function model = MRFC_CVlambda_learn(X, y,  method, lambdas)
% Learn a Markov Random Field classifier for Gaussian data X
% learn a separate multivariate Gaussain P(X|Y=i) distribution for each 
% value i of the discrete class variable Y
% % method - algorithm used to find MRF (e.g. 'covsel', 'projected_grad',
% 'our_covsel', 

[n,p] = size(X);
labels = unique(y);   n_labels = size(labels,1);

if nargin < 3
    method='our_covsel';
end
for j=1:n_labels
    prior(j)=size(find(y==labels(j)),1);
end
prior = prior/n;
        
for j=1:n_labels
    ind = find(y==labels(j));
    
    % choose the lambda that yields highest-likelihood  
    best_lambda(j)=0; best_loglik(j) = -Inf;
    for lambda=lambdas       
       cond = learn_mrf(X(ind,:),method, lambda);
       CC = inv(cond.C);
       CC = (CC + CC')/2;
       current_loglik =  sum(log(mvnpdf(X,cond.mu, CC ))); 
       if current_loglik > best_loglik(j)
           best_loglik(j) = current_loglik;
           best_lambda(j) = lambda;
           cond_model(j) = cond;         
       end
    end

end

 
model = struct('nvars',p,'nlabels',n_labels, 'labels', labels, 'class_cond',cond_model,'class_prior',prior,'class_lambda',best_lambda,'class_loglik',best_loglik);  
 