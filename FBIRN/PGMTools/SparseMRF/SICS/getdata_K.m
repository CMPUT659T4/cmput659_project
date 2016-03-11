function [A, Cstart, fstart, Wstart, EC, Sbase, lambda, K, EmpCov] = getdata_K(n, N, nzfrac, rho, rs)

rand('seed',rs); randn('seed',rs);

nnzr = nzfrac*n*n;          % Number of nonzero coefficients in A
% Generate A, the original inverse covariance, with random sparsity pattern...
A=eye(n);
for i=1:nnzr
    A(ceil(n*rand),ceil(n*rand))=sign(rand()-.5);
end
A=A*A'; % A is the gound truth inverse covariance matrix
B = inv(A); % B is the ground-truth covariance matrix
B=(B+B')/2; % B = B + 1e-6*eye(p); B = (B+B')/2;
data = mvnrnd(zeros(N,n),B);
EmpCov = (1/N)*data'*data;
Cstart=eye(n);
Wstart=eye(n);
%     N = 5*n;
EC=N*0.5*EmpCov;
fstart= - trace(EC*Cstart);
K=N*0.5;
Sbase=ones(n,n);
lambda = rho;
lambdaold=0;
% S=Sbase*lambda;
fstart=fstart -sum(sum(abs((lambda-lambdaold)*Sbase.*Cstart)));

end
