function [A, B, prec] = getdata_L(n, rs)

rand('seed',rs); randn('seed',rs);
ZERO = 1.0e-4; % small perturbation
signoise = .15; % Signal to noise ratio
% Algorithm ....
prec = 1e-1;          % Numerical precision
nnzr =(n*(n+1)*.5-n)*.01; % Number of nonzer coefficients in A
% Generate A, the original inverse covariance, with random sparsity pattern...
A=eye(n);
for i=1:nnzr
    A(ceil(n*rand),ceil(n*rand))=sign(rand()-.5);
end
A=A'*A;
B=2*(rand(n,n)-.5*ones(n,n)); % Add noise...
B=(B+B')*.5*signoise+inv(A+ZERO*eye(n));
B=B-min([min(eig(B))-ZERO,0])*eye(n); % (Set min eig > 0 if noise level too high)
B=(B+B')*.5;

end
