% Test script for Covariance selection code.
% Initialize parameters ****************
% Data ....
signoise=.15; % Signal to noise ratio
n=20;            % Dimension
nnz=8;          % Number of nonzer coefficients in A
% Algorithm ....
rho=0.6;            % Controls sparsity
prec=1e-1;         % Numerical precision
maxiter=2;         % Maximum Number of sweeps in coordinate descent
algot='nest';        % Algorith for solving BoxQP: 'sedumi' uses SEDUMI, 'nest' uses smooth minimization (usually faster)
maxnest=1000;   % Maximum number of iterations in the smooth BoxQP solver


% Form random test matrix 
e=ones(n,n);
rand('seed',20);
% Generate A, the original inverse covariance, with random sparsity pattern...
A=eye(n);
for i=1:nnz
    A(ceil(n*rand),ceil(n*rand))=sign(rand()-.5);
end
A=A'*A;
B=2*(rand(n,n)-.5*ones(n,n)); % Add noise...
B=(B+B')*.5*signoise+inv(A);
B=B-min([min(eig(B)),0])*eye(n); % (Set min eig > 0 if noise level too high)


% Test block COVSEL code
[Umat,X] = spmlcdvec(B,rho,maxiter,prec,maxnest,algot); % Inverse covariance is stored in Umat

% Plot results, comparing A, the solution Umat and the noisy inverse inv(B).
subplot(1,3,1);colorspy(A);xlabel('A');pbaspect('manual');set(gca,'XTickLabel',{' '});set(gca,'YTickLabel',{' '});
subplot(1,3,3);colorspy(inv(B));xlabel('Inverse');pbaspect('manual');set(gca,'XTickLabel',{' '});set(gca,'YTickLabel',{' '});
subplot(1,3,2);colorspy(Umat);xlabel('Solution');pbaspect('manual');set(gca,'XTickLabel',{' '});set(gca,'YTickLabel',{' '});