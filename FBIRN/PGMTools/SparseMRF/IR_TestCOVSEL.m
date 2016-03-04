% Test script for Covariance selection code.
% Initialize parameters ****************
% Data ....
% % % signoise=.15; % Signal to noise ratio
% % % 
% % % nnz=8;          % Number of nonzer coefficients in A
% % % % Algorithm ....
% % % rho=0.6;            % Controls sparsity
% % % prec=1e-20;         % Numerical precision
% % % maxiter=2;         % Maximum Number of sweeps in coordinate descent
% % % algot='nest';        % Algorith for solving BoxQP: 'sedumi' uses SEDUMI, 'nest' uses smooth minimization (usually faster)
% % % maxnest=1000;   % Maximum number of iterations in the smooth BoxQP solver

p=50;  % number of variables
N = 100; % number of samples
nzfrac = 0.05;
nnz=nzfrac*p*p;          % Number of nonzer coefficients in A


precision = 0.05;
lstart = 30;


% Form random test matrix 
e=ones(p,p);
rand('seed',20);
% Generate A, the original inverse covariance, with random sparsity pattern...
A=eye(p);
for i=1:nnz
    A(ceil(p*rand),ceil(p*rand))=sign(rand()-.5);
end
A=A'*A;  % A is the gound truth inverse covariance matrix

B = inv(A); % B is the ground-truth covariance matrix
p=100;  % number of variables
N = 30; % number of samples
nzfrac = 0.04;% fraction of nonzeros
nnz=nzfrac*p*p;          % Number of nonzer coefficients in A


precision = 0.05;
lstart = 200;


% Form random test matrix 
e=ones(p,p);
rand('seed',20);
% Generate A, the original inverse covariance, with random sparsity pattern...
A=eye(p);
for i=1:nnz
    A(ceil(p*rand),ceil(p*rand))=sign(rand()-.5);
end
A=A'*A;  % A is the gound truth inverse covariance matrix

B = inv(A); % B is the ground-truth covariance matrix
data = mvnrnd(zeros(N,p),B);
data = mvnrnd(zeros(N,p),B);

EmpCov = (1/N)*data'*data;
b= sum(sum(abs(EmpCov + 0.00001*eye(p,p))))/(p*p);

true_pos = [];
true_neg = [];
% total # of positves (1s) and negatives (zeros)
total_zeros = size(find(A == 0),1);
total_ones = p*p - total_zeros;

%%% old way - COVSEl's sampling
% % % B=2*(rand(p,p)-.5*ones(p,p)); % Add noise...
% % % B=(B+B')*.5*signoise+inv(A);
% % % B=B-min([min(eig(B)),0])*eye(n); % (Set min eig > 0 if noise level
% too high)

i=0;
%rho_range=[0.6:0.01:0.6];
lambda_range=[20:-0.5:5];
rho_range = lambda_range*2/N;

Cstart=eye(p);
Wstart=eye(p);
EC=N*0.5*EmpCov;
fstart= - trace(EC*Cstart);
K=N*0.5;
%Sbase=rand(p,p);
%Sbase=Sbase+Sbase';
Sbase=ones(p,p);
lambdaold=0;

solution_path = zeros(size(rho_range,2),p,p);

figure(1);
subplot(1,4,1);colorspy(A);xlabel('true inverse cov');pbaspect('manual');set(gca,'XTickLabel',{' '});set(gca,'YTickLabel',{' '});
subplot(1,4,3);colorspy(inv(EmpCov));xlabel('Inverse of EmpCov');pbaspect('manual');set(gca,'XTickLabel',{' '});set(gca,'YTickLabel',{' '});

% % %%%%%%%%%%% Test alternating minimization for \rho selection
 [C,our_lambda] = AltMin(EmpCov, N, b, precision, lstart);
 our_lambda
 
 for lambda=lambda_range
    i  = i+1
    % Test primal COVSEL code (katya's)
    S=Sbase*lambda;
    fstart=fstart -sum(sum(abs((lambda-lambdaold)*Sbase.*Cstart)));

    [Csol, Wsol, fsol]=primal_cd_invcov(Cstart, fstart, Wstart, EC, S, p, K);
    Umat = Csol;
    
    Cstart=Csol;
    fstart=fsol;
    fall(i)=fsol;
    Wstart=Wsol;
    lambdaold=lambda;
    
   %  COVSEL
   % [Umat,X] = spmlcdvec(EmpCov,rho,maxiter,prec,maxnest,algot); % Inverse covariance is stored in Umat
  
  % alternatively - Yuanqing's implementation of projected gradient
  % [Umat,W]=projGrad_SparseGMRF(EmpCov,rho,prec,maxiter);
%     figure(20);
%     hist(log(abs(Umat)))
% 
% 
    true_pos = size(find(Umat & A),1);
    true_neg = size(find(~Umat  & ~A),1);
 
    COVSEL_acc(i) =  (true_pos+true_neg)/(p*p);   
    COVSEL_true_pos(i) = true_pos/total_ones;
    COVSEL_true_neg(i) = true_neg/total_zeros;
    COVSEL_probab_correct(i) =  COVSEL_true_pos(i) * total_ones/(p*p)  +  COVSEL_true_neg(i) *  total_zeros/(p*p);
 

    
    rho = lambda*2/N;
    LogLik_objective(i) = 0.5*N*log(det(Umat)) - 0.5*N*trace(Umat*EmpCov) +p*p* log(lambda/2) - lambda*sum(sum(abs(Umat))) -lambda*b;
    %LogLik_objective(i) = log(det(Umat)) -(1/N)*trace(Umat*EmpCov)+(2/N)*p*p* log(lambda/2) - (2/N)lambda*sum(sum(abs(Umat)));
    %LogLik_objective(i) = log(det(Umat))-(1/N)*trace(Umat*EmpCov)+(2/N)*p*p* log(rho*N/4) - rho*sum(sum(abs(Umat)));
    %COVSEL_objective(i) =  0.5*N* (log(det(Umat)) - trace(Umat*EmpCov) - rho*sum(sum(abs(Umat))));
    %objective<-  - 0.5 * mm + dd * 0.5 * log ( det(a$wi)  )+ ss*ss* log (lambda/2)  - lambda * norm1 ;

    solution_path(i,:,:) = Umat;
    

end

% % %%%%%%%%%%% Test alternating minimization for \rho selection
 [our_C,our_lambda] = AltMin(EmpCov, N, b, precision, lstart);
 our_lambda
 our_LogLik_objective = 0.5*N*log(det(our_C)) - 0.5*N*trace(our_C*EmpCov) +p*p* log(lambda/2) - lambda*sum(sum(abs(our_C))) - lambda*b;

% subplot(1,4,4);colorspy(C);xlabel('AltMIn Solution');pbaspect('manual');set(gca,'XTickLabel',{' '});set(gca,'YTickLabel',{' '});
 
% COVSEL_objective_our_rho =  0.5*N* (  log(det(Umat)) - trace(Umat*EmpCov) - our_rho*sum(sum(abs(Umat))) );

our_true_pos = size(find(C & A),1);
our_true_neg = size(find(~C  & ~A),1);
%     
% COVSEL_acc_our_rho =  (true_pos+true_neg)/(p*p)  
% COVSEL_true_pos_our_rho = true_pos/total_ones
% COVSEL_true_neg_our_rho = true_neg/total_zeros
    

figure(1);
[a,ind] = max(COVSEL_true_pos+COVSEL_true_neg);
Umat = squeeze(solution_path(ind,:,:));
figure(1);
subplot(1,4,2);colorspy(Umat);xlabel(sprintf('COVSEL Solution (lambda = %.1f)',lambda_range(ind)));pbaspect('manual');set(gca,'XTickLabel',{' '});set(gca,'YTickLabel',{' '});
subplot(1,4,4);colorspy(our_C);xlabel(sprintf('AltMIn Solution (lambda = %.1f)',our_lambda));pbaspect('manual');set(gca,'XTickLabel',{' '});set(gca,'YTickLabel',{' '});
   
    
figure(2);
plot(lambda_range,COVSEL_true_pos , 'b*-',lambda_range, COVSEL_true_neg,'ro-', lambda_range, COVSEL_acc, 'k-');
xlabel('lambda','fontsize',16);
legend('true pos', 'true neg', 'prob(correct)');
[s,e]=sprintf('COVSEL true pos and neg (%d samples)', N); 
title(s,'fontsize',12);
        
figure(3);
plot(lambda_range,LogLik_objective, 'b*-',our_lambda,our_LogLik_objective,'ro');
xlabel('lambda','fontsize',16);
ylabel('log likelihood');
[s,e]=sprintf('Loglikelihood objective (%d samples)', N); 
title(s,'fontsize',12);
     
     
 
%     
     



