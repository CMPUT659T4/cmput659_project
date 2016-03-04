% Test script for Covariance selection code.
% Initialize parameters ****************
% Data ....
signoise=.15; % Signal to noise ratio
p=100;            % Dimension
nnz=8;          % Number of nonzer coefficients in A
% Algorithm ....
rho=0.6;            % Controls sparsity
prec=1e-1;         % Numerical precision
maxiter=2;         % Maximum Number of sweeps in coordinate descent
algot='nest';        % Algorith for solving BoxQP: 'sedumi' uses SEDUMI, 'nest' uses smooth minimization (usually faster)
maxnest=1000;   % Maximum number of iterations in the smooth BoxQP solver

N = 30; % number of samples

% 1. Generate a scale-free Gaussian Markov network
% 2. Sample N  instenaces

% use B-A Scale-Free Network Generation and Visualization package by by
% Mathew Neil George
% seed network
seed =[0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1; 
       1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 
       1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;  
       1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0;  
       1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;  
       1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0;  
       1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;  
       1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0;  
       1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;  
       1 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0;
       1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
       1 0 0 0 0 0 0 0 0 1 0 0 0 0 1 0;
       1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0;
       1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
       1 0 0 1 0 1 0 1 0 1 0 1 1 0 0 0;
       1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; 
% generate a scale-free network's adjacency matrix

% Generate A, the original inverse covariance, with random sparsity pattern...

Net = SFNG(p, 5, seed);
Net = Net | Net'; % make it symmetric

%maximum value for cov that does not violate pos-def constaint depends on p
const = 0.15; % for p=300; covariance value for each pair of nodes setting to 0.15 and larger violates positive-definiteness

% generate sparse inverse covariance matrix (capturing independence relationships
A = Net*const + eye(p);
%S = inv(S1);
%A = A'*A;
 
B = inv(A); % B is the ground-truth covariance matrix
data = mvnrnd(zeros(N,p),B);

EmpCov = (1/N)*data'*data;

true_pos = [];
true_neg = [];
% total # of positves (1s) and negatives (zeros)
total_zeros = size(find(A == 0),1);
total_ones = p*p - total_zeros;

%%% old way - COVSEl's sampling
% % % B=2*(rand(p,p)-.5*ones(p,p)); % Add noise...
% % % B=(B+B')*.5*signoise+inv(A);
% % % B=B-min([min(eig(B)),0])*eye(n); % (Set min eig > 0 if noise level too high)

 
i=0;
rho_range=[0.8:0.01:0.8];
lambda_range = rho_range*N/2;

for rho=rho_range
    i  = i+1;
    % Test block COVSEL code
   
       
  %  COVSEL
    [Umat,X] = spmlcdvec(EmpCov,rho,maxiter,prec,maxnest,algot); % Inverse covariance is stored in Umat
  
  % alternatively - Yuanqing's implementation of projected gradient
  % [Umat,W]=projGrad_SparseGMRF(EmpCov,rho,prec,maxiter);

    true_pos = size(find(Umat & A),1);
    true_neg = size(find(~Umat  & ~A),1);
 
    COVSEL_acc(i) =  (true_pos+true_neg)/(p*p);   
    COVSEL_true_pos(i) = true_pos/total_ones;
    COVSEL_true_neg(i) = true_neg/total_zeros;

    lambda = rho*N/2;
    
    LogLik_objective(i) = 0.5*N*log(det(Umat)) - 0.5*N*trace(Umat*EmpCov) +p*p* log(lambda/2) - lambda*sum(sum(abs(Umat)));
    %or, equivalently, 
    %LogLik_objective(i) = log(det(Umat)) - trace(Umat*EmpCov)+(2/N)*p*p* log(lambda/2) - (2/N)lambda*sum(sum(abs(Umat)));
    % alternative expression as a function of \rho to match COVSEL formulation
    % LogLik_objective(i) =
    %log(det(Umat))-trace(Umat*EmpCov)+(2/N)*p*p* log(rho*N/4) - rho*sum(sum(abs(Umat)));
    
    % 
    COVSEL_objective(i) =  log(det(Umat) - trace(Umat*EmpCov) - rho*sum(sum(abs(Umat))));

        
        
    %objective<-  - 0.5 * mm + dd * 0.5 * log ( det(a$wi)  )+ ss*ss* log (lambda/2)  - lambda * norm1 ;
% Plot results, comparing A, the solution Umat and the noisy inverse inv(B).
hist(log(abs(Umat)))

figure(2);
subplot(1,4,1);colorspy(A);xlabel('true inverse cov');pbaspect('manual');set(gca,'XTickLabel',{' '});set(gca,'YTickLabel',{' '});
subplot(1,4,3);colorspy(inv(EmpCov));xlabel('Inverse of EmpCov');pbaspect('manual');set(gca,'XTickLabel',{' '});set(gca,'YTickLabel',{' '});
subplot(1,4,2);colorspy(Umat);xlabel('Solution');pbaspect('manual');set(gca,'XTickLabel',{' '});set(gca,'YTickLabel',{' '});

end


% % 
% % %%%%%%%%%%% Test alternating minimization for \rho selection
 [C,our_rho] = AltMin(EmpCov, N, maxiter,prec,maxnest,algot);
 subplot(1,4,4);colorspy(C);xlabel('AltMIn Solution');pbaspect('manual');set(gca,'XTickLabel',{' '});set(gca,'YTickLabel',{' '});
 
COVSEL_objective_our_rho =  0.5*N* (  log(det(Umat)) - trace(Umat*EmpCov) - our_rho*sum(sum(abs(Umat))) );
our_rho
our_lambda = our_rho*N/2
true_pos = size(find(C & A),1);
true_neg = size(find(~C  & ~A),1);
    
COVSEL_acc_our_rho =  (true_pos+true_neg)/(p*p)  
COVSEL_true_pos_our_rho = true_pos/total_ones
COVSEL_true_neg_our_rho = true_neg/total_zeros
    
    
    
    figure(10);
    plot(rho_range,COVSEL_acc , 'b*-',rho_range,LogLik_objective, 'ro-');
    ylabel('objective','fontsize',16);
    xlabel('rho','fontsize',16);
    legend('accuracy', 'loglikelihood');
    [s,e]=sprintf('COVSEL accuracy (%d samples)', N); 
    title(s,'fontsize',12);
    
    figure(11);
    plot(lambda_range,COVSEL_acc , 'b*-',rho_range,LogLik_objective, 'ro-');
    ylabel('objective','fontsize',16);
    xlabel('lambda','fontsize',16);
    legend('accuracy', 'loglikelihood');
            [s,e]=sprintf('COVSEL true pos and neg (%d samples)', N); 
    title(s,'fontsize',12);
    
    
    figure(15);
    plot(rho_range,COVSEL_true_pos , 'b*-',rho_range, COVSEL_true_neg,'ro-');
    xlabel('rho','fontsize',16);
    legend('true pos', 'true neg');
        [s,e]=sprintf('COVSEL true pos and neg (%d samples)', N); 
    title(s,'fontsize',12);
        
    figure(20);
    plot(rho_range,LogLik_objective, 'b*-');
    xlabel('rho','fontsize',16);
     ylabel('log likelihood');
      [s,e]=sprintf('Loglikelihood objective (%d samples)', N); 
    title(s,'fontsize',12);
     
     
    figure(21);
    plot(lambda_range,LogLik_objective, 'b*-');
    xlabel('lambda','fontsize',16);
    ylabel('log likelihood');
    [s,e]=sprintf('Loglikelihood objective (%d samples)', N); 
    title(s,'fontsize',12);
    
     



