% Test script for Covariance selection code.
% Initialize parameters ****************
% Data ....


p=100;  % number of variables
N = 5000; % number of samples
nzfrac = 0.01;
nnzr=nzfrac*p*p;          % Number of nonzero coefficients in A
%nnzr=100;

precision = 0.005;
lstart = 10;


% Form random test matrix 
e=ones(p,p);
rand('seed',20);
% Generate A, the original inverse covariance, with random sparsity pattern...
A=eye(p);
for i=1:nnzr
    A(ceil(p*rand),ceil(p*rand))=sign(rand()-.5);
end
A=A*A';
% if want to use scale free, download matrix:
%A=load('SFNetworks/mat2net1sp21ml13.txt');

factor=1;
Atmp=zeros(factor*p, factor*p);
for i=1:factor
  Atmp((i-1)*p+1:i*p,(i-1)*p+1:i*p)=A;
end
A=Atmp;
N=(factor)*N;
p=factor*p;

% A=(A'+A)/2;  % A is the gound truth inverse covariance matrix
B = inv(A); % B is the ground-truth covariance matrix
B=(B+B')/2;
data = mvnrnd(zeros(N,p),B);

EmpCov = (1/N)*data'*data;

b= sum(sum(abs(inv(EmpCov + 0.001*eye(p,p)))))/(p*p-1);
%b=N*b/2
%b=27;
 true_pos = [];
 true_neg = [];
% % total # of positives (1s) and negatives (zeros)
 total_zeros = size(find(A == 0),1);
 total_ones = p*p - total_zeros;


lambda_range=[40:-1:1];
lambda_range=[300, 250, 200, 150, 100, lambda_range, 0.5, 0.1, 0.05, 0.01];
lambda_range=(factor)*lambda_range;
% lambda_range=[300];
rho_range = lambda_range*2/N;
tol=0.000001;
maxit=100000;
Cstart=eye(p);
Wstart=eye(p);
EC=N*0.5*EmpCov;
%load ECmat
fstart= - trace(EC*Cstart);
K=N*0.5;
%Sbase=rand(p,p);

%Sbase=Sbase+Sbase';
Sbase=ones(p,p);
%Sbase=Sbase-eye(p);
%for i=1:p
%  Sbase(i,i)=0.0001;
%end
lambdaold=0;

solution_path_l = zeros(size(lambda_range,2),p,p);
COVSEL_true_pos_l=zeros(1,length(lambda_range));
COVSEL_false_pos_l=zeros(1,length(lambda_range));
figure(1);
%subplot(1,4,1);colorspy(A);xlabel('true inverse cov');pbaspect('manual');set(gca,'XTickLabel',{' '});set(gca,'YTickLabel',{' '});
%subplot(1,4,3);colorspy(inv(EmpCov));xlabel('Inverse of EmpCov');pbaspect('manual');set(gca,'XTickLabel',{' '});set(gca,'YTickLabel',{' '});

% % %%%%%%%%%%% Test alternating minimization for \rho selection
% [C,our_lambda] = AltMin(EmpCov, N, Sbase, b, precision, tol, lstart);
% our_lambda
% keyboard
i=0;
t=cputime;
tic
 for lambda=lambda_range
    i  = i+1
    endi=i;
    % Test SINCO
    S=Sbase*lambda;
    fstart=fstart -(lambda-lambdaold)*sum(sum(abs(Sbase.*Cstart)));

%     [Csol, Wsol, fsol]=...
%        sinco(Cstart, fstart, Wstart, EC, Sbase, lambda, K, tol);
tic
[Csol, Wsol, fsol]=...
        sinco_matlab(Cstart, fstart, Wstart, EC, S, K, tol, maxit);
toc
    keyboard
    Umat = Csol;
    Cstart=Csol;
    fstart=fsol;
%     fall(i)=fsol;
%     lambdac(i)=lambda*sum(sum(abs(Sbase.*Csol)));
%     normC(i)=sum(sum(abs(Sbase.*Csol)));
%     tra(i)=trace(EC*Csol);
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
    solution_path_l(i, :, :)=Umat;
    true_pos = size(find(Umat & A),1);
    true_neg = size(find(~Umat  & ~A),1);
    false_pos=size(find(Umat  & ~A),1);
    tp(i)=true_pos;
 %   total_pos(i) = size(find(Umat),1);
 %   COVSEL_acc(i) =  (true_pos+true_neg)/(p*p);   
    COVSEL_true_pos_l(i) = true_pos/total_ones;
 %   COVSEL_true_neg(i) = true_neg/total_zeros;
    COVSEL_false_pos_l(i) = false_pos/total_zeros;
    
 %   COVSEL_probab_correct(i) =  COVSEL_true_pos_l(i) * total_ones/(p*p)  +  COVSEL_true_neg(i) *  total_zeros/(p*p);
 
    %if (true_pos/total_ones==1) 
    %  break
    %end
    rho = lambda*2/N;
%     LogLik_objective(i) = 0.5*N*log(det(Umat)) - 0.5*N*trace(Umat*EmpCov) ...
%                +p*p* log(lambda/2) - lambda*sum(sum(abs(Sbase.*Umat))) -...
%                 lambda*sum(abs(diag(Umat)));
% LogLik_objective(i) = 0.5*N*log(det(Umat)) - 0.5*N*trace(Umat*EmpCov) ...
%                +p*p* log(lambda/2) - lambda*sum(sum(abs(Sbase.*Umat))) -...
%                 lambda*b;
    %LogLik_objective(i) = log(det(Umat)) -(1/N)*trace(Umat*EmpCov)+(2/N)*p*p* log(lambda/2) - (2/N)lambda*sum(sum(abs(Umat)));
    %LogLik_objective(i) = log(det(Umat))-(1/N)*trace(Umat*EmpCov)+(2/N)*p*p* log(rho*N/4) - rho*sum(sum(abs(Umat)));
    %COVSEL_objective(i) =  0.5*N* (log(det(Umat)) - trace(Umat*EmpCov) - rho*sum(sum(abs(Umat))));
    %objective<-  - 0.5 * mm + dd * 0.5 * log ( det(a$wi)  )+ ss*ss* log (lambda/2)  - lambda * norm1 ;

%     solution_path(i,:,:) = Umat;
    

 end
time_lambda=cputime-t;
toc
%keyboard
% % %%%%%%%%%%% Test alternating minimization for \rho selection
% [our_C,our_lambda] = AltMin(EmpCov, N, b, precision, lstart);
% our_lambda
% our_LogLik_objective = 0.5*N*log(det(our_C)) - 0.5*N*trace(our_C*EmpCov) ...
%      +p*p* log(our_lambda) - our_lambda*sum(sum(abs(our_C))) - our_lambda*b;

% subplot(1,4,4);colorspy(C);xlabel('AltMIn Solution');pbaspect('manual');set(gca,'XTickLabel',{' '});set(gca,'YTickLabel',{' '});
 
% COVSEL_objective_our_rho =  0.5*N* (  log(det(Umat)) - trace(Umat*EmpCov) - our_rho*sum(sum(abs(Umat))) );

%our_true_pos = size(find(our_C & A),1);
%ur_true_neg = size(find(~our_C  & ~A),1);
%     
% COVSEL_acc_our_rho =  (true_pos+true_neg)/(p*p)  
% COVSEL_true_pos_our_rho = true_pos/total_ones
% COVSEL_true_neg_our_rho = true_neg/total_zeros
    

%figure(1);
%[a,ind] = max(COVSEL_true_pos+COVSEL_true_neg);
%Umat = squeeze(solution_path(ind,:,:));
%figure(1);
%subplot(1,4,2);colorspy(Umat);xlabel(sprintf('COVSEL Solution (lambda = %.1f)',lambda_range(ind)));pbaspect('manual');set(gca,'XTickLabel',{' '});set(gca,'YTickLabel',{' '});
%subplot(1,4,4);colorspy(our_C);xlabel(sprintf('AltMIn Solution (lambda = %.1f)',our_lambda));pbaspect('manual');set(gca,'XTickLabel',{' '});set(gca,'YTickLabel',{' '});
  plot(COVSEL_false_pos_l(1:endi), COVSEL_true_pos_l(1:endi)', 'r-x') 
  
    
% figure(2);
% plot(lambda_range,COVSEL_true_pos , 'b*-',lambda_range, COVSEL_false_pos,'ro-', lambda_range, COVSEL_acc, 'k-');
% xlabel('lambda','fontsize',16);
% legend('true pos', 'true neg', 'prob(correct)');
% [s,e]=sprintf('COVSEL true pos and neg (%d samples)', N); 
% title(s,'fontsize',12);
% 

% figure(3);
% plot(lambda_range,LogLik_objective, 'b*-',our_lambda,our_LogLik_objective,'ro');
% xlabel('lambda','fontsize',16);
% ylabel('log likelihood');
% [s,e]=sprintf('Loglikelihood objective (%d samples)', N); 
% title(s,'fontsize',12);
%      
     
 
%     
     



