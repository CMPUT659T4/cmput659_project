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

p=200;  % number of variables
N = 100; % number of samples
nzfrac = 0.01;
nnzr=nzfrac*p*p;          % Number of nonzero coefficients in A
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
% A=load('mat1net1sp21ml13.txt');
A=(A'+A)/2;  % A is the gound truth inverse covariance matrix
B = inv(A); % B is the ground-truth covariance matrix
B=(B+B')/2;
data = mvnrnd(zeros(N,p),B);

EmpCov = (1/N)*data'*data;

b= sum(sum(abs(inv(EmpCov + 0.001*eye(p,p)))))/(p*p-1);
%b=N*b/2
%b=27;
true_pos = [];
true_neg = [];
% total # of positives (1s) and negatives (zeros)
total_zeros = size(find(A == 0),1);
total_ones = p*p - total_zeros;

%%% old way - COVSEl's sampling
% % % B=2*(rand(p,p)-.5*ones(p,p)); % Add noise...
% % % B=(B+B')*.5*signoise+inv(A);
% % % B=B-min([min(eig(B)),0])*eye(n); % (Set min eig > 0 if noise level
% too high)


%rho_range=[0.6:0.01:0.6];
lambda_range=[300:-10:10];
nlambdas=length(lambda_range);
%lambda_range=[500:50:1000];
rho_range = lambda_range*2/N;
tol=0.000001;
Cstart=eye(p);
Wstart=eye(p);
Wstart_loo=zeros(p,p,N);
Cstart_loo=Wstart_loo;
Umat=zeros(p,p,N);
EC=N*0.5*EmpCov;
%load ECmat

% for loo_j=1:N
%   fstart_loo(loo_j)=-0.5*trace([data(1:loo_j-1,:);  data(loo_j+1:N,:)]'*...
%              [data(1:loo_j-1,:);  data(loo_j+1:N,:)]*Cstart);
%   Wstart_loo(:, :, loo_j)=eye(p);
%   Cstart_loo(:, :, loo_j)=eye(p);
% end
%Sbase=rand(p,p);

%Sbase=Sbase+Sbase';
Sbase=ones(p,p);
%Sbase=Sbase-eye(p);
%for i=1:p
%  Sbase(i,i)=0.0001;
%end
lambdaold=0;

solution_path = zeros(size(rho_range,2),p,p);
fall=zeros(nlambdas, p);
%subplot(1,4,1);colorspy(A);xlabel('true inverse cov');pbaspect('manual');set(gca,'XTickLabel',{' '});set(gca,'YTickLabel',{' '});
%subplot(1,4,3);colorspy(inv(EmpCov));xlabel('Inverse of EmpCov');pbaspect('manual');set(gca,'XTickLabel',{' '});set(gca,'YTickLabel',{' '});

% % %%%%%%%%%%% Test alternating minimization for \rho selection
%[C,our_lambda] = AltMin(EmpCov, N, Sbase, b, precision, tol, lstart);
% our_lambda
i=0;
fstart= -0.5*trace(data(1:N-1,:)'*data(1:N-1,:)*Cstart);
 for lambda=lambda_range
    i  = i+1
    S=Sbase*lambda;
    for loo_j=1:N
        loo_j
      EC_loo=0.5*[data(1:loo_j-1,:);  data(loo_j+1:N,:)]'*...
             [data(1:loo_j-1,:);  data(loo_j+1:N,:)];
    K=(N-1)/2;         
    % Call Sinco
    if loo_j==1 
      fstart=fstart -(lambda-lambdaold)*sum(sum(abs(Sbase.*Cstart)));
      fstart=fstart +0.5*trace(data(loo_j,:)'*data(loo_j,:)*Cstart)-...
          0.5*trace(data(N,:)'*data(N,:)*Cstart);
    else
      fstart=fstart +0.5*trace(data(loo_j,:)'*data(loo_j,:)*Cstart)-...
          0.5*trace(data(loo_j-1,:)'*data(loo_j-1,:)*Cstart);
    end
    [Csol, Wsol, fsol]=...
        sinco(Cstart, fstart, Wstart, ...
              EC_loo, Sbase, lambda, K, tol);
%     [Csol1, Wsol1, fsol1]=primal_cd_invcov(Cstart_loo(:,:,loo_j), fstart_loo(loo_j), Wstart_loo(:,:, loo_j), ...
%               EC_loo, S,  p, K);
    ftrue=K*log(det(Csol))-trace(EC_loo*Csol)-lambda*sum(sum(abs(Csol)));
%     if (abs(ftrue-fsol)>1.0e-6)
%         ftrue
%         fsol
%         keyboard
%     end
    Umat(:,:, loo_j)=Csol;
    Cstart=Csol;
    Wstart=Wsol;
    fstart=fsol;
    fall(i,loo_j)=fsol;
%     lambdac(i)=lambda*sum(sum(abs(Sbase.*Csol)));
%     normC(i)=sum(sum(abs(Sbase.*Csol)));
%     tra(i)=trace(EC*Csol);
%    Wstart_loo(:, :, loo_j)=Wsol;     
    end
    lambdaold=lambda; 
    Umatnnz=zeros(p,p);
    for loo_j=1:N
      Umatnnz=Umatnnz+(abs(Umat(:,:,loo_j))>1.0e-8);
    end
    Umatnnz=Umatnnz/N;
    Umatnnz=Umatnnz-diag(diag(Umatnnz));
    consist(i)=nnz(Umatnnz>0.9);
    total(i)=nnz(Umatnnz);
    robus(i)=consist(i)/total(i); 
end
plot(lambda_range, robus)
     



