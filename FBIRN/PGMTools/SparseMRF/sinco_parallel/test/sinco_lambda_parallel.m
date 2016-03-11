function [Csol, Wsol, fsol] = sinco_lambda_parallel(EmpCov, N, Sbase, lstart, ifvec, precision, tol)
%  call sequence to sinco with fixed lambda
% 
% Csol - empirical covariance matrix
% lambda - reg parameter found
 
%lambda = rho*N/2, i.e. rho = (2/N)*lambda where rho is used by Banerjee's formulation

lambda = lstart;
p = size(EmpCov,1);
 
%% Katya's COVSEL settings
Cstart=eye(p);
Wstart=eye(p);
EC=N*0.5*EmpCov;
fstart= - trace(EC*Cstart);
K=N*0.5;
if (ifvec)
 lambdaold=zeros(p,1);
 lambdaouter=zeros(p,1);
else
  lambdaold=0;
  lambdaouter=0;  
end
mult=1;
% solve the problem once
if (ifvec) 
 for i=1:p
  S(:,i)=Sbase(:,i)*lambda(i);
  fstart=fstart - sum(abs(lambda(i)*Sbase(:,i).*Cstart(:,i))-...
                       abs(lambdaold(i)*Sbase(:,i).*Cstart(:,i)));
  %fstart=fstart - sum(abs((lambda(i)-lambdaold(i))*Sbase(:,i).*Cstart(:,i)));
 end
else 
  S=Sbase*lambda;
  fstart=fstart - sum(sum(abs((lambda-lambdaold)*Sbase.*Cstart)));
end
% [Csol, Wsol, fsol]=primal_cd_invcov(Cstart, fstart, Wstart, EC, Sbase, p, K);
%[Csol, Wsol, fsol]=...
%        sinco(Cstart, fstart, Wstart, EC, S, mult, K, tol);
[Csol, Wsol, fsol]=...
        sinco_matlab_parallel(Cstart, fstart, Wstart, EC, S, mult, K, tol);

