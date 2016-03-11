function [Csol,lambda] = AltMin_vec(EmpCov, N, Sbase, ifvec, b, usediag, precision, tol, lstart)
%  alternating minimization algorithm for finding inv cov matrix AND lambda
% Input:
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
[Csol, Wsol, fsol]=...
        sinco(Cstart, fstart, Wstart, EC, S, mult, K, tol);
i=1;
Cstart=Csol;

fstart=fsol;
Wstart=Wsol;

%my_fsol = 0.5*N*log(det(Csol)) - 0.5*N*trace(Csol*EmpCov) - sum(sum(abs(S.*Csol)))
 fall(i)=fsol;

if (~usediag)
  if (ifvec)
     obj = fsol + sum(p*log(lambda/2)) - b*lambda;
  else
     obj  = fsol +p*p*log(lambda/2) - lambda*b;
  end
else
  if (ifvec)
   obj = fsol + sum(p*log(lambda/2)) - diag(Csol)'*lambda;
  else
   obj  = fsol +p*p*log(lambda/2) - lambda*sum(diag(Csol))
  end
end    
    
    
while norm(lambda - lambdaouter)/norm(lambda) > precision
       
       lambdaold = lambda;  % previous lambda for line search
       lambdaouter = lambda; % previous lambda for alternating minimization
       
       % fix Csol, solve for lambda; f(rho) = const + p^2*log rho - rho ||C||_1,
       % where p is the number of variables, i.e. size(1,C)
       % thus df/d lambda = p/lambda_i - ( ||C||_S_i + b_i) = 0 implies lambda_i =
       % p/(||C||_S_i + b_i)
       if (~usediag)   
         if (ifvec) 
           lambda  =  (size(Csol,1))./((sum(abs(Sbase.*Csol))) + b )';
         else 
            lambda= (size(Csol,1)^2)/( sum(sum(abs(Sbase.*Csol))) + b ); 
         end
       else
         if (ifvec) 
           lambda  =  (size(Csol,1))./((sum(abs(Sbase.*Csol))) + diag(Csol)' )';
         else 
           lambda = (size(Csol,1)^2)/( sum(sum(abs(Sbase.*Csol))) + sum(diag(Csol)) );
         end
       end
       % steps back if needed
       alpha = 1;
       while alpha > precision
       % fix lambda, solve for C      
       % [C,X] = spmlcdvec(CovEmp,rho,maxiter,prec,maxnest,algot); % Inverse covariance is stored in C
          if (ifvec)   
            for i=1:p
               S(:,i)=Sbase(:,i)*lambda(i);
               fstart=fstart - sum(abs(lambda(i)*Sbase(:,i).*Cstart(:,i))-...
                       abs(lambdaold(i)*Sbase(:,i).*Cstart(:,i)));
            end 
          else
            S=Sbase*lambda;
            fstart=fstart - sum(sum(abs((lambda-lambdaold)*Sbase.*Cstart)));
          end
%            my_fstart =  0.5*N*log(det(Cstart)) - 0.5*N*trace(Cstart*EmpCov) - lambda*sum(sum(abs(Cstart)))
            
          [Csol, Wsol, fsol]=sinco(Cstart, fstart, Wstart, EC, S, mult, K, tol);
          
%          my_fsol = 0.5*N*log(det(Csol)) - 0.5*N*trace(Csol*EmpCov) - sum(sum(abs(S.*Csol))) 
            
%              lambda
%              lambdaouter
%              lambdaold
 %            fsol
%             my_fsol = 0.5*N*log(det(Csol)) - 0.5*N*trace(Csol*EmpCov) - lambda*sum(sum(abs(Csol)))   

            Cstart=Csol;
            fstart=fsol;
            i = i+1;
            fall(i)=fsol;
            Wstart=Wsol;
            lambdaold = lambda;
            if (~usediag)
              if (ifvec)
                obj_new = fsol + sum(p*log(lambda/2)) - b*lambda;
              else
                obj_new  = fsol +p*p*log(lambda/2) - lambda*b;
              end
            else
              if (ifvec)
                 obj_new = fsol + sum(p*log(lambda/2)) - diag(Csol)'*lambda;
              else
                 obj_new  = fsol +p*p*log(lambda/2) - lambda*sum(diag(Csol))
              end
            end    
%            LogLik_obj = 0.5*N*log(det(Csol)) - 0.5*N*trace(Csol*EmpCov) - lambda*sum(sum(abs(Csol)))  +p*p* log(lambda/2)  - lambda*b

            if obj_new < obj + precision  
                alpha = alpha/2;
                lambda = (lambda+lambdaouter)/2;
            else
                obj = obj_new;   
                break;
            end
                
       end  
end
% for i=1:p
%    S(:,i)=Sbase(:,i)*lambda(i);
% %   fstart=fstart - sum(abs((lambda(i)-lambdaold(i))*Sbase(:,i).*Cstart(:,i)));
%    fstart=fstart - sum(abs(lambda(i)*Sbase(:,i).*Cstart(:,i))-...
%                        abs(lambdaold(i)*Sbase(:,i).*Cstart(:,i)));
% end    
% [Csol, Wsol, fsol]=sinco(Cstart, fstart, Wstart, EC, S, mult, K, tol);
        

 
