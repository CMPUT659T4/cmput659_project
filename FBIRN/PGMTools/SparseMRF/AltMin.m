function   [Csol,lambda] = AltMin(EmpCov, N, b, precision, lstart)
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
Sbase=ones(p,p);
lambdaold=0;
lambdaouter=0;

% solve the problem once
S=Sbase*lambda;
fstart=fstart - sum(sum(abs((lambda-lambdaold)*Sbase.*Cstart)));
[Csol, Wsol, fsol]=primal_cd_invcov(Cstart, fstart, Wstart, EC, S, p, K);
 
i=1;
Cstart=Csol;

fstart=fsol;
%my_fsol = 0.5*N*log(det(Csol)) - 0.5*N*trace(Csol*EmpCov) - lambda*sum(sum(abs(Csol)))   


fall(i)=fsol;
Wstart=Wsol;

       
obj = fsol +p*p*log(lambda/2) - lambda*b;

    
while abs(lambda - lambdaouter) > precision
       
       lambdaold = lambda;  % previous lambda for line search
       lambdaouter = lambda; % previous lambda for alternating minimization
       
       % fix Csol, solve for lambda; f(rho) = const + p^2*log rho - rho ||C||_1,
       % where p is the number of variables, i.e. size(1,C)
       % thus df/d lambda = p^2/lambda - ( ||C||_S + b) = 0 implies lambda =
       % p^2/(||C||_1 + b)
             
       lambda  =  (size(Csol,1)^2)/( sum(sum(abs(S.*Csol))) + b );    

       % steps back if needed
       alpha = 1;
       while alpha > precision
       % fix lambda, solve for C      
       % [C,X] = spmlcdvec(CovEmp,rho,maxiter,prec,maxnest,algot); % Inverse covariance is stored in C
            
            S=Sbase*lambda;
            fstart=fstart -(lambda-lambdaold)*sum(sum(abs(Sbase.*Cstart)));
%            my_fstart =  0.5*N*log(det(Cstart)) - 0.5*N*trace(Cstart*EmpCov) - lambda*sum(sum(abs(Cstart)))
            
            [Csol, Wsol, fsol]=primal_cd_invcov(Cstart, fstart, Wstart, EC, S, p, K);
            
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
            obj_new  = fsol +p*p*log(lambda/2) - lambda*b;
%            LogLik_obj = 0.5*N*log(det(Csol)) - 0.5*N*trace(Csol*EmpCov) - lambda*sum(sum(abs(Csol)))  +p*p* log(lambda/2)  - lambda*b

            if obj_new < obj + precision  
                alpha = alpha/2;
                lambda = (lambda+lambdaouter)/2;
            else
                break;
            end
            obj = obj_new;       
       end  
end

fstart=fstart -(lambda-lambdaold)*sum(sum(abs(Sbase.*Cstart)));
[Csol, Wsol, fsol]=primal_cd_invcov(Cstart, fstart, Wstart, EC, S, p, K);
        

 
