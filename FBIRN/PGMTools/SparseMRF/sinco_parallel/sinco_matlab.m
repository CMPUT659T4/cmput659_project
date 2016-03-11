function [Csol, Wsol, fsol]=sinco_matlab(Cstart, fstart, Wstart, A, S, K, tol, maxit);
% 
%  SINCO_MATLAB is the matlab implementation of SINCO algorithms, 
%  that applies greedy coordinate descent to solve
%  the following problem 
% 
% max _C[ K*log(det(C)-tr(AC)-\lambda ||S.*C||_1]
% 
% K and lambda are positive scalars, A is a given symmetric matrix, 
% S is a given nonnegative matrix
% (not necessarily symmetric), S.*C is a elementwhise product and ||*||_1 is the sum
% of all absolute values.
% 
% The method works by evaluating the best function value improvement that can be 
% obtained by changing only the (i,j) and (j,i) entry of matrix C' or C''. 
% For each pair a potential improving step is computed (by solving a
% quadratic equation) given the matrix W - the inverse of C=C'+C'. 
% Then he best (most improving) step is chosen. The C', C'' and W are
% updated accrdingly.
%
%
% Initialize the current solution, split it into positive (C') and negative (C'') 
% parts and compute the gradient with respect to C' and C'' - G' and G''.
%
%
p=size(A,1);
fsol=fstart;
Wsol=Wstart;
Cp=Cstart.*(Cstart>0);
Cpp=-Cstart.*(Cstart<0);
Gradp=-S+K*Wsol-A;
Gradpp=-S-K*Wsol+A;

%Initilaize a matrix that keeps 0/1 information that indicates if a step
%for particular pair (i,j) needs to be recomputed or that the past
%information is still valid on this iteration, becuase W_{ij} did not
%change. This may be useful when W is fairly sparse, but if not, we may
%ignore keeping and updating UpInd, to avoid the overhead. Right now it is
%used every time)
UpInd=ones(p,p);

% matrices that will store all steps (alphas), directions (updates=+1/-1) and
% all corresponding objective function updates (fchanges))
alphas=zeros(p,p);
fchange=zeros(p,p);
updates=zeros(p,p);
% initilize stopping triggers
stop=0;
iter=0;
% start the main loop
while (stop~=1 & iter<maxit )
  alphamax=0;
  fmax=fsol;
  iter=iter+1;
% For each pair (i,j) compute the step, if W(i,j), W(i,i) or W(j,j) were
% updated on the prior iteration as recorded by UpInd matrix.
  for i=1:p
      for j=1:i
        if (UpInd(i,j) || UpInd(i,i) || UpInd(j,j))
          alphas(i,j)=0;
 % If the G'(i,j)+G'(j,i)>0 then we want to increase C'(i,j), but only if
 % C''(i,j) is already zero. Find the step-lengths for pisitive step via solving the
 % quadratic equation.
          if (Gradp(i,j)+Gradp(j,i) > tol) & (Cpp(i,j)<=tol)
            updates(i,j)=1;
            alphas(i,j)=findposstep(K, Wsol, i, j, A, S, updates(i,j));  
 % If the G''(i,j)+G''(j,i)>0 then we want to increase C''(i,j), but only if
 % C'(i,j) is already zero. Find the step-lengths for positive step via solving the
 % quadratic equation.
          elseif (Gradpp(i,j)+Gradpp(j,i) > tol) & (Cp(i,j)<=tol)
            updates(i,j)=-1;
            alphas(i,j)=findposstep(K, Wsol, i, j, A, S, updates(i,j)); 
 % If the G'(i,j)+G'(j,i)<0 then we want to deccrease C'(i,j), but only if
 % C'(i,j) is not zero. Find the step-lengths for negative step via solving the
 % quadratic equation.
          elseif (Gradp(i,j)+Gradp(j,i)<-tol) & (Cp(i,j)> tol )
            updates(i,j)=1;
            alphas(i,j)=findnegstep(K, Wsol, i, j, A, S, Cp, Cpp, updates(i,j)); 
 % If the G''(i,j)+G''(j,i)<0 then we want to deccrease C''(i,j), but only if
 % C''(i,j) is not zero. Find the step-lengths for positive step via solving the
 % quadratic equation.
          elseif (Gradpp(i,j)+Gradpp(j,i)<-tol) & (Cpp(i,j)> tol) 
            updates(i,j)=-1;
            alphas(i,j)=findnegstep(K, Wsol, i, j, A, S, Cp, Cpp, updates(i,j)); 
          end
          if ( abs(alphas(i,j))> 10^(-6))
%  If the step is not too small, compute the change in the obective function
%  for each (i,j)
           fchange(i,j)=funvalue_update(alphas(i,j), K, Wsol, i, j, A, S, updates(i,j));
%          Correctness check for debugging purposes
%          fnewtest=K*log(det(Cptemp-Cpptemp))-sum(sum(A.*(Cptemp-Cpptemp)))-sum(sum(S.*(Cptemp+Cpptemp)));
%          if (abs(f+fchange(i,j)-fnewtest)/(1+abs(fnewtest)) >10^(-6))
%              display 'bad function update'
%             keyboard
%          end
          else
           fchange(i,j)=0;
          end
        end
%  Select the maximum of the objective function improvement computed so
%  far.
        fnew=fsol+fchange(i,j);
        if fnew > fmax
          fmax=fnew;
          imax=i;
          jmax=j;
          alphamax=alphas(i,j);
          updatemax=updates(i,j);
        end
      end
  end
 % Check is the maximum improvment is above the toleraqnce level 
  if alphamax==0 || max(max(fchange))/(abs(fsol)+1)<tol;
     display 'The improvements are below tolerance' 
     stop=1;
  end
  % Make the appropriate updates to C', C'', G' and G'', as well as to
  % W=C^{-1}.
  update=updatemax;
  alpha=alphamax;
  i=imax;
  j=jmax;
  if update==1
   Cp(i,j)=Cp(i,j)+alpha;
   Cp(j,i)=Cp(j,i)+alpha;
  else
   Cpp(i,j)=Cpp(i,j)+alpha;
   Cpp(j,i)=Cpp(j,i)+alpha;
  end
  fsol=fmax;
% Checking the function value computation for debigging purposes.
%  aaa=f-(K*log(det(Cp-Cpp))-sum(sum(A.*(Cp-Cpp)))-sum(sum(S.*(Cp+Cpp))));
%  if ( abs(aaa) > 10^(-8))
%    keyboard
%  end
  [Gradp, Gradpp, UpInd, Wsol]=invupdate(alpha, Gradp, Gradpp, Wsol, K, i, j, update);
end
% Recodr the current solution
Csol=Cp-Cpp;
% Control  for debigging puposes
% if fsol< fstart
%        display 'negative improvement'
%        keyboard
% end
