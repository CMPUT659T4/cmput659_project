function [Csol, Wsol, fsol]=sinco_matlab_parallel(Cstart, fstart, Wstart, A, S, K, tol, maxit)
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

  % Build up a sliced variable of job parameters for the parfor
  fprintf(1,'Preparing parfor job parameters\n');
  jobIdx = 1;
  for i=1:p
    for j=1:i
      currentParams.i = i;
      currentParams.j = j;
      currentParams.UpInds.ii = UpInd(i,i);
      currentParams.UpInds.jj = UpInd(j,j);
      currentParams.UpInds.ij = UpInd(i,j);
      currentParams.UpInds.ji = UpInd(j,i);
      currentParams.Gradps.ii = Gradp(i,i);
      currentParams.Gradps.jj = Gradp(j,j);
      currentParams.Gradps.ij = Gradp(i,j);
      currentParams.Gradps.ji = Gradp(j,i);
      currentParams.Gradpps.ii = Gradpp(i,i);
      currentParams.Gradpps.jj = Gradpp(j,j);
      currentParams.Gradpps.ij = Gradpp(i,j);
      currentParams.Gradpps.ji = Gradpp(j,i);
      currentParams.Cps.ii = Cp(i,i);
      currentParams.Cps.jj = Cp(j,j);
      currentParams.Cps.ij = Cp(i,j);
      currentParams.Cps.ji = Cp(j,i);
      currentParams.Cpps.ii = Cpp(i,i);
      currentParams.Cpps.jj = Cpp(j,j);
      currentParams.Cpps.ij = Cpp(i,j);
      currentParams.Cpps.ji = Cpp(j,i);
      currentParams.Ws.ii = Wsol(i,i);
      currentParams.Ws.jj = Wsol(j,j);
      currentParams.Ws.ij = Wsol(i,j);
      currentParams.Ws.ji = Wsol(j,i);
      currentParams.As.ii = A(i,i);
      currentParams.As.jj = A(j,j);
      currentParams.As.ij = A(i,j);
      currentParams.As.ji = A(j,i);
      currentParams.Ss.ii = S(i,i);
      currentParams.Ss.jj = S(j,j);
      currentParams.Ss.ij = S(i,j);
      currentParams.Ss.ji = S(j,i);
      currentParams.updates.ij = updates(i,j);
      currentParams.alphas.ij = alphas(i,j);
      currentParams.tol = tol;
      currentParams.K = K;

      job_params{jobIdx} = currentParams;
      jobIdx = jobIdx + 1;
    end
  end
  
  numJobs = length(job_params);

% For each pair (i,j) compute the step, if W(i,j), W(i,i) or W(j,j) were
% updated on the prior iteration as recorded by UpInd matrix.
  parfor jobNum = 1:numJobs
     
        params = job_params{jobNum};
        % move these outside if statements so we'll still have a value if conditions never set
        % these and parfor still needs a result
        currentUpdates = params.updates.ij;
        currentAlphas = params.alphas.ij;

        if (params.UpInds.ij || params.UpInds.ii || params.UpInds.jj)
          currentAlphas = 0;   %alphas(i,j)=0;
 % If the G'(i,j)+G'(j,i)>0 then we want to increase C'(i,j), but only if
 % C''(i,j) is already zero. Find the step-lengths for pisitive step via solving the
 % quadratic equation.
          if (params.Gradps.ij+params.Gradps.ji > params.tol) & (params.Cpps.ij<=params.tol)
            currentUpdates = 1;    %updates(i,j)=1;
            %alphas(i,j)=findposstep(K, Wsol, i, j, A, S, updates(i,j));
            currentAlphas = findposstep_parallel(params.K, params.Ws, params.i, params.j, params.As, params.Ss, currentUpdates);
  
 % If the G''(i,j)+G''(j,i)>0 then we want to increase C''(i,j), but only if
 % C'(i,j) is already zero. Find the step-lengths for positive step via solving the
 % quadratic equation.
          elseif (params.Gradpps.ij+params.Gradpps.ji > params.tol) & (params.Cps.ij<=params.tol)
             currentUpdates = 1; %updates(i,j)=-1;
             %alphas(i,j)=findposstep(K, Wsol, i, j, A, S, updates(i,j)); 
             currentAlphas = findposstep_parallel(params.K, params.Ws, params.i, params.j, params.As, params.Ss, currentUpdates);

 % If the G'(i,j)+G'(j,i)<0 then we want to deccrease C'(i,j), but only if
 % C'(i,j) is not zero. Find the step-lengths for negative step via solving the
 % quadratic equation.
          elseif (params.Gradps.ij+params.Gradps.ji<-params.tol) & (params.Cps.ij> params.tol )
            currentUpdates = 1;  %updates(i,j)=1;
            %alphas(i,j)=findnegstep(K, Wsol, i, j, A, S, Cp, Cpp, updates(i,j)); 
            currentAlphas=findnegstep_parallel(params.K, params.Ws, params.i, params.j, params.As, params.Ss, params.Cps, params.Cpps, currentUpdates);

 % If the G''(i,j)+G''(j,i)<0 then we want to deccrease C''(i,j), but only if
 % C''(i,j) is not zero. Find the step-lengths for positive step via solving the
 % quadratic equation.
          elseif (params.Gradpps.ij+params.Gradpps.ji<-params.tol) & (params.Cpps.ij> params.tol) 
            currentUpdates = 1; %updates(i,j)=-1;
            %alphas(i,j)=findnegstep(K, Wsol, i, j, A, S, Cp, Cpp, updates(i,j)); 
            currentAlphas=findnegstep_parallel(params.K, params.Ws, params.i, params.j, params.As, params.Ss, params.Cps, params.Cpps, currentUpdates);
          end

          %if ( abs(alphas(i,j))> 10^(-6))
          if ( abs(currentAlphas) > 10^(-6))
%  If the step is not too small, compute the change in the obective function
%  for each (i,j)
           %fchange(i,j)=funvalue_update(alphas(i,j), K, Wsol, i, j, A, S, updates(i,j));
           currentfchange=funvalue_update_parallel(currentAlphas, params.K, params.Ws, params.i, params.j, params.As, params.Ss, currentUpdates);

%          Correctness check for debugging purposes
%          fnewtest=K*log(det(Cptemp-Cpptemp))-sum(sum(A.*(Cptemp-Cpptemp)))-sum(sum(S.*(Cptemp+Cpptemp)));
%          if (abs(f+fchange(i,j)-fnewtest)/(1+abs(fnewtest)) >10^(-6))
%              display 'bad function update'
%             keyboard
%          end
          else
           currentfchange = 0;  %fchange(i,j)=0;
          end
        end


%  Select the maximum of the objective function improvement computed so far.

% this breaks the independence of the parfor so we need to do it another way
%        fnew=fsol+fchange(i,j);
%        if fnew > fmax
%          fmax=fnew;
%          imax=i;
%          jmax=j;
%          alphamax=alphas(i,j);
%          updatemax=updates(i,j);
%        end

     % Store the results we computed in this iteration so we can unroll it later
     job_results{jobNum}.updates = currentUpdates;
     job_results{jobNum}.alphas = currentAlphas;
     job_results{jobNum}.fchange = currentfchange;
      
  end  % end parfor

  % We finished the parfor, we need to unroll the results from sliced results variable
  jobIdx = 1;
  for i = 1:p
    for j = 1:i
      updates(i,j) = job_results{jobIdx}.updates;
      alphas(i,j) = job_results{jobIdx}.alphas;
      fchange(i,j) = job_results{jobIdx}.fchange;
  
      fnew=fsol+fchange(i,j);
      if fnew > fmax
        fmax=fnew;
        imax=i;
        jmax=j;
        alphamax=alphas(i,j);
        updatemax=updates(i,j);
      end

      jobIdx = jobIdx + 1;
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
