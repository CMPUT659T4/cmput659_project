function [betasparse indsbeta residue outgain] = larsen4(y, p, lambda2, stop, gain, trace)
%
%   Modified by Rahul Garg to carry out Granger causality determination
%   Assume a global UniverseData which stores time series in rows
%   Returns lasso/elastic nets coefficients in betasparse (row vector)
%   And the corresponding indices into the rows of UniverseData
%   Also, does not store betas corresponding to the whole regression
%   path. Saves memory just by storing last beta.
%
% LARSEN2  The LARSEN algorithm for elastic net regression using globally
%    available data UniverseData, which stores time series in rows.
%    BETA = LARSEN(X, Y) performs elastic net regression on the variables
%    in X to approximate the response Y.  Variables X are assumed to be
%    normalized (zero mean, unit length), the response Y is assumed to be
%    centered. The ridge term coefficient, lambda2, has a default value of
%    1e-6. This keeps the ridge influence low while making p > n possible.
%    BETA = LARSEN(X, Y, LAMBDA2) adds a user-specified LAMBDA2. LAMBDA2 =
%    0 produces the lasso solution.
%    BETA = LARSEN(X, Y, LAMBDA2, STOP) with nonzero STOP will perform
%    elastic net regression with early stopping. If STOP is negative, its 
%    absolute value corresponds to the desired number of variables. If STOP
%    is positive, it corresponds to an upper bound on the L1-norm of the
%    BETA coefficients.
%    BETA = LARSEN(X, Y, LAMBDA2, STOP, TRACE) with nonzero TRACE will
%    print the adding and subtracting of variables as all elastic net
%    solutions are found.
%    Returns BETA where each row contains the predictor coefficients of
%    one iteration. A suitable row is chosen using e.g. cross-validation,
%    possibly including interpolation to achieve sub-iteration accuracy.
%
% Author: Karl Skoglund, IMM, DTU, kas@imm.dtu.dk
% Reference: 'Regularization and Variable Selection via the Elastic Net' by
% Hui Zou and Trevor Hastie, 2005.

global UniverseData;

%% Input checking
if nargin < 6
  trace = 0;
end
if nargin < 5
  gain = 1e10;
end
if nargin < 4
  stop = 0;
end
if nargin < 3
  lambda2 = 1e-6;
end

%% Elastic net variable setups
n = size(UniverseData, 2);
maxk = 8*(n+p); % Maximum number of iterations

if lambda2 < eps
  nvars = min(n-1,p); %Pure LASSO
else
  nvars = p; % Elastic net
end
if stop > 0,
  stop = stop/sqrt(1 + lambda2);
end


if stop == 0
%  beta = zeros(2*nvars, p);
elseif stop < 0
%  beta = zeros(2*round(-stop), p);
else
%  beta = zeros(100, p);
end
% TBD -- remove the last 7 lines
beta = zeros(1, p);
last_beta = zeros(1, p);


mu = zeros(n, 1); % current "position" as LARS-EN travels towards lsq solution
I = 1:p; % inactive set
A = []; % active set

R = []; % Cholesky factorization R'R = X'X where R is upper triangular

lassocond = 0; % Set to 1 if LASSO condition is met
stopcond = 0; % Set to 1 if early stopping condition is met
k = 0; % Algorithm step count
vars = 0; % Current number of variables

d1 = sqrt(lambda2); % Convenience variables d1 and d2
d2 = 1/sqrt(1 + lambda2); 

if trace
  disp(sprintf('Step\tAdded\tDropped\t\tActive set size'));
end

smallerr = 0;
n1 = 0; n2 = 0; e1 = 0; e2 = 0;

%% Elastic net main loop
while vars < nvars && ~stopcond && k < maxk && ~smallerr
  k = k + 1;
  c = UniverseData*(y - mu)*d2;
  [C j] = max(abs(c(I)));
  j = I(j);

  if ~lassocond % if a variable has been dropped, do one iteration with this configuration (don't add new one right away)
    tmp = UniverseData(j, :);
    tmp2 = UniverseData(A, :);
    R = cholinsert(R, tmp', tmp2', lambda2);
    A = [A j];
    I(I == j) = [];
    vars = vars + 1;
    if trace
      disp(sprintf('%d\t\t%d\t\t\t\t\t%d', k, j, vars));
    end
  end

  s = sign(c(A)); % get the signs of the correlations

  GA1 = R\(R'\s);  
  AA = 1/sqrt(sum(GA1.*s));
  w = AA*GA1;
  u1 = d2 * (w' * UniverseData(A, :))' ; % equiangular direction (unit vector) part 1
  u2 = zeros(p, 1); u2(A) = d1*d2*w; % part 2
  if vars == nvars % if all variables active, go all the way to the lsq solution
    gamma = C/AA;
  else
    a = (UniverseData*u1 + d1*u2)*d2; % correlation between each variable and eqiangular vector
    temp = [(C - c(I))./(AA - a(I)); (C + c(I))./(AA + a(I))];
    gamma = min([temp(temp > 0); C/AA]);
  end

  % LASSO modification
  lassocond = 0;
  temp = -last_beta(1, A)./w';
  [gamma_tilde] = min([temp(temp > 0) gamma]);
  j = find(temp == gamma_tilde);
  if gamma_tilde < gamma,
    gamma = gamma_tilde;
    lassocond = 1;
  end

%  e1 = (y - mu)' * (y - mu);
  mu = mu + gamma*u1;
  e2 = (y - mu)' * (y - mu);
  residue = e2;
  if (e2) <= eps
    smallerr = 1;
  end;

%  n1 = sum(abs(beta(k,:)));
%  beta(k+1,A) = beta(k,A) + gamma*w';
   beta(1, A) = last_beta(1, A) + gamma*w';
   n2 = sum(abs(beta(1,:)));

  if ((gamma > eps))
    outgain = (n2 - n1) / (e1 - e2);
    if trace
        disp(sprintf('(%f %f %e) (%f %f %e) %f\n', e1, e2, e1 - e2, n1, n2, n2 - n1, outgain));
    end; 
    if outgain >= gain
        smallerr = 1;
    end;
  end;  

  % Early stopping at specified bound on L1 norm of beta
  if stop > 0
    t2 = sum(abs(beta(1,:)));
    if t2 >= stop
      t1 = sum(abs(last_beta(1,:)));
      s = (stop - t1)/(t2 - t1); % interpolation factor 0 < s < 1
      beta(1,:) = last_beta(1,:) + s*(beta(1,:) - last_beta(1,:));
      stopcond = 1;
    end
  end

  % If LASSO condition satisfied, drop variable from active set
  if lassocond == 1
    R = choldelete(R,j);
    I = [I A(j)];
    A(j) = [];
    vars = vars - 1;
    if trace
      disp(sprintf('%d\t\t\t\t%d\t\t\t%d', k, j, vars));
    end
  end

  % Early stopping at specified number of variables
  if stop < 0
    stopcond = vars >= -stop;
  end

  last_beta = beta;
  n1 = n2; e1 = e2;
  
end

% trim beta
%if size(beta,1) > k+1
%  beta(k+2:end, :) = [];
%end

% divide by d2 to avoid double shrinkage
beta = beta/d2;


if k == maxk
  disp('LARS-EN warning: Forced exit. Maximum number of iteration reached.');
end

indsbeta = find(beta(1, :));
betasparse = beta(1, indsbeta);



%% Fast Cholesky insert and remove functions
% Updates R in a Cholesky factorization R'R = X'X of a data matrix X. R is
% the current R matrix to be updated. x is a column vector representing the
% variable to be added and X is the data matrix containing the currently
% active variables (not including x).
function R = cholinsert(R, x, UniverseData, lambda)
diag_k = (x'*x + lambda)/(1 + lambda); % diagonal element k in UniverseData'UniverseData matrix
if isempty(R)
  R = sqrt(diag_k);
else
  col_k = x'*UniverseData/(1 + lambda); % elements of column k in X'X matrix
  R_k = R'\col_k'; % R'R_k = (UniverseData'UniverseData)_k, solve for R_k
  R_kk = sqrt(diag_k - R_k'*R_k); % norm(x'x) = norm(R'*R), find last element by exclusion
  R = [R R_k; [zeros(1,size(R,2)) R_kk]]; % update R
end

% Deletes a variable from the UniverseData'UniverseData matrix in a Cholesky factorisation R'R =
% X'X. Returns the downdated R. This function is just a stripped version of
% Matlab's qrdelete.
function R = choldelete(R,j)
R(:,j) = []; % remove column j
n = size(R,2);
for k = j:n
  p = k:k+1;
  [G,R(p,k)] = planerot(R(p,k)); % remove extra element in column
  if k < n
    R(p,k+1:n) = G*R(p,k+1:n); % adjust rest of row
  end
end
R(end,:) = []; % remove zero'ed out row



    
