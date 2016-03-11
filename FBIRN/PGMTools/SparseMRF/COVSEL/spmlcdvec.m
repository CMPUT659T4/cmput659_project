function [S,E,gvals,cputimes] = spmlcdvec(C,R,maxiters,prec,maxnest,algo)
% function [S,E] = spmlcdvec(C,R,maxiters)
% SPMLCDVEC solves the weighted one-norm penalized maximum 
% likelihood problem
% maximize     log det(S) - trace(C*S) - sum(sum(R.*abs(S))
% s.t.         S positive-definite
% in the dual form
% min_E           -log det(C+E) : abs(E) <= R, C+E p.d. ,
% using a block-coordinate descent method, updating one colum 
% (and corresponding row) at a time.
% inputs:
% C         nxn symmetric, positive-definite matrix
% R         symmetric componentwise non-negative nxn matrix 
%           (only the upper triangular part of R is considered)
% maxiters  maximum number of outer loops of the block-coordinate 
%           descent scheme (i.e. how many times we update all elements)
%           default: 1.
% prec  target precision
% maxnest   max number of iterations in the BOX-QP smooth minimization code
% algo  'sedumi' uses SEDUMI to solve the box QP, anything else means smooth minimization is used instead
%
% outputs:
% E     optimal nxn positive-definite matrix
% S     corresponding estimate of inverse of covariance matrix
%
% L. El Ghaoui, A d'Aspremont, March 2006

%%% preliminaries
% problem size
freq=5;disp(' ');
n = max(size(C));period=ceil(n/freq);
tic;cputimes=[];gvals=[];
% defaults & error checks
if nargin <= 2,
    maxiters = 1;
end
if nargin <= 3,
    prec = 0;
end
if nargin <= 4,
    maxnest=1000;
end
if nargin <= 5,
    algo='nest';
end
if max(size(R)) == 1, R = R.*ones(n,n); end
if min(min(R)) < 0, 
    error('Weighting matrix R must be componentwise non-negative.'); 
end
% make sure R is symmetric
R = triu(R) + tril(R',-1);

%%% do coordinate descent
E = diag(diag(R)); % this solves the diagonal part
W = C+E;
S=inv(W);
gap=trace(C*S)+trace(R*abs(S))-n;
gvals=[gvals;gap];cputimes=[cputimes;toc];objval=0;
% start loop
for k = 1:maxiters,
    % go through all rows/columns, updating the whole row and column at a time, except for the diagonal element
    for i = 1:n,
        indi = 1:n; indi(indi==i) = [];
        % form the inverse of W(indi,indi)
        invWi = inv(W(indi,indi));
        % set up the QP
        % bounds on variables
        prob.bux = C(indi,i)+R(indi,i); % upper bound
        prob.blx = C(indi,i)-R(indi,i); % lower bound
        if strcmp(algo,'sedumi')
            y=sedumiQP(invWi,prob.bux,prob.blx);
        else
            y=BoxQP(invWi, prob.blx, prob.bux, maxnest,maxnest,0);
       end
        % update row and column of W matrix
        W(indi,i) = y;
        W(i,indi) = y';
        if mod(i,period)==0
            S=inv(W);
            gap=trace(C*S)+trace(R*abs(S))-n;
            objval=log(det(W))+n;
            gvals=[gvals;gap];cputimes=[cputimes;toc];
            disp(['CD Iter: ', num2str(k),'   Primal: ',num2str(objval,'%.4e'),'   Gap: ',num2str(gap,'%.4e')]);
        end
        if gap<=prec
            break;    
        end
    end
    
    % Record function values
    Wt = tril(W')+triu(W,1);
    S = inv(Wt);
    S = tril(S')+triu(S,1);
    if gap<=prec % Test is gap is less than prec target
        break;
    end
end

% set output variables
W = tril(W')+triu(W,1);
E = W-C;
S = inv(W);
S = tril(S')+triu(S,1);