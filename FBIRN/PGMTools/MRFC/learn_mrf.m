function model = learn_mrf(data,method,rho)
% Learns a Gaussian Merkov Random Field from data
% data - n x p data matrix
% method - algorithm used to find MRF (e.g. 'covsel', 'projected_grad',
% 'katya',
% model - structure with fields C (inverse covariance matrix) and mu (mean)

[n,p] = size(data);
model = struct('C',zeros(p,p),'mu',mean(data));  % initialize the model

data = data - ones(n,1)*mean(data); % center the data before learning an MRF
EmpCov = (1/n)*data'*data; % empirical covariance matrix

tic;

%%%MArio ==> Added these codes (not clear what these does )
rho=0.7;             % Controls sparsity
prec=1e-2;
maxiter=1000;
switch lower(method)
    case 'projected_gradient'
        % Yuanqing's implementation of projected gradient
        [C,W]=projGrad_SparseGMRF(EmpCov,rho,prec,maxiter);
    case 'sinco'
        tol=0.000001;
        N=n;
        Cstart=eye(p);
        Wstart=eye(p);
        EC=N*0.5*EmpCov;
        fstart= - trace(EC*Cstart);
        K=N*0.5;
        precision=0.01;
        Sbase=ones(p,p);
        ifvec=0; % vector formulation(1) or scalar formulation(0)
        usediag=0; % 0 - regularize diagonal, 1 - "regularized likelihood"
        if (usediag)
            Sbase=Sbase-eye(p);
        end
        
        lstart = rho*(n/2);
        large_lambdas = [50 10];
        
        %     C = sinco_lambda(EmpCov, N, Sbase, lstart, ifvec, precision, tol);
        
        % if lambda is too low, SINCO will take long time to converge; it's
        % better to start with a 'warm-start' solution for a sparse problem
        
        [Csol_set,C] = sinco_lambda_range(EmpCov, N, [large_lambdas lstart], Cstart, Wstart, fstart, Sbase, ifvec, precision, tol);
        %C corresponds to the last lambda
    case 'alm'
        %% Call ALM to solve the problem
        % the parameters can be tuned
        opts.mxitr = 500; % max iteration number
        opts.mu0 = n; opts.mu0 = 1e-1;  % initial mu
        opts.muf = 1e-3; % final mu
        opts.rmu = 1/4; % ratio of decreasing mu
        opts.tol_gap = 1e-2; % tolerance for duality gap
        opts.tol_frel = 1e-7; % tolerance for relative change of obj value
        opts.tol_Xrel = 1e-10; % tolerance for relative change of X
        opts.tol_Yrel = 1e-10; % tolerance for relative change of Y
        % opts.tol_pinf = 1e-3; % tolerance for infeasibility
        opts.numDG = 10; % every numDG iterations, we compute duality gap since it's expensive
        opts.record = 1; % print stats
        opts.sigma = 1e-10; % sigma is the smoothness parameter
        
        tic; out = SICS_ALM(EmpCov,rho,opts); solveALM = toc;
        C = out.Y;
        
    case  'varsel_mrf'
        %   EmpCov: sample covariance
        %   rho: regularization parameter for the L1 norm (off-diagonal elements)
        %   rho_diag: regularization parameter for the L1 norm (diagonal elements)
        %   tau: regularization parameter for the L1/L2 or L1/Linf norm
        %   group_norm: 2 for L1/L2, inf for L1/Linf
        %   T: number of iterations (10 is usually enough)
        rho_diag = rho;
        tau = 100*rho;
        group_norm = 2;
        T = 10;
        
        [C obj_func] = learn_featsel_ggm(EmpCov, rho, rho_diag, tau, group_norm, T);
        
    case  'no_varsel_mrf'
        %   EmpCov: sample covariance
        %   rho: regularization parameter for the L1 norm (off-diagonal elements)
        %   rho_diag: regularization parameter for the L1 norm (diagonal elements)
        %   tau: regularization parameter for the L1/L2 or L1/Linf norm
        %   group_norm: 2 for L1/L2, inf for L1/Linf
        %   T: number of iterations (10 is usually enough)
        rho_diag = rho;
        tau = 0;
        group_norm = 2;
        T = 10;
        
        [C obj_func] = learn_featsel_ggm(EmpCov, rho, rho_diag, tau, group_norm, T);
        
    case 'covsel'
    otherwise
        %  COVSEL
        %  set Algorithm's parameters ....
        %rho=0.6;             % Controls sparsity
        prec=1e-1;           % Numerical precision
        maxiter=2;           % Maximum Number of sweeps in coordinate descent
        algot='nest';        % Algorith for solving BoxQP: 'sedumi' uses SEDUMI, 'nest' uses smooth minimization (usually faster)
        maxnest=1000;        % Maximum number of iterations in the smooth BoxQP solver
        
        [C,X] = spmlcdvec(EmpCov,rho,maxiter,prec,maxnest,algot); % Inverse covariance is stored in Umat
end

model.C = C;

C1 = C - diag(diag(C));
nonzeros = nnz(sum(C1))
toc;



