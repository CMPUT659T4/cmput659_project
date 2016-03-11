
rho = 50; % weighting parameter
n = 200;  % dimension of the matrix
rs = 20; % random seed

getdata = 'L';  %getdata = 'K'; % the way to get the data
if strcmp(getdata,'L') == 1
    [A, B, prec] = getdata_L(n, rs);
elseif strcmp(getdata,'K') == 1
    N = 5*n; nzfrac = 0.01; prec = 0.1;
    [A, Cstart, fstart, Wstart, EC, Sbase, lambda, K, EmpCov] = getdata_K(n, N, nzfrac, rho, rs);
    B = EmpCov;
end


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

tic; out = SICS_ALM(B,rho,opts); solveALM = toc;
% compute the dual gap corresponding to X
[V,D] = mexeig(out.X); d = diag(D); Ux = V*diag(1./d)*V' - B; Ux = min(rho,max(-rho,Ux));
fX = -sum(log(d))+sum(sum(B.*out.X))+rho*sum(sum(abs(out.X)));
[V,D] = mexeig(B+Ux); dUx = sum(log(diag(D))) + n;
gapX = fX - dUx;
% compute the dual gap corresponding to Y
[V,D] = mexeig(out.Y); d = diag(D);
if min(d) <= 0
    gapY = inf;
else
    Uy = V*diag(1./d)*V' - B; Uy = min(rho,max(-rho,Uy));
    fY = -sum(log(d))+sum(sum(B.*out.Y))+rho*sum(sum(abs(out.Y)));
    [V,D] = mexeig(B+Uy); dUy = sum(log(diag(D))) + n;
    gapY = fY - dUy;
end
%%% print output
fprintf('n:%d,rho:%3.2f,iter:%d,gap:%3.1e,gapX:%3.2e,gapY:%3.2e,time:%3.2f\n',n,rho,out.iter,out.gap,gapX,gapY,solveALM);
