% GGM with feature selection
%   [X_allvars obj_func] = learn_featsel_ggm(S_allvars, rho, rho_diag, tau, group_norm, T);
% Input:
%   S_allvars: sample covariance
%   rho: regularization parameter for the L1 norm (off-diagonal elements)
%   rho_diag: regularization parameter for the L1 norm (diagonal elements)
%   tau: regularization parameter for the L1/L2 or L1/Linf norm
%   group_norm: 2 for L1/L2, inf for L1/Linf
%   T: number of iterations (10 is usually enough)
% Output:
%   X_allvars: precision matrix (estimated inverse covariance)
%   obj_func: objective function value after each iteration
% Author: Jean Honorio
% E-mail: jean_honorio@hotmail.com
function [X_allvars obj_func] = learn_featsel_ggm(S_allvars, rho, rho_diag, tau, group_norm, T);

% Check input
if group_norm ~= 2 & group_norm ~= inf
    error('group_norm must be either 2 or inf');
end

% Remove variables from optimization
N = size(S_allvars,1);
if group_norm == 2
    threshold = rho + tau/(2*sqrt(N-1));
else
    threshold = rho + tau/(2*(N-1));
end
idx_zero = find(max(abs(S_allvars-diag(diag(S_allvars))),[],1) <= threshold);
idx_opti = setdiff(1:N,idx_zero);
N = length(idx_opti);
fprintf('%d variables for optimization, %d variables removed...\n',N,length(idx_zero));
S = S_allvars(idx_opti,idx_opti);
X = diag(1 ./ diag(S));

% Iterate T times
obj_func = [];
if N > 0
    first = true;
    for t = 1:T
        fprintf('Iteration %d variable ', t);
        for n = 1:N
            if n > 1
                for rep = 1:fix(log10(n-1)+1)
                    fprintf('\b');
                end
            end
            fprintf('%d', n);
            idx = [1:n-1 n+1:N];
            W = X(idx,idx);
            u = S(idx,n);
            v = S(n,n);
            y = X(idx,n);
            z = X(n,n);
            if first

                % At the beginning, H is inverse of diagonal terms of W
                H = diag(1 ./ diag(W));
                first = false;
            else

                % Compute inverse of W by Woodbury formula twice in order
                % to update previous inverse H
                if n == 1
                    idx_prev = [2:N-1 1];
                    H_prev = H_prev(idx_prev,idx_prev);
                    y_prev = y_prev(idx_prev);
                end
                delta = y_prev - y;
                if n == 1
                    pos = N-1;
                else
                    pos = n-1;
                end
                delta(pos) = (z_prev - z) / 2;
                H_prev_U = [H_prev * delta, H_prev(:,pos)];
                Vt_H_prev = [H_prev(pos,:); delta.' * H_prev];
                Vt_H_prev_U = [Vt_H_prev * delta, Vt_H_prev(:,pos)];
                H = H_prev - H_prev_U * inv(eye(2) + Vt_H_prev_U) * Vt_H_prev;
            end

            % Solve in closed form
            if group_norm == 2
                m = sum((W - diag(diag(W))).^2);
            else
                m = max(abs(W - diag(diag(W))));
            end
            y = solve(v + rho_diag,H,u,rho,tau,group_norm,m,y);
            z = (1 / (v + rho_diag)) + y.' * H * y;

            % Update matrices
            X(idx,n) = y;
            X(n,idx) = y.';
            X(n,n) = z;
            H_prev = H;
            y_prev = y;
            z_prev = z;
        end

        % Compute error
        [U p] = chol(X);
        if p ~= 0
            error('Matrix is not positive definite at iteration %d',t);
        end
        Xd = X - diag(diag(X));
        if group_norm == 2
            E = 2*sum(log(diag(U))) - sum(S(:).*X(:)) - rho*sum(abs(Xd(:))) - rho_diag*(sum(diag(X))) - tau*sum(sqrt(sum(Xd.^2)));
        else
            E = 2*sum(log(diag(U))) - sum(S(:).*X(:)) - rho*sum(abs(Xd(:))) - rho_diag*(sum(diag(X))) - tau*sum(max(abs(Xd)));
        end
        fprintf(' obj.fun. %0.12f\n',E);
        obj_func = [obj_func E];
    end
end

% Return solution
N = size(S_allvars,1);
X_allvars = zeros(N);
X_allvars(idx_opti,idx_opti) = X;
for n = idx_zero
    X_allvars(n,n) = 1 / (S_allvars(n,n) + rho_diag);
end

% Solve for one specific variable at a time
function y = solve(v,H,u,rho,tau,group_norm,m,y0);

N = size(u,1);
if N == 1
    y = -sign(u) * max(abs(u) - rho - tau,0) / (v * H);
    return
end
for n = 1:N
    idx = [1:n-1 n+1:N];
    q = v * H(n,n);
    c = -(v * (H(n,idx) * y0(idx)) + u(n));
    if group_norm == 2
        a = sum(y0(idx).^2);
    else
        a = max(abs(y0(idx)));
    end
    b = m(n);
    if group_norm == 2

        % Solve L1/L2 lasso subproblem by Newton's method
        Tinner = 10;
        epsilon = 1e-4;
        a = a + epsilon;
        b = b + epsilon;
        dl = -c - rho;
        dr = -c + rho;
        if dl <= 0 & dr >= 0
            x = 0;
        elseif dl > 0
            x = 0;
            for tinner = 1:Tinner
                x2 = x*x;
                sqrtx2a = sqrt(x2 + a);
                sqrtx2b = sqrt(x2 + b);
                dl = q*x - c - rho + tau*x/(2*sqrtx2a) ...
                                   + tau*x/(2*sqrtx2b);
                hl = q + tau*a/(2*sqrtx2a*sqrtx2a*sqrtx2a) ...
                       + tau*b/(2*sqrtx2b*sqrtx2b*sqrtx2b);
                if dl <= 0
                    break
                end
                x = x - dl/hl;
            end
        else
            x = 0;
            for tinner = 1:Tinner
                x2 = x*x;
                sqrtx2a = sqrt(x2 + a);
                sqrtx2b = sqrt(x2 + b);
                dr = q*x - c + rho + tau*x/(2*sqrtx2a) ...
                                   + tau*x/(2*sqrtx2b);
                hr = q + tau*a/(2*sqrtx2a*sqrtx2a*sqrtx2a) ...
                       + tau*b/(2*sqrtx2b*sqrtx2b*sqrtx2b);
                if dr >= 0
                    break
                end
                x = x - dr/hr;
            end
        end
        y0(n) = x;
    else

        % Solve L1/Linf lasso subproblem
        d1 = min(a,b);
        d2 = max(a,b);
        if d2 == 0
            w0 = rho + tau;
            w1 = nan;
            w2 = nan;
            d1 = inf;
            d2 = inf;
            range = 3:4;
        elseif d1 == 0
            w0 = rho + tau/2;
            w1 = tau/2;
            w2 = nan;
            d1 = d2;
            d2 = inf;
            range = 2:5;
        else
            w0 = rho;
            w1 = tau/2;
            w2 = tau/2;
            range = 1:6;
        end
        best_x = [];
        best_fx = inf;
        for i = range
            switch i
                case 1
                    x = min((c + w0 + w1 + w2) / q, -d2);
                case 2
                    x = min(max((c + w0 + w1) / q, -d1), -d2);
                case 3
                    x = min(max((c + w0) / q, -d1), 0);
                case 4
                    x = min(max((c - w0) / q, 0), d1);
                case 5
                    x = max(max((c - w0 - w1) / q, d1), d2);
                case 6
                    x = max((c - w0 - w1 - w2) / q, d2);
            end
            fx = 0.5*q*x^2 - c*x + w0*abs(x);
            if ~isnan(w1)
                fx = fx + w1*max(abs(x),d1);
            end
            if ~isnan(w2)
                fx = fx + w2*max(abs(x),d2);
            end
            if fx < best_fx
                best_x = x;
                best_fx = fx;
            end
        end
        y0(n) = best_x;
    end
end
y = y0;
