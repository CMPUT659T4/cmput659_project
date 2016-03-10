function [X,W]=projGrad_SparseGMRF(S,lamda,tolrel,maxIter)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% solve optimization min_X -logdet(X)+tr(S*X)+lamda*sum_{i~=j}|X_{ij}|  via dual  %%%%%%
%%% S           MxM matrix      empirical covariance matrix                         %%%%%%
%%% lamda       scalor          L1-norm regularization parameter                    %%%%%%
%%% tolrel      scalor          tolerance of relative duality gap (stopping criterion) %%%
%%% X           MxM matrix      estimated inversion covariance matrix                %%%%%
%%% W           MxM matrix      estimated covariance matrix                          %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%      Yuanqing Lin     linyuanq@seas.upenn.edu    12/22/2007  %%%%%%%%%%%%%%%%%%%%%%
%%%%%%      revised on 05/14/2008  %%%%%%%%%%%%%%%%%%%%%%
tic
if nargin<3 tolrel=10^(-13); end
if nargin<4 maxIter=10000; end

M=length(S(:,1));
% initialization
W=(S-diag(diag(S)))*(1-10^(-5))+diag(diag(S)); %S.*(1-10^(-8)*mean(diag(S)))+eye(M)*10^(-8)*mean(diag(S));
[LL,p]=chol((W+W')/2);
objf_d=2*sum(log(diag(LL)))+M;%dual objective function;
X=inv(W);
[LL,p]=chol((X+X')/2);
objf_p=-2*sum(log(diag(LL)))+sum(sum(S.*X))+lamda*(sum(sum(abs(X)))-sum(abs(diag(X)))); % primal objective function 
gap=objf_p-objf_d;
gaprel=gap/abs(objf_p);
i=1;
U=S+lamda;
L=S-lamda;
I_diag=((0:M-1)*M)+(1:M);

%%%% line search parameters
BETA=0.5;
ALPHA=0.01;
t0=1;

%%%%%% projected gradient ascent
while(gaprel>=tolrel)&(i<=maxIter)
    % computing the gradient
    %dW=inv(W);
    dW=X;
    dW(I_diag)=0;
    
    %% line search, decreasing stepsize
    t=t0;
    doneF=0;
    while(~doneF)
        Wp=W+t*dW;
        %projection
        I=Wp>=U; Wp(I)=U(I);
        I=Wp<=L; Wp(I)=L(I);
        %check positiveness using chol
        [LL,p]=chol((Wp+Wp')/2);
        if (p==0) %matrix is positive definite
            objf_d_new=2*sum(log(diag(LL)))+M;
            if (objf_d_new-objf_d>ALPHA*sum(sum(dW.*(Wp-W))))&(t>=t0*10^(-20))
                doneF=1;
            else
                t=BETA*t;
            end
            if t<=10^(-20)*t0
                break;
                fprintf('Too small stepsize... stop.\n')
            end
        else
            t=BETA*t;
        end
    end
    %% line search, increasing stepsize
    if t==t0
        doneF=0;
        t=t/BETA;
        while(~doneF)
            Wp=W+t*dW;
            %projection
            I=Wp>=U; Wp(I)=U(I);
            I=Wp<=L; Wp(I)=L(I);
            %check positiveness using chol
            [LL,p]=chol((Wp+Wp')/2);
            if (p==0) %matrix is positive definite
                objf_d_new=2*sum(log(diag(LL)))+M;
                if ((objf_d_new-objf_d<ALPHA*sum(sum(dW.*(Wp-W))))|(t>10^(10)*t0)) %% if not significant increase
                    doneF=1;
                else
                    t=t/BETA;
                end
            else
                doneF=1;
            end
        end
        t=t*BETA; %% one step back
    end
    
    %update the results
    t0=t;
    W=W+t*dW;
    %projection
    I=W>=U; W(I)=U(I);
    I=W<=L; W(I)=L(I);
    %compute objective functions
    [LL,p]=chol((W+W')/2);
    objf_d=2*sum(log(diag(LL)))+M;
    X=inv(W);
    I1=W==U; I1(((0:M-1)*M)+(1:M))=1;
    I2=W==L;
    X_p=X;  %% primal variable
    X_p(~(I1+I2))=0;   %% if W_ij is strictly inside the bound, X_ij must be zero. 
    [LL,p]=chol((X_p+X_p')/2);
    if p==0   %% it is possible that X_p is not positive definite
        objf_p=-2*sum(log(diag(LL)))+sum(sum(S.*X_p))+lamda*(sum(sum(abs(X_p)))-sum(abs(diag(X_p))));
    else
        objf_p=inf;
    end
    gap=objf_p-objf_d;
    gaprel=gap/abs(objf_d);
    gaps(i)=gap;
    gaprels(i)=gaprel;
    figure(10); subplot(2,1,1); semilogy(gaps,'r:.'); title('Duality gap')
    subplot(2,1,2); semilogy(gaprels,'r:.'); title('Relative duality gap'); xlabel('Iterations')
    if (i-1)/10==floor((i-1)/10)
        fprintf('iter %5d, t=%5.3e, gaprel=%8.4e,\n',i,t,gaprel);
    end
    %drawnow;
    i=i+1;
end
if i==maxIter
    drawnow;
    statusF='Maximum iteration reached';
else
    drawnow;
    statusF='solved';
    fprintf('Solved! iter %d, t=%5.3e, gaprel=%8.4e,\n',i,t,gaprel);
end
toc



