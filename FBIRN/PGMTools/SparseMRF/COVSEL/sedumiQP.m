function res=sedumiQP(A,Up,Low)
% Solves he following QP using SEDUMI
% Min x'Ax s.t. Low_i <= x_i <= Up_i

n=max(size(A));
Amat=[];bvec=[];
Amat=[Amat;eye(n),zeros(n,1)];bvec=[bvec;Up];
Amat=[Amat;-eye(n),zeros(n,1)];bvec=[bvec;-Low];
Amat=[Amat;[zeros(1,n),-1];[chol(A),zeros(n,1)]];bvec=[bvec;zeros(n+1,1)];
%Amat=[Amat;[zeros(1,n),-1];[eye(n),zeros(n,1)]];bvec=[bvec;zeros(n+1,1)];
objvec=-[zeros(n,1);1];
K.l=2*n;
K.q=n+1;
pars.fid=0;
[xsol,ysol]=sedumi(Amat,objvec,bvec,K,pars);
res=ysol(1:n);