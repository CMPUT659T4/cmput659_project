% Test script for Covariance selection code.
% Initialize parameters ****************
% Data ....

p=200;            % Dimension

% Generate A, the original inverse covariance, with random sparsity pattern...
N=1;

nmat=20;

nzrange = [1 2 3];
 
% random seeds
s=[ 65 64 61 60 59  57 56 55 54 49  45 41 40 32 29 28 25 21 20 19   18 17 16 14 13 9 7 4 2 0 ; % p=100, nz=1
     67 66 65 63 62 61 59 58 56 55  54 49 45 41 40 32 28 25 21 20   19 18 16 14 13 9 7 4 2 0; % p=100, nz=2
     31 30 29 28 27 25 24 23 22 21  20 19 18 17 16 15 14 13 12 11   10  9 8  7  6 5 4 3 2 1      % p=100, nz=3
    ]; 
     
     
for j = 1:length(nzrange)
nz = nzrange(j);
nzfrac=nz/100;
nnzr = nzfrac*p*p;


for nm = 1:nmat
    sprintf('s=%d',s(j,nm))
    
    rand('seed',s(j,nm));
% Generate A, the original inverse covariance, with random sparsity pattern...
    A=eye(p);
    for i=1:nnzr
        A(ceil(p*rand),ceil(p*rand))=sign(rand()-.5);
    end
    A=A*A';  % A is the gound truth inverse covariance matrix
    A=(A'+A)/2;
    size(find(A))/(p*p)
    B = inv(A); % B is the ground-truth covariance matrix
    B=(B+B')/2;
    data = mvnrnd(zeros(N,p),B);
    EmpCov = (1/N)*data'*data;

    A = double(A);
    fnm=sprintf('Random/p%dmat%dpar%d.txt',p,nm,nz);
    save(fnm,'A','-ascii', '-tabs');
end

end
 

 

