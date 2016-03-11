% Test script for Covariance selection code.
% Initialize parameters ****************
% Data ....

p=100;            % Dimension
 
% 1. Generate a scale-free Gaussian Markov network
% 2. Sample N  instenaces

% use B-A Scale-Free Network Generation and Visualization package by by
% Mathew Neil George
% seed network
   
 seed_id = 1;
 seed =[0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1; 
        1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 
        1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;  
        1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 
        1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;        
        1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;    
        1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 
        1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;         
        1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 
        1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;         
        1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;        
        1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;         
        1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;         
        1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;         
        1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;       
        1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;        
        1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;         
        1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;         
        1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; ]        
        
        
        % generate a scale-free network's adjacency matrix

% Generate A, the original inverse covariance, with random sparsity pattern...
N=1;
networks = 10;
nmat = 10;
mlinks=3 %sp=5% %mlinks = 13; %sp=21% %mlinks = 19; %sp=30%   

%maximum value for cov that does not violate pos-def constaint depends on p
const = 0.1; % for p=300; covariance value for each pair of nodes setting to 0.15 and larger violates positive-definiteness

rand('state',0);

for mlinks = [3 13 19]
mlinks
    for ns = 1:networks
    Net = SFNG(p, mlinks, seed);
    spar = floor(100*size(find(Net ~=0),1)/(p*p));%sparsity - % of nonzeros
    for nm = 1:nmat
        % generate sparse inverse covariance matrix (capturing independence relationships
        A =  eye(p) + Net.*const; %(const*(rand(p)));   
%         w = floor(2*rand(p));
% %         w = sign(w+w');
%         A =  eye(p) + Net.*(const*w);       
 
        B = inv(A); % B is the ground-truth covariance matrix
        B = (B+B')/2;
        %chol(B);
        data = mvnrnd(zeros(N,p),B);
        EmpCov = (1/N)*data'*data;
        
                 
%figure(1);hist(sum(Net));
% figure(2);hist(sum(abs(A)));
% figure(3);hist(sum(abs(B)));
        %save(sprintf('net%dmat%d',ns,i),'A');
        A = double(A);
        fnm=sprintf('SFNetworks/mat%dnet%dsp%dml%d.txt',nm,ns,spar,mlinks);
        save(fnm,'A','-ascii', '-tabs');
 
    end
end
end 

 
%%%%%%
% % % seed_id = 1;
% % % seed =[0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1; 
% % %        1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 
% % %        1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;  
% % %        1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0;  
% % %        1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;  
% % %        1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0;  
% % %        1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;  
% % %        1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0;  
% % %        1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;  
% % %        1 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0;
% % %        1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
% % %        1 0 0 0 0 0 0 0 0 1 0 0 0 0 1 0;
% % %        1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0;
% % %        1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
% % %        1 0 0 1 0 1 0 1 0 1 0 1 1 0 0 0;
% % %        1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; 

