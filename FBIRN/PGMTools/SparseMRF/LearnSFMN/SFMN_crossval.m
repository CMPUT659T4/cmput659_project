function LearnSFMN
% 1. Generate a scale-free Gaussian Markov network
% 2. Sample N  instenaces

% use B-A Scale-Free Network Generation and Visualization package by
% by Mathew Neil George

% seed network
seed =[0 1 0 0 1;1 0 0 1 0;0 0 0 1 0;0 1 1 0 0;1 0 0 0 0]
% generate a scale-free network's adjacency matrix

p=100; % number of nodes/variables

Net = SFNG(p, 2, seed);

%maximum value for cov that does not violate pos-def constaint depends on p
const = 0.2; % for p=100; covariance value for each pair of nodes setting to 0.15 and larger violates positive-definiteness
%const = 0.13; % for p=300; covariance value for each pair of nodes setting to 0.15 and larger violates positive-definiteness

% generate covariance matrix C = L'L
S = Net*const + eye(p);
%take  Cholesly factorization of S=C'C 
C=chol(S);
    
true_pos = [];
true_neg = [];
% total # of positves (1s) and negatives (zeros)
total_zeros = size(find(Net == 0),1);
total_ones = p*p - total_zeros;
    
% for different number of samples
ii = 0;
inds = [50:50:250];
for n = inds 
    Our_graph=[];
    neighborhood = [];
    beta = [];
    data=[];
    ii = ii +1;
    % generate Gaussian entries with zero mean and variance 1 in a nxp matrix Z
    Z = random('norm',0,1,[n,p]);
    ZZ = normalize(Z);

    %and compute data = ZC as your data
    data= Z*C;
    %you'll still get zero means for each column i in X and the correponding
    %variable X_i will be a mixture of independent Gaussian variables corresponding to non-zero entries in the i-th column of S;
    %thus we get dependencies among columns of X that correspond to non-zero-pattern inS;
    %then we can normalize X to get unit variance.

    % perform Lasso (using LARS w/ lasso modification) for each node,
    % obtain all solutions
    
 
    for i=1:p
        y = data(:,i); % i-th variable is now the response variable 
        if i==1
            X = data(:,[2:p]);
        elseif i==p
            X = data(:,[1:(p-1)]);
        else
            X = [data(:,[1:(i-1)]) data(:,[(i+1):p])];
        end
  
        tt = floor(0.7*n);
        rand('state',0);
        
        best_beta = [];
        best_ind=[];
        maxrep = 10;
        for rep=1:maxrep
            % randomize the data and 
            % split data into 70% training and 30% cross-validation set
            indrnd = randperm(n);
            X = X(indrnd, :);
            y = y(indrnd);
 
            X_train = X(1:tt,:);
            y_train = y(1:tt,:);
            X_test = X((tt+1):n,:);
            y_test = y((tt+1):n);
        
            % to run Lasso, normalize X (so the variables have unit variance) and center y
            X_train = normalize(X_train);
            y_train=center(y_train);

            beta = lars(X, y, 'lasso', 0, 0, [], 1);

            % choose the least-squares solution over all iterations 
            L2err=[];
            X_test = normalize(X_test);
            y_test=center(y_test);
            for j=1:size(beta,1)
                er = y_test-X_test*beta(j,:)';
                L2err(j) = er'*er;
            end
            [err,ind]=min(L2err);

            % this solution defines the neighborhood of i-th node using adaptive method (ours)
            %create row in the new adjacency matrix (include the node itself with 0
            if i==1
                sol = [0 beta(ind,:)];
            elseif i==p
                sol = [beta(ind,:) 0];
            else
                sol = [beta(ind,[1:(i-1)]) 0 beta(ind,[i:(p-1)])];
            end
            best_beta = [best_beta; sol]; 
            best_ind = [best_ind; ind];
        end
        possible_neighborhoods = (best_beta ~= 0);
        ANDneighborhood =  (sum(possible_neighborhoods) == maxrep); % all folds agree on this edge
        neighborhood = [neighborhood;  ANDneighborhood];
            
    end

    Our_graph = neighborhood' & neighborhood; % choose edge if both nodes agree on it
    
    % compute true postive (precision) and true negative (recall) ratios
    true_pos(ii) = size(find(Our_graph & Net),1)/total_ones;
    true_neg(ii) = size(find(~Our_graph  & ~Net),1)/total_zeros;
end

 
figure(1);
plot(inds,true_pos, 'bx-',inds, true_neg, 'ro-');
ylabel('Accuracy','fontsize',16);
xlabel('number of samples','fontsize',16);
[s,e]=sprintf('AND-rule: graphs with %d nodes', p); 
title(s,'fontsize',12);
legend('true positives (precision)','true negatives (recall)');

figure(2);
plot(inds,true_neg);
ylabel('True positives \% (precision)','fontsize',16);
xlabel('number of samples','fontsize',16);
[s,e]=sprintf('AND-rule: graphs with %d nodes', p); 
title(s,'fontsize',12);

% map each constraint t to corresponding value of lambda (knowing opt solution)

% choose best solution for each node (our method)

% choose fixed-lambda solution where lambda is suggested by  MB(meinhausen&bulman)  

% for each method, use OR and AND rules to reconstruct edges

% compute precision and recall (true positive and true negative errs)

 
