function LearnSFMN
% 1. Generate a scale-free Gaussian Markov network
% 2. Sample N  instenaces

% use B-A Scale-Free Network Generation and Visualization package by
% by Mathew Neil George

% seed network
seed =[0 1 0 0 1;1 0 0 1 0;0 0 0 1 0;0 1 1 0 0;1 0 0 0 0]
% generate a scale-free network's adjacency matrix

p=100; % number of nodes/variables
max_trials = 20;
inds = [100:100:1000];
th_inds=[0.05:0.02:0.15];

corr_true_pos = zeros(length(inds),length(th_inds),max_trials);
corr_true_neg = zeros(length(inds),length(th_inds),max_trials); 
L1_true_pos = zeros(length(inds),length(th_inds),max_trials);
L1_true_neg = zeros(length(inds),length(th_inds),max_trials); 
 
for trial = 1:max_trials
Net = SFNG(p, 2, seed);

%maximum value for cov that does not violate pos-def constaint depends on p
const = 0.15; % for p=100; covariance value for each pair of nodes setting to 0.15 and larger violates positive-definiteness
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


for n = inds 
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

    % compute adjacency matrix based on thresholded correlations
    corr_mat = corrcoef(data) - eye(p,p);
    % compute true postive (precision) and true negative (recall) ratios
    i1 = 0;

    for threshold = th_inds
        i1 = i1 + 1;
        Corr_graph = (abs(corr_mat) > threshold);
        corr_true_pos(ii,i1,trial) = size(find(Corr_graph & Net),1)/total_ones;
        corr_true_neg(ii,i1,trial) = size(find(~Corr_graph  & ~Net),1)/total_zeros;
    
        % plot the dependence of reconstruction accuracy based on thresholding
        % the (empirical) correlation matrix
   
        % perform Lasso (using LARS w/ lasso modification) for each node,
        % obtain all solutions
        best_beta=[];
        best_ind=[];
        for i=1:p
            y = data(:,i); % i-th variable is now the response variable 
            if i==1
                X = data(:,[2:p]);
            elseif i==p
                X = data(:,[1:(p-1)]);
            else
                X = [data(:,[1:(i-1)]) data(:,[(i+1):p])];
            end

            num = sum(Corr_graph(i,:),2); % bound on the number of neighbors according to the thresholded corr matrix
            beta = lars(X, y, 'lasso', -num, 0, [], 1);

            % choose the least-squares solution over all iterations 
            L2err=[];
            X = normalize(X);
            y=center(y);
            for j=1:size(beta,1)
                er = y-X*beta(j,:)';
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
        neighborhood = (best_beta ~= 0);
      
        Our_graph = neighborhood' & neighborhood; % choose edge if both nodes agree on it
    
        % compute true postive (precision) and true negative (recall) ratios
        L1_true_pos(ii,i1,trial) = size(find(Our_graph & Net),1)/total_ones;
        L1_true_neg(ii,i1,trial) = size(find(~Our_graph  & ~Net),1)/total_zeros;
    end
end
end

total_corr_true_pos = sum(corr_true_pos,3)/max_trials;
total_corr_true_neg = sum(corr_true_neg,3)/max_trials;
total_L1_true_pos = sum(L1_true_pos,3)/max_trials;
total_L1_true_neg = sum(L1_true_neg,3)/max_trials;

 i1 = 0;
 for threshold = th_inds
     i1 = i1 + 1;
     figure(i1);
    plot(inds,total_corr_true_pos(:,i1), 'bx-',inds, total_corr_true_neg(:,i1), 'ro-');
    ylabel('Accuracy','fontsize',16);
    xlabel('number of samples','fontsize',16);
    [s,e]=sprintf('Thresholded-Correlation: threshold = %.2f', threshold); 
    title(s,'fontsize',12);
    legend('true positives (precision)','true negatives (recall)');
    
    figure(100+i1);
plot(inds,total_L1_true_pos(:,i1), 'bx-',inds, total_L1_true_neg(:,i1), 'ro-');
ylabel('Accuracy','fontsize',16);
xlabel('number of samples','fontsize',16);
[s,e]=sprintf('Lasso with Correlation-Prior(AND-rule): threshold = %.2f', threshold); 
title(s,'fontsize',12);
legend('true positives (precision)','true negatives (recall)');
 end
% map each constraint t to corresponding value of lambda (knowing opt solution)

% choose best solution for each node (our method)

% choose fixed-lambda solution where lambda is suggested by  MB(meinhausen&bulman)  

% for each method, use OR and AND rules to reconstruct edges

% compute precision and recall (true positive and true negative errs)

 
