% Test script for SINCO - sparse inverse covariance selection
% Scale-free networks

ifvec=0; % vector formulation(1) or scalar formulation(0)
usediag=0; % 0 - regularize diagonal, 1 - "regularized likelihood"


if (ifvec)
    sumtabf = 'res_sum_table_exp_vector_l_fixed_SF.txt';
    tabf = 'res_table_exp_vector_l_fixed_SF.txt';
    fid_table = fopen(tabf,'a');
    fprintf(fid_table,'%% N density LLtest TP FP\n');   
    fid_sum_table = fopen(sumtabf,'a');
    fprintf(fid_sum_table,'%% N density LLtest TP FP \n');      
else
    sumtabf = 'res_sum_table_exp_scalar_l_fixed_SF.txt';
    tabf = 'res_table_exp_scalar_l_fixed_SF.txt';  
    fid_sum_table = fopen(sumtabf,'a');
    fprintf(fid_sum_table,'%% N density b b_old_avg lambda_avg LLtest TP FP \n');   
    fid_table = fopen(tabf,'a');
    fprintf(fid_table,'%% N density b b_old lambda LLtest TP FP\n');
    
    
end

%fid = fopen('res_exp_scalar_SF.txt','a');


p=100;  % number of variables



Nvec = [floor(p/2) p 2*p 5*p 10*p 50*p 100*p 1000*p];
spvec = [5 21 30];
mlvec = [3 13 19];
%bvec = [0.01 0.1 1 10 100 500 1000 5000 10000];
bvec=[1];

rand('state',0);       

lambda_avg = zeros(length(Nvec),length(spvec)); %include b extimated from data
b_old_avg = zeros(length(Nvec),length(spvec));
test_ll_avg = zeros(length(Nvec),length(spvec));
TP_avg = zeros(length(Nvec),length(spvec));
FP_avg = zeros(length(Nvec),length(spvec));

for  ni = 1:length(Nvec) % number of samples
  N= Nvec(ni);
  for si=1:3
        runs1 = 0;  % total number of runs to average over  
        
        for mats=1:5
            for nets=1:5
                runs1 = runs1+1;
                A=load(sprintf('../SFNetworks/mat%dnet%dsp%dml%d.txt',mats,nets,spvec(si),mlvec(si)));

                % total # of positives (1s) and negatives (zeros)
                total_zeros = size(find(A == 0),1);
                total_ones = p*p - total_zeros - p; % don't count the diagonal elements

                B = inv(A); % B is the ground-truth covariance matrix
                %B=(B+B')/2;
                data = mvnrnd(zeros(N,p),B);
                test_data = mvnrnd(zeros(N,p),B);

                EmpCov = (1/N)*data'*data;

                if (ifvec) 
                    b_old = sum(abs(inv(EmpCov + 0.001*eye(p,p))))/(p);
                else
                    b_old = sum(sum(abs(inv(EmpCov + 0.001*eye(p,p)))))/(p^2-1)
                end

                %b= b_old; %%%%%%%%%%%%%%%%%%%%%%%%%%%% temporarily
                %commented out
                
                b = 1;
                
                tol=0.000001;
                Cstart=eye(p);
                Wstart=eye(p);
                EC=N*0.5*EmpCov;
                %load ECmat
                fstart= - trace(EC*Cstart);
                K=N*0.5;
                %Sbase=rand(p,p);
                precision=0.01;
                %Sbase=Sbase+Sbase';
                Sbase=ones(p,p);
                if (usediag) 
                    Sbase=Sbase-eye(p);
                end
                %for i=1:p
                %  Sbase(i,i)=0.0001;
                %end
                %lambdaold=0;

                %solution_path = zeros(size(rho_range,2),p,p);

                %figure(1);
                %subplot(1,4,1);colorspy(A);xlabel('true inverse cov');pbaspect('manual');set(gca,'XTickLabel',{' '});set(gca,'YTickLabel',{' '});
                %subplot(1,4,3);colorspy(inv(EmpCov));xlabel('Inverse of EmpCov');pbaspect('manual');set(gca,'XTickLabel',{' '});set(gca,'YTickLabel',{' '});

                % % %%%%%%%%%%% Test alternating minimization for \rho selection
                if (ifvec)
                    lstart=100*(ones(p,1));
                else
                    lstart=100;
                end
                %[C,our_lambda] = AltMin_vec(EmpCov, N, Sbase, ifvec, b, usediag, precision, tol, lstart);
                
                % try fixed lambda
                
                our_lambda = 1/b_old;
                
                C = sinco_lambda(EmpCov, N, Sbase, ifvec, precision, tol, our_lambda);
                
                true_pos = size(find(C & A),1)-p; % don't count diagonal elements
                true_neg = size(find(~C  & ~A),1);
 
                TP = true_pos/total_ones
                TN  = true_neg/total_zeros
                FP = 1-TN
                
                TP_avg(ni,si) = TP_avg(ni,si) +TP;
                FP_avg(ni,si) = FP_avg(ni,si) +FP;
                
                test_ll = 0; % temporary
                test_ll_avg(ni,si)  = test_ll_avg(ni,si) + test_ll;
                density = total_ones/(p*p - p); % off-diagonal density 
                
                if (ifvec)            
                    fprintf(fid_table,'%d %.2f %f %f %f\n', N, density, test_ll, TP, FP);     
                else                   
                    lambda_avg(ni,si) = lambda_avg(ni,si) + our_lambda;            
                    b_old_avg(ni,si) = b_old_avg(ni,si) + b_old;
                    fprintf(fid_table,'%d %.2f %f %f %f %f %f %f\n', N, density, b, b_old, our_lambda, test_ll, TP, FP);
                end
                %keyboard
            end %nets
        end % mats                
        TP_avg(ni,si) = TP_avg(ni,si)/runs1;
        FP_avg(ni,si) = FP_avg(ni,si)/runs1;        
        test_ll_avg(ni,si)  = test_ll_avg(ni,si)/runs1;
        
        if (ifvec)            
            fprintf(fid_sum_table,'%d %.2f %f %f %f\n', N, density, test_ll_avg(ni,si), TP_avg(ni,si), FP_avg(ni,si));     
        else
            lambda_avg(ni,si) = lambda_avg(ni,si) /runs1;
            b_old_avg(ni,si) = b_old_avg(ni,si)/runs1; 
            fprintf(fid_sum_table,'%d %.2f %f %f %f %f %f %f\n', N, density, b, b_old_avg(ni,si), lambda_avg(ni,si), test_ll_avg(ni,si), TP_avg(ni,si), FP_avg(ni,si));
        end
        %        keyboard
  end %si=1:3
end % ni = 1:length(Nvec)

fclose(fid_sum_table);
fclose(fid_table);
%fclose(fid);

keyboard
 

        
                
 
