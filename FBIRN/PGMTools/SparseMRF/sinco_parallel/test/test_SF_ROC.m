% Test script for SINCO - sparse inverse covariance selection
% Scale-free networks

ifvec=0; % vector formulation(1) or scalar formulation(0)
usediag=0; % 0 - regularize diagonal, 1 - "regularized likelihood"


%     sumtabf = 'res_sum_table_ROC_SF_p100sp1_SINCO_10.txt';
%     tabf = 'res_table_ROC_SF_p100sp1_SINCO_10.txt';  
%     fid_sum_table = fopen(sumtabf,'a');
%     fprintf(fid_sum_table,'%% N density lambda LLtest TP FP\n'); 
%     fid_table = fopen(tabf,'a');
%     fprintf(fid_table,'%% N density lambda LLtest TP FP\n');

    sumtabf = 'res_sum_table_ROC_SFblock_SINCO_10.txt';
    tabf = 'res_table_ROC_SFblock_SINCO_10.txt';  
    fid_sum_table = fopen(sumtabf,'a');
    fprintf(fid_sum_table,'%% n p density lambda LLtest TP FP\n'); 
    fid_table = fopen(tabf,'a');
    fprintf(fid_table,'%% n p density lambda LLtest TP FP\n');
    
    lambda_range<-c(300,250,200,150,100,40,39,38,37,36,35,34,33,32,31,30,29,28,27,26,25,24,23,22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0.5,0.1,0.05,0.01)

 
p=100;  % number of variables

%Nvec = [floor(p/2) p 2*p 5*p 10*p 50*p 100*p 1000*p];
Nvec = [50 100 200 500 1000 5000 10000];
spvec = [5 21 30];
mlvec = [3 13 19];
 
%bvec = [0.01 0.1 1 10 100 500 1000 5000 10000];
%lambda_range = [500 450 400 350 300 250 200 150 100 50 40 30 20 10 5 1];
 
 lambda_range = [300 250 200 150 100 40 39 38 37 36 35 34 33 32 31 30 29 28 27 26 25 24 23 22 21 20 19 18 17 16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1 0.5 0.1 0.05 0.01];

 
test_ll_avg = zeros(length(Nvec),length(spvec),length(lambda_range));
TP_avg = zeros(length(Nvec),length(spvec),length(lambda_range));
FP_avg = zeros(length(Nvec),length(spvec),length(lambda_range));

for  ni = 1:length(Nvec) % number of samples
  N= Nvec(ni);
  for si=1:1
      runs1 = 0;  % total number of runs to average over  
      for mats=1:10
        for nets=1:10 
          runs1 = runs1+1;
                A=load(sprintf('../SFNetworks/mat%dnet%dsp%dml%d.txt',mats,nets,spvec(si),mlvec(si)));

          % total # of positives (1s) and negatives (zeros)
          total_zeros = size(find(A == 0),1);
          total_ones = p*p - total_zeros - p; % don't count the diagonal elements

          B = inv(A); % B is the ground-truth covariance matrix
          B=(B+B')/2;
          data = mvnrnd(zeros(N,p),B);
          test_data = mvnrnd(zeros(N,p),B);

          EmpCov = (1/N)*data'*data;
                
          tol=0.000001;
          Cstart=eye(p);
          Wstart=eye(p);
          EC=N*0.5*EmpCov;
          fstart= - trace(EC*Cstart);
          K=N*0.5;
          precision=0.01;
          Sbase=ones(p,p);
          if (usediag) 
                    Sbase=Sbase-eye(p);
          end
 
          % try lambda range

          i=0;
          lambdaold = 0;
          for lambda=lambda_range
                i  = i+1
                 
                S=Sbase*lambda;
                fstart=fstart -sum(sum(abs((lambda-lambdaold)*Sbase.*Cstart)));

                [Csol, Wsol, fsol]=...
                    sinco(Cstart, fstart, Wstart, EC, Sbase, lambda, K, tol);
                
                Umat = Csol;
                Cstart=Csol;
                fstart=fsol;               
%                 fall(i)=fsol;
%                 lambdac(i)=lambda*sum(sum(abs(Sbase.*Csol)));
%                 normC(i)=sum(sum(abs(Sbase.*Csol)));
%                 tra(i)=trace(EC*Csol);
                Wstart=Wsol;
                lambdaold=lambda;

                true_pos(i) = size(find(Csol & A),1)-p; % don't count diagonal elements
                true_neg(i) = size(find(~Csol  & ~A),1);
 
                TP(i) = true_pos(i)/total_ones
                TN(i)  = true_neg(i)/total_zeros
                FP(i) = 1-TN(i)
                
                TP_avg(ni,si,i) = TP_avg(ni,si,i) +TP(i);
                FP_avg(ni,si,i) = FP_avg(ni,si,i) +FP(i);
 
    
                test_ll = 0; % temporary
                test_ll_avg(ni,si,i)  = test_ll_avg(ni,si,i) + test_ll;
                density = total_ones/(p*p - p); % off-diagonal density  
 
                fprintf(fid_table,'%d %.2f %f  %f %f %f\n', N, density, lambda, test_ll, TP(i), FP(i));
          end % for lambda
         
        end % nets 
      end %mats
      
      i=0;
      for lambda=lambda_range
        i  = i+1

        TP_avg(ni,si,i) = TP_avg(ni,si,i)/runs1;
        FP_avg(ni,si,i) = FP_avg(ni,si,i)/runs1;        
        test_ll_avg(ni,si,i)  = test_ll_avg(ni,si,i)/runs1;
    
        fprintf(fid_sum_table,'%d %.2f %f   %f %f %f \n', N, density, lambda,  test_ll_avg(ni,si,i), TP_avg(ni,si,i), FP_avg(ni,si,i));
      end % lambda
        %        keyboard
  end %si=1:3
end % ni = 1:length(Nvec)

fclose(fid_sum_table);
fclose(fid_table);
%fclose(fid);

%keyboard
 

        
                
 
