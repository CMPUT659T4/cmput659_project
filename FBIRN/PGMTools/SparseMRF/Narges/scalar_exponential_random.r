library(MASS)
library(glasso)


## Initialization: 

file_name<-"results_exponential_scalar_random.txt"
fnm_table<-"res_table_exp_scalar_random.txt"
fnm_sum_table<-"res_sum_table_exp_scalar_random.txt"

write("%n density b b_old  lambda obj_iter lse TP FP ",fnm_table,append=TRUE);
write("%n density_avg b b_old_avg  lambda_avg obj_iter_avg lse  TP FP",fnm_sum_table,append=TRUE);


seed_r<-c(12,32,42,52,62,72,82,92,102,100,200,300,500,555,546,223,999,597,97,99,888,832,732,777,735,733,723,713,703,603)




N_vector<-c(30,50,500,1000)

%N_vector<-c(30,50,500,1000,400,2000,100)
b_vector<-c(1,1/0.1,1/0.01,1/0.005,1/0.001); #,1/0.0001);

vv<-1
for(vv in 1 : 4){
	for (b in b_vector) {
		b_old_avg<-0;
		lambda_avg<-0;
		glcalls_avg<-0;
		P_avg<-0
		FP_avg<-0
		TP_avg<-0
		density_avg<-0

		runs1<-0;
		runs<-1;

		for(runs in 1:5){
			error_flag<-0
			set.seed(seed_r[runs])
			ss<-100

			# generate random matrix
			nnz<-4*ss*runif(1)

			A<-diag(ss)
			for(i in 1:nnz){
				A[ss*runif(1),ss*runif(1)]=sign(runif(1)-.5);
			}

			A1<-t(A)%*%A

			Sig=solve(A1) # Sig is the ground-truth covariance matrix
                  #Sig = (Sig + t(Sig))/2; # to make sure Sig is symmetric (up to high precision)
 

			mtest<-100;

			x_test<-mvrnorm(mtest, mu= c(rep(0,ss) ), Sig )

			dd<-N_vector[vv];
			data<-mvrnorm(dd,rep(0,ss),Sig)
			x<-matrix(data,ncol=ss)
			s<- var(x)


			# compute b
			eps<-0.001
			Sinv<-solve(s + eps*diag(1,ss))
			b_old<-0
			b_old<-sum(abs(Sinv))
			b_old<-b_old/(ss*ss-1) # changed by IR from ss*ss

			#b_old<-b_old*(ss/dd) # alternatively, scale b

			A1_nnz<-0
			for(i in 1:ss)
				for(j in 1:ss)
					if(A1[i,j]!=0) A1_nnz<-A1_nnz+1
			print(A1_nnz)

			N<-dd
			m<-1# number of iterations
			obj_opt<-matrix(rep(0,m))
			obj_lambda_opt<-matrix(rep(0,m),ncol=m)
			obj_cov_opt<-matrix(rep(0,ss*ss*m) ,ncol=ss*ss)
			obj_iteration<-matrix(rep(0,m))
			obj_nnz<-matrix(rep(0,m))
			obj_trace<-matrix(rep(0,m*500),ncol=500)

			glcalls<-0; # added by Irina
			for(kk in 1:m){
				print("******************************ITERATION: ******************** ")
				print(kk)
				done<-0
				objective<--10000000
				objective_last<--10000000
				lambda<- kk*5
				last_lambda<- kk*5
				new_lambda<- kk*5
				it<-0


				### SEARCH:	
				while(done==0){
					it<-it+1
					objective_last<-objective
					a<-glasso(s, rho=(2/N) *lambda ,thr=1.0e-3 ,penalize.diagonal=TRUE)
					glcalls = glcalls +1; # added by Irina
					norm1<-0
					for (i in 1:ss){
						for(j in 1:ss){
							norm1<-norm1+ abs(a$wi[i,j])
						}
					}
					mm<-0
					for(i in 1:dd){
						mm<-mm+( x[i,] %*% a$wi %*% x[i,] )
					}

		
					ll<- ss*ss*log(lambda/2)- lambda* norm1-b*lambda
	

					objective<-  - 0.5 * mm + dd * 0.5 * log ( det(a$wi)  )  + ll
					if(det(a$wi)<10^-320)
						{error_flag<-1}
					else if( (objective > objective_last) || it==1){
						print("good") 
						print(it)
						new_lambda<-(ss*ss)/(norm1 + b)   ##( (ss*ss/last_lambda) - norm11 )		
					}
					else{
						print(" **************************backtrachking **************************************")
						alpha<-1  
						while(objective < objective_last && alpha>0.001){	      
							##print(alpha)
							new_lambda<-last_lambda + alpha * (new_lambda - last_lambda)    ##( (ss*ss/last_lambda) - norm11 )		
			
							alpha<-alpha/2;
							a<-glasso(s, rho=(2/N) *new_lambda ,thr=1.0e-3 ,penalize.diagonal=TRUE)
							glcalls = glcalls +1; # added by Irina

							norm1<-0
							for (i in 1:ss){
								for(j in 1:ss){
									norm1<-norm1+ abs(a$wi[i,j])
								}	
							}
							mm<-0
							for(i in 1:dd){
								mm<-mm+( x[i,] %*% a$wi %*% x[i,] )

							}
							ll<- ss*ss*log(new_lambda/2)- new_lambda * norm1 -b*new_lambda

							objective<-  - 0.5 * mm + dd * 0.5 * log ( det(a$wi)  )  + ll
						} #end while
					}#end if
					last_lambda<-lambda
					lambda<-new_lambda
					obj_trace[kk,it]<-objective
					if(error_flag==1){
						done<-1
					}
					else if( objective - objective_last  < 0.01){
		 				done<-1
					}
	
					print("objective:")
					print(objective)
				} #end while


				for (g in 1:ss){
					for(gg in 1:ss){
						if(a$wi[g,gg]!=0) obj_nnz[kk]<-obj_nnz[kk]+1
					}
				}
				obj_opt[kk]<-objective
				obj_lambda_opt[kk]<-lambda
	
				obj_cov_opt[kk,]<-a$wi
				obj_iteration[kk]<- it
				glcalls_avg = glcalls + glcalls_avg; # added by Irina
			} # end for kk

			if(error_flag==0){
				runs1<-runs1+1

				A1_nnz<-0
				for(i in 1:ss)
					for(j in 1:ss)
						if(A1[i,j]!=0) A1_nnz<-A1_nnz+1
 				A1_nnz

 				nnz<-0
				p2<-ss*ss
 				for(i in 1:p2){
					##for(j in 1:ss)
					if(abs(obj_cov_opt[m,i])>0.0){ nnz<-nnz+1}
				}
 				nnz

				fp<-0
				for(i in 1:ss)
					for(j in 1:ss)
						if(A1[i,j]==0 && obj_cov_opt[m,(i-1)*ss+j]!=0) fp<-fp+1

				fp

				tp<-nnz-fp

				tp<- (tp-ss) / (A1_nnz - ss)

				fp<- fp / (ss*ss - A1_nnz)

				lse<-0
				error<-0
				n<-ss
				for(f in 1:n){
					sigma12<-matrix(rep(0,n-1), nrow=1)
					sigma22<- matrix(rep(0,(n-1)*(n-1) ), nrow=n-1, ncol= n-1)

					sigma12<-a$w[f,-f]; ## the covariance of variable f with the rest of variables
					sigma22<-a$w[-f,-f];
					sigma22inv<-solve(sigma22)
					sigma12<-as.matrix(sigma12, nrow=n-1, ncol=1)
					sigma22inv<-as.matrix(sigma22inv, nrow=n-1, ncol=n-1)
					y_1<-matrix(rep(0,(n-1)*(n-1)), nrow=n-1, ncol=n-1)
					y_test<-0
					for (i in 1:mtest){
						xval<-as.matrix(t(x_test[i,-f]), nrow=n-1, ncol=1)
						y1<-t(sigma12) %*% sigma22inv
						y_test<- y1 %*% t(xval)
						error = error + ( y_test - x_test [i,f] )^2
					}
				}
				lse<-error/ (mtest * n)

				b_old_avg<-b_old_avg + b_old

				lambda_avg<-lambda_avg + obj_lambda_opt

				P_avg<-P_avg + lse

				TP_avg<- TP_avg + tp

				FP_avg<- FP_avg + fp


				### write the results in the file:


 				write("----------------------------------------------------------",file_name,append=TRUE)
 				write("P=100  NOT Scaled  and no normalization of the data ",file_name,append=TRUE)
 
				write("density= ",file_name,append=TRUE)

 				write(A1_nnz/(ss*ss),file_name,append=TRUE)

 				write("N= ",file_name,append=TRUE)

 				write(dd,file_name,append=TRUE)

 				write("b= ",file_name,append=TRUE)

 				write(b,file_name,append=TRUE)

 				write("learned lambda",file_name,append=TRUE)
 				write(obj_lambda_opt,file_name,append=TRUE)

 				write("Prediction Mean Squared Error of the learned model ",file_name,append=TRUE)
 				write(lse,file_name,append=TRUE)

				write("TP of the learned C ",file_name,append=TRUE)
 				write(tp,file_name,append=TRUE)

				write("FP of the learned C ",file_name,append=TRUE)
 				write(fp,file_name,append=TRUE)

				# added by Irina: write a table
				density=A1_nnz/(ss*ss);
				density_avg<- density_avg + density

				res<-c(dd,density,b,b_old,obj_lambda_opt,glcalls,lse,tp,fp);
				write(res,fnm_table,9,append=TRUE);

			} #end  of if(error_flag==0)
		} #end for runs
		write("----------------------------------------------------------",file_name,append=TRUE)

 		write("lambda_avg",file_name,append=TRUE)
 		write(lambda_avg/runs1,file_name,append=TRUE)
 
 		write("lse_avg",file_name,append=TRUE)
 		write(P_avg/runs1,file_name,append=TRUE)

 		write("TP_avg",file_name,append=TRUE)
 		write(TP_avg/runs1,file_name,append=TRUE)

 		write("FP_avg",file_name,append=TRUE)
 		write(FP_avg/runs1,file_name,append=TRUE)


		# added by Irina: write a summary table
      	glcalls_avg = glcalls_avg/runs1
		b_old_avg = b_old_avg/runs1;
		l_avg = lambda_avg/runs1;
		lse_avg = P_avg/runs1;
		tp_avg = TP_avg/runs1;
		fp_avg = FP_avg/runs1;
                density_avg = density_avg/runs1;
		res<-c(dd,density_avg,b,b_old_avg,l_avg,glcalls_avg,lse_avg,tp_avg,fp_avg);
		write(res,fnm_sum_table,9,append=TRUE);
	} # for b
} # for vv



































































































	








