library(MASS)
library(glasso)


## Initialization: 

file_name<-"res_CV_SF.txt"
fnm_table<-"res_table_CV_SF.txt"
fnm_sum_table<-"res_sum_table_CV_SF.txt"

write("%n density lambda_avg lse TP FP",fnm_sum_table,append=TRUE);
write("%n density lambda lse TP FP ",fnm_table,append=TRUE);

sp<-c(5,21,30);
ml<-c(3,13,19);

N_vector<- c(30,50,500,1000)

vv<-1
for(vv in 1 : 4){  ## over 4 values of N
	for(sl in 1:3){
		lambda_avg<-0
		P_avg<-0
		FP_avg<-0
		TP_avg<-0

		runs1<-0;
		runs<-1;

	      ss<-100; # the number of variables (p)
################### average over networks (nets) and matrices(mats), for each fixed sparsity level sl and number of samples vv
write("%n density lambda lse TP FP ",fnm_table,append=TRUE);

         	for(nets in 1:5){
			for(mats in 1:5){
				error_flag<-0

				name<-paste("SFNetworks/mat",mats,"net",nets,"sp",sp[sl],"ml",ml[sl],".txt", sep="")
				A1 <- matrix(scan(name, 0), ncol=ss, byrow=TRUE)

				Sig=solve(A1)

				mtest<-100;

				x_test<-mvrnorm(mtest, mu= c(rep(0,ss) ), Sig )

				dd<-N_vector[vv];
				data<-mvrnorm(dd,rep(0,ss),Sig)
				x<-matrix(data,ncol=ss)
				s<- var(x)

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

				lambda_range<-c(0.1,0.2,0.5,1,5,10,20,30,50,100,200,300,400,500,1000,2000)
				obj_cv<-rep(0,16)
				for (index in 1:16){
					gg<-lambda_range[index]
					lambda<-gg

##############################             Cross validation:         ##########################
## assume k-fold cross validation where k=5
#print("CV")
					k<-5
					k_d<-dd/k
					avg_likelihood<-0
					test_indicator<-matrix(rep(0,dd));
					train<-dd-k_d
					train_data<-matrix(rep(0,train*ss),ncol=ss);
					mc<-0
					for(f in 1:k){
						test_indicator<-matrix(rep(0,dd));
						for(ff in 1:k_d){
							test_indicator[ff + mc] <-1
						}
						t<-1;
						for(i in 1:dd){
							if(test_indicator[i]==0){
								train_data[t,]<-x[i,];
								t<-t+1;
							}
						}
						sc<- var(train_data)
						aa<-glasso(sc, rho= lambda *( 2/train ),penalize.diagonal=TRUE)
						mm<-0
						for(i in 1:dd){
							if(test_indicator[i]!=0){
								mm<-mm+( x[i,] %*% aa$wi %*% x[i,] )
							}
						}
						norm1c<-0
						for (ii in 1:ss){
							for(jj in 1:ss){
								norm1c<-norm1c+ abs(aa$wi[ii,jj])
							}
						}
						avg_likelihood<- avg_likelihood - 0.5 * mm + k_d* 0.5 * log ( det(aa$wi)  ) - lambda * norm1c
						mc<-mc+k_d
					}
					avg_likelihood<- avg_likelihood/k
					obj_cv[index]<-avg_likelihood

###########################################################################################
				}  
		 
				lambda_cv<-0
				max_cv<-max(obj_cv)
				for(i in 1:16){
					if(obj_cv[i]==max_cv){
						lambda_cv<-lambda_range[i]
						print(lambda_cv)
					}
				}

				ag<-glasso(s, rho= lambda_cv *( 2/dd ),penalize.diagonal=TRUE)

				if(error_flag==0){
					runs1<-runs1+1

					A1_nnz<-0
					for(i in 1:ss)
						for(j in 1:ss)
							if(A1[i,j]!=0) A1_nnz<-A1_nnz+1
 					A1_nnz

					nnz<-0
					for(i in 1:ss){
						for(j in 1:ss)
							if(abs(ag$wi[i,j])>0.0){ nnz<-nnz+1}
					}
 					nnz

					fp<-0
					for(i in 1:ss)
						for(j in 1:ss)
							if(A1[i,j]==0 && ag$wi[i,j]!=0) fp<-fp+1

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

						sigma12<-ag$w[f,-f]; ## the covariance of variable f with the rest of variables
						sigma22<-ag$w[-f,-f];
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

					lambda_avg<-lambda_avg + lambda_cv

					P_avg<-P_avg + lse

					TP_avg<- TP_avg + tp

					FP_avg<- FP_avg + fp

					### write the results in the file:

 					write("----------------------------------------------------------",file_name,append=TRUE)
 					write("P=100  with Banerjee ",file_name,append=TRUE)
					write("density= ",file_name,append=TRUE)
					write(A1_nnz/(ss*ss),file_name,append=TRUE)
					write("N= ",file_name,append=TRUE)
					write(dd,file_name,append=TRUE)
 					write("cv lambda",file_name,append=TRUE)
 					write(lambda_cv,file_name,append=TRUE)

 					write("Prediction Mean Squared Error of the learned model ",file_name,append=TRUE)
 					write(lse,file_name,append=TRUE)

					write("TP of the learned C ",file_name,append=TRUE)
 					write(tp,file_name,append=TRUE)

					write("FP of the learned C ",file_name,append=TRUE)
 					write(fp,file_name,append=TRUE)

					# added by Irina: write a table
					density=A1_nnz/(ss*ss);
					res<-c(dd,density,lambda_cv,lse,tp,fp);
					write(res,fnm_table,6,append=TRUE);

				} # end of if(error_flag==0) 
			} # end for mats
		} #end for nets
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
		l_avg = lambda_avg/runs1;
		lse_avg = P_avg/runs1;
		tp_avg = TP_avg/runs1;
		fp_avg = FP_avg/runs1;
		res<-c(dd,density,l_avg,lse_avg,tp_avg,fp_avg);
		write(res,fnm_sum_table,6,append=TRUE);

	} # end of for sl

} # end of for vv



































































































	








