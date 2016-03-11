library(MASS)
library(glasso)


## Initialization: 

file_name<-"res_exp_vector_random.txt"
fnm_table<-"res_table_exp_vector_random.txt"
fnm_sum_table<-"res_sum_table_exp_vector_random.txt"
	

write("%n density  glcalls lse TP FP ",fnm_table,append=TRUE);
write("%n density  glcalls_avg lse  TP FP",fnm_sum_table,append=TRUE);

seed_r<-c(12,32,42,52,62,72,82,92,102,100,200,300,500,555,546,223,999,597,97,99,888,832,732,777,735,733,723,713,703,603)

vv<-1
for(vv in 1 : 4){
lambda_avg<-0
glcalls_avg<-0
P_avg<-0
FP_avg<-0
TP_avg<-0


runs<-1;
runs1<-0;

for(runs in 1:5){
error_flag<-0
set.seed(seed_r[runs])
ss<-100

nnz<-4*ss*runif(1)

A<-diag(ss)
for(i in 1:nnz){
 A[ss*runif(1),ss*runif(1)]=sign(runif(1)-.5);
}

A1<-t(A)%*%A

Sig=solve(A1)

N_vector<- c(30,50,500,1000)


mtest<-100;

x_test<-mvrnorm(mtest, mu= c(rep(0,ss) ), Sig )






dd<-N_vector[vv];
data<-mvrnorm(dd,rep(0,ss),Sig)
x<-matrix(data,ncol=ss)
s<- var(x)

eps<-0.001
Sinv<-solve(s + eps*diag(1,ss)) ##solve(s,tol = 10^(-22))
b<-rep(0,ss)
for(i in 1:ss){
b[i]<-sum(abs(Sinv[i,]))
}
b<-b/(ss)
#b<-b*(ss/dd)

A1_nnz<-0
for(i in 1:ss)
for(j in 1:ss)
if(A1[i,j]!=0) A1_nnz<-A1_nnz+1
print(A1_nnz)

N<-dd
m<-1# number of iterations
obj_opt<-matrix(rep(0,m))
obj_lambda_opt<-matrix(rep(0,m*ss),ncol=ss)
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
	lambda<- rep(kk*5, ss)
	last_lambda<-rep(kk*5,ss) 
	new_lambda<-rep(kk*5,ss) 
	it<-0


### SEARCH:	
while(done==0){
	it<-it+1
	objective_last<-objective
	a<-glasso(s, rho=(2/N) *lambda,thr=1.0e-3  ,penalize.diagonal=TRUE)
	glcalls = glcalls +1; # added by Irina
	norm1<-matrix(rep(0,ss))
	for (i in 1:ss){
		for(j in 1:ss){
			norm1[i]<-norm1[i]+ abs(a$wi[i,j])
		}
	}
	mm<-0
	for(i in 1:dd){
		mm<-mm+( x[i,] %*% a$wi %*% x[i,] )
	}

	ll<-0
	for(j in 1:ss){
		ll<-ll + ss*log(lambda[j]/2)- lambda[j] * norm1[j] -b[j]*lambda[j]
	}

	objective<-  - 0.5 * mm + dd * 0.5 * log ( det(a$wi)  )  + ll
	if(det(a$wi)<10^-320)
	{error_flag<-1}
	else if( (objective > objective_last) || it==1){
		print("good") 
		print(it)
		for(i in 1:ss){
			new_lambda[i]<-ss/(norm1[i] + b[i])   ##( (ss*ss/last_lambda) - norm11 )		
		}
	}
	else{
		print(" **************************backtrachking **************************************")
		alpha<-1  
		while(objective < objective_last && alpha>0.001){	      
			##print(alpha)
			for(i in 1:ss){
				new_lambda[i]<-last_lambda[i] + alpha * (new_lambda[i] - last_lambda[i])    ##( (ss*ss/last_lambda) - norm11 )		
			}
			alpha<-alpha/2;
			a<-glasso(s, rho=(2/N) *new_lambda,thr=1.0e-3  ,penalize.diagonal=TRUE)
			glcalls = glcalls +1; # added by Irina
			norm1<-matrix(rep(0,ss))
			for (i in 1:ss){
				for(j in 1:ss){
					norm1[i]<-norm1[i]+ abs(a$wi[i,j])
				}	
			}
			mm<-0
			for(i in 1:dd){
				mm<-mm+( x[i,] %*% a$wi %*% x[i,] )
			}
			ll<-0
			for(j in 1:ss){
				ll<-ll + ss*log(new_lambda[j]/2)- new_lambda[j] * norm1[j] -b[j]*new_lambda[j]

			}
			objective<-  - 0.5 * mm + dd * 0.5 * log ( det(a$wi)  )  + ll
		}
	}
	last_lambda<-lambda
	lambda<-new_lambda
	obj_trace[kk,it]<-objective
	if( error_flag==1 )
	{
		done<-1
	}
	else if (objective - objective_last  < 0.01 ){
		 done<-1
	}
	
	print("objective:")
	print(objective)
}


	for (g in 1:ss){
		for(gg in 1:ss){
			if(a$wi[g,gg]!=0) obj_nnz[kk]<-obj_nnz[kk]+1
		}
	}
	obj_opt[kk]<-objective
	for(g in 1:ss){
		obj_lambda_opt[kk,g]<-lambda[g]
	}
	obj_cov_opt[kk,]<-a$wi
	obj_iteration[kk]<- it
glcalls_avg = glcalls + glcalls_avg; # added by Irina

}



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
		res<-c(dd,density,glcalls,lse,tp,fp);
		write(res,fnm_table,6,append=TRUE);


}

}
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
		
		l_avg = lambda_avg/runs1;
		lse_avg = P_avg/runs1;
		tp_avg = TP_avg/runs1;
		fp_avg = FP_avg/runs1;
		res<-c(dd,density,glcalls_avg,lse_avg,tp_avg,fp_avg);
		write(res,fnm_sum_table,6,append=TRUE);

}



































































































	








