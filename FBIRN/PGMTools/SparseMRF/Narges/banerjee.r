library(MASS)
library(glasso)


## Initialization: 

file_name<-"results_banerjee_random.txt"
seed_r<-c(12,32,42,52,62,72,82,92,102,100,200,300,500,555,546,223,999,597,97,99,888,832,732,777,735,733,723,713,703,603)

vv<-1
for(vv in 1 : 4){  ## over 4 values of N
lambda_avg<-0
P_avg<-0
FP_avg<-0
TP_avg<-0

runs1<-0;
runs<-1;

for(runs in 1:30){
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



var_x<-matrix(rep(0,ss))
for(ii in 1:ss){var_x[ii]<-var(x[,ii])}
sm<-0
sm2<-0
for(ii in 1:ss){
	jj<-ii+1
	while(jj<=ss){			sm2<-var_x[ii] * var_x[jj]
		if (sm2 > sm){
			 sm<-sm2
		}	
		jj<-jj+1
	}
}

h<- abs(qt(0.001/(ss*ss*2),df=dd-2 ) ) 
gg<- (dd/2) * sqrt(sm) * (h)/sqrt(dd-2 + h^2)
print(gg)
ag<-glasso(s, rho= gg *( 2/dd ),penalize.diagonal=TRUE)


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


lambda_avg<-lambda_avg + gg

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

 write("banerjee lambda",file_name,append=TRUE)
 write(gg,file_name,append=TRUE)

 write("Prediction Mean Squared Error of the learned model ",file_name,append=TRUE)
 write(lse,file_name,append=TRUE)

write("TP of the learned C ",file_name,append=TRUE)
 write(tp,file_name,append=TRUE)

write("FP of the learned C ",file_name,append=TRUE)
 write(fp,file_name,append=TRUE)

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



}



































































































	








