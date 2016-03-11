#experiments with the fMRI data with Gaussian prior
# you need to first generate the input files -  this is done by process.m 

library(MASS)
library(glasso)

file_name<-"fmri_results_sub14_k3_p200_prior.txt"


## Initialization:

set.seed(400)
ss<-201
dd<-704
mtest<-704

x <- matrix(scan("sub1_k15_200.txt", 0), ncol=ss, byrow=TRUE)
x_test<- matrix(scan("sub1_k15_200_test.txt", 0), ncol=ss, byrow=TRUE)

x<-0.1 * x
x_test <- 0.1* x_test
s<- var(x)

N<-dd

m<-1 # number of iterations
obj_opt<-matrix(rep(0,m))
obj_lambda_opt<-matrix(rep(0,m))
obj_cov_opt<-matrix(rep(0,ss*ss*m) ,ncol=ss*ss)
obj_iteration<-matrix(rep(0,m))
obj_trace<-matrix(rep(0,1000*m),ncol=1000)
obj_lse<-matrix(rep(0,m))

obj_det<-matrix(rep(0,1000))

kk<-1
done<-0
objective<-0
objective_last<-0
lambda<-50  ## just an initial point
lambda_last<-lambda

first<-1
it<-1


eps<-0.001
Sinv<-solve(s + eps*diag(1,ss)) ##solve(s,tol = 10^(-22))
b<-sum(abs(Sinv))
b<-b/(ss*ss)
b<-1/b
#b<-b*2*(ss/dd)

print(b)




while(done==0){

print("******************************ITERATION: ******************** ")
print(it)
print(lambda)

objective_last<-objective

it<-it+1



lambdav<-rep(lambda,ss)
lambdav[1]<-0.01
a<-glasso(s,thr=1.0e-3, rho=lambdav * 2/N, trace=TRUE, penalize.diagonal=TRUE)


norm1<-0
for (i in 2:ss){
	for(j in 2:ss){
		if(i!=1) norm1<-norm1+ abs(a$wi[i,j])
	}
}
mm<-0

for(i in 1:dd){
mm<-mm+( x[i,] %*% a$wi %*% x[i,] )
}

objective<-  - 0.5 * mm + N * 0.5 * log ( det(a$wi)  )  +  (ss-1)*(ss-1)* log (lambda/2) - lambda * norm1 -0.5*(lambda-b)^2

bb<- norm1  -b
new_lambda<- 0.5*( -bb + sqrt(bb^2 + 4*(ss-1)^2) ) ##(ss-1) *( ss-1) / (norm1 + b) 


if( (objective > objective_last) || first==1){
	
	print("good")
}

else if(objective <= objective_last){
print("not good")
}
last_lambda<-lambda
lambda<-new_lambda

if(it<1000){ obj_trace[kk,it]<-objective}

if( ( (objective - objective_last)  < 1) && (first==0)){
 done<-1
}

first<-0;

print("lambda:")
print(lambda)
print("det of inverse covariance:")
print(det(a$wi))
print("norm")
print(norm1)
print("objective:")
print(objective)
}

obj_opt[kk]<-objective

obj_lambda_opt[kk]<-lambda

obj_cov_opt[kk,]<-a$wi

obj_iteration[kk]<-it

obj_det[it]<-det(a$wi)

###########################     Prediction Error  ##########################################
#print("LSE")
n<-ss
corr<-0
sigma12<-matrix(rep(0,n-1), nrow=1)
sigma22<- matrix(rep(0,(n-1)*(n-1) ), nrow=n-1, ncol= n-1)
sigma12<-a$w[1,2:n];
sigma22<-a$w[2:n,2:n];
sigma22inv<-solve(sigma22)
sigma12<-as.matrix(sigma12, nrow=n-1, ncol=1)
sigma22inv<-as.matrix(sigma22inv, nrow=n-1, ncol=n-1)
error=0
y_1<-matrix(rep(0,(n-1)*(n-1)), nrow=n-1, ncol=n-1)
y_test<-matrix(rep(0,mtest))
for (i in 1:mtest){
xval<-as.matrix(t(x_test[i,2:n]), nrow=n-1, ncol=1)
y1<-t(sigma12) %*% sigma22inv
y_test[i]<- y1 %*% t(xval)
#print(y_test)
error = error + ( y_test[i] - x_test [i,1] )^2

##corr<-corr+ y_test[i]*x_test[i,1]

}
lse=error/mtest
obj_lse[kk]<-lse

corr<- cov(y_test,x_test[,1]) / sqrt(var(x_test[,1]) * var(y_test) )
print(corr)


########################## Calculate Banerjees Lambda:
var_x<-matrix(rep(0,ss))
for(ii in 1:ss){
	var_x[ii]<-var(x[,ii])
}
sm<-0
sm2<-0
for(ii in 1:ss){
jj<-ii+1
	while(jj<=ss){
	sm2<-var_x[ii] * var_x[jj]
	if (sm2 > sm){ sm<-sm2}
	jj<-jj+1
	}
}

h<- abs(qt(0.01/(ss*ss*2),df=dd-2 ) )
gg<- (dd/2) * sqrt(sm) * (h)/sqrt(dd-2 + h^2)


###########################     Prediction Error for Banerjee's  ##########################################


a_g<-glasso(s, rho=gg * 2/N, penalize.diagonal=TRUE)
corrg<-0
n<-ss
sigma12<-matrix(rep(0,n-1), nrow=1)
sigma22<- matrix(rep(0,(n-1)*(n-1) ), nrow=n-1, ncol= n-1)
sigma12<-a_g$w[1,2:n];
sigma22<-a_g$w[2:n,2:n];
sigma22inv<-solve(sigma22)
sigma12<-as.matrix(sigma12, nrow=n-1, ncol=1)
sigma22inv<-as.matrix(sigma22inv, nrow=n-1, ncol=n-1)
error=0
y_1<-matrix(rep(0,(n-1)*(n-1)), nrow=n-1, ncol=n-1)
y_test<-matrix(rep(0,mtest))
for (i in 1:mtest){
xval<-as.matrix(t(x_test[i,2:n]), nrow=n-1, ncol=1)
y1<-t(sigma12) %*% sigma22inv
y_test[i]<- y1 %*% t(xval)
#print(y_test)
#corrg<-corrg+ (y_test*x_test[i,1])

error = error + ( y_test[i] - x_test [i,1] )^2
}
lse_g=error/mtest
corrg<- cov(y_test,x_test[,1]) / sqrt(var(x_test[,1]) * var(y_test) )

################   report the results:


p2<-ss*ss

 nnz<-0
 for(i in 1:p2){
##for(j in 1:ss)
if(abs(obj_cov_opt[kk,i])>0.0){ nnz<-nnz+1}
}

 nnzg<-0
 for(i in 1:ss){
for(j in 1:ss)
if(abs(a_g$wi[i,j])>0.0){ nnzg<-nnzg+1}
}


 write("----------------------------",file_name,append=TRUE)

write("N",file_name,append=TRUE)
write(dd,file_name,append=TRUE)
write("P",file_name,append=TRUE)
write(ss,file_name,append=TRUE)



 write("nnz of learned C ",file_name,append=TRUE)
 write(nnz,file_name,append=TRUE)

 write("learned lambda",file_name,append=TRUE)
 write(obj_lambda_opt[kk],file_name,append=TRUE)

 write("Prediction Mean Squared Error of the learned model ",file_name,append=TRUE)
 write(obj_lse[kk],file_name,append=TRUE)

 write("Iterations of the line search ",file_name,append=TRUE)
 write(obj_iteration[kk],file_name,append=TRUE)


 write("Lambda learned by Banerjee ",file_name,append=TRUE)
 write(gg,file_name,append=TRUE)


 write("Prediction Mean Squared Error of the Banerjee Lambda ",file_name,append=TRUE)
 write(lse_g,file_name,append=TRUE)

 write("nnz of Banerjee  C ",file_name,append=TRUE)
 write(nnzg,file_name,append=TRUE)

 write("b ",file_name,append=TRUE)
 write(b,file_name,append=TRUE)

 write("corr ",file_name,append=TRUE)
 write(corr,file_name,append=TRUE)


 write("corrg",file_name,append=TRUE)
 write(corrg,file_name,append=TRUE)






































