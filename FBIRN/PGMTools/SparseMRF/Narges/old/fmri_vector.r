library(MASS)
library(glasso)

file_name<-"fmri_results_sub1_k3_p200_vector.txt"


## Initialization:

set.seed(400)
ss<-201
dd<-704
mtest<-704

x <- matrix(scan("sub1_k3_200.txt", 0), ncol=ss, byrow=TRUE)
x_test<- matrix(scan("sub1_k3_200_test.txt", 0), ncol=ss, byrow=TRUE)

x<-0.1 * x
x_test <- 0.1* x_test
s<- var(x)

N<-dd

m<-1 # number of iterations
obj_opt<-matrix(rep(0,m))
obj_lambda_opt<-matrix(rep(0,m*ss),ncol=ss)
obj_cov_opt<-matrix(rep(0,ss*ss*m) ,ncol=ss*ss)
obj_iteration<-matrix(rep(0,m))
obj_trace<-matrix(rep(0,1000*m),ncol=1000)
obj_lse<-matrix(rep(0,m))

obj_det<-matrix(rep(0,1000))

kk<-1
done<-0
objective<-0
objective_last<-0
lambda<- rep(kk*5, ss)
last_lambda<-rep(kk*5,ss) 
new_lambda<-rep(kk*5,ss) 
first<-1
it<-0


eps<-0.001
Sinv<-solve(s + eps*diag(1,ss)) ##solve(s,tol = 10^(-22))
b<-rep(0,ss)
for(i in 1:ss){
b[i]<-ss/sum(abs(Sinv[i,]))
}



while(done==0){
	it<-it+1
	objective_last<-objective
	a<-glasso(s, rho=(2/N) *lambda  ,thr=1.0e-3,penalize.diagonal=TRUE)
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
		ll<-ll + ss*log(lambda[j]/2)- lambda[j] * norm1[j] - 0.5*(lambda[j]-b[j])^2
	}

	objective<-  - 0.5 * mm + dd * 0.5 * log ( det(a$wi)  )  + ll
	if( (objective > objective_last) || it==1){
		print("good") 
		print(it)
		for(i in 1:ss){
			bb<- norm1[i]  -b[i]
			new_lambda[i]<- 0.5*( -bb + sqrt(bb^2 + 4*ss) ) 		
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
			a<-glasso(s, rho=(2/N) *new_lambda ,thr=1.0e-3 ,penalize.diagonal=TRUE)
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
				ll<-ll + ss*log(new_lambda[j]/2)- new_lambda[j] * norm1[j] - 0.5*(new_lambda[j]-b[j])^2


			}
			objective<-  - 0.5 * mm + dd * 0.5 * log ( det(a$wi)  )  + ll
		}
	}
	last_lambda<-lambda
	lambda<-new_lambda
	obj_trace[kk,it]<-objective
	if( abs(objective - objective_last)  < 0.01){
		 done<-1
	}
	
	print("objective:")
	print(objective)
}


	obj_opt[kk]<-objective
	for(g in 1:ss){
		obj_lambda_opt[kk,g]<-lambda[g]
	}
	obj_cov_opt[kk,]<-a$wi
	obj_iteration[kk]<- it


















first<-0;

print("Sum lambda:")
print(sum(lambda))
print("det of inverse covariance:")
print(det(a$wi))

#print("norm")
#print(norm1)
print("objective:")
print(objective)



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



################   report the results:


p2<-ss*ss

 nnz<-0
 for(i in 1:p2){
##for(j in 1:ss)
if(abs(obj_cov_opt[kk,i])>0.0){ nnz<-nnz+1}
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

 write("b ",file_name,append=TRUE)
 write(b,file_name,append=TRUE)

 write("corr ",file_name,append=TRUE)
 write(corr,file_name,append=TRUE)




































