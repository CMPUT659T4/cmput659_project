# experiments on synthetic data with exponential prior
 
 
#setwd("C:\\Tools\\MyMatlabTools\\SparseMRF\\Narges")

library(MASS)
library(glasso)
##library(corpcor)
## Initialization: 



file_name<-"report_p_wp_d50_n_Final.txt"
b<-0
set.seed(32)
ss<-100

nnz<-900


A<-diag(ss)
for(i in 1:nnz){
 A[ss*runif(1),ss*runif(1)]=sign(runif(1)-.5);
}
A1<-t(A)%*%A

A1_nnz<-0
for(ii in 1:ss)
for(jj in 1:ss)
if(A1[ii,jj]!=0) A1_nnz<-A1_nnz+1
print(A1_nnz)

Sig=solve(A1)
N_vector<- c(30,50,100,200,500,1000)
vv<-1

mtest<-100;

x_test<-mvrnorm(mtest, mu= c(rep(0,ss) ), Sig )
xv_test<-x_test
for(z in 1:ss){
		x_test[,z]<-x_test[,z]-mean(x_test[,z])
}
for(z in 1:ss){
		x_test[,z]<-x_test[,z]/sqrt(var(x_test[,z]))
}



vv<-1
for(vv in 1 : 6){

dd<-N_vector[vv];
data<-mvrnorm(dd,rep(0,ss),Sig)
x<-matrix(data,ncol=ss)
xv<-x

for(z in 1:ss){
		x[,z]<-x[,z]-mean(x[,z])
}
for(z in 1:ss){
		x[,z]<-x[,z]/sqrt(var(x[,z]))
}


s<- var(xv)  #########################################3 NOT NORMALIZED

eps<-0.001
Sinv<-solve(s + eps*diag(1,ss)) ##solve(s,tol = 10^(-22))
b<-sum(abs(Sinv))
b<-b/(ss*ss)
b<-b*2*(ss/dd)


m<-300# number of iterations
obj<-matrix(rep(0,m))
obj_lambda<-matrix(rep(0,m))
#obj_cov<-matrix(rep(0,ss*ss*m) ,ncol=ss*ss)
obj_norm<-matrix(rep(0,m))
obj_nnz<-matrix(rep(0,m))
obj_lse<-matrix(rep(0,m))
obj_fp<-matrix(rep(0,m))
#obj_fp_avg<-matrix(rep(0,m))
obj_fn<-matrix(rep(0,m))
#obj_fn_avg<-matrix(rep(0,m))
#obj_nnz_avg<-matrix(rep(0,m))
obj_det<-matrix(rep(0,m))
accuracy<-matrix(rep(ss*ss,m),nrow=m)
#accuracy_avg<-matrix(rep(ss*ss,m),nrow=m)
obj_cv<-matrix(rep(0,m))
obj_lik<-matrix(rep(0,m))
#diff<-matrix(rep(0,m))
#p1<-matrix(rep(0,m))
#p2<-matrix(rep(0,m))


LV<-matrix(rep(0,m))
LV[1:5]<-c(0.1,0.2,0.4,0.6,0.8)
i=6;
k<-1
for(i in 6:100){
	LV[i]<-k
	k<-k+1	
}
i<-101
k<-51
for(i in 101:m){
	LV[i]<-2*k
	k<-k+1
}

########################################### LOOP  #######################333

for(kk in 1:m){

lambda<-LV[kk]
a<-glasso(s, rho= lambda *( 2/dd ),penalize.diagonal=TRUE)

print(kk)

norm1<-0
for (i in 1:ss){
	for(j in 1:ss){
		norm1<-norm1+ abs(a$wi[i,j])	
		if(a$wi[i,j]!=0){ obj_nnz[kk]<-obj_nnz[kk]+1}
	}
}

obj_norm[kk]<-norm1

mm<-0
for(i in 1:dd){
	mm<-mm+( x[i,] %*% a$wi %*% x[i,] )
}


objective<-  - 0.5 * mm + dd * 0.5 * log ( det(a$wi)  )+ (ss*ss)* log (lambda/2)  - lambda * norm1 -b*lambda

obj_lik[kk]<- mm
obj_det[kk]<-dd * 0.5 * log ( det(a$wi)  )
obj_norm[kk]<- norm1

#print("lambda:")
#print(lambda)
#print("inverse covariance:")
#print(a$wi)
#print("objective:")
#print(objective)
obj[kk]<-objective
obj_lambda[kk]<-lambda
#obj_cov[kk,]<-a$wi
 

###########################     Prediction Error  ##########################################
#print("LSE")
n<-ss
sigma12<-matrix(rep(0,n-1), nrow=1)
sigma22<- matrix(rep(0,(n-1)*(n-1) ), nrow=n-1, ncol= n-1)
sigma12<-a$w[1,2:n];
sigma22<-a$w[2:n,2:n];
sigma22inv<-solve(sigma22)
sigma12<-as.matrix(sigma12, nrow=n-1, ncol=1)
sigma22inv<-as.matrix(sigma22inv, nrow=n-1, ncol=n-1)
error=0
y_1<-matrix(rep(0,(n-1)*(n-1)), nrow=n-1, ncol=n-1)
y_test<-0
for (i in 1:mtest){
xval<-as.matrix(t(xv_test[i,2:n]), nrow=n-1, ncol=1)
y1<-t(sigma12) %*% sigma22inv
y_test<- y1 %*% t(xval)
#print(y_test)
error = error + ( y_test - xv_test [i,1] )^2
}
lse=error/mtest
obj_lse[kk]<-lse

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
			train_data[t,]<-xv[i,];
			t<-t+1;
		}
	}
	sc<- var(train_data)
	aa<-glasso(sc, rho= lambda *( 2/dd ),penalize.diagonal=TRUE)
	mm<-0
	for(i in 1:dd){
		if(test_indicator[i]!=0){
			mm<-mm+( xv[i,] %*% aa$wi %*% xv[i,] )
		}
	}
	norm1c<-0
	for (ii in 1:ss){
		for(jj in 1:ss){
			norm1c<-norm1c+ abs(aa$wi[ii,jj])
		}
	}
	avg_likelihood<- avg_likelihood - 0.5 * mm + train * 0.5 * log ( det(aa$wi)  ) - lambda * norm1c
	mc<-mc+k_d
}
avg_likelihood<- avg_likelihood/k
obj_cv[kk]<-avg_likelihood

###########################################################################################
#avg_nnz<-0
#nnz_el<-0
for(h in 1:ss){
	for(hh in 1:ss){
#		if (a$wi[h,hh]!=0 && h!=hh){
#			avg_nnz<-avg_nnz + abs (a$wi[h,hh])
			#nnz_el<-nnz_el+1;
		#}
		if(  ( (A1[h,hh]!=0) && (a$wi[h,hh]==0) ) ||  (A1[h,hh]==0 && a$wi[h,hh]!=0)  ){
			accuracy[kk]<-accuracy[kk]-1
		}
		if(  (A1[h,hh]!=0) && (a$wi[h,hh]==0) ) {
			obj_fn[kk]<-obj_fn[kk]+1
		}
		if(  (A1[h,hh]==0) && (a$wi[h,hh]!=0) ) {
			obj_fp[kk]<-obj_fp[kk]+1
		}
	}
}

}






tp<-obj_nnz-obj_fp

##plot(tp)
## ROC CURVE:
##plot(obj_fp,tp)


Lambda_obj<-0
max_obj<-max(obj)
for(i in 1:m){
	if(obj[i]==max_obj){

		obj_lambda1<-obj_lambda[i]
		print(obj_lambda)
		obj_tp<-tp[i]
		obj_fp1<-obj_fp[i]
		obj_lse1<-obj_lse[i]
		print(tp[i])
		print(obj_fp[i])
		print(obj_lse[i])
	}
}

Lambda_cv<-0
max_cv<-max(obj_cv)
for(i in 1:m){
	if(obj_cv[i]==max_cv){
		cv_lambda<-obj_lambda[i]
		print(cv_lambda)
		cv_tp<-tp[i]
		cv_fp<-obj_fp[i]
		cv_lse<-obj_lse[i]
		print(tp[i])
		print(obj_fp[i])
		print(obj_lse[i])
	}
}


Lambda_p<-0
max_p<-min(obj_lse)
for(i in 1:m){
	if(obj_lse[i]==max_p){
		p_lambda<-obj_lambda[i]
		print(p_lambda)
		p_tp<-tp[i]
		p_fp<-obj_fp[i]
		p_lse<-obj_lse[i]
		print(tp[i])
		print(obj_fp[i])
		print(obj_lse[i])
	}
}



######### Calculate Banerjees Lambda:

var_x<-matrix(rep(0,ss))
for(ii in 1:ss){
var_x[ii]<-var(xv[,ii])
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

h<- abs(qt(0.001/(ss*ss*2),df=dd-2 ))
h
gg<- (dd/2) * sqrt(sm) * (h)/sqrt(dd-2 + h^2)

print(gg)
#sv<-var(xv)
ag<-glasso(s, rho= gg *( 2/dd ),penalize.diagonal=TRUE)
###########################     Prediction Error  ##########################################
#print("LSE")
n<-ss
sigma12<-matrix(rep(0,n-1), nrow=1)
sigma22<- matrix(rep(0,(n-1)*(n-1) ), nrow=n-1, ncol= n-1)
sigma12<-ag$w[1,2:n];
sigma22<-ag$w[2:n,2:n];
sigma22inv<-solve(sigma22)
sigma12<-as.matrix(sigma12, nrow=n-1, ncol=1)
sigma22inv<-as.matrix(sigma22inv, nrow=n-1, ncol=n-1)
error=0
y_1<-matrix(rep(0,(n-1)*(n-1)), nrow=n-1, ncol=n-1)
y_test<-0
for (i in 1:mtest){
xval<-as.matrix(t(xv_test[i,2:n]), nrow=n-1, ncol=1)
y1<-t(sigma12) %*% sigma22inv
y_test<- y1 %*% t(xval)
#print(y_test)
error = error + ( y_test - xv_test [i,1] )^2
}
lse_ag=error/mtest

tp_ag<-0
fp_ag<-0
nnz_ag<-0
for(h in 1:ss){
	for(hh in 1:ss){
		if(ag$wi[h,hh]!=0) {nnz_ag<-nnz_ag+1}
		if(  (A1[h,hh]==0) && (ag$wi[h,hh]!=0) ) {
			fp_ag<-fp_ag+1
		}
	}
}


tp_ag<-nnz_ag - fp_ag
print(tp_ag)
print(fp_ag)
print(lse_ag)




### write the results in the file:


 write("-------***---------------------------------------------------",file_name,append=TRUE)
 write("density=52% ,P=100  scaled: 2(p/N)  and no normalization of the data Final Simulation! ",file_name,append=TRUE)

 write("N= ",file_name,append=TRUE)

 write(dd,file_name,append=TRUE)

 write("b= ",file_name,append=TRUE)

 write(b,file_name,append=TRUE)

 write("learned lambda",file_name,append=TRUE)
 write(obj_lambda1,file_name,append=TRUE)

 write("Prediction Mean Squared Error of the learned model ",file_name,append=TRUE)
 write(obj_lse1,file_name,append=TRUE)

write("TP of the learned C ",file_name,append=TRUE)
 write(obj_tp,file_name,append=TRUE)

write("FP of the learned C ",file_name,append=TRUE)
 write(obj_fp1,file_name,append=TRUE)

 write("Lambda learned by Banerjee ",file_name,append=TRUE)
 write(gg,file_name,append=TRUE)

 write("TP  Banerjee  C ",file_name,append=TRUE)
 write(tp_ag,file_name,append=TRUE)

 write("FP  Banerjee  C ",file_name,append=TRUE)
 write(fp_ag,file_name,append=TRUE)

 write(" Prediction Error  Banerjee  C ",file_name,append=TRUE)
 write(lse_ag,file_name,append=TRUE)



 write("CV learned lambda",file_name,append=TRUE)
 write(cv_lambda,file_name,append=TRUE)

 write("Prediction Mean Squared Error of the learned model ",file_name,append=TRUE)
 write(cv_lse,file_name,append=TRUE)

write("TP of the learned C ",file_name,append=TRUE)
 write(cv_tp,file_name,append=TRUE)

write("FP of the learned C ",file_name,append=TRUE)
 write(cv_fp,file_name,append=TRUE)





 write("prediction lambda",file_name,append=TRUE)
 write(p_lambda,file_name,append=TRUE)

 write("PE ",file_name,append=TRUE)
 write(p_lse,file_name,append=TRUE)

write("TP of the learned C ",file_name,append=TRUE)
 write(p_tp,file_name,append=TRUE)

write("FP of the learned C ",file_name,append=TRUE)
 write(p_fp,file_name,append=TRUE)






}
