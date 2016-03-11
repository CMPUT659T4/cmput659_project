
file_name<-"RandomInverseCovariances.txt"
seed_r<-c(12,32,42,52,62,72,82,92,102,100,200,300,500,555,546,223,999,597,97,99,888,832,732,777,735,733,723,713,703,603)

runs<-0;
ss<-100

for(runs in 1:30){
	set.seed(seed_r[runs])
	nnz<-4*ss*runif(1)
	A<-diag(ss)
	for(i in 1:nnz){
	 A[ss*runif(1),ss*runif(1)]=sign(runif(1)-.5);
	}
	A1<-t(A)%*%A
	Sig=solve(A1)
	write("random matrix",file_name,append=TRUE)
	write(runs,file_name,append=TRUE)
	write(A1,file_name,ncolumns=100,append=TRUE)
}
