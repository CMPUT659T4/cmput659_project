#include "BoxQP.h"

// BLAS wrapper for win32 version
#ifdef WIN32
void cblas_dscal( mwSignedIndex N, double alpha, double *X, mwSignedIndex incX)
{
	dscal(&N,&alpha,X,&incX);
}

void cblas_dcopy(mwSignedIndex N,double *X,mwSignedIndex incX,double *Y,mwSignedIndex incY)
{
	dcopy(&N,X,&incX,Y,&incY);
}

void cblas_dgemm(enum CBLAS_ORDER Order,enum CBLAS_TRANSPOSE transA, enum CBLAS_TRANSPOSE transB,
                 mwSignedIndex M, mwSignedIndex N, mwSignedIndex K, double alpha, double *A, mwSignedIndex lda,
                 double *B, mwSignedIndex ldb, double beta, double *C, mwSignedIndex ldc)
{
	char ta[1],tb[1];
	if (transA==111)
	{
		*ta='N';
	}
	else
	{
		*ta='T';
	};
	if (transB==111)
	{
		*tb='N';
	}
	else
	{
		*tb='T';
	};
	dgemm(ta,tb,&M,&N,&K,&alpha,A,&lda,B,&ldb,&beta,C,&ldc);
}

void cblas_dgemv(enum CBLAS_ORDER Order,enum CBLAS_TRANSPOSE transA,
                 mwSignedIndex M, mwSignedIndex N, double alpha, double *A, mwSignedIndex lda,
                 double *B, mwSignedIndex incB, double beta, double *C, mwSignedIndex incC)
{
	char ta[1];
	if (transA==111)
	{
		*ta='N';
	}
	else
	{
		*ta='T';
	};
	dgemv(ta,&M,&N,&alpha,A,&lda,B,&incB,&beta,C,&incC);
}

void cblas_daxpy(mwSignedIndex N,double alpha,double *X,mwSignedIndex incX,double *Y,mwSignedIndex incY)
{
	daxpy(&N,&alpha,X,&incX,Y,&incY);
}

#endif WIN32

// BLAS wrapper for the linux version
#ifdef linuxp
void cblas_dscal( mwSignedIndex N, double alpha, double *X, mwSignedIndex incX)
{
	dscal(&N,&alpha,X,&incX);
}

void cblas_dcopy(mwSignedIndex N,double *X,mwSignedIndex incX,double *Y,mwSignedIndex incY)
{
	dcopy(&N,X,&incX,Y,&incY);
}

void cblas_dgemm(enum CBLAS_ORDER Order,enum CBLAS_TRANSPOSE transA, enum CBLAS_TRANSPOSE transB,
                 mwSignedIndex M, mwSignedIndex N, mwSignedIndex K, double alpha, double *A, mwSignedIndex lda,
                 double *B, mwSignedIndex ldb, double beta, double *C, mwSignedIndex ldc)
{
	char ta[1],tb[1];
	if (transA==111)
	{
		*ta='N';
	}
	else
	{
		*ta='T';
	};
	if (transB==111)
	{
		*tb='N';
	}
	else
	{
		*tb='T';
	};
	dgemm(ta,tb,&M,&N,&K,&alpha,A,&lda,B,&ldb,&beta,C,&ldc);
}

void cblas_dgemv(enum CBLAS_ORDER Order,enum CBLAS_TRANSPOSE transA,
                 mwSignedIndex M, mwSignedIndex N, double alpha, double *A, mwSignedIndex lda,
                 double *B, mwSignedIndex incB, double beta, double *C, mwSignedIndex incC)
{
	char ta[1];
	if (transA==111)
	{
		*ta='N';
	}
	else
	{
		*ta='T';
	};
	dgemv(ta,&M,&N,&alpha,A,&lda,B,&incB,&beta,C,&incC);
}

void cblas_daxpy(mwSignedIndex N,double alpha,double *X,mwSignedIndex incX,double *Y,mwSignedIndex incY)
{
	daxpy(&N,&alpha,X,&incX,Y,&incY);
}
#endif


// Some useful functions ...
double doubsum(double *xmat, mwSignedIndex n)
{
	mwSignedIndex i;
	double res=0.0;
	for (i=0;i<n;i++){res+=xmat[i];};
	return res;
}


double doubdot(double *xvec, double *yvec, mwSignedIndex n)
{
	mwSignedIndex i;
	double res=0.0;
	for (i=0;i<n;i++){res+=xvec[i]*yvec[i];};
	return res;
}

mwSignedIndex idxmax(double *xmat, mwSignedIndex n)
{
	mwSignedIndex i;
	mwSignedIndex res=0;
	for (i=0;i<n;i++)
	{
		if (xmat[i]>xmat[res]) {res=i;}
	}
	return res;
}


double doubasum(double *xmat, mwSignedIndex n)
{
	mwSignedIndex i;
	double res=0.0;
	for (i=0;i<n;i++){res+=dabsf(xmat[i]);};
	return res;
}

double doubnorm2(double *xmat, mwSignedIndex n)
{
	mwSignedIndex i;
	double res=0.0;
	for (i=0;i<n;i++){res+=xmat[i]*xmat[i];};
	return sqrt(res);
}

double dsignf(double x)
{
	if (x>=0)
		return 1.0;
	else
		return -1.0;
}

double dminif(double x, double y)
{
	if (x>=y)
		return y;
	else
		return x;
}

mwSignedIndex imaxf(mwSignedIndex x, mwSignedIndex y)
{
	if (x>=y)
		return x;
	else
		return y;
}

double dabsf(double x)
{
	if (x>=0)
		return x;
	else
		return -x;
}

void dispmat(double *xmat, mwSignedIndex n, mwSignedIndex m)
{
	mwSignedIndex i,j;

	for (i=0; i<n; i++)
	{
		for (j=0;j<m;j++)
		{
			prmwSignedIndexf("%+.4f ",xmat[j*n+i]);
		}
		prmwSignedIndexf("\n");
	}
	prmwSignedIndexf("\n");
}

double maxeig(double *xmat, mwSignedIndex n) 
{
	// xmat is symmetric n x n matrix
	mwSignedIndex incx=1,indmax,maxloop=10000,k=0;
	double alpha, beta, dmax, dmax_temp,dmax_tol;
	double *bufveca=(double *) calloc(n,sizeof(double));
	double *bufvecb=(double *) calloc(n,sizeof(double));

	dmax_tol=.001;
	// do power series to get approximation to max eigenvalue of A+X
	alpha=0.0;cblas_dscal(n,alpha,bufveca,incx);bufveca[0]=1.0; // x_0 = [1,0,0,...] do something better later
	beta=0.0;dmax=1.0;dmax_temp=0.0;
	while ((dabsf(dmax-dmax_temp)>dmax_tol)&&(k<=maxloop)){
		dmax_temp=dmax;
		alpha=1.0;
		cblas_dgemv(CblasColMajor,CblasNoTrans,n,n,alpha,xmat,n,bufveca,incx,beta,bufvecb,incx);
		indmax=idxmax(bufvecb,n);dmax=bufvecb[indmax];
		alpha=1.0/dmax;cblas_dscal(n,alpha,bufvecb,incx);
		cblas_dcopy(n,bufvecb,incx,bufveca,incx);
		k++;
	}
	alpha=1.0;	
	// compute Rayleigh Quotient to approximate max eigenvalue of A+X
	cblas_dgemv(CblasColMajor,CblasNoTrans,n,n,alpha,xmat,n,bufvecb,incx,beta,bufveca,incx);
	dmax=doubdot(bufvecb,bufveca,n);alpha=doubnorm2(bufvecb,n);dmax=dmax/alpha/alpha;
	free(bufveca);free(bufvecb);
	return dmax;
}
