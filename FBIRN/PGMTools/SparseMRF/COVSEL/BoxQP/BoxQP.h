#ifdef MACVERSION 
#define mac 
#endif

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <mex.h>
#include <matrix.h>

#ifdef WIN32
#include <malloc.h>
#endif

#ifdef mac
#include "/System/Library/Frameworks/vecLib.framework/Headers/cblas.h"
#include "/System/Library/Frameworks/vecLib.framework/Headers/clapack.h"
#endif

// Calling LAPACK on the mac
#ifdef mac
#define dsyev dsyev_ 
#define dgesvd dgesvd_
#endif

// Calling LAPACK on linux 
#ifdef linuxp
#define dsyev dsyev_
#define dgesvd dgesvd_
// #define dgemv dgemv_ *** Some systems require underscore in calling the fortran BLAS ***
#endif

#ifdef WIN32
enum CBLAS_ORDER 	{CblasRowMajor=101, CblasColMajor=102};
enum CBLAS_TRANSPOSE 	{CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113};
enum CBLAS_UPLO		{CblasUpper=121, CblasLower=122};
enum CBLAS_DIAG		{CblasNonUnit=131, CblasUnit=132};
enum CBLAS_SIDE		{CblasLeft=141, CblasRight=142};
void cblas_dscal(mwSignedIndex N,double alpha, double *X,mwSignedIndex incX);
void cblas_dcopy(mwSignedIndex N,double *X,mwSignedIndex incX,double *Y,mwSignedIndex incY);
void cblas_dgemm(enum CBLAS_ORDER Order,enum CBLAS_TRANSPOSE transA, enum CBLAS_TRANSPOSE transB, mwSignedIndex M, mwSignedIndex N, mwSignedIndex K, double alpha, double *A, mwSignedIndex lda, double *B, mwSignedIndex ldb, double beta, double *C, mwSignedIndex ldc);
void cblas_dgemv(enum CBLAS_ORDER Order,enum CBLAS_TRANSPOSE transA, mwSignedIndex M, mwSignedIndex N, double alpha, double *A, mwSignedIndex lda, double *B, mwSignedIndex incB, double beta, double *C, mwSignedIndex incC);
void cblas_daxpy(mwSignedIndex N,double alpha,double *X,mwSignedIndex incX,double *Y,mwSignedIndex incY);
double cblas_ddot(mwSignedIndex n, double *x, mwSignedIndex incx, double *y, mwSignedIndex incy); 
#endif

// CBLAS declarations on linux 
#ifdef linuxp
enum CBLAS_ORDER 	{CblasRowMajor=101, CblasColMajor=102};
enum CBLAS_TRANSPOSE 	{CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113};
enum CBLAS_UPLO		{CblasUpper=121, CblasLower=122};
enum CBLAS_DIAG		{CblasNonUnit=131, CblasUnit=132};
enum CBLAS_SIDE		{CblasLeft=141, CblasRight=142};
void cblas_dscal(mwSignedIndex N,double alpha, double *X,mwSignedIndex incX);
void cblas_dcopy(mwSignedIndex N,double *X,mwSignedIndex incX,double *Y,mwSignedIndex incY);
void cblas_dgemm(enum CBLAS_ORDER Order,enum CBLAS_TRANSPOSE transA, enum CBLAS_TRANSPOSE transB, mwSignedIndex M, mwSignedIndex N, mwSignedIndex K, double alpha, double *A, mwSignedIndex lda, double *B, mwSignedIndex ldb, double beta, double *C, mwSignedIndex ldc);
void cblas_dgemv(enum CBLAS_ORDER Order,enum CBLAS_TRANSPOSE transA, mwSignedIndex M, mwSignedIndex N, double alpha, double *A, mwSignedIndex lda, double *B, mwSignedIndex incB, double beta, double *C, mwSignedIndex incC);
void cblas_daxpy(mwSignedIndex N,double alpha,double *X,mwSignedIndex incX,double *Y,mwSignedIndex incY);
#endif


// Local functions
double doubsum(double *xmat, mwSignedIndex n);

double doubdot(double *xvec, double *yvec, mwSignedIndex n);

double doubasum(double *xmat, mwSignedIndex n);

double doubnorm2(double *xmat, mwSignedIndex n);

mwSignedIndex idxmax(double *xmat, mwSignedIndex n);

void BoxQP(double *Amat, mwSignedIndex n, double *lvec, double *uvec, mwSignedIndex MaxIter, mwSignedIndex info, mwSignedIndex verbose, double *xvecout);

double dsignf(double x);

double dminif(double x, double y);

mwSignedIndex imaxf(mwSignedIndex x, mwSignedIndex y);

double dabsf(double x);


void dispmat(double *xmat, mwSignedIndex n, mwSignedIndex m);

double maxeig(double *xmat, mwSignedIndex n); 
