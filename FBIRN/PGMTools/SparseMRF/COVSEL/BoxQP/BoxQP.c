/* Main function
Solves the follwoinng problem

  min x'A'x  s.t. l <= x <= u

This code implements Nesterov's smooth minimization algorithm. 
See: Y. Nesterov "Smooth Minimization of NonSmooth Functions", CORE DP 2003/12. 

Last Modified: A. d'Aspremont, March 2006.
http://www.carva.org/alexandre.daspremont
*/

#include "BoxQP.h"

//BoxQP(Amat,n,lvec,uvec,(mwSignedIndex)(MaxIter),(mwSignedIndex)(info),(mwSignedIndex)(verbose),xvecout);
void BoxQP(double *Amat, mwSignedIndex n, double *lvec, double *uvec, mwSignedIndex MaxIter, mwSignedIndex info, mwSignedIndex verbose, double *xvec)
{
	// Hard parameters
	mwSignedIndex Nperiod=imaxf(1,info);
	// Working variables
	double sig,L;
	double fval;
	mwSignedIndex incx=1,precision_flag=0,iteration_flag=0;
	mwSignedIndex k=0,i,n2=n*n;
	double alpha,beta,lambda,sumalpha,kalpha;
	double cputime,last_time=(double)clock();double start_time=(double)clock();mwSignedIndex left_h=0,left_m=0,left_s=0;
	double *fvec=(double *) calloc(n,sizeof(double));
	double *gvec=(double *) calloc(n,sizeof(double));
	double *bufveca=(double *) calloc(n,sizeof(double));
	double *bufvecb=(double *) calloc(n,sizeof(double));
	
	// Start...
	if (verbose>=1)
	{
		mexPrmwSignedIndexf("BoxQP starting ... \n");
		mexEvalString("drawnow;");
	}
	// Test malloc results
	if ((fvec==NULL) || (gvec==NULL) || (bufveca==NULL) || (bufvecb==NULL))
	{
		mexPrmwSignedIndexf("BoxQP: memory allocation failed ... \n");
		mexEvalString("drawnow;");return;
	}
	// First, compute some local params
	sig=1.0;kalpha=.5;
	sumalpha=kalpha;
	//L=dabsf(maxeig(Amat,n));
	L=sqrt(doubdot(Amat,Amat,n2));
	cputime=start_time;
	// mwSignedIndexial pomwSignedIndex
	for (i=0;i<n;i++) {xvec[i]=(lvec[i]+uvec[i])/2.0;}
	while ((precision_flag+iteration_flag)==0)
	{
		// Compute gradient and obj
		alpha=2.0;beta=0.0;
		cblas_dgemv(CblasColMajor,CblasNoTrans,n,n,alpha,Amat,n,xvec,incx,beta,gvec,incx);
		fval=0.5*doubdot(gvec,xvec,n);
		// update gradient's weighted average 
		alpha=kalpha;
		cblas_daxpy(n,alpha,gvec,incx,fvec,incx);
		// find a projection of x-Gmu/L on feasible set 
		cblas_dcopy(n,xvec,incx,bufveca,incx);
		alpha=-1.0/L;
		cblas_daxpy(n,alpha,gvec,incx,bufveca,incx);
		for (i=0;i<n;i++){bufveca[i]=dminif(-dminif(-bufveca[i],-lvec[i]),uvec[i]);}
		// project the weigthed gradient z= max([l';min([u';((l+u)/2-w)'])]);
		alpha=-(sig/L);
		cblas_dcopy(n,fvec,incx,bufvecb,incx);
		cblas_dscal(n,alpha,bufvecb,incx);
		for (i=0;i<n;i++){bufvecb[i]=dminif(-dminif(-0.5*(lvec[i]+uvec[i])-bufvecb[i],-lvec[i]),uvec[i]);}
		// update X
		kalpha=exp(0.1*log((double)(k)+2.0))/2.0;
		sumalpha+=kalpha;lambda=alpha/sumalpha;
		for (i=0;i<n;i++){xvec[i]=lambda*bufvecb[i]+(1.0-lambda)*bufveca[i];}
		// Report status
		cputime=((double)clock()-start_time)/CLOCKS_PER_SEC;
		if ((k%Nperiod==0)||(((double)(clock())/CLOCKS_PER_SEC-last_time)>=900))
		{
			last_time=(double)(clock())/CLOCKS_PER_SEC;
			if (k>=MaxIter) iteration_flag=1;
			// report iteration, gap and time left
			if (verbose>=1.0){
				if (k>0){
					left_h=(mwSignedIndex)floor(cputime/3600);left_m=(mwSignedIndex)floor(cputime/60-left_h*60);left_s=(mwSignedIndex)floor(cputime-left_h*3600-left_m*60);
				mexPrmwSignedIndexf("Iter: %.3e   Obj: %.4e   CPU Time: %2dh %2dm %2ds\n",(double)(k),fval,left_h,left_m,left_s);
				mexEvalString("drawnow;");}
			}
		}
		k++;
	}
	// Return primal solution
	for (i=0;i<n;i++){xvec[i]=bufveca[i];}
	// Free everything
	free(fvec);
	free(gvec);
	free(bufveca);
	free(bufvecb);
}


