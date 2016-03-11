// MEX wrapper for the sparse_rank_one function in sparsesvd.c
// function [X,U,u,F,k] = sparse_rank_one_mex(A,rho,tol,MaxIter,Info,[F],[X0],[k0]) 

// Mex wrapper for Box QP

// Calling sequence: BoxQP(A, l, u, maxiters,info,verbose)

// Last Modified: A. d'Aspremont, Laurent El Ghaoui, Dec. 2004.
// http://www.carva.org/alexandre.daspremont


#include "BoxQP.h" 

void mexFunction(mwSignedIndex nlhs, mxArray *plhs[], mwSignedIndex nrhs, const mxArray *prhs[])
{
  double *Amat, *lvec, *uvec, *xvecout;
  double MaxIter,info,verbose;
  mwSignedIndex n;

  Amat = mxGetPr(prhs[0]);
  n = mxGetN(prhs[0]);
  lvec = mxGetPr(prhs[1]);
  uvec = mxGetPr(prhs[2]);
  MaxIter = mxGetScalar(prhs[3]);
  info = mxGetScalar(prhs[4]);
  verbose = mxGetScalar(prhs[5]);

  plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
  xvecout = mxGetPr(plhs[0]);

  BoxQP(Amat,n,lvec,uvec,(mwSignedIndex)(MaxIter),(mwSignedIndex)(info),(mwSignedIndex)(verbose),xvecout);
}

