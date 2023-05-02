#include "mex.h"
#include <math.h>


/* 
Coded by Ezequiel López-Rubio, April 2012.

In order to compile this function, type the following at the Matlab prompt:
>> mex distanceMEX.c

This function computes the Euclidean distance matrix
E = distanceMEX(A,B)

Inputs:
	A is a (DxM) matrix 
	B is a (DxN) matrix

Output:
    E is a (MxN) matrix with the Euclidean distances between vectors in A and B

*/

void mexFunction(int nlhs, mxArray* plhs[],
                 int nrhs, const mxArray* prhs[])
{   
	int D,M,N,NdxA,NdxB,NdxDim;
	double *ptrA,*ptrB,*ptrE;
	mxArray *E;
	register double result,diff;

    /* Obtain working variables */
    ptrA=mxGetPr(prhs[0]);
	ptrB=mxGetPr(prhs[1]);
	D=mxGetM(prhs[0]);
	M=mxGetN(prhs[0]);
	N=mxGetN(prhs[1]);

	/* Create output mxArray */
	plhs[0]=mxCreateDoubleMatrix(M,N,mxREAL);
	ptrE=mxGetPr(plhs[0]);

	/* Compute Euclidean distance */
	for(NdxA=0;NdxA<M;NdxA++)
	{
		for(NdxB=0;NdxB<N;NdxB++)
		{
			result=0.0;
			for(NdxDim=0;NdxDim<D;NdxDim++)
			{
				diff=ptrA[NdxA*D+NdxDim]-ptrB[NdxB*D+NdxDim];
				result+=diff*diff;
			}
			ptrE[NdxA+NdxB*M]=sqrt(result);
		}
	}

}    
