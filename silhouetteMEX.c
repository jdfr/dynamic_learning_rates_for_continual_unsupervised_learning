#include "mex.h"
#include <math.h>
#include <float.h>


/* 
Coded by Ezequiel López-Rubio, September 2013.

In order to compile this function, type the following at the Matlab prompt:
>> mex silhouetteMEX.c

This function computes the silhouette coefficients
Silhouettes = silhouetteMEX(Samples,Labels,NumClusters);

Samples is a D x NumSamples matrix
Labels is a 1 x NumSamples matrix which contains integer labels in the range 0...NumClusters-1
Silhouettes is a 1 x NumSamples matrix with the silhouette coefficients in the range -1...1

*/


/*
This function computes the squared Euclidean distance matrix
E = PairwiseDistance(A,B)

Inputs:
	A is a (DxM) matrix 
	B is a (DxN) matrix

Output:
    E is a (MxN) matrix with the squared Euclidean distances between vectors in A and B

*/

void PairwiseDistance(double *ptrA,double *ptrB,double *ptrE,
                 int D,int M,int N)
{   
	int NdxA,NdxB,NdxDim;
	register double result,diff;

	/* Compute Euclidean distance */
	for(NdxA=0;NdxA<M;NdxA++)
	{
		for(NdxB=NdxA+1;NdxB<N;NdxB++)
		{
			result=0.0;
			for(NdxDim=0;NdxDim<D;NdxDim++)
			{
				diff=ptrA[NdxA*D+NdxDim]-ptrB[NdxB*D+NdxDim];
				result+=diff*diff;
			}
			/*ptrE[NdxA+NdxB*M]=sqrt(result);*/
			ptrE[NdxA+NdxB*M]=result;
		}
	}

}    

/*
This function computes the squared Euclidean distance 
E = EuclideanDistance(A,B)

Inputs:
	A is a (Dx1) vector 
	B is a (Dx1) vector

Output:
    A scalar with the squared Euclidean distances between A and B

*/

double EuclideanDistance(double *ptrA,double *ptrB,
                 int D)
{   
	int NdxDim;
	register double result,diff;

	/* Compute Euclidean distance */

	result=0.0;
	for(NdxDim=0;NdxDim<D;NdxDim++)
	{
		diff=ptrA[NdxDim]-ptrB[NdxDim];
		result+=diff*diff;
	}

	return result;

} 


void mexFunction(int nlhs, mxArray* plhs[],
                 int nrhs, const mxArray* prhs[])
{   
	int NumClusters,NdxSample1,NdxSample2,NdxCluster;
	int D,NumSamples,ClosestCluster;
	double *ptrSamples,*ptrLabels,*ptrSilhouettes;
	double *ptrDissimilarities,*ptrNumClusters,*ptrClusterSizes;
	double ClosestClusterDissimilarity,ThisClusterDissimilarity;
	double ThisDistance;

    /* Obtain working variables */
    ptrSamples=mxGetPr(prhs[0]);
	ptrLabels=mxGetPr(prhs[1]);
	NumClusters=(int)mxGetScalar(prhs[2]);
	D=(int)mxGetM(prhs[0]);
	NumSamples=(int)mxGetN(prhs[0]);

	/* Create output mxArray */
	plhs[0]=mxCreateDoubleMatrix(1,NumSamples,mxREAL);
	ptrSilhouettes=mxGetPr(plhs[0]);

	/* Create auxiliary matrices */
	ptrDissimilarities=mxMalloc(NumClusters*sizeof(double));
	ptrClusterSizes=mxMalloc(NumClusters*sizeof(double));

	/* Compute the cluster sizes */
	memset(ptrClusterSizes,0,NumClusters*sizeof(double));
	for(NdxSample1=0;NdxSample1<NumSamples;NdxSample1++)
	{
		ptrClusterSizes[(int)ptrLabels[NdxSample1]]++;
	}

	/* For each sample compute its silhouette coefficient */
	for(NdxSample1=0;NdxSample1<NumSamples;NdxSample1++)
	{
		/* For each other sample, accumulate the dissimilarity with
		    this sample */
		memset(ptrDissimilarities,0,NumClusters*sizeof(double));
		for(NdxSample2=0;NdxSample2<NumSamples;NdxSample2++)
		{
			ptrDissimilarities[(int)ptrLabels[NdxSample2]]+=
				EuclideanDistance(ptrSamples+NdxSample1*D,
					ptrSamples+NdxSample2*D,D);
		}

		/* Compute the average dissimilarity per cluster, and obtain
		the cluster with the smallest average dissimilarity excluding the cluster
		that the sample belongs to */
		ClosestClusterDissimilarity=DBL_MAX;
		ClosestCluster=0;
		for(NdxCluster=0;NdxCluster<NumClusters;NdxCluster++)
		{
			if (NdxCluster!=ptrLabels[NdxSample1])
			{
				ptrDissimilarities[NdxCluster]/=ptrClusterSizes[NdxCluster];
			}
			else
			{
				ptrDissimilarities[NdxCluster]/=(ptrClusterSizes[NdxCluster]-1.0);
			}
			if ((NdxCluster!=ptrLabels[NdxSample1]) && 
				(ptrDissimilarities[NdxCluster]<ClosestClusterDissimilarity))
			{
				ClosestCluster=NdxCluster;
				ClosestClusterDissimilarity=ptrDissimilarities[NdxCluster];
			}
		}

		/* Compute the silhouette coefficient for this sample */
		ThisClusterDissimilarity=ptrDissimilarities[(int)ptrLabels[NdxSample1]];
		if (ThisClusterDissimilarity>ClosestClusterDissimilarity)
		{
			ptrSilhouettes[NdxSample1]=(ClosestClusterDissimilarity-ThisClusterDissimilarity)/
				ThisClusterDissimilarity;
		}
		else
		{
			ptrSilhouettes[NdxSample1]=(ClosestClusterDissimilarity-ThisClusterDissimilarity)/
				ClosestClusterDissimilarity;
		}
	}

	/* Release dynamic memory */
	mxFree(ptrDissimilarities);
	mxFree(ptrClusterSizes);
}    
