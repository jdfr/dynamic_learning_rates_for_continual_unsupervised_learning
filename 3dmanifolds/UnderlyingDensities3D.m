function [Densities]=UnderlyingDensities3D(NdxDataset,TestSamples,StdNoise)
% Estimate the underlying probability density function for some 3D datasets

% Values of the Precision parameter to get more than 1000 kernel centers
Precisions=[40 40 40 1000 1000 1500 1500];

% Obtain kernel centers for the underlying manifold
KernelCenters=GenerateManifolds3D(NdxDataset,Precisions(NdxDataset),0);

% Compute the squared distances from the kernel centers to the test samples
KernelSampleDists=sqrdistanceMEX(KernelCenters,TestSamples);

% Compute the spherical Gaussian pdf with the specified standard deviation
GaussianDensities=((2*pi)^(-3/2))*(StdNoise^(-3))*exp(-0.5*KernelSampleDists(:)/(StdNoise^2));

% The underlying density is the mean of the Gaussian densities
Densities=mean(reshape(GaussianDensities,size(KernelSampleDists)),1);
