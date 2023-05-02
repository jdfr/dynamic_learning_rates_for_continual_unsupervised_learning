function [Densities]=ExactDensitySimple3D(DatasetName,StdNoise,Samples)
% Generate random samples from some simple 3D datasets
% Inputs:
%   DatasetName=Name of the dataset to use
%   StdNoise=Standard deviation of the Gaussian noise to use
%   NumSamples=Number of samples to generate
% Outputs:
%   Samples=Random samples from the 3D distribution
%   SampleProjections=Projections of Samples on the underlying manifold
%   Densities=Exact probability densities of the samples

switch DatasetName
    case 'UnitCircle'
        % Circle of unit radius on the XY plane, with center in the coordinate origin
        Densities=(1/pi)*...
            (1/(StdNoise*sqrt(2*pi)))*exp(-(0.5/(StdNoise^2))*Samples(3,:).^2);
        SqrDists=sum(Samples(1:2,:).^2,1);
        Densities(SqrDists>1)=0;
    case 'UnitSquare'
        % Square of unit side length on the XY plane, with vertices in
        % (0,0,0), (1,0,0), (1,1,0) and (0,1,0)
        Densities=(1/(StdNoise*sqrt(2*pi)))*exp(-(0.5/(StdNoise^2))*Samples(3,:).^2); 
        Densities((Samples(1,:)<0) | (Samples(1,:)>1))=0;
        Densities((Samples(2,:)<0) | (Samples(2,:)>1))=0;
    case 'UnitSegment'
        % Segment of unit length in the X axis, with vertices in
        % (0,0,0) and (1,0,0)
        Densities=(1/(2*pi*(StdNoise^2)))*exp(-(0.5/(StdNoise^2))*...
            sum(Samples(2:3,:).^2,1));
        Densities((Samples(1,:)<0) | (Samples(1,:)>1))=0;
    case 'UnitSphere'
        % Gaussian data inside the unit sphere
        % Note that 1/chi2cdf(1,3)=5.031496081211187, where chi2cdf(1,3)
        % is the probability that a sample from the unit spherical Gaussian
        % lies within the unit sphere
        Densities=5.031496081211187*((2*pi)^(-3/2))*exp(-0.5*sum(Samples.^2,1));   
        SqrDists=sum(Samples.^2,1);
        Densities(SqrDists>1)=0;

end