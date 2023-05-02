function [Samples,SampleProjections,Densities]=GenerateSimple3D(DatasetName,StdNoise,NumSamples)
% Generate random samples from some simple 3D datasets
% Inputs:
%   DatasetName=Name of the dataset to use
%   StdNoise=Standard deviation of the Gaussian noise to use
%   NumSamples=Number of samples to generate
% Outputs:
%   Samples=Random samples from the 3D distribution
%   SampleProjections=Projections of Samples on the underlying manifold
%   Densities=Exact probability densities of the samples

%Samples=zeros(3,NumSamples);
SampleProjections=zeros(3,NumSamples);
%Densities=zeros(1,NumSamples);
switch DatasetName
    case 'UnitCircle'
        % Circle of unit radius on the XY plane, with center in the coordinate origin
        Angle=2*pi*rand(1,NumSamples);
        Radius=sqrt(rand(1,NumSamples));
        SampleProjections(1,:)=Radius.*cos(Angle);
        SampleProjections(2,:)=Radius.*sin(Angle);
        Samples=SampleProjections;
        Samples(3,:)=Samples(3,:)+StdNoise*randn(1,NumSamples);
        Densities=(1/pi)*...
            (1/(StdNoise*sqrt(2*pi)))*exp(-(0.5/(StdNoise^2))*Samples(3,:).^2);
    case 'UnitSquare'
        % Square of unit side length on the XY plane, with vertices in
        % (0,0,0), (1,0,0), (1,1,0) and (0,1,0)
        SampleProjections(1,:)=rand(1,NumSamples);
        SampleProjections(2,:)=rand(1,NumSamples);
        Samples=SampleProjections;
        Samples(3,:)=Samples(3,:)+StdNoise*randn(1,NumSamples);
        Densities=(1/(StdNoise*sqrt(2*pi)))*exp(-(0.5/(StdNoise^2))*Samples(3,:).^2);        
    case 'UnitSegment'
        % Segment of unit length in the X axis, with vertices in
        % (0,0,0) and (1,0,0)
        SampleProjections(1,:)=rand(1,NumSamples);
        Samples=SampleProjections;
        Samples(2:3,:)=Samples(2:3,:)+StdNoise*randn(2,NumSamples);
        Densities=(1/(2*pi*(StdNoise^2)))*exp(-(0.5/(StdNoise^2))*...
            sum(Samples(2:3,:).^2,1));
    case 'UnitSphere'
        % Gaussian data inside the unit sphere
        Samples=zeros(3,NumSamples);
        NumSamplesRemaining=NumSamples;
        while NumSamplesRemaining>0
            PreSamples=randn(3,NumSamplesRemaining);
            PreSamplesDists=sum(PreSamples.^2,1);
            NdxNewSamples=find(PreSamplesDists<1);
            NumNewSamples=numel(NdxNewSamples);
            Samples(:,NumSamples-NumSamplesRemaining+1:NumSamples-NumSamplesRemaining+NumNewSamples)=...
                PreSamples(:,NdxNewSamples);
            NumSamplesRemaining=NumSamplesRemaining-NumNewSamples;
        end
        % Note that 1/chi2cdf(1,3)=5.031496081211187, where chi2cdf(1,3)
        % is the probability that a sample from the unit spherical Gaussian
        % lies within the unit sphere
        Densities=5.031496081211187*((2*pi)^(-3/2))*exp(-0.5*sum(Samples.^2,1));        

end