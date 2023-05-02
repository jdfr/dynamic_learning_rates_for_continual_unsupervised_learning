function [Winners,Errors] = TestGNG(Model,Samples)
% Simulate any Neural Network Model depending on the distance
% Inputs:
%   Model = GNG model
%   Samples = Matrix of input data (columns are samples)

NumSamples=size(Samples,2);
Winners = zeros(NumSamples,1);
Errors = zeros(NumSamples,1);
Prototypes = Model.Means;
NumNeuro=size(Prototypes,2);

% Main loop
for NdxSample=1:NumSamples,            
    SquaredDistances=sum((repmat(Samples(:,NdxSample),1,NumNeuro)-Prototypes).^2,1);
    [Minimum,NdxWinner]=min(SquaredDistances);
    Winners(NdxSample)=NdxWinner;
    Errors(NdxSample)=Minimum;    
end