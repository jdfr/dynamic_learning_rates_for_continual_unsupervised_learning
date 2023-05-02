function [TE] = TopographicErrorSOM(Model,Samples)
% Compute the topographic error for a GHBSOM model
% Input:
%   Model = SOM Model 
%	Samples = Test samples (one sample per column)
% Output:
%   TE = The topographic error (real number in the interval [0,1])

Distances = GetDistancesMapGHBSOM(Model.Prototypes,Samples,1);
NumNeurons = size(Model.TopolDist,1);
[Minima,Winners] = min(Distances);
WinnerLinearIndices = sub2ind(size(Distances),Winners,1:size(Samples,2));
Distances(WinnerLinearIndices) = inf;
[Minima,SecondWinners] = min(Distances);

TE=0;
for NdxNeuro=1:NumNeurons,    
    MySamples=(Winners==NdxNeuro); 
        TopolDistances = Model.TopolDist{NdxNeuro}(SecondWinners(MySamples));
        TE = TE + sum(TopolDistances>1);
end
TE = TE/size(Samples,2);

