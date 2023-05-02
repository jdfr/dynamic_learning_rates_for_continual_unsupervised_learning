function [NewCentroids] = GetCentroidsGHNF(Model)
% Get recursively the centroids of the GHNF model
% E.J. Palomo
% Inputs:
%   Model=GHNF model
% Output:
%   NewCentroids=Planar prototypes of the GHNF model

NdxValidNeurons = find(isfinite(Model.Means(1,:)));
NewCentroids = [];

for NdxNeuro=NdxValidNeurons,        
    if ~isempty(Model.Child{NdxNeuro}),        
        ChildCentroids = GetCentroidsGHNF(Model.Child{NdxNeuro}); 
        NewCentroids = [NewCentroids ChildCentroids];
    else
        NewCentroids = [NewCentroids Model.Means(:,NdxNeuro)];
    end    
end
