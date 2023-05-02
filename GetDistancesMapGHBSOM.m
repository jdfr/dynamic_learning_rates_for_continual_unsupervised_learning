function Distances = GetDistancesMapGHBSOM(Prototypes,Samples,NdxType)
% Get the distances from the neurons of a map of the GHBSOM model to the 
% input samples of the map
% Inputs:
%   Prototypes = Neurons' weight vectors of the map
%	Samples = Test samples (one sample per column)
%   NdxType = Bregman divergence type (see valid types below), as discussed in:
%       Banerjee, A.; Merugu, S.; Dhillon, I. S.; Ghosh, J. (2005). Clustering 
%       with Bregman divergences. Journal of Machine Learning Research 6,
%       1705-1749.
% Output:
%   Distances = Distances from each neuron of the map to input samples

if isempty(Samples),
    Distances = [];
else       
    NumNeurons = size(Prototypes,2)*size(Prototypes,3);
    NumSamples = size(Samples,2);
    Distances = [];
%     Distances = zeros(NumNeurons,NumSamples);
    
    for NdxNeuro=1:NumNeurons,
%         epsilon = 10^-4;
        [i,j] = ind2sub([size(Prototypes,2) size(Prototypes,3)],NdxNeuro);
        WeightVector = Prototypes(:,i,j);
        RepW = repmat(WeightVector,1,NumSamples);
        Quotient = Samples./RepW;
    %     Quotient(Quotient==0) = epsilon;

        switch NdxType
            case 1
                % Squared Euclidean distance (standard Kohonen's SOFM)
                Distances(NdxNeuro,:) = sum((Samples-RepW).^2,1);            			
            case 2
                % Generalized I-divergence, strictly positive data only
                Distances(NdxNeuro,:) = sum(Samples.*log(Quotient)-Samples+RepW,1);
            case 3
                % Itakura-Saito distance, strictly positive data only
                Distances(NdxNeuro,:) = sum(Quotient-log(Quotient)-1,1);
            case 4
                % Exponential loss, strictly positive data only
                Distances(NdxNeuro,:) = sum(exp(Samples)-exp(RepW)-(Samples-RepW).*exp(RepW),1);
            case 5
                % Logistic loss, data inside the [0,1] interval only
                Distances(NdxNeuro,:) = sum(Samples.*log(Quotient)+(1-Samples).*log( (1-Samples)./(1-RepW)),1);
            case 6
                % Kullback-Leibler
                Distances(NdxNeuro,:) = sum(Samples.*log(Quotient),1);
        end
    end
end
