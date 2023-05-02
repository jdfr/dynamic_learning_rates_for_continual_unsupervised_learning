function Data = NormalizeData(Data)
% *************************************************************************
% Description: function that performs a standard normalization of the input
%              data in [0,1]
% Input:
%   Data: input data matrix where rows are features and columns are samples
% Output:
%   Data: normalized input data matrix
% *************************************************************************

for NdxFeat=1:size(Data,1),
    
    % Get the maximum and the minimum value of the feature
    MaxValue = max(Data(NdxFeat,:));
    MinValue = min(Data(NdxFeat,:));
    
    if (MaxValue > MinValue),
        Data(NdxFeat,:) = (Data(NdxFeat,:) - MinValue)/(MaxValue - MinValue);        
    else
        Data(NdxFeat,:) = Data(NdxFeat,:) - MinValue;
    end
end
Data=0.9998*Data+0.0001;