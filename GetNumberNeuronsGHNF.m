function [NumNeuronas,NumNeuronasHoja] = GetNumberNeuronsGHNF(Model)
% Get the numer of leaf and non-leaf neurons of the GHNF model.
% Input:
%   Model = GHNG Model 
% Output:
%   NumNeuronas = Number of neurons of the GHNG model
%   NumNeuronasHoja = Number of leaf neurons of the GHNG model
    
NdxValidNeurons = find(isfinite(Model.Means(1,:)));
NumNeuronas = numel(NdxValidNeurons);
NumNeuronasHoja = NumNeuronas;

for NdxNeuro=NdxValidNeurons
    if isfield(Model,'Child') && ~isempty(Model.Child{NdxNeuro})
        [NumNeuronasHijo,NumNeuronasHojaHijo] = GetNumberNeuronsGHNF(Model.Child{NdxNeuro});
        NumNeuronas = NumNeuronas + NumNeuronasHijo;
        NumNeuronasHoja = NumNeuronasHoja + NumNeuronasHojaHijo - 1;
    end    
end

