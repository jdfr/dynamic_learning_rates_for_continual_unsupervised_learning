function [TE] = TopographicErrorGHNG(Model,Samples)
% Compute the topographic error for a GHNG model
% Input:
%   Model = GHNG Model 
%	Samples = Test samples (one sample per column)
% Output:
%   TE = The topographic error (real number in the interval [0,1])

TE = 0;
if ~isempty(Samples)
    Distances = GetDistancesMapGHBSOM(Model.Means,Samples,1);
    ValidNeurons = find(isfinite(Model.Means(1,:)));
    NumNeurons = length(ValidNeurons);
    [Minima,Winners] = min(Distances);
    WinnerLinearIndices = sub2ind(size(Distances),Winners,1:size(Samples,2));
    Distances(WinnerLinearIndices) = inf;
    [Minima,SecondWinners] = min(Distances);
    
    for NdxNeuro=ValidNeurons,    
        MySamples = (Winners==NdxNeuro);         
        if isfield(Model,'Child') && ~isempty(Model.Child{NdxNeuro})        
            ThisTE = TopographicErrorGHNG(Model.Child{NdxNeuro},Samples(:,MySamples)); 
            TE = TE + ThisTE*nnz(MySamples);
        else    
            Connections = full(Model.Connections);
            MyConnections = Connections(NdxNeuro,SecondWinners(MySamples));
            TE = TE + nnz(MyConnections==0);    % Aqui habia un error porque se sumaban los resultados de find que son índices
        end            
    end
    TE = TE/size(Samples,2);
end

