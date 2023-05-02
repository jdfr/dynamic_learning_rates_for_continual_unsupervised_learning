function Connections = makeSOMConnections(Model)
Connections = zeros(Model.NumRowsMap*Model.NumColsMap);
for NdxNeuro=1:size(Connections,1)
    Connections(NdxNeuro,:) = (Model.TopolDist{NdxNeuro} == 1);
end

