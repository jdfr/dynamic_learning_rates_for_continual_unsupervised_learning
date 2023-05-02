function [Handle]=PlotGNG(MaxUnits, Means, Connections)

%Handle=figure;
hold on
plot(Means(1,:),Means(2,:),'or')

for NdxUnit=1:MaxUnits
    if isfinite(Means(1,NdxUnit))
        NdxNeighbors=find(Connections(NdxUnit,:));
        for NdxMyNeigh=1:numel(NdxNeighbors)
            line([Means(1,NdxUnit) Means(1,NdxNeighbors(NdxMyNeigh))],...
                [Means(2,NdxUnit) Means(2,NdxNeighbors(NdxMyNeigh))]);
        end
    end
end

%hold off

