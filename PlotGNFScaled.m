function [Handle]=PlotGNFScaled(Model,MyAxis)
% Plot a GNF model in 2D
% E.J. Palomo 2015
% Inputs:
%   Model = GNF model
%   MyAxis = Vector of four components to scale the plot

Handle = [];
% Handle=figure;
hold on
plot(Model.Means(1,:)*MyAxis(2)+MyAxis(1),Model.Means(2,:)*MyAxis(4)+MyAxis(3),'or','LineWidth',2,'MarkerFaceColor',[1 0 0],'MarkerSize',5)

for NdxUnit=1:Model.MaxUnits
    if isfinite(Model.Means(1,NdxUnit))
        NdxNeighbors=find(Model.SpanningTree(NdxUnit,:));
        for NdxMyNeigh=1:numel(NdxNeighbors)
            line([Model.Means(1,NdxUnit)*MyAxis(2)+MyAxis(1) Model.Means(1,NdxNeighbors(NdxMyNeigh))*MyAxis(2)+MyAxis(1)],...
                [Model.Means(2,NdxUnit)*MyAxis(4)+MyAxis(3) Model.Means(2,NdxNeighbors(NdxMyNeigh))*MyAxis(4)+MyAxis(3)],'Color',[1 1 0]);
        end
    end
end

hold off
