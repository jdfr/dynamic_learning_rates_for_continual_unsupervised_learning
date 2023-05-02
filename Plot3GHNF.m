function [Handle]=Plot3GHNF(Model)
% Plot a GHNF model in 3D
% E.J. Palomo
% Inputs:
%   Model = GHNF model

% Handle=figure;
Handle=[];
hold on
Plot3GHNFRec(Model,1);

function Plot3GHNFRec(Model,Level)

if ~isempty(Model) && (Level <= 4),
    
    Colors = [0.4 0.2 0;1 0 0;1 0.65 0;0 1 1];
    MyColor = Colors(Level,:);
    
    % Draw the neurons
    plot3(Model.Means(1,:),Model.Means(2,:),Model.Means(3,:),'or','LineWidth',1,'MarkerFaceColor',MyColor,'MarkerSize',15*(1/Level));

    % Draw the vertical connections
    NumNeurons = size(Model.Means,2);
    for NdxUnit=1:NumNeurons,
        if isfinite(Model.Means(1,NdxUnit))
%             NdxNeighbors = find(Model.Connections(NdxUnit,:));
            NdxNeighbors=find(Model.SpanningTree(NdxUnit,:));
            for NdxMyNeigh=1:numel(NdxNeighbors)
                line([Model.Means(1,NdxUnit) Model.Means(1,NdxNeighbors(NdxMyNeigh))],...
                    [Model.Means(2,NdxUnit) Model.Means(2,NdxNeighbors(NdxMyNeigh))],...
                    [Model.Means(3,NdxUnit) Model.Means(3,NdxNeighbors(NdxMyNeigh))],'Color',MyColor);
            end
        end
    end

    % Draw the child maps
    for NdxNeuro=1:numel(Model.Child),    
        Plot3GHNFRec(Model.Child{NdxNeuro},Level+1);    
    end
end