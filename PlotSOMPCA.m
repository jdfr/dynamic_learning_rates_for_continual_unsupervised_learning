function PlotSOMPCA(NumRowsMap, NumColsMap, Prototypes)

[~,score,~,~,~,~] = pca(Prototypes(:,:)');

score = score';

hold on;
plot(score(1,:),score(2,:),'or','LineWidth',2,'MarkerFaceColor',[1 0 0],'MarkerSize',7);

score = reshape(score, size(Prototypes));

% Draw the horizontal connections
for NdxRow=1:NumRowsMap
    for NdxCol=1:(NumColsMap-1)
        LineX=[score(1,NdxRow,NdxCol) score(1,NdxRow,NdxCol+1)];
        LineY=[score(2,NdxRow,NdxCol) score(2,NdxRow,NdxCol+1)];        
        plot(LineX,LineY,'-','LineWidth',2,'Color',[0 0 1]);
    end
end

% Draw the vertical connections
for NdxRow=1:(NumRowsMap-1)
    for NdxCol=1:NumColsMap
        LineX=[score(1,NdxRow,NdxCol) score(1,NdxRow+1,NdxCol)];
        LineY=[score(2,NdxRow,NdxCol) score(2,NdxRow+1,NdxCol)];        
        plot(LineX,LineY,'-','LineWidth',2,'Color',[0 0 1]);
    end
end   


% Version to plot the two first principal components of the Means after
% performing a PCA to the Means.
function [Handle]=PlotSOMPCA_old(Modelos)

Handle=figure;

for i=1:length(Modelos)
subplot(3,4,i);
Model = Modelos{i};

[~,score,~,~,~,~] = pca(Model.Means');

NumRowsMap=Model.NumRowsMap;
NumColsMap=Model.NumColsMap;

% Draw the horizontal connections
for NdxRow=1:NumRowsMap
    for NdxCol=1:(NumColsMap-1)
        NdxNeuro = sub2ind([NumRowsMap,NumColsMap],NdxRow,NdxCol);
        NdxNeuro2 = sub2ind([NumRowsMap,NumColsMap],NdxRow,NdxCol+1);
        LineX=[score(NdxNeuro,1) score(NdxNeuro2,1)];
        LineY=[score(NdxNeuro,2) score(NdxNeuro2,2)];        
        plot(LineX,LineY,'-','LineWidth',2,'Color',[1 1 0]);
        hold on
    end
end

% Draw the vertical connections
for NdxRow=1:(NumRowsMap-1)
    for NdxCol=1:NumColsMap
        NdxNeuro = sub2ind([NumRowsMap,NumColsMap],NdxRow,NdxCol);
        NdxNeuro2 = sub2ind([NumRowsMap,NumColsMap],NdxRow+1,NdxCol);
        LineX=[score(NdxNeuro,1) score(NdxNeuro2,1)];
        LineY=[score(NdxNeuro,2) score(NdxNeuro2,2)];
        plot(LineX,LineY,'-','LineWidth',2,'Color',[1 1 0]);
    end
end   

% Draw the neurons
for NdxRow=1:NumRowsMap
    for NdxCol=1:NumColsMap
        NdxNeuro = sub2ind([NumRowsMap,NumColsMap],NdxRow,NdxCol);
        CircleX=score(NdxNeuro,1);
        CircleY=score(NdxNeuro,2);        
        plot(CircleX,CircleY,'or','LineWidth',2,'MarkerFaceColor',[1 0 0],'MarkerSize',5);
    end
end

end

hold off


% [~,score,~,~,~,~] = pca(Model.Means');
% plot(score(:,1),score(:,2),'or','LineWidth',2,'MarkerFaceColor',[1 0 0],'MarkerSize',5)
% hold on
% 
% NumNeuro = Model.NumRowsMap*Model.NumColsMap;
% for NdxUnit=1:NumNeuro
%     if isfinite(Model.Means(1,NdxUnit))
%         NdxNeighbors=find(Model.Connections(NdxUnit,:));
%         for NdxMyNeigh=1:numel(NdxNeighbors)
%             line([score(NdxUnit,1) score(NdxNeighbors(NdxMyNeigh),1)],...
%                 [score(NdxUnit,2) score(NdxNeighbors(NdxMyNeigh),2)],'Color',[1 1 0],'LineWidth',2);
%         end
%     end
% end
% 
% hold off

