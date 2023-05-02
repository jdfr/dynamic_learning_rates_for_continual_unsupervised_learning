function [Handle]=PlotGNGPCA(MaxUnits, Means, Connections)


finiteMeans = isfinite(Means(1,:);

[~,score,~,~,~,~] = pca(Means(:,finiteMeans)');

score = score';

fullscore = zeros(size(Means));
fullscore(:,finiteMeans) = score;

%Handle=figure;
hold on
plot(score(1,:),score(2,:),'or')

for NdxUnit=1:MaxUnits
    if finiteMeans(NdxUnit)
        NdxNeighbors=find(Connections(NdxUnit,:));
        for NdxMyNeigh=1:numel(NdxNeighbors)
            line([fullscore(1,NdxUnit) fullscore(1,NdxNeighbors(NdxMyNeigh))],...
                [fullscore(2,NdxUnit) fullscore(2,NdxNeighbors(NdxMyNeigh))]);
        end
    end
end

%hold off

