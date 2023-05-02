function PlotCompetitivePCA(Model)

[~,score,~,~,~,~] = pca(Model.Prototypes');

plot(score(:,1),score(:,2),'or','LineWidth',2,'MarkerFaceColor',[1 0 0],'MarkerSize',7);

%% Version to plot the two first principal components of the Means after
%% performing a PCA to the Means.
%function [Handle]=PlotCompetitivePCA(Modelos)
%
%Handle=figure;
%
%for i=1:length(Modelos)
%subplot(3,4,i);
%Model = Modelos{i};
%Prototypes=Model.Prototypes;

%[~,score,~,~,~,~] = pca(Model.Means');


%% Draw the neurons
%for NdxNeuron=1:size(Prototypes,2)
%    CircleX=score(NdxNeuron,1);
%    CircleY=score(NdxNeuron,2);
%    plot(CircleX,CircleY,'or','LineWidth',2,'MarkerFaceColor',[1 0 0],'MarkerSize',7);
%    %xlim([-0.6 0.15]);
%    %ylim([-0.4 0.4]);
%    hold on;
%end
%end
%%hold off


