function [Handle]=PlotCompetitive(prots)%Model)%Modelos)

plot(prots(1,:),prots(2,:),'or','LineWidth',2,'MarkerFaceColor',[1 0 0],'MarkerSize',7);

%plot(Model.Prototypes(1,:),Model.Prototypes(2,:),'or','LineWidth',2,'MarkerFaceColor',[1 0 0],'MarkerSize',7);


% Handle=figure;
% 
% %for i=1:length(Modelos)
% %subplot(3,4,i);
% %Model = Modelos{i};
% Prototypes=Model.Prototypes;
% 
% Draw the neurons
% for NdxNeuron=1:size(Prototypes,2)
%     CircleX=Prototypes(1,NdxNeuron);
%     CircleY=Prototypes(2,NdxNeuron);
%     plot(CircleX,CircleY,'or','LineWidth',2,'MarkerFaceColor',[1 0 0],'MarkerSize',7);
%     %xlim([0 0.5]);
%     %ylim([0.22 0.45]);
%     hold on;
% end
% %end
% 
% %hold off
