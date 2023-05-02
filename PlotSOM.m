function [Handle]=PlotSOM(NumRowsMap, NumColsMap, Prototypes)

hold on;
plot(Prototypes(1,:),Prototypes(2,:),'or','LineWidth',2,'MarkerFaceColor',[1 0 0],'MarkerSize',7);


% Draw the horizontal connections
for NdxRow=1:NumRowsMap
    for NdxCol=1:(NumColsMap-1)
        LineX=[Prototypes(1,NdxRow,NdxCol) Prototypes(1,NdxRow,NdxCol+1)];
        LineY=[Prototypes(2,NdxRow,NdxCol) Prototypes(2,NdxRow,NdxCol+1)];        
        plot(LineX,LineY,'-','LineWidth',2,'Color',[0 0 1]);
    end
end

% Draw the vertical connections
for NdxRow=1:(NumRowsMap-1)
    for NdxCol=1:NumColsMap
        LineX=[Prototypes(1,NdxRow,NdxCol) Prototypes(1,NdxRow+1,NdxCol)];
        LineY=[Prototypes(2,NdxRow,NdxCol) Prototypes(2,NdxRow+1,NdxCol)];        
        plot(LineX,LineY,'-','LineWidth',2,'Color',[0 0 1]);
    end
end   



if false

Handle=figure;

for i=1:length(Modelos)
subplot(3,4,i);
Model = Modelos{i};

NumRowsMap=Model.NumRowsMap;
NumColsMap=Model.NumColsMap;
Prototypes=Model.Prototypes;

% Draw the horizontal connections
for NdxRow=1:NumRowsMap
    for NdxCol=1:(NumColsMap-1)
        LineX=[Prototypes(1,NdxRow,NdxCol) Prototypes(1,NdxRow,NdxCol+1)];
        LineY=[Prototypes(2,NdxRow,NdxCol) Prototypes(2,NdxRow,NdxCol+1)];        
        plot(LineX,LineY,'-','LineWidth',2,'Color',[0 0 1]);
        hold on
    end
end

% Draw the vertical connections
for NdxRow=1:(NumRowsMap-1)
    for NdxCol=1:NumColsMap
        LineX=[Prototypes(1,NdxRow,NdxCol) Prototypes(1,NdxRow+1,NdxCol)];
        LineY=[Prototypes(2,NdxRow,NdxCol) Prototypes(2,NdxRow+1,NdxCol)];        
        plot(LineX,LineY,'-','LineWidth',2,'Color',[0 0 1]);
    end
end   

% Draw the neurons
for NdxRow=1:NumRowsMap
    for NdxCol=1:NumColsMap
        CircleX=Prototypes(1,NdxRow,NdxCol);
        CircleY=Prototypes(2,NdxRow,NdxCol);        
        plot(CircleX,CircleY,'or','LineWidth',2,'MarkerFaceColor',[1 0 0],'MarkerSize',7);
    end
end

end

%hold off

end
