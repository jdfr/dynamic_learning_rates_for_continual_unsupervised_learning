function [Model]=TrainKohonenSOM(Samples,NumRowsMap,NumColsMap,NumSteps)
% Train a Self-Organizing Map (SOM)
% Inputs:
%	Samples=Training samples (one sample per column)
%	NumSteps=Number of training steps
% Output:
%   Model=Trained SOM model

[Dimension,NumSamples]=size(Samples);
NumNeuro=NumRowsMap*NumColsMap;  
Model.NumColsMap=NumColsMap;
Model.NumRowsMap=NumRowsMap;

% Prototypes initialization
NumPatIni=max([Dimension+1,ceil(NumSamples/(NumRowsMap*NumColsMap))]);
Model.Prototypes=zeros(Dimension,NumRowsMap,NumColsMap);
for NdxRow=1:NumRowsMap
    for NdxCol=1:NumColsMap
        MySamples=Samples(:,ceil(NumSamples*rand(1,NumPatIni)));
        Model.Prototypes(:,NdxRow,NdxCol)=mean(MySamples,2);         
    end
end        

% Precompute topological distances
[AllXCoords AllYCoords]=ind2sub([NumRowsMap NumColsMap],1:NumNeuro);
AllCoords(1,:)=AllXCoords;
AllCoords(2,:)=AllYCoords;
TopolDist=cell(NumNeuro,1);
for NdxNeuro=1:NumNeuro    
    TopolDist{NdxNeuro}=sum((repmat(AllCoords(:,NdxNeuro),1,NumNeuro)-AllCoords).^2,1);
end
Model.TopolDist = TopolDist;

% Training
Distances=zeros(NumNeuro,NumSamples);
MaxRadius=(NumRowsMap+NumColsMap)/8;
NdxPermSamples = randperm(NumSteps);

for NdxStep=1:NumSteps    
    %NdxSample = ceil(NumSamples*rand(1));
    NdxSample = mod(NdxPermSamples(NdxStep)-1,NumSamples)+1;    
    MySample=Samples(:,NdxSample);
    if NdxStep<0.5*NumSteps   
        % Ordering phase: linear decay
        LearningRate=0.4*(1-NdxStep/NumSteps);
        MyRadius=MaxRadius*(1-(NdxStep-1)/NumSteps);
    else
        % Convergence phase: constant
        LearningRate=0.01;
        MyRadius=0.1;
    end
    
    % Euclidean distance (standard Kohonen's SOFM)
    RepMySample=repmat(MySample,1,NumNeuro);        
    MyDistances=sqrt(sum((RepMySample-Model.Prototypes(:,:)).^2,1));            			    
    Distances(:,NdxSample) = MyDistances';        
   
    % Update the neurons
    [~,NdxWinner]=min(MyDistances);
    Model.Winners(NdxSample)=NdxWinner;
    Coef=repmat(LearningRate*exp(-TopolDist{NdxWinner}/(MyRadius^2)),Dimension,1);
    Model.Prototypes(:,:)=Coef.*repmat(MySample,1,NumNeuro)+...
        (1-Coef).*Model.Prototypes(:,:);    
end

% Store the means (same than prototypes but arranged in a row)
Model.Means=zeros(Dimension,NumNeuro);
for NdxRow=1:NumRowsMap
    for NdxCol=1:NumColsMap                
        NdxNeuro = sub2ind([NumRowsMap,NumColsMap],NdxRow,NdxCol);
        Model.Means(:,NdxNeuro)=Model.Prototypes(:,NdxRow,NdxCol);
    end
end

% Store the connections
Connections = zeros(NumNeuro,NumNeuro);
for NdxNeuro=1:NumNeuro
    Connections(NdxNeuro,:) = (TopolDist{NdxNeuro} == 1);
end
Model.Connections = sparse(Connections);
Model.Samples=Samples;  

    
    
        
