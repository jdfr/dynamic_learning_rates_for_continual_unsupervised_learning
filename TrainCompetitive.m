function [Model]=TrainCompetitive(Samples,NumNeuro,NumSteps)
% Train a Competitive Neural Network
% Inputs:
%	Samples=Training samples (one sample per column)
%	NumSteps=Number of training steps
% Output:
%   Model=Trained Competitive model

[Dimension,NumSamples]=size(Samples);

% Prototypes initialization
NumPatIni=max([Dimension+1,ceil(NumSamples/NumNeuro)]);
Model.Prototypes=zeros(Dimension,NumNeuro);
for NdxNeuro=1:NumNeuro   
    MySamples=Samples(:,ceil(NumSamples*rand(1,NumPatIni)));
    Model.Prototypes(:,NdxNeuro)=mean(MySamples,2);             
end        

% Training
Distances=zeros(NumNeuro,NumSamples);
NdxPermSamples = randperm(NumSteps);

for NdxStep=1:NumSteps    
    %NdxSample = ceil(NumSamples*rand(1));
    NdxSample = mod(NdxPermSamples(NdxStep)-1,NumSamples)+1;    
    MySample=Samples(:,NdxSample);
    if NdxStep<0.5*NumSteps   
        % Ordering phase: linear decay
        LearningRate=0.4*(1-NdxStep/NumSteps);
    else
        % Convergence phase: constant
        LearningRate=0.01;        
    end
    
    % Euclidean distance
    RepMySample=repmat(MySample,1,NumNeuro);        
    MyDistances=sqrt(sum((RepMySample-Model.Prototypes(:,:)).^2,1));            			    
    Distances(:,NdxSample) = MyDistances';        
   
    % Update the neurons
    [~,NdxWinner]=min(MyDistances);
    Model.Winners(NdxSample)=NdxWinner;
    Coef=repmat(LearningRate,Dimension,1);
    Model.Prototypes(:,NdxWinner)=Coef.*MySample+...
        (1-Coef).*Model.Prototypes(:,NdxWinner);    
end

% Store the means
Model.Means = Model.Prototypes;
Model.Samples = Samples;


    
    
        
