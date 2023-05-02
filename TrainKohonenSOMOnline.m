function [Model, quantization_errors, consensus_quantization_errors, learning_rates]=TrainKohonenSOMOnline(Samples,NumNeuro,NumSteps, HistorySize, FunctionType, FunctionConstant, MaxLearningRate, MinLearningRate)
% Train a Self-Organizing Map (SOM)
% Inputs:
%	Samples=Training samples (one sample per column)
%	NumSteps=Number of training steps
% Output:
%   Model=Trained SOM model

[Dimension,NumSamples]=size(Samples);
Model.NumRowsMap=NumNeuro(1);
NumRowsMap=NumNeuro(1);
Model.NumColsMap=NumNeuro(2);
NumColsMap=NumNeuro(2);
NumNeuro=NumRowsMap*NumColsMap;  

quantization_error_history    = zeros(HistorySize, 1);
assigned_neuron_history       = zeros(HistorySize, 1);
running_total_by_neuron       = zeros(NumNeuro, 1);
num_asigned_by_neuron         = zeros(NumNeuro, 1);
history_ind                   = 1;
quantization_errors           = zeros(NumSteps, 1);
consensus_quantization_errors = zeros(NumSteps, 1);
learning_rates                = zeros(NumSteps, 1);


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
%Distances=zeros(NumNeuro,NumSamples);
MaxRadius=(NumRowsMap+NumColsMap)/8;
%NdxPermSamples = randperm(NumSteps);

for NdxStep=1:NumSteps    
    %NdxSample = ceil(NumSamples*rand(1));
    %NdxSample = mod(NdxPermSamples(NdxStep)-1,NumSamples)+1;    
    NdxSample = mod(NdxStep-1,NumSamples)+1;    
    MySample=Samples(:,NdxSample);

    %if NdxStep<0.5*NumSteps   
    %    % Ordering phase: linear decay
    %    LearningRate=0.4*(1-NdxStep/NumSteps);
    %    MyRadius=MaxRadius*(1-(NdxStep-1)/NumSteps);
    %else
    %    % Convergence phase: constant
    %    LearningRate=0.01;
    %    MyRadius=0.1;
    %end
    
    % Euclidean distance (standard Kohonen's SOFM)
    RepMySample=repmat(MySample,1,NumNeuro);        
    MyDistances=sqrt(sum((RepMySample-Model.Prototypes(:,:)).^2,1));            			    
    %Distances(:,NdxSample) = MyDistances';        
   
    % Update the neurons
    [~,NdxWinner]=min(MyDistances);
    Model.Winners(NdxSample)=NdxWinner;

    %dstvec = MySample-Model.Prototypes(:,NdxWinner);
    %dst = sqrt(sum(dstvec.*dstvec));
    %quantization_error_history(history_ind) = dst;
    if NdxStep>HistorySize
      running_total_by_neuron(assigned_neuron_history(history_ind)) = running_total_by_neuron(assigned_neuron_history(history_ind)) - quantization_error_history(history_ind);
      num_asigned_by_neuron(  assigned_neuron_history(history_ind)) = num_asigned_by_neuron(  assigned_neuron_history(history_ind)) - 1;
    end
    quantization_error_history(history_ind) = MyDistances(NdxWinner);
    assigned_neuron_history(history_ind)    = NdxWinner;
    running_total_by_neuron(NdxWinner)      = running_total_by_neuron(assigned_neuron_history(history_ind)) + MyDistances(NdxWinner);
    num_asigned_by_neuron(NdxWinner)        = num_asigned_by_neuron(NdxWinner) + 1;
    history_ind                             = mod(history_ind, HistorySize)+1;
    %if NdxStep<HistorySize
    %  consensus_error = median(quantization_error_history(1:NdxStep));
    %else
    %  consensus_error = median(quantization_error_history);
    %end
    nonzero         = num_asigned_by_neuron~=0;
    means           = running_total_by_neuron(nonzero)./num_asigned_by_neuron(nonzero);
    consensus_error = median(means);
    switch FunctionType
      case 'l'
        LearningRate = (consensus_error)*FunctionConstant;
      case 'q'
        LearningRate = (consensus_error^2)*FunctionConstant;
      case 'c'
        LearningRate = (consensus_error^3)*FunctionConstant;
      case 'e'
        LearningRate = exp(consensus_error)*FunctionConstant;
      case 'o'
        if NdxStep<0.5*NumSteps   
            % Ordering phase: linear decay
            LearningRate=0.4*(1-NdxStep/NumSteps);
            MyRadius=MaxRadius*(1-(NdxStep-1)/NumSteps);
        else
            % Convergence phase: constant
            LearningRate=0.01;        
            MyRadius=0.1;
        end
    end
    LearningRate = max(MinLearningRate, min(MaxLearningRate, LearningRate));
    if FunctionType~='o'
      MyRadius = LearningRate/MaxLearningRate*MaxRadius;
    end
    quantization_errors(NdxStep)           = MyDistances(NdxWinner);
    consensus_quantization_errors(NdxStep) = consensus_error;
    learning_rates(NdxStep)                = LearningRate;


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

    
    
        
