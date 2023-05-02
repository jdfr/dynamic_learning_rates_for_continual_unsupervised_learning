function [Model, quantization_errors, consensus_quantization_errors, learning_rates]=TrainCompetitiveOnline(Samples,NumNeuro,NumSteps, HistorySize, FunctionType, FunctionConstant, mqe, MaxLearningRate, MinLearningRate)
% Train a Competitive Neural Network
% Inputs:
%	Samples=Training samples (one sample per column)
%	NumSteps=Number of training steps
% Output:
%   Model=Trained Competitive model

[Dimension,NumSamples]=size(Samples);

quantization_error_history    = zeros(HistorySize, 1);
assigned_neuron_history       = zeros(HistorySize, 1);
running_total_by_neuron       = zeros(NumNeuro, 1);
num_asigned_by_neuron         = zeros(NumNeuro, 1);
history_ind                   = 1;
quantization_errors           = zeros(NumSteps, 1);
consensus_quantization_errors = zeros(NumSteps, 1);
learning_rates                = zeros(NumSteps, 1);

% Prototypes initialization
NumPatIni=max([Dimension+1,ceil(NumSamples/NumNeuro)]);
Model.Prototypes=zeros(Dimension,NumNeuro);
for NdxNeuro=1:NumNeuro   
    MySamples=Samples(:,ceil(NumSamples*rand(1,NumPatIni)));
    Model.Prototypes(:,NdxNeuro)=mean(MySamples,2);             
end        

% Training
%Distances=zeros(NumNeuro,NumSamples);
%NdxPermSamples = randperm(NumSteps);

if FunctionType=='e'
  FunctionConstant = FunctionConstant*0.4/exp(0.5);
end

consensus_error = initial_consensus_error;
LearningRate = initial_LearningRate;

for NdxStep=1:NumSteps
    %NdxSample = ceil(NumSamples*rand(1));
    %NdxSample = mod(NdxPermSamples(NdxStep)-1,NumSamples)+1;    
    NdxSample = mod(NdxStep-1,NumSamples)+1;    
    MySample=Samples(:,NdxSample);
    
    % Euclidean distance
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
    running_total_by_neuron(NdxWinner)      = running_total_by_neuron(NdxWinner) + MyDistances(NdxWinner);
    num_asigned_by_neuron(NdxWinner)        =   num_asigned_by_neuron(NdxWinner) + 1;

    if slidingWindow || history_ind==HistorySize
      switch mqeMode
        # with sliding window: median of all errors (without binning by assigned neuron) from last N steps
        case 1
          if NdxStep<HistorySize
            consensus_error = median(quantization_error_history(1:NdxStep));
          else
            consensus_error = median(quantization_error_history);
          end
        # with sliding window: mean   of all errors (without binning by assigned neuron) from last N steps
        case 2
          if NdxStep<HistorySize
            consensus_error = mean(quantization_error_history(1:NdxStep));
          else
            consensus_error = mean(quantization_error_history);
          end
        # with sliding window: bin errors by neuron, discard neurons with 0 errors, compute mean error by neuron, compute median of these means
        case 3
          nonzero         = num_asigned_by_neuron~=0;
          values          = running_total_by_neuron(nonzero))./num_asigned_by_neuron(nonzero);
          consensus_error = median(values);
        # with sliding window: bin errors by neuron, discard neurons with 0 errors, compute mean error by neuron, compute mean   of these means
        case 4
          nonzero         = num_asigned_by_neuron~=0;
          values          = running_total_by_neuron(nonzero))./num_asigned_by_neuron(nonzero);
          consensus_error = median(values);
        # with sliding window: bin errors by neuron, compute mean error by neuron, compute median of these means
        case 5
          values          = running_total_by_neuron./num_asigned_by_neuron;
          consensus_error = median(values);
        # with sliding window: bin errors by neuron, compute mean error by neuron, compute mean   of these means
        case 6
          values          = running_total_by_neuron./num_asigned_by_neuron;
          consensus_error = mean(values);
      end    
      
      switch FunctionType
        case 1
          LearningRate = (consensus_error)*FunctionConstant;
        case 2
          LearningRate = (consensus_error.^2)*FunctionConstant;
        case 3
          LearningRate = (consensus_error.^3)*FunctionConstant;
        case 4
          LearningRate = exp(consensus_error)*FunctionConstant;
        case 5
          if NdxStep<0.5*NumSteps   
              % Ordering phase: linear decay
              LearningRate=0.4*(1-NdxStep/NumSteps);
          else
              % Convergence phase: constant
              LearningRate=0.01;        
          end
      end
      LearningRate = max(MinLearningRate, min(MaxLearningRate, LearningRate));
    end

    history_ind                             = mod(history_ind, HistorySize)+1;
    %if NdxStep<HistorySize
    %  consensus_error = median(quantization_error_history(1:NdxStep));
    %else
    %  running_total_by_neuron(assigned_neuron_history(history_ind_prev)) = running_total_by_neuron(assigned_neuron_history(history_ind_prev)) + quantization_error_history(history_ind_prev);
    %  consensus_error = median(quantization_error_history);
    %end

    quantization_errors(NdxStep)           = MyDistances(NdxWinner);
    consensus_quantization_errors(NdxStep) = consensus_error;
    learning_rates(NdxStep)                = LearningRate;

    Coef=repmat(LearningRate,Dimension,1);
    Model.Prototypes(:,NdxWinner)=Coef.*MySample+...
        (1-Coef).*Model.Prototypes(:,NdxWinner);    
end

% Store the means
Model.Means = Model.Prototypes;
Model.Samples = Samples;

