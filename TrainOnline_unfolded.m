function [varargout]=TrainOnline(trainMode, varargin)

if strcmp(trainMode, 'singleLR')
  [varargout{1:nargout}] = TrainOnlineSingle(   varargin{:});
else
  [varargout{1:nargout}] = TrainOnlineMultiRate(varargin{:});
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Model, quantization_errors, consensus_quantization_errors, learning_rates, meanQuantizationErrors]=TrainOnlineSingle(Samples,netType, params,NumSteps, HistorySize, FunctionType, FunctionConstant, mqeMode, slidingWindow, initial_consensus_error, initial_LearningRate, MaxLearningRate, MinLearningRate, VideoInfo, SamplesForInit, intermediateTimesToRecord, SamplesCell)
% Train a Competitive Neural Network
% Inputs:
%	Samples=Training samples (one sample per column)
%	NumSteps=Number of training steps
% Output:
%   Model=Trained Competitive model

isCompetitive = strcmp(netType, 'competitive');
isSOM         = strcmp(netType, 'som');
isGNG         = strcmp(netType, 'gng');
meanQuantizationErrors = {};

doVideo = numel(VideoInfo)>0;
if doVideo
  if iscell(VideoInfo)
    videoname     = VideoInfo{1};
    fps           = VideoInfo{2};
    skipframes    = VideoInfo{3};
    sample_limits = VideoInfo{4};
    vw = mp4_video(videoname, fps);
    h = figure;
    set(h, 'Position', get(0, 'Screensize'));
  else
    Model.all_prototypes = cell(size(Samples, 2), 1);
  end
end

if isCompetitive
  NumNeuro = params.NumNeurons;
  [Dimension,NumSamples]=size(Samples);

  % Prototypes initialization
  NumPatIni=max([Dimension+1,ceil(NumSamples/NumNeuro)]);
  Model.Prototypes=zeros(Dimension,NumNeuro);
  for NdxNeuro=1:NumNeuro   
      %MySamples=Samples(:,ceil(NumSamples*rand(1,NumPatIni)));
      MySamples=Samples(:,ceil(size(SamplesForInit, 2)*rand(1,NumPatIni)));
      Model.Prototypes(:,NdxNeuro)=mean(MySamples,2);             
  end
elseif isSOM
  NumNeuro = params.NumNeurons;
  [Dimension,NumSamples]=size(Samples);
  Model.NumRowsMap=NumNeuro(1);
  NumRowsMap=NumNeuro(1);
  Model.NumColsMap=NumNeuro(2);
  NumColsMap=NumNeuro(2);
  NumNeuro=NumRowsMap*NumColsMap;  

  % Prototypes initialization
  NumPatIni=max([Dimension+1,ceil(NumSamples/NumNeuro)]);
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

  MaxRadius=(NumRowsMap+NumColsMap)/8;
elseif isGNG
  NumNeuro = params.MaxUnits;
  MaxUnits=params.MaxUnits;
  Lambda=params.Lambda;
  Alpha=params.Alpha;
  AMax=params.AMax;
  D=params.D;
  ratio_epsilon_N_B = params.ratio_epsilon_N_B; %0.03:
  Model.MaxUnits=MaxUnits;
  Model.Lambda=Lambda;
  Model.Alpha=Alpha;
  Model.AMax=AMax;
  Model.D=D;
  Model.NumSteps=NumSteps;
  Model.Winners=zeros(1,size(Samples,2));


  [Dimension,NumSamples]=size(Samples);

  % Prototype vectors
  Model.Prototypes=nan*ones(Dimension,MaxUnits);

  % Accumulated errors
  Model.Errors=zeros(1,MaxUnits);

  % Matrix of connections. An absent connection is represented by a zero value.
  % A present connection is represented by a positive value which is the
  % number of steps remaining until connection removal, i.e. its time to live
  Model.Connections=sparse(MaxUnits,MaxUnits);

  % Initialization (two units and a connection between them)
  Model.Prototypes(:,1:2)=Samples(:,ceil(rand(2,1)*NumSamples));
  Model.Connections(1,2)=AMax;
  Model.Connections(2,1)=AMax;
end


quantization_error_history    = zeros(HistorySize, 1);
assigned_neuron_history       = zeros(HistorySize, 1);
running_total_by_neuron       = zeros(NumNeuro, 1);
num_asigned_by_neuron         = zeros(NumNeuro, 1);
history_ind                   = 1;
quantization_errors           = zeros(NumSteps, 1);
consensus_quantization_errors = zeros(NumSteps, 1);
learning_rates                = zeros(NumSteps, 1);

if FunctionType=='e'
  FunctionConstant = FunctionConstant/exp(1);%*0.4/exp(0.5);
end

consensus_error = initial_consensus_error;
LearningRate = initial_LearningRate;

% Training
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
    if isGNG
    % Determine the first (S1) and second (S2) best matching units
      otherDistances = MyDistances;
      otherDistances(NdxWinner) = nan;
      S1 = NdxWinner;
      [~,S2]=min(otherDistances);
    end

    if NdxStep>HistorySize
      running_total_by_neuron(assigned_neuron_history(history_ind)) = running_total_by_neuron(assigned_neuron_history(history_ind)) - quantization_error_history(history_ind);
      num_asigned_by_neuron(  assigned_neuron_history(history_ind)) = num_asigned_by_neuron(  assigned_neuron_history(history_ind)) - 1;
    end
    %dstvec = MySample-Model.Prototypes(:,NdxWinner);
    %dst = sqrt(sum(dstvec.*dstvec));
    %quantization_error_history(history_ind) = dst;
    quantization_error_history(history_ind) = MyDistances(NdxWinner);
    assigned_neuron_history(history_ind)    = NdxWinner;
    running_total_by_neuron(NdxWinner)      = running_total_by_neuron(NdxWinner) + MyDistances(NdxWinner);
    num_asigned_by_neuron(NdxWinner)        =   num_asigned_by_neuron(NdxWinner) + 1;

    if slidingWindow || history_ind==HistorySize
      switch mqeMode
        % median of all errors (without binning by assigned neuron) from last N steps
        case 0
          if NdxStep<HistorySize
            consensus_error = median(quantization_error_history(1:NdxStep));
          else
            consensus_error = median(quantization_error_history);
          end
        % mean   of all errors (without binning by assigned neuron) from last N steps
        case 1
          if NdxStep<HistorySize
            consensus_error = mean(quantization_error_history(1:NdxStep));
          else
            consensus_error = mean(quantization_error_history);
          end
        % bin errors by neuron, discard neurons with 0 errors, compute mean error by neuron, compute median of these means
        case 2
          nonzero         = num_asigned_by_neuron~=0;
          values          = running_total_by_neuron(nonzero)./num_asigned_by_neuron(nonzero);
          consensus_error = median(values);
        % bin errors by neuron, discard neurons with 0 errors, compute mean error by neuron, compute mean   of these means
        case 3
          nonzero         = num_asigned_by_neuron~=0;
          values          = running_total_by_neuron(nonzero)./num_asigned_by_neuron(nonzero);
          consensus_error = mean(values);
        % bin errors by neuron, discard neurons with 0 errors, add up errors by neuron, compute median of these summations
        case 4
          nonzero         = num_asigned_by_neuron~=0;
          values          = running_total_by_neuron(nonzero);
          consensus_error = median(values);
        % bin errors by neuron, discard neurons with 0 errors, add up errors by neuron, compute mean   of these summations
        case 5
          nonzero         = num_asigned_by_neuron~=0;
          values          = running_total_by_neuron(nonzero);
          consensus_error = mean(values);
        % bin errors by neuron, add up errors by neuron, compute median of these summations
        case 6
          values          = running_total_by_neuron;
          consensus_error = median(values);
        % bin errors by neuron, add up errors by neuron, compute mean   of these summations
        case 7
          values          = running_total_by_neuron;
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
          LearningRate = FunctionConstant;
        case 6
          LearningRate = FunctionConstant/consensus_error;
        case 111
          if NdxStep<0.5*NumSteps   
              % Ordering phase: linear decay
              LearningRate=0.4*(1-NdxStep/NumSteps);
              if isSOM
                MyRadius=MaxRadius*(1-(NdxStep-1)/NumSteps);
              end
          else
              % Convergence phase: constant
              LearningRate=0.01;        
              if isSOM
                MyRadius=0.1;
              end
          end
      end
      LearningRate = max(MinLearningRate, min(MaxLearningRate, LearningRate));
      if isSOM% && FunctionType~=5
        MyRadius = LearningRate/MaxLearningRate*MaxRadius;
      end
      if isGNG
        EpsilonB = LearningRate;
        EpsilonN = LearningRate*ratio_epsilon_N_B;
      end
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

    if isCompetitive
      Coef=repmat(LearningRate,Dimension,1);
      Model.Prototypes(:,NdxWinner)= Coef.*MySample+(1-Coef).*Model.Prototypes(:,NdxWinner);
    elseif isSOM
      Coef=repmat(LearningRate*exp(-TopolDist{NdxWinner}/(MyRadius^2)),Dimension,1);
      Model.Prototypes(:,:)        = Coef.*repmat(MySample,1,NumNeuro)+(1-Coef).*Model.Prototypes(:,:);
    elseif isGNG
      % Decrease the time to live of all edges emanating from S1
      Model.Connections(S1,:)=max(0,Model.Connections(S1,:)-1);
      Model.Connections(:,S1)=max(0,Model.Connections(:,S1)-1);
      % Add the squared distance of S1 to the input sample to the error
      % counter of S1
      Model.Errors(S1)=Model.Errors(S1)+MyDistances(NdxWinner);
      % Move S1 and its topological neighbors towards the input sample
      Model.Prototypes(:,S1)=(1-EpsilonB)*Model.Prototypes(:,S1)+EpsilonB*MySample;
      Neighbors=find(Model.Connections(S1,:));
      Model.Prototypes(:,Neighbors)=(1-EpsilonN)*Model.Prototypes(:,Neighbors)+...
          EpsilonN*repmat(MySample,1,numel(Neighbors));
      % Create or refresh the connection between S1 and S2
      Model.Connections(S1,S2)=AMax;
      Model.Connections(S2,S1)=AMax;
      % Remove units with no emanating edges
      NdxNoEdges=find(sum(Model.Connections>0,1)==0);
      Model.Prototypes(:,NdxNoEdges)=nan;
  %     fprintf('Antes de la creación de unidades\n');
      % Unit creation
      if mod(NdxStep,Lambda)==0
          % Find the unit with the largest error
          [Maximum NdxMaxError]=max(Model.Errors);
          % Find its neighbor with the largest error
          [Maximum NdxNeighbor]=max(Model.Errors.*(Model.Connections(NdxMaxError,:)>0));
          % Create the new unit, if possible. Otherwise, finish
          NdxNewUnit=find(isnan(Model.Prototypes(1,:)),1,'first');
          if ~isempty(NdxNewUnit)        
              % Set the new prototype vector
              Model.Prototypes(:,NdxNewUnit)=0.5*(Model.Prototypes(:,NdxMaxError)+...
                  Model.Prototypes(:,NdxNeighbor));
              % Remove the connection between the two old units
              Model.Connections(NdxMaxError,NdxNeighbor)=0;
              Model.Connections(NdxNeighbor,NdxMaxError)=0;    
              % Create connections among the new unit and the two old ones
              Model.Connections(NdxNewUnit,[NdxMaxError NdxNeighbor])=AMax;
              Model.Connections([NdxMaxError NdxNeighbor],NdxNewUnit)=AMax;
              % Decrease the errors of the old units and set the error of the new
              % one
              Model.Errors([NdxMaxError NdxNeighbor])=Alpha*...
                  Model.Errors([NdxMaxError NdxNeighbor]);
              Model.Errors(NdxNewUnit)=Model.Errors(NdxMaxError);
          end
      end
      % Decrease all error variables by multiplying them by D
      Model.Errors=D*Model.Errors;
    end

    if doVideo
      if iscell(VideoInfo) && (mod(NdxStep-1, skipframes)==0)
        minidx = sample_limits(find(NdxStep<sample_limits, 1)-1);
        %suptitle(sprintf('%06d/%06d', NdxStep, NumSamples));
        subplot(2,1,1);
        t=1:NdxStep;
        plot(t, quantization_errors(1:NdxStep), 'b', t, consensus_quantization_errors(1:NdxStep), 'r', t, learning_rates(1:NdxStep), 'y');
        xlim([0 size(Samples, 2)]);
        ylim([0 1]);
        subplot(2,1,2);
        plot(Samples(1,minidx:NdxStep), Samples(2,minidx:NdxStep), '*y', Model.Prototypes(1,:),Model.Prototypes(2,:),'or','LineWidth',2,'MarkerFaceColor',[1 0 0],'MarkerSize',7);
        xlim([0 1]);  
        ylim([0 1]);
        axis square;  
        vw.addFrame(h);
      else
        Model.all_prototypes{NdxStep} = Model.Prototypes;
      end
    end
    if any(NdxStep==intermediateTimesToRecord)
      idxIntermediate=find(NdxStep==intermediateTimesToRecord);
      if isfield(Model, 'firstPrototypes')
        meanQuantizationErrors{end+1}   = compute_consensus_quantization_error(mqeMode, Model.Prototypes, SamplesCell{idxIntermediate});
        Model.firstPrototypes{end+1}    = Model.Prototypes;
        if isGNG
          Model.firstConnections{end+1} = Model.Connections;
        end
      else
        meanQuantizationErrors   = {compute_consensus_quantization_error(mqeMode, Model.Prototypes, SamplesCell{idxIntermediate})};
        Model.firstPrototypes    = {Model.Prototypes};
        if isGNG
          Model.firstConnections = {Model.Connections};
        end
      end
    end
end

if isCompetitive
  % Store the means
  Model.Means = Model.Prototypes;
  Model.Samples = Samples;
elseif isSOM
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
elseif isGNG
  Model.Means = Model.Prototypes;
  Model.Samples = Samples;
end

if doVideo && iscell(VideoInfo)
  vw.save();
  close(h);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Model, quantization_errors, consensus_quantization_errors, learning_rates, meanQuantizationErrors]=TrainOnlineMultiRate(Samples,netType, params,NumSteps, HistorySize, FunctionType, FunctionConstant, mqeMode, slidingWindow, initial_consensus_error, initial_LearningRate, MaxLearningRate, MinLearningRate, VideoInfo, SamplesForInit, intermediateTimesToRecord, SamplesCell)
% Train a Competitive Neural Network
% Inputs:
%	Samples=Training samples (one sample per column)
%	NumSteps=Number of training steps
% Output:
%   Model=Trained Competitive model

isCompetitive = strcmp(netType, 'competitive');
isSOM         = strcmp(netType, 'som');
isGNG         = strcmp(netType, 'gng');
meanQuantizationErrors = {};

doVideo = numel(VideoInfo)>0;
if doVideo
  if iscell(VideoInfo)
    videoname     = VideoInfo{1};
    fps           = VideoInfo{2};
    skipframes    = VideoInfo{3};
    sample_limits = VideoInfo{4};
    vw = mp4_video(videoname, fps);
    h = figure;
    set(h, 'Position', get(0, 'Screensize'));
  else
    Model.all_prototypes = cell(size(Samples, 2), 1);
  end
end

if isCompetitive
  NumNeuro = params.NumNeurons;
  [Dimension,NumSamples]=size(Samples);

  % Prototypes initialization
  NumPatIni=max([Dimension+1,ceil(NumSamples/NumNeuro)]);
  Model.Prototypes=zeros(Dimension,NumNeuro);
  for NdxNeuro=1:NumNeuro   
      %MySamples=Samples(:,ceil(NumSamples*rand(1,NumPatIni)));
      MySamples=Samples(:,ceil(size(SamplesForInit, 2)*rand(1,NumPatIni)));
      Model.Prototypes(:,NdxNeuro)=mean(MySamples,2);             
  end
elseif isSOM
  NumNeuro = params.NumNeurons;
  [Dimension,NumSamples]=size(Samples);
  Model.NumRowsMap=NumNeuro(1);
  NumRowsMap=NumNeuro(1);
  Model.NumColsMap=NumNeuro(2);
  NumColsMap=NumNeuro(2);
  NumNeuro=NumRowsMap*NumColsMap;  

  % Prototypes initialization
  NumPatIni=max([Dimension+1,ceil(NumSamples/NumNeuro)]);
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

  MaxRadius=(NumRowsMap+NumColsMap)/8;
elseif isGNG
  NumNeuro = params.MaxUnits;
  MaxUnits=params.MaxUnits;
  Lambda=params.Lambda;
  Alpha=params.Alpha;
  AMax=params.AMax;
  D=params.D;
  ratio_epsilon_N_B = params.ratio_epsilon_N_B; %0.03:
  Model.MaxUnits=MaxUnits;
  Model.Lambda=Lambda;
  Model.Alpha=Alpha;
  Model.AMax=AMax;
  Model.D=D;
  Model.NumSteps=NumSteps;
  Model.Winners=zeros(1,size(Samples,2));


  [Dimension,NumSamples]=size(Samples);

  % Prototype vectors
  Model.Prototypes=nan*ones(Dimension,MaxUnits);

  % Accumulated errors
  Model.Errors=zeros(1,MaxUnits);

  % Matrix of connections. An absent connection is represented by a zero value.
  % A present connection is represented by a positive value which is the
  % number of steps remaining until connection removal, i.e. its time to live
  Model.Connections=sparse(MaxUnits,MaxUnits);

  % Initialization (two units and a connection between them)
  Model.Prototypes(:,1:2)=Samples(:,ceil(rand(2,1)*NumSamples));
  Model.Connections(1,2)=AMax;
  Model.Connections(2,1)=AMax;
end


history_size_per_neuron       = zeros(NumNeuro, 1);
quantization_error_history    = zeros(HistorySize, NumNeuro);
history_ind                   = ones(NumNeuro, 1);;
quantization_errors           = zeros(NumSteps, 1);
consensus_quantization_errors = zeros(NumSteps, 1);
learning_rates                = zeros(NumSteps, 1);

if FunctionType=='e'
  FunctionConstant = FunctionConstant/exp(1);%*0.4/exp(0.5);
end

consensus_error = zeros(NumNeuro,1)+initial_consensus_error;
LearningRate    = zeros(NumNeuro,1)+initial_LearningRate;

% Training
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
    if isGNG
    % Determine the first (S1) and second (S2) best matching units
      otherDistances = MyDistances;
      otherDistances(NdxWinner) = nan;
      S1 = NdxWinner;
      [~,S2]=min(otherDistances);
    end

    quantization_error_history(history_ind(NdxWinner),NdxWinner) = MyDistances(NdxWinner);
    history_size_per_neuron(NdxWinner)                           = history_size_per_neuron(NdxWinner)+1;

    if slidingWindow || history_ind(NdxWinner)==HistorySize
      switch mqeMode
        % median of all errors (without binning by assigned neuron) from last N steps
        case 0
          if history_size_per_neuron(NdxWinner)<HistorySize
            consensus_error(NdxWinner) = median(quantization_error_history(1:history_size_per_neuron(NdxWinner),NdxWinner));
          else
            consensus_error(NdxWinner) = median(quantization_error_history(:,NdxWinner));
          end
        % mean   of all errors (without binning by assigned neuron) from last N steps
        case 1
          if history_size_per_neuron(NdxWinner)<HistorySize
            consensus_error(NdxWinner) = mean(quantization_error_history(1:history_size_per_neuron(NdxWinner),NdxWinner));
          else
            consensus_error(NdxWinner) = mean(quantization_error_history(:,NdxWinner));
          end
        otherwise
          error(sprintf('mqeMODE %d does not make sense when having a learning rate per neuron!', mqeMode));
      end    
      
      switch FunctionType
        case 1
          LearningRate(NdxWinner) = (consensus_error(NdxWinner))*FunctionConstant;
        case 2
          LearningRate(NdxWinner) = (consensus_error(NdxWinner).^2)*FunctionConstant;
        case 3
          LearningRate(NdxWinner) = (consensus_error(NdxWinner).^3)*FunctionConstant;
        case 4
          LearningRate(NdxWinner) = exp(consensus_error(NdxWinner))*FunctionConstant;
        case 5
          LearningRate(NdxWinner) = FunctionConstant;
        case 6
          LearningRate(NdxWinner) = FunctionConstant/consensus_error(NdxWinner);
        case 111
          if NdxStep<0.5*NumSteps   
              % Ordering phase: linear decay
              LearningRate(NdxWinner)=0.4*(1-NdxStep/NumSteps);
              if isSOM
                MyRadius=MaxRadius*(1-(NdxStep-1)/NumSteps);
              end
          else
              % Convergence phase: constant
              LearningRate(NdxWinner)=0.01;        
              if isSOM
                MyRadius=0.1;
              end
          end
      end
      LearningRate(NdxWinner) = max(MinLearningRate, min(MaxLearningRate, LearningRate(NdxWinner)));
      if isSOM% && FunctionType~=5
        MyRadius = LearningRate(NdxWinner)/MaxLearningRate*MaxRadius;
      end
      if isGNG
        EpsilonB = LearningRate(NdxWinner);
        EpsilonN = LearningRate(NdxWinner)*ratio_epsilon_N_B;
      end
    end

    history_ind(NdxWinner)                 = mod(history_ind(NdxWinner), HistorySize)+1;

    quantization_errors(NdxStep)           = MyDistances(NdxWinner);
    consensus_quantization_errors(NdxStep) = consensus_error(NdxWinner);
    learning_rates(NdxStep)                = LearningRate(NdxWinner);

    if isCompetitive
      Coef=repmat(LearningRate(NdxWinner),Dimension,1);
      Model.Prototypes(:,NdxWinner)= Coef.*MySample+(1-Coef).*Model.Prototypes(:,NdxWinner);
    elseif isSOM
      Coef=repmat(LearningRate(NdxWinner)*exp(-TopolDist{NdxWinner}/(MyRadius^2)),Dimension,1);
      Model.Prototypes(:,:)        = Coef.*repmat(MySample,1,NumNeuro)+(1-Coef).*Model.Prototypes(:,:);
    elseif isGNG
      % Decrease the time to live of all edges emanating from S1
      Model.Connections(S1,:)=max(0,Model.Connections(S1,:)-1);
      Model.Connections(:,S1)=max(0,Model.Connections(:,S1)-1);
      % Add the squared distance of S1 to the input sample to the error
      % counter of S1
      Model.Errors(S1)=Model.Errors(S1)+MyDistances(NdxWinner);
      % Move S1 and its topological neighbors towards the input sample
      Model.Prototypes(:,S1)=(1-EpsilonB)*Model.Prototypes(:,S1)+EpsilonB*MySample;
      Neighbors=find(Model.Connections(S1,:));
      Model.Prototypes(:,Neighbors)=(1-EpsilonN)*Model.Prototypes(:,Neighbors)+...
          EpsilonN*repmat(MySample,1,numel(Neighbors));
      % Create or refresh the connection between S1 and S2
      Model.Connections(S1,S2)=AMax;
      Model.Connections(S2,S1)=AMax;
      % Remove units with no emanating edges
      NdxNoEdges=find(sum(Model.Connections>0,1)==0);
      Model.Prototypes(:,NdxNoEdges)=nan;
      history_ind(NdxNoEdges) = 1;
      history_size_per_neuron(NdxNoEdges) = 0;
  %     fprintf('Antes de la creación de unidades\n');
      % Unit creation
      if mod(NdxStep,Lambda)==0
          % Find the unit with the largest error
          [Maximum NdxMaxError]=max(Model.Errors);
          % Find its neighbor with the largest error
          [Maximum NdxNeighbor]=max(Model.Errors.*(Model.Connections(NdxMaxError,:)>0));
          % Create the new unit, if possible. Otherwise, finish
          NdxNewUnit=find(isnan(Model.Prototypes(1,:)),1,'first');
          if ~isempty(NdxNewUnit)        
              % Set the new prototype vector
              Model.Prototypes(:,NdxNewUnit)=0.5*(Model.Prototypes(:,NdxMaxError)+...
                  Model.Prototypes(:,NdxNeighbor));
              % Remove the connection between the two old units
              Model.Connections(NdxMaxError,NdxNeighbor)=0;
              Model.Connections(NdxNeighbor,NdxMaxError)=0;    
              % Create connections among the new unit and the two old ones
              Model.Connections(NdxNewUnit,[NdxMaxError NdxNeighbor])=AMax;
              Model.Connections([NdxMaxError NdxNeighbor],NdxNewUnit)=AMax;
              % Decrease the errors of the old units and set the error of the new
              % one
              Model.Errors([NdxMaxError NdxNeighbor])=Alpha*...
                  Model.Errors([NdxMaxError NdxNeighbor]);
              Model.Errors(NdxNewUnit)=Model.Errors(NdxMaxError);
          end
      end
      % Decrease all error variables by multiplying them by D
      Model.Errors=D*Model.Errors;
    end

    if doVideo
      if iscell(VideoInfo) && (mod(NdxStep-1, skipframes)==0)
        minidx = sample_limits(find(NdxStep<sample_limits, 1)-1);
        %suptitle(sprintf('%06d/%06d', NdxStep, NumSamples));
        subplot(2,1,1);
        t=1:NdxStep;
        plot(t, quantization_errors(1:NdxStep), 'b', t, consensus_quantization_errors(1:NdxStep), 'r', t, learning_rates(1:NdxStep), 'y');
        xlim([0 size(Samples, 2)]);
        ylim([0 1]);
        subplot(2,1,2);
        plot(Samples(1,minidx:NdxStep), Samples(2,minidx:NdxStep), '*y', Model.Prototypes(1,:),Model.Prototypes(2,:),'or','LineWidth',2,'MarkerFaceColor',[1 0 0],'MarkerSize',7);
        xlim([0 1]);  
        ylim([0 1]);
        axis square;  
        vw.addFrame(h);
      else
        Model.all_prototypes{NdxStep} = Model.Prototypes;
      end
    end
    if any(NdxStep==intermediateTimesToRecord)
      idxIntermediate=find(NdxStep==intermediateTimesToRecord);
      if isfield(Model, 'firstPrototypes')
        meanQuantizationErrors{end+1}   = compute_consensus_quantization_error(mqeMode, Model.Prototypes, SamplesCell{idxIntermediate});
        Model.firstPrototypes{end+1}    = Model.Prototypes;
        if isGNG
          Model.firstConnections{end+1} = Model.Connections;
        end
      else
        meanQuantizationErrors   = {compute_consensus_quantization_error(mqeMode, Model.Prototypes, SamplesCell{idxIntermediate})};
        Model.firstPrototypes    = {Model.Prototypes};
        if isGNG
          Model.firstConnections = {Model.Connections};
        end
      end
    end
end

if isCompetitive
  % Store the means
  Model.Means = Model.Prototypes;
  Model.Samples = Samples;
elseif isSOM
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
elseif isGNG
  Model.Means = Model.Prototypes;
  Model.Samples = Samples;
end

if doVideo && iscell(VideoInfo)
  vw.save();
  close(h);
end
