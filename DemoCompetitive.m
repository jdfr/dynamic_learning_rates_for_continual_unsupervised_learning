function DemoCompetitive(prefixpath, FileNames, do_load, do_show, trainingType, netType, mqeMode, slidingWindow, initial_LearningRate, MaxLearningRate, NumSamplesPerImage);

if ~(strcmp(netType, 'competitive') || strcmp(netType, 'som') || strcmp(netType, 'gng'))
  error(sprintf('netType must be som, competitive or gng'));
end
isCompetitive = strcmp(netType, 'competitive');
isSOM         = strcmp(netType, 'som');
isGNG         = strcmp(netType, 'gng');

trainMode = 'singleLR';
if ~(strcmp(trainMode, 'singleLR') || strcmp(trainMode, 'multiLR'))
  error(sprintf('trainMode must be singleLR, or multiLR'));
end

%FileName='Irregular.bmp';
%NumSamples = 20000;
NumEpochs = 2;
NumStepsPerImage = NumSamplesPerImage; %NumSamples*NumEpochs;
if isCompetitive
  params.NumNeurons = 25;
elseif isSOM
  NumRowsMap = 5;
  NumColsMap = 5;
  params.NumNeurons = [NumRowsMap NumColsMap];
elseif isGNG
  params.MaxUnits = 25; %50;
  params.Lambda=100;
  params.Alpha=0.5;
  params.AMax=50;
  params.D=0.995;
  params.ratio_epsilon_N_B = 0.03; %so that if episilonB is 0.2, epsilonN=epsilonB*ratio=0.006
end
HistorySizes  = [10, 50, 100]; %10:10:100;
%LinearFactors = 0.1:0.1:1.2;
%QuadraticFactors = 1:10:70;%0.2:0.2:2.4;
%CubicFactors = 1:50:700;%0.4:0.4:4.8;
%ExponentialFactors = 0.1:0.1:1.2;
%Factors = CubicFactors; %ExponentialFactors; %CubicFactors; %QuadraticFactors; %LinearFactors;
%MaxLearningRate = 1;%0.4;
MinLearningRate = 0;%0.001;

[Samples, NumSteps, FusedFileNames] = generateSamplesMultiImgs(prefixpath, FileNames, NumSamplesPerImage);
UniqueSamples = cell(size(Samples));
for k=1:numel(Samples)
  UniqueSamples{k} = unique(Samples{k}', 'rows')';
end
%Samples = GenerateSamplesImg(FileName,NumSamples);

if true && strcmp(trainingType, 'online')
  a=0;
  a=a+1;
  experiments{a}.type    = 1;
  experiments{a}.name    = 'linear';
  experiments{a}.eq      = 'linear';
  experiments{a}.factors = logspace(-3, 3, 300);%[0.1:0.1:1.2 1.5:0.5:10];

  a=a+1;
  experiments{a}.type    = 2;
  experiments{a}.name    = 'quadratic';
  experiments{a}.eq      = 'quadratic';
  experiments{a}.factors = logspace(-3, 3, 300);%[0.2:0.2:2.4 3:10 20:10:70];

if false
  a=a+1;
  experiments{a}.type    = 3;
  experiments{a}.name    = 'cubic';
  experiments{a}.eq      = 'cubic';
  experiments{a}.factors = logspace(-1, 3, 30);%[0.4:0.4:4.8 6:10 100:100:700];
end

  a=a+1;
  experiments{a}.type    = 6;
  experiments{a}.name    = 'inverse';
  experiments{a}.eq      = 'inverse';
  experiments{a}.factors = logspace(-3, 3, 300);%[0.4:0.4:4.8 6:10 100:100:700];

%  a=a+1;
%  experiments{a}.type    = 4;
%  experiments{a}.name    = 'exponential';
%  experiments{a}.eq      = 'LR=min(F*exp(MQE)*0.4/exp(0.5),0.4)';
%  experiments{a}.factors = logspace(-1, 1, 10);%[0.1:0.1:1.2 1.5:0.5:10];

  a=a+1;
  experiments{a}.type    = 5;
  experiments{a}.name    = 'constant';
  experiments{a}.eq      = 'constant';
  experiments{a}.factors = logspace(-3, 0, 100);%[0.1:0.1:1.2 1.5:0.5:10];

  for a=1:length(experiments)
    experiments{a}.FileNames          = FileNames;
    experiments{a}.NumSamplesPerImage = NumSamplesPerImage;
    experiments{a}.NumStepsPerImage   = NumStepsPerImage;
    experiments{a}.params             = params;
    experiments{a}.HistorySizes       = HistorySizes;
    experiments{a}.MaxLearningRate    = MaxLearningRate;
    experiments{a}.MinLearningRate    = MinLearningRate;
    experiments{a}.mqeMode            = mqeMode;
    experiments{a}.slidingWindow      = slidingWindow;
    experiments{a}.initial_LearningRate=initial_LearningRate;
    if strcmp(netType, 'som')
      experiments{a}.NumRowsMap      = NumRowsMap;
      experiments{a}.NumColsMap      = NumColsMap;
    end
    experiments{a}.netType           = netType;
    experiments{a}.meanQuantizationErrors1 = zeros(length(experiments{a}.HistorySizes), length(experiments{a}.factors));
    experiments{a}.meanQuantizationErrors = zeros(length(experiments{a}.HistorySizes), length(experiments{a}.factors));
  end

  if do_load
    load([FusedFileNames 'experiments.mat']);
  else
    best_run.quantization_errors           =  cell(length(experiments),1);
    best_run.consensus_quantization_errors =  cell(length(experiments),1);
    best_run.learning_rates                =  cell(length(experiments),1);
    best_run.mqe                           = zeros(length(experiments),1)+1e6;
    best_run.Model                         =  cell(length(experiments),1);
    best_run.parameters                    = zeros(length(experiments),2);
    for a=1:length(experiments)
      for i=1:size(experiments{a}.meanQuantizationErrors,1)
        for j=1:size(experiments{a}.meanQuantizationErrors,2)
          disp(sprintf('%s, experiment %d: T=%d, H=%d, F=%f', FusedFileNames, a, experiments{a}.type, experiments{a}.HistorySizes(i), experiments{a}.factors(j)));
          for k=1:length(FileNames)
            Samples{k} = Samples{k}(:,randperm(experiments{a}.NumSamplesPerImage));
          end
          [Model, quantization_errors, consensus_quantization_errors, learning_rates, meanQuantizationErrors1]=TrainOnline(trainMode, cat(2, Samples{:}),netType, experiments{a}.params,NumSteps, experiments{a}.HistorySizes(i), experiments{a}.type, experiments{a}.factors(j), mqeMode, slidingWindow, 0, initial_LearningRate, experiments{a}.MaxLearningRate, experiments{a}.MinLearningRate, [], Samples{1}, size(Samples{1}, 2), Samples);%UniqueSamples);

          experiments{a}.meanQuantizationErrors1(i,j) = meanQuantizationErrors1{1};
          experiments{a}.meanQuantizationErrors(i,j) = compute_consensus_quantization_error(mqeMode, Model.Prototypes, Samples{end}); %UniqueSamples{end});
          if experiments{a}.meanQuantizationErrors(i,j) < best_run.mqe(a)
            best_run.mqe(a)                           = experiments{a}.meanQuantizationErrors(i,j);
            best_run.quantization_errors{a}           = quantization_errors;
            best_run.consensus_quantization_errors{a} = consensus_quantization_errors;
            best_run.learning_rates{a}                = learning_rates;
            best_run.Model{a}                         = Model;
            best_run.parameters(a,:)                  = [experiments{a}.HistorySizes(i), experiments{a}.factors(j)];
          end
       end
      end 
    end
    save([FusedFileNames 'experiments.mat'], 'experiments', 'best_run');
  end

  if do_show && false
    h = figure('Name', ['experiments for ' FusedFileNames]);% ' (Y axes: global MQE after training)']);
    nrows = 2;%1;
    ncols = 2;%numel(experiments);
    colors = {'k', 'r', 'b'};
    for a=1:numel(experiments)
      subplot(nrows, ncols, a); hold on;
      x = experiments{a}.meanQuantizationErrors1;
      y = experiments{a}.meanQuantizationErrors;
      labs = cell(size(experiments{a}.meanQuantizationErrors));
      for i=1:size(experiments{a}.meanQuantizationErrors,1)
        for j=1:size(experiments{a}.meanQuantizationErrors,2)
          labs{i,j} = sprintf('a=%0.1f', experiments{a}.factors(j));
        end
        plot(x(i,:),y(i,:),['.' colors{i}]);
      end
      text(x(:),y(:),labs(:),'VerticalAlignment','bottom','HorizontalAlignment','right');
      xlabel(['MQE first shape. ' experiments{a}.eq]);
      ylabel('MQE second shape');
      axis equal;
      legend('N=10', 'N=50', 'N=100');
    end
  end

  if do_show
    exp_show = numel(experiments);
    mn = min(experiments{1}.meanQuantizationErrors, [], 'all');
    mx = max(experiments{1}.meanQuantizationErrors, [], 'all');
    mn1 = min(experiments{1}.meanQuantizationErrors1, [], 'all');
    mx1 = max(experiments{1}.meanQuantizationErrors1, [], 'all');
    mn2 = min([min(best_run.quantization_errors{1}, [], 'all'), min(best_run.consensus_quantization_errors{1}, [], 'all'), min(best_run.learning_rates{1}, [], 'all')]);
    mx2 = max([max(best_run.quantization_errors{1}, [], 'all'), max(best_run.consensus_quantization_errors{1}, [], 'all'), max(best_run.learning_rates{1}, [], 'all')]);
    for a=2:exp_show
      mn = min(min(experiments{a}.meanQuantizationErrors, [], 'all'), mn);
      mx = max(max(experiments{a}.meanQuantizationErrors, [], 'all'), mx);
      mn1 = min(min(experiments{a}.meanQuantizationErrors1, [], 'all'), mn1);
      mx1 = max(max(experiments{a}.meanQuantizationErrors1, [], 'all'), mx1);
      mn2 = min(min([min(best_run.quantization_errors{a}, [], 'all'), min(best_run.consensus_quantization_errors{a}, [], 'all'), min(best_run.learning_rates{a}, [], 'all')]), mn2);
      mx2 = max(max([max(best_run.quantization_errors{a}, [], 'all'), max(best_run.consensus_quantization_errors{a}, [], 'all'), max(best_run.learning_rates{a}, [], 'all')]), mx2);
    end
    t = 1:experiments{exp_show}.NumSamplesPerImage*length(FileNames);
    h = figure('Name', ['experiments for ' FusedFileNames]);% ' (Y axes: global MQE after training)']);
    if isCompetitive
      m = 'competitive';
    elseif isSOM
      m = 'SOM';
    elseif isGNG
      m = 'GNG';
    end
    switch mqeMode
    case 0
      mqeMode_s = '(median of all errors)';
    case 1
      mqeMode_s = '(mean of all errors)';
    case 2
      mqeMode_s = '(for neurons with samples, bin errors by neuron, compute mean error by neuron, compute median of means)';
    case 3
      mqeMode_s = '(for neurons with samples, bin errors by neuron, compute mean error by neuron, compute mean of means)';
    case 4
      mqeMode_s = '(for neurons with samples, bin errors by neuron, add up errors by neuron, compute median of sums)';
    case 5
      mqeMode_s = '(for neurons with samples, bin errors by neuron, add up errors by neuron, compute mean of sums)';
    case 6
      mqeMode_s = '(for all neurons, bin errors by neuron, add up errors by neuron, compute median of sums)';
    case 7
      mqeMode_s = '(for all neurons, bin errors by neuron, add up errors by neuron, compute mean of sums)';
    end
    %suptitle
    sgtitle(sprintf('model=%s, slidingWindow=%d, MqeMode=%s, initLR=%f', m, slidingWindow, mqeMode_s, initial_LearningRate));
    pidx = 1;
    nrows = exp_show;%numel(experiments);
    alsoFirstMQE=false;
    alsoFirstImage=true;
    ncols = 3;
    if alsoFirstMQE
      ncols=ncols+1;
    end
    if alsoFirstImage
      ncols=ncols+1;
      img2 = imread(FileNames{1});
      if ndims(img2)==3
        img2 = mean(img2, 3);
      end
      img2 = logical(img2);
    end
    img = imread(FileNames{end});
    if ndims(img)==3
      img = mean(img, 3);
    end
    img = logical(img);
    for a=1:nrows
      if alsoFirstMQE
      subplot(nrows, ncols, pidx); pidx=pidx+1;
      plot(experiments{a}.factors, experiments{a}.meanQuantizationErrors1(1,:), 'g', experiments{a}.factors, experiments{a}.meanQuantizationErrors1(2,:), 'b', experiments{a}.factors, experiments{a}.meanQuantizationErrors(3,:), 'k', 'LineWidth',1);
      legend({sprintf('$N=%d$', experiments{a}.HistorySizes(1)), sprintf('$N=%d$', experiments{a}.HistorySizes(2)), sprintf('$N=%d$', experiments{a}.HistorySizes(3))},'Interpreter','latex')
      if mn1<mx1
        ylim([mn1 mx1]);
      end
      end

      subplot(nrows, ncols, pidx); pidx=pidx+1;
      %experiments{a}.factors
      %experiments{a}.meanQuantizationErrors(1,:)
      %experiments{a}.meanQuantizationErrors(2,:)
      %experiments{a}.meanQuantizationErrors(3,:)
      plot(experiments{a}.factors, experiments{a}.meanQuantizationErrors(1,:), 'g', experiments{a}.factors, experiments{a}.meanQuantizationErrors(2,:), 'b', experiments{a}.factors, experiments{a}.meanQuantizationErrors(3,:), 'k', 'LineWidth',1);
      legend({sprintf('$N=%d$', experiments{a}.HistorySizes(1)), sprintf('$N=%d$', experiments{a}.HistorySizes(2)), sprintf('$N=%d$', experiments{a}.HistorySizes(3))},'Interpreter','latex')
      %xlabel(['F in ' experiments{a}.eq]);
      xlabel([experiments{a}.eq sprintf(' scaling factor. Minimal MQE achieved with $N=%d$ and factor %g', best_run.parameters(a,1), best_run.parameters(a,2))],'Interpreter','latex');
      if mn<mx
        ylim([mn mx]);
      end
      %ylabel('global MQE after training');
      ylabel('MQE');
      if alsoFirstImage
        subplot(nrows, ncols, pidx); pidx=pidx+1;
        hold on;
        %imshow(img2);
        colormap([0 0 0; 1 1 1]);
        imagesc([0 1], [1 0], img2);
        %image([0 1], [0 1], img2);
        %plot(Samples(1,:),Samples(2,:),'*y');
        if isCompetitive
          plot(best_run.Model{a}.firstPrototypes{1}(1,:),best_run.Model{a}.firstPrototypes{1}(2,:),'or','LineWidth',2,'MarkerFaceColor',[1 0 0],'MarkerSize',7);
          %PlotCompetitive(best_run.Model{a}.Prototypes);
        elseif isSOM
          PlotSOM(best_run.Model{a}.NumRowsMap, best_run.Model{a}.NumColsMap, best_run.Model{a}.firstPrototypes{1});
        elseif isGNG
          PlotGNG(best_run.Model{a}.MaxUnits, best_run.Model{a}.firstPrototypes{1}, best_run.Model{a}.firstConnections{1});
        end
        axis equal;
        axis off;
      end
      subplot(nrows, ncols, pidx); pidx=pidx+1;
      hold on;
      %imshow(img);
      colormap([0 0 0; 1 1 1]);
      imagesc([0 1], [1 0], img);
      %image([0 1], [0 1], img);
      %plot(Samples(1,:),Samples(2,:),'*y');
      if isCompetitive
        PlotCompetitive(best_run.Model{a}.Prototypes);
      elseif isSOM
        PlotSOM(best_run.Model{a}.NumRowsMap, best_run.Model{a}.NumColsMap, best_run.Model{a}.Prototypes);
      elseif isGNG
        PlotGNG(best_run.Model{a}.MaxUnits, best_run.Model{a}.Prototypes, best_run.Model{a}.Connections);
      end
      axis equal;
      axis off;
      %set(gca,'YDir','normal');
      xlabel(sprintf('best run: neurons for N=%d and factor=%f', best_run.parameters(a,1), best_run.parameters(a,2)));
      subplot(nrows, ncols, pidx); pidx=pidx+1;
      plot(t, best_run.quantization_errors{a}, 'c', t, best_run.consensus_quantization_errors{a}, 'b', t, best_run.learning_rates{a}, 'k', 'LineWidth',1);
      %plot(t, best_run.quantization_errors{a}, 'b', t, best_run.consensus_quantization_errors{a}, 'r', t, best_run.learning_rates{a}, 'y');
      %legend({'step-by-step QE', 'consensus QE', 'learning rate'});
      legend({'$\mathrm{qe}(t)$', '$C(t)$', '$\alpha(t)$'},'Interpreter','latex');
      %xlabel('best run: training steps');
      xlabel('$t$','Interpreter','latex');
      if mn2<mx2
        ylim([mn2 mx2]);
      end
    end
    
    savefig(h, [FusedFileNames 'experiments.fig'], 'compact');
    close(h);
  end
end

if true && strcmp(trainingType, 'traditional')
  [Model, quantization_errors, consensus_quantization_errors, learning_rates]=TrainOnline(trainMode, Samples{1},netType, params,NumSamplesPerImage*NumEpochs, 10, 111, 1, 1, True, 0, 0, MaxLearningRate, MinLearningRate, [], [], {});
  mqe = compute_consensus_quantization_error(mqeMode, Model.Prototypes, Samples{1});
  save([FileNames{1} '.baseline.mat'], 'Model', 'quantization_errors', 'consensus_quantization_errors', 'learning_rates', 'mqe');
  h = figure('Name', ['baseline experiment (time-decaying LR) for ' FileNames{1}]);
  img = imread(FileNames{1});
  if ndims(img)==3
    img = mean(img, 3);
  end
  img = logical(img);
  t = 1:NumSamplesPerImage*NumEpochs;
  subplot(1, 2, 1);
  hold on;
  imagesc([0 1], [1 0], img);
  if isCompetitive
    PlotCompetitive(Model.Prototypes);
  elseif isSOM
    PlotSOM(Model);
  elseif isGNG
    PlotGNG(Model);
  end
  xlabel('neurons for baseline training');
  subplot(1, 2, 2);
  plot(t, quantization_errors, 'b', t, consensus_quantization_errors, 'r', t, learning_rates, 'y');
  legend({'online QE', 'online MQE', 'learning rate'});
  xlabel('baseline training: training steps');  
  savefig(h, [FileNames{1} '.baseline.fig'], 'compact');
  close(h);
end

