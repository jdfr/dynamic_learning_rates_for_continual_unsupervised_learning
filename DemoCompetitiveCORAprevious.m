function DemoCompetitiveCORA(cora, prefixpath, numRepeats, numfactors, do_load, do_show, trainingType, netType, trainMode, mqeMode, slidingWindow, initial_LearningRate, MaxLearningRate, NumSamples);

if prefixpath(end)~='/'
  prefixpath = [prefixpath '/'];
end

if ~(strcmp(netType, 'competitive') || strcmp(netType, 'som') || strcmp(netType, 'gng'))
  error(sprintf('netType must be som, competitive or gng'));
end
isCompetitive = strcmp(netType, 'competitive');
isSOM         = strcmp(netType, 'som');
isGNG         = strcmp(netType, 'gng');

if ~(strcmp(trainMode, 'singleLR') || strcmp(trainMode, 'multiLR'))
  error(sprintf('trainMode must be singleLR, or multiLR'));
end

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


[Samples, originalSamples] = GenerateSamplesCORA(cora, NumSamples);
dissimMatrix = squareform(pdist(cora.paperWords));
originalSamplesTranposed = originalSamples';

if true && strcmp(trainingType, 'online')
  a=0;

  a=a+1;
  experiments{a}.type    = 1;
  experiments{a}.name    = 'linear';
  experiments{a}.eq      = 'linear';
  experiments{a}.factors = numfactors{a};%logspace(-3, log10(5), 100);%[0.1:0.1:1.2 1.5:0.5:10];

  a=a+1;
  experiments{a}.type    = 2;
  experiments{a}.name    = 'quadratic';
  experiments{a}.eq      = 'quadratic';
  experiments{a}.factors = numfactors{a};%logspace(-3, log10(5), 100);%[0.2:0.2:2.4 3:10 20:10:70];

if false
  a=a+1;
  experiments{a}.type    = 3;
  experiments{a}.name    = 'cubic';
  experiments{a}.eq      = 'cubic';
  experiments{a}.factors = numfactors{a};%logspace(-3, log10(5), 100);%[0.4:0.4:4.8 6:10 100:100:700];
end

%  a=a+1;
%  experiments{a}.type    = 4;
%  experiments{a}.name    = 'exponential';
%  experiments{a}.eq      = 'LR=min(F*exp(MQE)*0.4/exp(0.5),0.4)';
%  experiments{a}.factors = factors{a};%logspace(-1, 1, 10);%[0.1:0.1:1.2 1.5:0.5:10];

  a=a+1;
  experiments{a}.type    = 6;
  experiments{a}.name    = 'inverse';
  experiments{a}.eq      = 'inverse';
  experiments{a}.factors = numfactors{a};%logspace(-3, log10(5), 100);%[0.1:0.1:1.2 1.5:0.5:10];

  a=a+1;
  experiments{a}.type    = 5;
  experiments{a}.name    = 'constant';
  experiments{a}.eq      = 'constant';
  %experiments{a}.factors = linspace(0.005, 1, 200);%[0.1:0.1:1.2 1.5:0.5:10];
  experiments{a}.factors = numfactors{a};%logspace(-3, 0, 80);%[0.1:0.1:1.2 1.5:0.5:10];

  for a=1:length(experiments)
    experiments{a}.numRepeats         = numRepeats;
    experiments{a}.prefixpath         = prefixpath;
    experiments{a}.NumSamples         = NumSamples;
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
    experiments{a}.meanQuantizationErrors = zeros(length(experiments{a}.HistorySizes), length(experiments{a}.factors), numRepeats);
    experiments{a}.CalinskiHarabasz       = zeros(length(experiments{a}.HistorySizes), length(experiments{a}.factors), numRepeats);
    experiments{a}.DaviesBouldin          = zeros(length(experiments{a}.HistorySizes), length(experiments{a}.factors), numRepeats);
    experiments{a}.Silhouette             = zeros(length(experiments{a}.HistorySizes), length(experiments{a}.factors), numRepeats);
    experiments{a}.TopographicError       = zeros(length(experiments{a}.HistorySizes), length(experiments{a}.factors), numRepeats);
    experiments{a}.DunnIndex              = zeros(length(experiments{a}.HistorySizes), length(experiments{a}.factors), numRepeats);
    experiments{a}.receptiveFieldSizes = cell(length(experiments{a}.HistorySizes), length(experiments{a}.factors), numRepeats);
  end

  if do_load
    load([prefixpath 'experiments.mat']);
  else
    best_run.quantization_errors           =  cell(length(experiments),1);
    best_run.receptiveFieldSizes           =  cell(length(experiments),1);
    best_run.consensus_quantization_errors =  cell(length(experiments),1);
    best_run.learning_rates                =  cell(length(experiments),1);
    best_run.mqe                           = zeros(length(experiments),1)+1e6;
    best_run.Model                         =  cell(length(experiments),1);
    best_run.parameters                    = zeros(length(experiments),2);
    for a=1:length(experiments)
      for i=1:size(experiments{a}.meanQuantizationErrors,1)
        for j=1:size(experiments{a}.meanQuantizationErrors,2)
          %disp(sprintf('CORA dataset, experiment %d, T=%d, H=%d, F=%f', a, experiments{a}.type, experiments{a}.HistorySizes(i), experiments{a}.factors(j)));
          runs = cell(1,numRepeats);
          meanQuantizationErrors = zeros(numRepeats,1);
          receptiveFieldSizes = cell(numRepeats,1);
          parfor b=1:numRepeats
            disp(sprintf('CORA dataset, experiment %d, repeat %d: T=%d, H=%d, F=%f', a, b, experiments{a}.type, experiments{a}.HistorySizes(i), experiments{a}.factors(j)));
            TheseSamples = Samples(:,randperm(experiments{a}.NumSamples));
            [Model, quantization_errors, consensus_quantization_errors, learning_rates, meanQuantizationErrors1]=TrainOnline(trainMode,TheseSamples,netType, experiments{a}.params,NumSamples, experiments{a}.HistorySizes(i), experiments{a}.type, experiments{a}.factors(j), mqeMode, slidingWindow, 0, initial_LearningRate, experiments{a}.MaxLearningRate, experiments{a}.MinLearningRate, [], Samples, [], {});
            if isCompetitive
              NumClusters = params.NumNeurons;
              centroids = Model.Prototypes';
            elseif isSOM
              NumClusters = prod(params.NumNeurons);
              centroids = Model.Prototypes(:,:)';
            elseif isGNG
              NumClusters = find(isfinite(Model.Means(1,:)));
              validNeurons = isfinite(Model.Prototypes(1,:));
              centroids = Model.Prototypes(:,validNeurons)';
            end
            [meanQuantizationErrors(b), receptiveFieldSizes{b}, winners] = compute_consensus_quantization_error(mqeMode, Model.Prototypes, originalSamples);
            metrics = computeMetrics(Model, NumClusters, winners, centroids, originalSamples, originalSamplesTranposed, dissimMatrix, isCompetitive, isSOM, isGNG);
            runs{b} = {Model, quantization_errors, consensus_quantization_errors, learning_rates, metrics};

            %[experiments{a}.meanQuantizationErrors(i,j,b), experiments{a}.receptiveFieldSizes{i,j,b}] = compute_consensus_quantization_error(mqeMode, Model.Prototypes, originalSamples);
            %if experiments{a}.meanQuantizationErrors(i,j,b) < best_run.mqe(a)
            %  best_run.mqe(a)                           = experiments{a}.meanQuantizationErrors(i,j,b);
            %  best_run.receptiveFieldSizes{a}           = experiments{a}.receptiveFieldSizes{i,j,b};
            %  best_run.quantization_errors{a}           = quantization_errors;
            %  best_run.consensus_quantization_errors{a} = consensus_quantization_errors;
            %  best_run.learning_rates{a}                = learning_rates;
            %  best_run.Model{a}                         = Model;
            %  best_run.parameters(a,:)                  = [experiments{a}.HistorySizes(i), experiments{a}.factors(j)];
            %end
          end
          experiments{a}.meanQuantizationErrors(i,j,:) = meanQuantizationErrors;
          experiments{a}.receptiveFieldSizes(i,j,:) = receptiveFieldSizes;

          for b=1:numRepeats
            experiments{a}.CalinskiHarabasz(i,j,b) = runs{b}{5}(1);
            experiments{a}.DaviesBouldin(i,j,b)    = runs{b}{5}(2);
            experiments{a}.Silhouette(i,j,b)       = runs{b}{5}(3);
            experiments{a}.TopographicError(i,j,b) = runs{b}{5}(4);
            experiments{a}.DunnIndex(i,j,b)        = runs{b}{5}(5);
            if experiments{a}.meanQuantizationErrors(i,j,b) < best_run.mqe(a)
              best_run.mqe(a)                           = experiments{a}.meanQuantizationErrors(i,j,b);
              best_run.receptiveFieldSizes{a}           = experiments{a}.receptiveFieldSizes{i,j,b};
              best_run.quantization_errors{a}           = runs{b}{2};
              best_run.consensus_quantization_errors{a} = runs{b}{3};
              best_run.learning_rates{a}                = runs{b}{4};
              best_run.Model{a}                         = runs{b}{1};
              best_run.parameters(a,:)                  = [experiments{a}.HistorySizes(i), experiments{a}.factors(j)];
            end
          end
          clear runs;
       end
      end 
    end
    save([prefixpath 'experiments.mat'], 'experiments', 'best_run');
  end

  if do_show
    exp_show = numel(experiments);
    mn = min(experiments{1}.meanQuantizationErrors, [], 'all');
    mx = max(experiments{1}.meanQuantizationErrors, [], 'all');
    mn2 = min([min(best_run.quantization_errors{1}, [], 'all'), min(best_run.consensus_quantization_errors{1}, [], 'all'), min(best_run.learning_rates{1}, [], 'all')]);
    mx2 = max([max(best_run.quantization_errors{1}, [], 'all'), max(best_run.consensus_quantization_errors{1}, [], 'all'), max(best_run.learning_rates{1}, [], 'all')]);
    for a=2:exp_show
      mn = min(min(experiments{a}.meanQuantizationErrors, [], 'all'), mn);
      mx = max(max(experiments{a}.meanQuantizationErrors, [], 'all'), mx);
      mn2 = min(min([min(best_run.quantization_errors{a}, [], 'all'), min(best_run.consensus_quantization_errors{a}, [], 'all'), min(best_run.learning_rates{a}, [], 'all')]), mn2);
      mx2 = max(max([max(best_run.quantization_errors{a}, [], 'all'), max(best_run.consensus_quantization_errors{a}, [], 'all'), max(best_run.learning_rates{a}, [], 'all')]), mx2);
    end
    t = 1:experiments{exp_show}.NumSamples;
    h = figure('Name', ['experiments for CORA']);% ' (Y axes: global MQE after training)']);
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
    ncols = 3;
    for a=1:nrows
      subplot(nrows, ncols, pidx); pidx=pidx+1;
      meanValues = mean(experiments{a}.meanQuantizationErrors, 3);
      plot(experiments{a}.factors, meanValues(1,:), 'g', experiments{a}.factors, meanValues(2,:), 'b', experiments{a}.factors, meanValues(3,:), 'k', 'LineWidth',1);
      legend({sprintf('$N=%d$', experiments{a}.HistorySizes(1)), sprintf('$N=%d$', experiments{a}.HistorySizes(2)), sprintf('$N=%d$', experiments{a}.HistorySizes(3))},'Interpreter','latex')
      xlabel([experiments{a}.eq sprintf(' scaling. Minimal MQE %g for $N=%d$ and factor %g', best_run.mqe(a), best_run.parameters(a,1), best_run.parameters(a,2))],'Interpreter','latex');
      if mn<mx
        ylim([mn mx]);
      end
      ylabel('MQE');

      subplot(nrows, ncols, pidx); pidx=pidx+1;
      hold on;
      [coeffs,score,~,~,explanation,mu] = pca(double(cora.paperWords));
      plot(score(:,1),score(:,2),'.b');
      prots = best_run.Model{a}.Prototypes;
      if isSOM
        sizeSOM = size(prots);
        prots = prots(:,:);
      end
      prots = ((prots'-mu)*coeffs)';
      if isSOM
        prots = reshape(prots, sizeSOM);
      end
      if isCompetitive
        PlotCompetitive(prots);
      elseif isSOM
        PlotSOM(best_run.Model{a}.NumRowsMap, best_run.Model{a}.NumColsMap, prots);
      elseif isGNG
        PlotGNG(best_run.Model{a}.MaxUnits, prots, best_run.Model{a}.Connections);
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
    
    savefig(h, [prefixpath 'experiments.fig'], 'compact');
    close(h);
  end
end

if true && strcmp(trainingType, 'traditional')
  NumEpochs = 2;
  [Model, quantization_errors, consensus_quantization_errors, learning_rates]=TrainOnline(trainMode,Samples{1},netType, params,NumSamplesPerImage*NumEpochs, 10, 111, 1, 1, True, 0, 0, MaxLearningRate, MinLearningRate, [], [], {});
  mqe = compute_consensus_quantization_error(mqeMode, Model.Prototypes, Samples{1});
  save([prefixpath '.baseline.mat'], 'Model', 'quantization_errors', 'consensus_quantization_errors', 'learning_rates', 'mqe');
  h = figure('Name', ['baseline experiment (time-decaying LR) for CORA']);
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function metrics = computeMetrics(Model, NumClusters, winners, centroids, originalSamples, originalSamplesTranposed, dissimMatrix, isCompetitive, isSOM, isGNG);

[winners, NumClusters_adjusted] = normalize_winners(winners);

CalinskiHarabasz = evalclusters(originalSamplesTranposed, winners, 'CalinskiHarabasz');
DaviesBouldin    = evalclusters(originalSamplesTranposed, winners, 'DaviesBouldin');
silhouette       = evalclusters(originalSamplesTranposed, winners, 'silhouette');
if isCompetitive
  TopographicError = 0;
elseif isSOM
  TopographicError = TopographicErrorSOM(Model, originalSamples);
elseif isGNG
  TopographicError = TopographicErrorGNG(Model, originalSamples);
end
dunnindex        = dunnIndex(NumClusters_adjusted, dissimMatrix, winners);

metrics = [CalinskiHarabasz.CriterionValues DaviesBouldin.CriterionValues silhouette.CriterionValues TopographicError dunnindex];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [normalized, numclusters] = normalize_winners(winners);

normalized = winners;

maxK = max(winners);

for k=maxK-1:-1:1
  areK = winners==k;
  if sum(areK)==0
    normalized(normalized>k) = normalized(normalized>k) - 1;
  end
end

numclusters = max(normalized);
