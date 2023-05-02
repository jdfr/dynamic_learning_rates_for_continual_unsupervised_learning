function DemoCompetitiveCORABatches(cora, prefixpath, numRepeats, numfactors, pre_load, do_load, do_show, trainingType, netType, trainMode, mqeMode, slidingWindow, initial_LearningRate, MaxLearningRate, NumEpochs, NumBatches);

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

paperWordsTransposed = double(cora.paperWords');
originalSamples = double(cora.paperWords');
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
    experiments{a}.NumBatches         = NumBatches;
    experiments{a}.NumEpochs          = NumEpochs;
    experiments{a}.params             = params;
    experiments{a}.HistorySizes       = HistorySizes;
    experiments{a}.MaxLearningRate    = MaxLearningRate;
    experiments{a}.MinLearningRate    = MinLearningRate;
    experiments{a}.mqeMode            = mqeMode;
    experiments{a}.slidingWindow      = slidingWindow;
    experiments{a}.initial_LearningRate=initial_LearningRate;
    if isSOM
      experiments{a}.NumRowsMap      = NumRowsMap;
      experiments{a}.NumColsMap      = NumColsMap;
    end
    experiments{a}.netType           = netType;
    experiments{a}.meanQuantizationErrors = zeros(length(experiments{a}.HistorySizes), length(experiments{a}.factors), NumBatches, numRepeats);
    experiments{a}.CalinskiHarabasz       = zeros(length(experiments{a}.HistorySizes), length(experiments{a}.factors), NumBatches, numRepeats);
    experiments{a}.DaviesBouldin          = zeros(length(experiments{a}.HistorySizes), length(experiments{a}.factors), NumBatches, numRepeats);
    experiments{a}.Silhouette             = zeros(length(experiments{a}.HistorySizes), length(experiments{a}.factors), NumBatches, numRepeats);
    experiments{a}.TopographicError       = zeros(length(experiments{a}.HistorySizes), length(experiments{a}.factors), NumBatches, numRepeats);
    experiments{a}.DunnIndex              = zeros(length(experiments{a}.HistorySizes), length(experiments{a}.factors), NumBatches, numRepeats);
    experiments{a}.accuracy               = zeros(length(experiments{a}.HistorySizes), length(experiments{a}.factors), NumBatches, numRepeats);
    %experiments{a}.receptiveFieldSizes    = cell( length(experiments{a}.HistorySizes), length(experiments{a}.factors), NumBatches, numRepeats);
    %experiments{a}.prototypes             = cell( length(experiments{a}.HistorySizes), length(experiments{a}.factors), NumBatches, numRepeats);
    %if isGNG
    %experiments{a}.connections            = cell( length(experiments{a}.HistorySizes), length(experiments{a}.factors), NumBatches, numRepeats);
    %end
    experiments{a}.alldataset_meanQuantizationErrors = zeros(length(experiments{a}.HistorySizes), length(experiments{a}.factors), numRepeats);
    experiments{a}.alldataset_CalinskiHarabasz       = zeros(length(experiments{a}.HistorySizes), length(experiments{a}.factors), numRepeats);
    experiments{a}.alldataset_DaviesBouldin          = zeros(length(experiments{a}.HistorySizes), length(experiments{a}.factors), numRepeats);
    experiments{a}.alldataset_Silhouette             = zeros(length(experiments{a}.HistorySizes), length(experiments{a}.factors), numRepeats);
    experiments{a}.alldataset_TopographicError       = zeros(length(experiments{a}.HistorySizes), length(experiments{a}.factors), numRepeats);
    experiments{a}.alldataset_DunnIndex              = zeros(length(experiments{a}.HistorySizes), length(experiments{a}.factors), numRepeats);
    experiments{a}.alldataset_accuracy               = zeros(length(experiments{a}.HistorySizes), length(experiments{a}.factors), numRepeats);
  end

  %experiments{1}.current.Experiment  = 1;
  %experiments{1}.current.HistorySize = 1;
  %experiments{1}.current.Factor      = 1;

  if do_load
    load([prefixpath 'experiments.mat']);
  else
    if pre_load && isfile([prefixpath 'experiments.mat'])
      load([prefixpath 'experiments.mat']);
    else
      best_run.quantization_errors           =  cell(length(experiments),1);
      %best_run.receptiveFieldSizes           =  cell(length(experiments),1);
      best_run.consensus_quantization_errors =  cell(length(experiments),1);
      best_run.learning_rates                =  cell(length(experiments),1);
      best_run.mqe                           = zeros(length(experiments),1)+1e6;
      best_run.Model                         =  cell(length(experiments),1);
      best_run.parameters                    = zeros(length(experiments),2);
    end
    %for a=experiments{1}.current.Experiment:length(experiments)
    %  for i=experiments{1}.current.HistorySize:size(experiments{a}.meanQuantizationErrors,1)
    %    for j=experiments{1}.current.Factor:size(experiments{a}.meanQuantizationErrors,2)
    for a=1:length(experiments)
      for i=1:size(experiments{a}.meanQuantizationErrors,1)
        for j=1:size(experiments{a}.meanQuantizationErrors,2)
          %disp(sprintf('CORA dataset, experiment %d, T=%d, H=%d, F=%f', a, experiments{a}.type, experiments{a}.HistorySizes(i), experiments{a}.factors(j)));
          runs = cell(NumBatches,numRepeats);
          runs2 = cell(numRepeats,1);
          meanQuantizationErrors = zeros(NumBatches, numRepeats);
          %receptiveFieldSizes = cell(NumBatches, numRepeats);
          %prototypes  = cell(NumBatches, numRepeats);
          %connections = cell(NumBatches, numRepeats);
          parfor b=1:numRepeats
          try
          %for b=1:numRepeats
            [batches, classes_per_batch, batches_expanded] = GenerateSamplesCORABatches(paperWordsTransposed, cora.paperClass, NumEpochs, NumBatches);
            %disp(sprintf('CORA dataset, experiment %d, repeat %d: T=%d, H=%d, F=%f', a, b, experiments{a}.type, experiments{a}.HistorySizes(i), experiments{a}.factors(j)));
            Model = [];
            for c=1:NumBatches
              disp(sprintf('CORA dataset, experiment %d, repeat %d, batch %d: T=%d, H=%d, F=%f', a, b, c, experiments{a}.type, experiments{a}.HistorySizes(i), experiments{a}.factors(j)));
              [Model, quantization_errors, consensus_quantization_errors, learning_rates, meanQuantizationErrors1]=TrainOnline(trainMode,batches_expanded{c},netType, experiments{a}.params,size(batches_expanded{c},2), experiments{a}.HistorySizes(i), experiments{a}.type, experiments{a}.factors(j), mqeMode, slidingWindow, 0, initial_LearningRate, experiments{a}.MaxLearningRate, experiments{a}.MinLearningRate, [], batches{c}, [], {}, Model);
              if isCompetitive
                NumClusters = params.NumNeurons;
                centroids = Model.Prototypes';
              elseif isSOM
                NumClusters = prod(params.NumNeurons);
                centroids = Model.Prototypes(:,:)';
              elseif isGNG
                validNeurons = isfinite(Model.Prototypes(1,:));
                NumClusters = sum(validNeurons);%sum(isfinite(Model.Means(1,:)));
                centroids = Model.Prototypes(:,validNeurons)';
              end
              %[meanQuantizationErrors(c,b), receptiveFieldSizes{c,b}, winners] = compute_consensus_quantization_error(mqeMode, Model.Prototypes, batches{c});
              [meanQuantizationErrors(c,b), ~, winners] = compute_consensus_quantization_error(mqeMode, Model.Prototypes, batches{c});
              %prototypes{c,b} = Model.Prototypes;
              %if isGNG
              %  connections{c,b} = Model.Connections;
              %end
              metrics = computeMetrics(classes_per_batch{c}, Model, NumClusters, winners, centroids, batches{c}, [], [], isCompetitive, isSOM, isGNG);
              runs{c,b} = {Model, quantization_errors, consensus_quantization_errors, learning_rates, metrics};
            end
            [meanQuantizationErrors_alldata, ~, winners] = compute_consensus_quantization_error(mqeMode, Model.Prototypes, originalSamples);
            runs2{b} = computeMetrics(cora.paperClass, Model, NumClusters, winners, centroids, originalSamples, originalSamplesTranposed, dissimMatrix, isCompetitive, isSOM, isGNG);
            runs2{b}(end+1) = meanQuantizationErrors_alldata;
          catch ME
            msgText = getReport(ME,'extended','hyperlinks','off');
            disp(sprintf(['\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n' msgText '\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n']));
          end
          end
          experiments{a}.meanQuantizationErrors(i,j,:,:) = meanQuantizationErrors;
          %experiments{a}.receptiveFieldSizes(i,j,:,:) = receptiveFieldSizes;
          %experiments{a}.prototypes(i,j,:,:) = prototypes;
          %if isGNG
          %  experiments{a}.connections(i,j,:,:) = connections;
          %end

          for b=1:numRepeats
            experiments{a}.alldataset_CalinskiHarabasz(i,j,b) = runs2{b}(1);
            experiments{a}.alldataset_DaviesBouldin(i,j,b)    = runs2{b}(2);
            experiments{a}.alldataset_Silhouette(i,j,b)       = runs2{b}(3);
            experiments{a}.alldataset_TopographicError(i,j,b) = runs2{b}(4);
            experiments{a}.alldataset_DunnIndex(i,j,b)        = runs2{b}(5);
            experiments{a}.alldataset_accuracy(i,j,b)         = runs2{b}(6);
            experiments{a}.alldataset_meanQuantizationErrors(i,j,b) = runs2{b}(7);
            for c=1:NumBatches
              experiments{a}.CalinskiHarabasz(i,j,c,b) = runs{c,b}{5}(1);
              experiments{a}.DaviesBouldin(i,j,c,b)    = runs{c,b}{5}(2);
              experiments{a}.Silhouette(i,j,c,b)       = runs{c,b}{5}(3);
              experiments{a}.TopographicError(i,j,c,b) = runs{c,b}{5}(4);
              experiments{a}.DunnIndex(i,j,c,b)        = runs{c,b}{5}(5);
              experiments{a}.accuracy(i,j,c,b)         = runs{c,b}{5}(6);
            end
            if experiments{a}.meanQuantizationErrors(i,j,end,b) < best_run.mqe(a)
              best_run.mqe(a)                           = experiments{a}.meanQuantizationErrors(i,j,end,b);
              %best_run.receptiveFieldSizes{a}           = experiments{a}.receptiveFieldSizes{i,j,end,b};
              best_run.quantization_errors{a}           = runs{end,b}{2};
              best_run.consensus_quantization_errors{a} = runs{end,b}{3};
              best_run.learning_rates{a}                = runs{end,b}{4};
              best_run.Model{a}                         = runs{end,b}{1};
              best_run.parameters(a,:)                  = [experiments{a}.HistorySizes(i), experiments{a}.factors(j)];
            end
          end
          clear runs, runs2;
          %experiments{1}.current.Factor      = j+1;
        end
        %experiments{1}.current.HistorySize = i+1;
        %if experiments{1}.current.HistorySize<=size(experiments{a}.meanQuantizationErrors,1)
        %  save([prefixpath 'experiments.mat'], 'experiments', 'best_run');
        %end
        save([prefixpath sprintf('experiments_E%d_HS%d.mat', a, i)], 'experiments', 'best_run');
      end 
      %experiments{1}.current.Experiment  = a+1;
      %if experiments{1}.current.Experiment<=length(experiments)
      %  save([prefixpath 'experiments.mat'], 'experiments', 'best_run');
      %end
    end
    save([prefixpath 'experiments.mat'], 'experiments', 'best_run');
  end

  showmode = 0;
  showmode = 1;

  if do_show && showmode==0
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
    t = 1:(experiments{exp_show}.NumSamplesPerBatch*experiments{exp_show}.NumBatches);
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
    sgtitle(sprintf('model=%s, N=%d, multi=%d, slidingWindow=%d, MqeMode=%s, initLR=%f', m, numRepeats, strcmp(trainMode, 'multiLR'), slidingWindow, mqeMode_s, initial_LearningRate));
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

  if do_show && showmode==1
    exp_show = numel(experiments);
    mn = [min(experiments{1}.meanQuantizationErrors, [], 'all'); ...
          min(experiments{1}.CalinskiHarabasz, [], 'all'); ...
          min(experiments{1}.DaviesBouldin, [], 'all'); ...
          min(experiments{1}.Silhouette, [], 'all'); ...
          min(experiments{1}.TopographicError, [], 'all'); ...
          min(experiments{1}.DunnIndex, [], 'all')];
    mx = [max(experiments{1}.meanQuantizationErrors, [], 'all'); ...
          max(experiments{1}.CalinskiHarabasz, [], 'all'); ...
          max(experiments{1}.DaviesBouldin, [], 'all'); ...
          max(experiments{1}.Silhouette, [], 'all'); ...
          max(experiments{1}.TopographicError, [], 'all'); ...
          max(experiments{1}.DunnIndex, [], 'all')];
    for a=2:exp_show
      mn = min(mn, ...
         [min(experiments{a}.meanQuantizationErrors, [], 'all'); ...
          min(experiments{a}.CalinskiHarabasz, [], 'all'); ...
          min(experiments{a}.DaviesBouldin, [], 'all'); ...
          min(experiments{a}.Silhouette, [], 'all'); ...
          min(experiments{a}.TopographicError, [], 'all'); ...
          min(experiments{a}.DunnIndex, [], 'all')]);
      mx = max(mx, ...
         [max(experiments{a}.meanQuantizationErrors, [], 'all'); ...
          max(experiments{a}.CalinskiHarabasz, [], 'all'); ...
          max(experiments{a}.DaviesBouldin, [], 'all'); ...
          max(experiments{a}.Silhouette, [], 'all'); ...
          max(experiments{a}.TopographicError, [], 'all'); ...
          max(experiments{a}.DunnIndex, [], 'all')]);
    end
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
    nrows = exp_show;%numel(experiments);
    ncols = 3;
    h = figure('Name', sprintf('model=%s, NumExperiments=%d, multiple learning rates=%d', m, numRepeats, strcmp(trainMode, 'multiLR')));% ' (Y axes: global MQE after training)']);
    pidx = 1;
    for a=1:nrows
      colidx = 1;
      doPlot(nrows, ncols, pidx, colidx, mn, mx, experiments{a}.factors, experiments{a}.meanQuantizationErrors, 'MQE');
      if pidx==1
        legend({sprintf('$N=%d$', experiments{a}.HistorySizes(1)), sprintf('$N=%d$', experiments{a}.HistorySizes(2)), sprintf('$N=%d$', experiments{a}.HistorySizes(3))},'Interpreter','latex');
      end
      pidx=pidx+1; colidx=colidx+1;
      hold on;
      doPlot(nrows, ncols, pidx, colidx, mn, mx, experiments{a}.factors, experiments{a}.CalinskiHarabasz,       'Cal-Har'); pidx=pidx+1; colidx=colidx+1;
      doPlot(nrows, ncols, pidx, colidx, mn, mx, experiments{a}.factors, experiments{a}.DaviesBouldin,          'Dav-Bou'); pidx=pidx+1; colidx=colidx+1;
    end
    savefig(h, [prefixpath '../../' sprintf('%s_NumExp%d_multipleLR%d_measures_MQE_CH_DB.fig', m, numRepeats, strcmp(trainMode, 'multiLR'))], 'compact');
    close(h);
    h = figure('Name', sprintf('model=%s, NumExperiments=%d, multiple learning rates=%d', m, numRepeats, strcmp(trainMode, 'multiLR')));% ' (Y axes: global MQE after training)']);
    pidx = 1;
    for a=1:nrows
      colidx = 4;
      doPlot(nrows, ncols, pidx, colidx, mn, mx, experiments{a}.factors, experiments{a}.Silhouette,             'Sil');
      if pidx==1
        legend({sprintf('$N=%d$', experiments{a}.HistorySizes(1)), sprintf('$N=%d$', experiments{a}.HistorySizes(2)), sprintf('$N=%d$', experiments{a}.HistorySizes(3))},'Interpreter','latex');
      end
      pidx=pidx+1; colidx=colidx+1;
      hold on;
      doPlot(nrows, ncols, pidx, colidx, mn, mx, experiments{a}.factors, experiments{a}.TopographicError,       'TE'); pidx=pidx+1; colidx=colidx+1;
      doPlot(nrows, ncols, pidx, colidx, mn, mx, experiments{a}.factors, experiments{a}.DunnIndex,              'Dunn'); pidx=pidx+1; colidx=colidx+1;
    end
    savefig(h, [prefixpath '../../' sprintf('%s_NumExp%d_multipleLR%d_measures_S_TE_DI.fig', m, numRepeats, strcmp(trainMode, 'multiLR'))], 'compact');
    close(h);
  end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function doPlot(nrows, ncols, pidx, colidx, mn, mx, xvalues, values, label)

subplot(nrows, ncols, pidx);
meanValues = mean(values, 3);
plot(xvalues, meanValues(1,:), 'g', xvalues, meanValues(2,:), 'b', xvalues, meanValues(3,:), 'k', 'LineWidth',1);
if mn(colidx)<mx(colidx)
  ylim([mn(colidx) mx(colidx)]);
end
ylabel(label);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function metrics = computeMetrics(paperClassPerSample, Model, NumClusters, winners, centroids, originalSamples, originalSamplesTranposed, dissimMatrix, isCompetitive, isSOM, isGNG);

if numel(originalSamplesTranposed)==0
  originalSamplesTranposed = originalSamples';
end

if numel(dissimMatrix)==0
  dissimMatrix = squareform(pdist(originalSamplesTranposed));
  %dissimMatrix = squareform(pdist(cora.paperWords));
end

classesForCentroids = zeros(NumClusters,1);
correctlyClassified = 0;

for k=1:NumClusters
  samplesForCentroid = winners==k;
  classesForCentroid = paperClassPerSample(samplesForCentroid);
  if numel(classesForCentroid)==0
    continue
  end
  [classesForCentroids(k),freq] = mode(classesForCentroid);
  correctlyClassified = correctlyClassified + freq;%sum(classesForCentroid==classesForCentroids(k));
end
accuracy = correctlyClassified/numel(winners);

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

metrics = {CalinskiHarabasz.CriterionValues, DaviesBouldin.CriterionValues, silhouette.CriterionValues, TopographicError, dunnindex, accuracy};
for k=1:numel(metrics)
  if numel(metrics{k})==0
    metrics{k} = nan;
  elseif numel(metrics{k})>1
    metrics{k} = metrics{k}(1);
    error(sprintf('One of the metrics was multidimensional!!!!!! CalinskiHarabasz==%s, DaviesBouldin=%s, silhouette=%s, TopographicError=%s, dunnindex=%s, accuracy=%s', mat2str(CalinskiHarabasz.CriterionValues), mat2str(DaviesBouldin.CriterionValues), mat2str(silhouette.CriterionValues), mat2str(TopographicError), mat2str(dunnindex), mat2str(accuracy)));
  end
end
metrics = [metrics{:}];


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
