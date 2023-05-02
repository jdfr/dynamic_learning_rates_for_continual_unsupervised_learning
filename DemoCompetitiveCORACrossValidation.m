function DemoCompetitiveCORACrossValidation(cora, prefixpath, numRepeats, numfactors, pre_load, do_load, do_show, trainingType, netType, trainMode, mqeMode, slidingWindow, initial_LearningRate, MaxLearningRate, NumSamples, kfold);

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

originalfullCORA = double(cora.paperWords');

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
    experiments{a}.NumSamples         = NumSamples*size(cora.paperWords,1);
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
    experiments{a}.accuracy               = zeros(length(experiments{a}.HistorySizes), length(experiments{a}.factors), numRepeats);
    %experiments{a}.receptiveFieldSizes = cell(length(experiments{a}.HistorySizes), length(experiments{a}.factors), numRepeats);
    %experiments{a}.prototypes             = cell( length(experiments{a}.HistorySizes), length(experiments{a}.factors), numRepeats);
    %if isGNG
    %experiments{a}.connections            = cell( length(experiments{a}.HistorySizes), length(experiments{a}.factors), numRepeats);
    %end
  end

  %experiments{1}.current.Experiment  = 1;
  %experiments{1}.current.HistorySize = 1;
  %experiments{1}.current.Factor      = 1;

  if do_load
    load([prefixpath 'experiments.mat']);
  else
    for a=1:length(experiments)
      for i=1:size(experiments{a}.meanQuantizationErrors,1)
        for j=1:size(experiments{a}.meanQuantizationErrors,2)
          %disp(sprintf('CORA dataset, experiment %d, T=%d, H=%d, F=%f', a, experiments{a}.type, experiments{a}.HistorySizes(i), experiments{a}.factors(j)));
          runs = cell(1,numRepeats);
          meanQuantizationErrors = zeros(numRepeats,1);
          parfor b=1:numRepeats
          %for b=1:numRepeats
            cv=cvpartition(size(cora.paperWords,1),'KFold', kfold);
            allfolds_metrics = zeros(7,kfold);
            for kf=1:kfold
              disp(sprintf('CORA dataset, experiment %d, repeat %d, fold %d: T=%d, H=%d, F=%f', a, b, kf, experiments{a}.type, experiments{a}.HistorySizes(i), experiments{a}.factors(j)));
              forTest            = cv.test(kf);
              forTraining        = ~forTest;
              NumSamplesTraining = NumSamples*sum(forTraining);
              SamplesForTraining = originalfullCORA(:,forTraining);
              SamplesForTest     = originalfullCORA(:,forTest);
              [Samples, ~]       = GenerateSamplesCORA(SamplesForTraining, NumSamplesTraining);
              TheseSamples       = Samples(:,randperm(NumSamplesTraining));
              [Model, quantization_errors, consensus_quantization_errors, learning_rates, meanQuantizationErrors1]=TrainOnline(trainMode,TheseSamples,netType, experiments{a}.params,NumSamplesTraining, experiments{a}.HistorySizes(i), experiments{a}.type, experiments{a}.factors(j), mqeMode, slidingWindow, 0, initial_LearningRate, experiments{a}.MaxLearningRate, experiments{a}.MinLearningRate, [], Samples, [], {});
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
              dissimMatrix = squareform(pdist(SamplesForTest'));
              [allfolds_metrics(7,kf), ~, winners] = compute_consensus_quantization_error(mqeMode, Model.Prototypes, SamplesForTest);
              allfolds_metrics(1:6,kf) = computeMetrics(cora.paperClass, Model, NumClusters, winners, centroids, SamplesForTest, SamplesForTest', dissimMatrix, isCompetitive, isSOM, isGNG);
            end
            meanQuantizationErrors(b) = nanmean(allfolds_metrics(7,:));
            runs{b} = nanmean(allfolds_metrics(1:6,:),2);
          end
          experiments{a}.meanQuantizationErrors(i,j,:) = meanQuantizationErrors;
          %experiments{a}.receptiveFieldSizes(i,j,:) = receptiveFieldSizes;
          %experiments{a}.prototypes(i,j,:) = prototypes;
          %if isGNG
          %  experiments{a}.connections(i,j,:) = connections;
          %end

          for b=1:numRepeats
            experiments{a}.CalinskiHarabasz(i,j,b) = runs{b}(1);
            experiments{a}.DaviesBouldin(i,j,b)    = runs{b}(2);
            experiments{a}.Silhouette(i,j,b)       = runs{b}(3);
            experiments{a}.TopographicError(i,j,b) = runs{b}(4);
            experiments{a}.DunnIndex(i,j,b)        = runs{b}(5);
            experiments{a}.accuracy(i,j,b)         = runs{b}(6);
          end
          clear runs;
          %experiments{1}.current.Factor      = j+1;
        end
        %experiments{1}.current.HistorySize = i+1;
        %if experiments{1}.current.HistorySize<=size(experiments{a}.meanQuantizationErrors,1)
        %  save([prefixpath 'experiments.mat'], 'experiments');
        %end
        save([prefixpath sprintf('experiments_E%d_HS%d.mat', a, i)], 'experiments');
      end
      %experiments{1}.current.Experiment  = a+1;
      %if experiments{1}.current.Experiment<=length(experiments)
      %  save([prefixpath 'experiments.mat'], 'experiments');
      %end
    end
    save([prefixpath 'experiments.mat'], 'experiments');
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
