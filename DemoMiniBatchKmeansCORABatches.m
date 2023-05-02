function [metrics1, meanmetrics1, metrics2, meanmetrics2, MiniBatchKmeans, all_centroids, all_winners] = DemoMiniBatchKmeansCORABatches(cora, numRepeats, NumEpochs, NumBatches, NumMetaBatches, NumClusters);

mqeMode=0;

all_centroids = cell(numRepeats,1);
all_winners = cell(NumBatches, numRepeats);

paperWordsTransposed = double(cora.paperWords');
originalSamples = double(cora.paperWords');
dissimMatrix = squareform(pdist(cora.paperWords));
originalSamplesTranposed = originalSamples';

runs = cell(NumBatches,numRepeats);
runs2 = cell(numRepeats,1);
meanQuantizationErrors = zeros(NumBatches, numRepeats);
%parfor b=1:numRepeats

Model = {};

try
for b=1:numRepeats
  fprintf('experiment %d/%d\n', b, numRepeats);
  [batches, classes_per_batch, batches_expanded] = GenerateSamplesCORABatches(paperWordsTransposed, cora.paperClass, NumEpochs, NumBatches);
  rnd = round(rand*1000000);
  inputpath  = [sprintf('tmp_in_%d.mat', rnd)];
  outputpath = [sprintf('tmp_out_%d.mat', rnd)];
  save(inputpath, 'batches_expanded');
  command = sprintf('python run.py --num_clusters %d --num_meta_repeats %d --num_repeats %d --input %s --output %s', NumClusters, NumMetaBatches, NumBatches, inputpath, outputpath);
  decodertimer=tic;
  system(command, '-echo');
  decodetime = toc(decodertimer);
  batch_centroids = load(outputpath);
  batch_centroids = batch_centroids.all_centroids;
  all_centroids{b} = batch_centroids;
  delete(inputpath);
  delete(outputpath);
  for c=1:NumBatches
    Prototypes = batch_centroids{c}';
    centroids = batch_centroids{c};
    %[meanQuantizationErrors(c,b), receptiveFieldSizes{c,b}, winners] = compute_consensus_quantization_error(mqeMode, Model.Prototypes, batches{c});
    [meanQuantizationErrors(c,b), ~, winners] = compute_consensus_quantization_error(mqeMode, Prototypes, batches{c});
    all_winners{c,b} = winners;
    metrics = computeMetrics(classes_per_batch{c}, Model, NumClusters, winners, centroids, batches{c}, [], [], true, false, false);
    runs{c,b} = metrics;
  end
  [meanQuantizationErrors_alldata, ~, winners] = compute_consensus_quantization_error(mqeMode, Prototypes, originalSamples);
  runs2{b} = computeMetrics(cora.paperClass, Model, NumClusters, winners, centroids, originalSamples, originalSamplesTranposed, dissimMatrix, true, false, false);
  runs2{b}(end+1) = meanQuantizationErrors_alldata;
end
catch ME
  msgText = getReport(ME,'extended','hyperlinks','off');
  disp(sprintf(['\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n' msgText '\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n']));
end
MiniBatchKmeans.meanQuantizationErrors = meanQuantizationErrors;
MiniBatchKmeans.alldataset_CalinskiHarabasz = zeros(numRepeats,1);
MiniBatchKmeans.alldataset_DaviesBouldin = zeros(numRepeats,1);
MiniBatchKmeans.alldataset_Silhouette = zeros(numRepeats,1);
MiniBatchKmeans.alldataset_TopographicError = zeros(numRepeats,1);
MiniBatchKmeans.alldataset_DunnIndex = zeros(numRepeats,1);
MiniBatchKmeans.alldataset_accuracy = zeros(numRepeats,1);
MiniBatchKmeans.alldataset_meanQuantizationErrors = zeros(numRepeats,1);
MiniBatchKmeans.CalinskiHarabasz = zeros(NumBatches,numRepeats,1);
MiniBatchKmeans.DaviesBouldin = zeros(NumBatches,numRepeats,1);
MiniBatchKmeans.Silhouette = zeros(NumBatches,numRepeats,1);
MiniBatchKmeans.TopographicError = zeros(NumBatches,numRepeats,1);
MiniBatchKmeans.DunnIndex = zeros(NumBatches,numRepeats,1);
MiniBatchKmeans.accuracy = zeros(NumBatches,numRepeats,1);

for b=1:numRepeats
   MiniBatchKmeans.alldataset_CalinskiHarabasz(b) = runs2{b}(1);
   MiniBatchKmeans.alldataset_DaviesBouldin(b)    = runs2{b}(2);
   MiniBatchKmeans.alldataset_Silhouette(b)       = runs2{b}(3);
   MiniBatchKmeans.alldataset_TopographicError(b) = runs2{b}(4);
   MiniBatchKmeans.alldataset_DunnIndex(b)        = runs2{b}(5);
   MiniBatchKmeans.alldataset_accuracy(b)         = runs2{b}(6);
   MiniBatchKmeans.alldataset_meanQuantizationErrors(b) = runs2{b}(7);
  for c=1:NumBatches
     MiniBatchKmeans.CalinskiHarabasz(c,b) = runs{c,b}(1);
     MiniBatchKmeans.DaviesBouldin(c,b)    = runs{c,b}(2);
     MiniBatchKmeans.Silhouette(c,b)       = runs{c,b}(3);
     MiniBatchKmeans.TopographicError(c,b) = runs{c,b}(4);
     MiniBatchKmeans.DunnIndex(c,b)        = runs{c,b}(5);
     MiniBatchKmeans.accuracy(c,b)         = runs{c,b}(6);
  end
end

%metrics = {CalinskiHarabasz.CriterionValues, DaviesBouldin.CriterionValues, silhouette.CriterionValues, TopographicError, dunnindex, accuracy};
metrics1 = zeros(7, NumBatches,numRepeats);
metrics1(1,:,:) = MiniBatchKmeans.CalinskiHarabasz;
metrics1(2,:,:) = MiniBatchKmeans.DaviesBouldin;
metrics1(3,:,:) = MiniBatchKmeans.Silhouette;
metrics1(4,:,:) = MiniBatchKmeans.TopographicError;
metrics1(5,:,:) = MiniBatchKmeans.DunnIndex;
metrics1(6,:,:) = MiniBatchKmeans.accuracy;
metrics1(7,:,:) = MiniBatchKmeans.meanQuantizationErrors;
metrics1 = metrics1([6,1,3,5,7,2,4],:,:);

meanmetrics1 = mean(reshape(metrics1,7,[]),2);

metrics2 = zeros(7, numRepeats);
metrics2(1,:) = MiniBatchKmeans.alldataset_CalinskiHarabasz;
metrics2(2,:) = MiniBatchKmeans.alldataset_DaviesBouldin;
metrics2(3,:) = MiniBatchKmeans.alldataset_Silhouette;
metrics2(4,:) = MiniBatchKmeans.alldataset_TopographicError;
metrics2(5,:) = MiniBatchKmeans.alldataset_DunnIndex;
metrics2(6,:) = MiniBatchKmeans.alldataset_accuracy;
metrics2(7,:) = MiniBatchKmeans.alldataset_meanQuantizationErrors;
metrics2 = metrics2([6,1,3,5,7,2,4],:,:);
meanmetrics2 = mean(metrics2,2);

clear runs, runs2;
%save([prefixpath 'MiniBatchKmeans.mat'], 'MiniBatchKmeans');


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
