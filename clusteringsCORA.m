function out = clusteringsCORA(varargin)

%out = kmeansCORAExploration(varargin{:});
%out = DBScanCORA(varargin{:});
out = DBScanCORAPerturbation(varargin{:});
%out = competitiveCORA(varargin{:});
%out = kmeansCORACanon(varargin{:});
%out = spectralCORA(varargin{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ret = spectralCORA(cora, numk, ntimes)

cora.dissimMatrix = squareform(pdist(cora.paperWords));
cora.originalSamplesTranposed = double(cora.paperWords);

metrics = zeros(7, ntimes);
%algorithm = @kmeans;
%algorithm = @kmedoids;

for k=1:ntimes
  metrics(:,k) = oCORASingle(@spectralcluster, cora, 0, numk);
end

meanmetrics = mean(metrics,2);

ret = {metrics, meanmetrics};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ret] = kmeansCORACanon(algorithm, cora,numk, ntimes)
%function err = kmeansCORAMany(cora,numk, ntimes)

cora.dissimMatrix = squareform(pdist(cora.paperWords));
cora.originalSamplesTranposed = double(cora.paperWords);

metrics = zeros(6, ntimes);
%algorithm = @kmeans;
%algorithm = @kmedoids;

for k=1:ntimes
  metrics(:,k) = kCORASingle(algorithm, cora, numk);
end

meanmetrics = mean(metrics,2);

ret = {func2str(algorithm), metrics, meanmetrics};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ret] = competitiveCORA(cora, netType, MaxLearningRate, ntimes)

cora.dissimMatrix = squareform(pdist(cora.paperWords));
cora.originalSamplesTranposed = double(cora.paperWords);

metrics = zeros(7, ntimes);
metrics(:) = nan;

for k=1:ntimes
  metrics(:,k) = competitiveCORASingle(cora, netType, MaxLearningRate);
end

meanmetrics = mean(metrics,2);

ret = {netType, MaxLearningRate, metrics, meanmetrics};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ret] = DBScanCORA(cora, epsilons, minPts)

%con DBSCAN: epsion=1.5, minPts=2 => N==25

cora.dissimMatrix = squareform(pdist(cora.paperWords));
cora.originalSamplesTranposed = double(cora.paperWords);

allmetrics=cell(numel(minPts),1);
for g=1:numel(allmetrics)
  metrics = zeros(9, numel(epsilons));
  for k=1:numel(epsilons)
    metrics(1:7,k) = oCORASingle(@dbscan, cora, 0, epsilons(k), minPts(g));
    metrics(8,k) = epsilons(k);
    metrics(9,k) = minPts(g);
    meanmetrics = mean(metrics(1:7,:),2);
  end
  allmetrics{g} = {metrics meanmetrics};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ret] = DBScanCORAPerturbation(cora, epsilon, minPts, shuffle, ntimes)

%con DBSCAN: epsion=1.5, minPts=2 => N==25

cora.dissimMatrix = squareform(pdist(cora.paperWords));
cora.originalSamplesTranposed = double(cora.paperWords);

metrics = zeros(7, ntimes);
metrics(:) = nan;

for k=1:ntimes
  fprintf('Doing %d/%d\n', k, ntimes);
  metrics(:,k) = oCORASingle(@dbscan, cora, shuffle, epsilon, minPts);
end

meanmetrics = mean(metrics,2);

%ret = {epsilon, minPts, perturbationSize, metrics, meanmetrics};
ret = {epsilon, minPts, metrics, meanmetrics};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%this is for dbscan & kmedoids
function metrics_def = oCORASingle(algorithm, cora, shuffle, varargin)

paperWords = double(cora.paperWords);
%if perturbationSize~=0
%  paperWords = paperWords + (rand(size(paperWords))-0.5)*perturbationSize;
%end
if shuffle
  paperWords = paperWords(randperm(size(paperWords,1)),:);
end
numPapers = size(paperWords,1);
[idxs] = algorithm(paperWords,varargin{:});
uniqueidxs = unique(idxs);
winners = zeros(numel(idxs),1);
numk = numel(uniqueidxs);
dims = size(cora.paperWords,2);
centroids = zeros(numk,dims);
quantization_errors = zeros(numPapers,1);
for k=1:numk
  idx = uniqueidxs(k);
  selected = idxs==idx;
  winners(selected) = k;
  points = double(cora.paperWords(selected,:));
  centroids(k,:) = mean(points,1);
end

for k=1:numPapers
  point = double(cora.paperWords(k,:));
  quantization_errors(k) = sqrt(sum((point-centroids(winners(k),:)).^2, 2));
end


err = sum(quantization_errors)/numk;
metrics = computeMetrics_bis(cora.paperClass, numk, winners, cora.originalSamplesTranposed, cora.dissimMatrix);
%metrics = {CalinskiHarabasz.CriterionValues, DaviesBouldin.CriterionValues, silhouette.CriterionValues, TopographicError, dunnindex, accuracy};
metrics_def = [metrics(6), metrics(1), metrics(3), metrics(5), err, metrics(2), numk];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%this is for dbscan & kmedoids
function metrics_def = competitiveCORASingle(cora, netType, MaxLearningRate)

NumSamples = numel(cora.paperIds)*4;

isCompetitive = strcmp(netType, 'competitive');
isSOM         = strcmp(netType, 'som');
isGNG         = strcmp(netType, 'gng');

[Samples, originalSamples] = GenerateSamplesCORA(cora, NumSamples);

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
trainMode = 'singleLR';
HistorySize = 10;
FunctionType = 111;
FunctionConstant = [];
mqeMode=0;
slidingWindow=true;
initial_consensus_error=0;
initial_LearningRate=0.4;
%MaxLearningRate = 0.4;
MinLearningRate = 0.01;
VideoInfo = [];
SamplesForInit = Samples;
intermediateTimesToRecord = [];
SamplesCell = {};

TheseSamples = Samples(:,randperm(NumSamples));

[Model, ~, ~, ~, ~]=TrainOnline(trainMode, TheseSamples,netType, params,NumSamples, HistorySize, FunctionType, FunctionConstant, mqeMode, slidingWindow, initial_consensus_error, initial_LearningRate, MaxLearningRate, MinLearningRate, VideoInfo, SamplesForInit, intermediateTimesToRecord, SamplesCell);

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
[meanQuantizationError, ~, winners] = compute_consensus_quantization_error(mqeMode, Model.Prototypes, originalSamples);

metrics = computeMetrics(cora.paperClass, Model, NumClusters, winners, centroids, originalSamples, originalSamples', cora.dissimMatrix, isCompetitive, isSOM, isGNG);

%metrics = {CalinskiHarabasz.CriterionValues, DaviesBouldin.CriterionValues, silhouette.CriterionValues, TopographicError, dunnindex, accuracy};
metrics_def = [metrics(6), metrics(1), metrics(3), metrics(5), meanQuantizationError, metrics(2), metrics(4)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%this is for kmeans & kmedoids
function metrics_def = kCORASingle(algorithm, cora,numk)

paperWords = cora.paperWords;
numPapers = size(paperWords,1);
[winners, centroids,  sumdists, alldists] = algorithm(paperWords,numk);
quantization_errors = zeros(numPapers,1);
%winners             = zeros(numPapers, 1);
alldists = sqrt(alldists);
for k=1:numPapers
    %dstvec                    = paperWords(k,:)-centroids(winners(k),:);
    %dst                       = sqrt(sum(dstvec.*dstvec));
    quantization_errors(k)    = alldists(k, winners(k));
end
err = sum(quantization_errors)/numk;
metrics = computeMetrics_bis(cora.paperClass, numk, winners, cora.originalSamplesTranposed, cora.dissimMatrix);
%metrics = {CalinskiHarabasz.CriterionValues, DaviesBouldin.CriterionValues, silhouette.CriterionValues, TopographicError, dunnindex, accuracy};
metrics_def = [metrics(6), metrics(1), metrics(3), metrics(5), err, metrics(2)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function err = kmeansCORAExploration(cora,numk)

paperWords = full(cora.paperWords);

[idxs, centroids,  sumdists, alldists] = kmeans(paperWords,numk);

numPapers = size(cora.paperWords,1);

quantization_errors = zeros(numPapers,1);

for k=1:numPapers
    dstvec                    = paperWords(k,:)-centroids(idxs(k),:);
    dst                       = sqrt(sum(dstvec.*dstvec));
    quantization_errors(k)    = dst;
end

err = sum(quantization_errors)/numk;
prototypes = centroids';
samples = paperWords';
err2 = compute_consensus_quantization_error(0, prototypes, samples);
[err3, receptive_field_sizes, winners] = compute_consensus_quantization_error_bis(prototypes, samples, idxs, paperWords, centroids,alldists);
err4 = sum(sumdists)/numk;

quantization_errors5 = zeros(numPapers,1);
alldists = sqrt(alldists);
for k=1:numPapers
    dstvec                    = paperWords(k,:)-centroids(idxs(k),:);
    dst                       = sqrt(sum(dstvec.*dstvec));
    quantization_errors5(k)    = alldists(k, idxs(k));
end
err5 = sum(quantization_errors5)/numk;

err = {err, err2, err3,err4,err5, receptive_field_sizes}';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function metrics = computeMetrics_bis(paperClassPerSample, NumClusters, winners, originalSamplesTranposed, dissimMatrix);

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
dunnindex        = dunnIndex(NumClusters_adjusted, dissimMatrix, winners);
TopographicError=nan;

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [err, receptive_field_sizes, winners]=compute_consensus_quantization_error_bis(prototypes, samples, idxs, paperWords, centroids,alldists)

[nevermind,nump]    = size(prototypes);
quantization_errors = zeros(size(samples, 2),1);
quantization_errors2 = zeros(size(samples, 2),1);
quantization_errors5 = zeros(size(samples, 2),1);
%assigned_neurons    = zeros(size(samples, 2),1);
%running_totals      = zeros(nump, 1);
%num_samples         = zeros(nump, 1);
receptive_field_sizes = zeros(nump,1);
winners               = zeros(nump, 1);
for k=1:size(samples, 2)
    RepMySample               = repmat(samples(:,k),1,nump);        
    MyDistances               = sqrt(sum((RepMySample-prototypes(:,:)).^2,1));            			    
    [~,NdxWinner]             = min(MyDistances);
    winners(k)                = NdxWinner;
    dstvec1                   = samples(:,k)-prototypes(:,NdxWinner);
    dst1                      = sqrt(sum(dstvec1.*dstvec1));
    dstvec                    = samples(:,k)-prototypes(:,idxs(k));
    dst                       = sqrt(sum(dstvec.*dstvec));
    receptive_field_sizes(NdxWinner) = receptive_field_sizes(NdxWinner)+1;
    %running_totals(NdxWinner) = running_totals(NdxWinner) + dst;
    %num_samples(NdxWinner)    = num_samples(NdxWinner) + 1;
    quantization_errors(k) = dst;
    %assigned_neurons(k)    = NdxWinner;
    dstvec2                    = paperWords(k,:)-centroids(idxs(k),:);
    dst2                       = sqrt(sum(dstvec.*dstvec));
    quantization_errors5(k) = alldists(k, idxs(k));
    quantization_errors2(k)    = dst2;    
end

err = sum(quantization_errors)/nump;
err2 = sum(quantization_errors2)/nump;
err5 = sum(quantization_errors5)/nump;
err3=err2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
