function DemoCORA(directory, pre_load, runExperiments, generateGraphs, NumBatches, trainMode, netType);

if false
DemoCORA(true,true);
DemoCORA('experiments_10batches_100exps', true, true, false, 10, {'singleLR', 'multiLR'}, {'competitive', 'som'});

%VER DONDE PONER
DemoCORA('experiments_10batches_100exps', true, true, false, 10, {'multiLR'}, {'competitive'});

%en ejecucion en serv4
DemoCORA('experiments_nobatches_100exps', true, true, false, [], {'singleLR'}, {'competitive'});
DemoCORA('experiments_nobatches_100exps', true, true, false, [], {'multiLR'}, {'competitive'});

DemoCORA('experiments_nobatches_100exps', true, true, false, [], {'singleLR'}, {'som'});
DemoCORA('experiments_10batches_100exps', true, true, false, 10, {'multiLR'}, {'som'});

DemoCORA('experiments_nobatches_100exps', true, true, false, [], {'multiLR'}, {'som'});
DemoCORA('experiments_10batches_100exps', true, true, false, 10, {'singleLR'}, {'som'});

%en ejecucion en serv1
DemoCORA('experiments_nobatches_100exps', true, true, false, [], {'singleLR'}, {'gng'});
DemoCORA('experiments_10batches_100exps', true, true, false, 10, {'singleLR'}, {'competitive'});

%en ejecucion en serv3
DemoCORA('experiments_10batches_100exps', true, true, false, 10, {'singleLR', 'multiLR'}, {'gng'});
DemoCORA('experiments_nobatches_100exps', true, true, false, [], {'multiLR'}, {'gng'});

%ya hecho

DemoCORA('experiments_10batches_100exps', true, true, false, 10, {'singleLR'}, {'som', 'gng'});
DemoCORA('experiments_10batches_100exps', true, true, false, 10, {'multiLR'}, {'som', 'gng'});
DemoCORA('experiments_10batches_100exps', true, true, false, 10, {'singleLR'}, {'gng'});
DemoCORA('experiments_10batches_100exps', true, true, false, 10, {'multiLR'}, {'gng'});
DemoCORA('experiments_10batches_100exps', true, true, false, 10, {'singleLR', 'multiLR'}, {'gng'});
DemoCORA('experiments_10batches_100exps', true, true, false, 10, {'singleLR', 'multiLR'}, {'som'});
DemoCORA('experiments_10batches_100exps', true, true, false, 10, {'singleLR', 'multiLR'}, {'competitive'});
DemoCORA('experiments_nobatches_100exps', true, true, false, [], {'singleLR'}, {'competitive', 'som', 'gng'});
DemoCORA('experiments_nobatches_100exps', true, true, false, [], {'multiLR'}, {'competitive', 'som', 'gng'});
DemoCORA('experiments_nobatches_100exps', true, true, false, [], {'singleLR', 'multiLR'}, {'gng'});
DemoCORA('experiments_nobatches_100exps', true, true, false, [], {'multiLR'}, {'gng'});
DemoCORA('experiments_nobatches_100exps', true, true, false, [], {'multiLR'}, {'som'});


DemoCORA('experiments_nobatches_100exps', false, false, true, [], {'singleLR'}, {'competitive'});

DemoCORA('experiments_nobatches_100exps', false, false, true, [], {'singleLR', 'multiLR'}, {'competitive', 'som', 'gng'});
DemoCORA('experiments_nobatches_100exps', false, false, true, 10, {'singleLR', 'multiLR'}, {'competitive', 'som', 'gng'});


end

cora = slurpCORA();

%netType={'competitive', 'som', 'gng'};
%netType={'competitive', 'gng'};
mqeMode=0;
slidingWindow=true;
experimentName='exploration';
MaxLearningRate=1;
%trainMode = {'singleLR', 'multiLR'};
%trainMode = {'singleLR'};
if numel(NumBatches)==0
  NumSamples = {numel(cora.paperIds)*1, numel(cora.paperIds)*2, numel(cora.paperIds)*4, 40000};
  NumSamples = {numel(cora.paperIds)*4};
else
  %NumSamples is really NumEpochs!!!!
  NumSamples = {4};
end
numRepeats = 100;%[10 100];

[numSamplesIds netTypeIds, trainModeIds, numRepeatsIds]=ndgrid(1:numel(NumSamples), 1:numel(netType), 1:numel(trainMode), 1:numel(numRepeats));
numSamplesIds = numSamplesIds(:);
netTypeIds    = netTypeIds(:);
trainModeIds  = trainModeIds(:);
numRepeatsIds = numRepeatsIds(:);

if runExperiments
  for k=1:numel(numSamplesIds)
  %parfor k=1:numel(numSamplesIds)
    Demo1(directory, cora, numRepeats(numRepeatsIds(k)), pre_load, false, experimentName, netType{netTypeIds(k)}, trainMode{trainModeIds(k)}, mqeMode, slidingWindow, MaxLearningRate, NumSamples{numSamplesIds(k)}, NumBatches);
  end
end
if generateGraphs
  for k=1:numel(numSamplesIds)
    Demo1(directory, cora, numRepeats(numRepeatsIds(k)), pre_load, true, experimentName, netType{netTypeIds(k)}, trainMode{trainModeIds(k)}, mqeMode, slidingWindow, MaxLearningRate, NumSamples{numSamplesIds(k)}, NumBatches);
  end
end

function Demo1(directory, cora, numRepeats, pre_load, makeResults, experimentName, netType, trainMode, mqeMode, slidingWindow, MaxLearningRate, NumSamples, NumBatches);

numfactors = {logspace(-3, log10(5), 100), logspace(-3, log10(5), 100), logspace(-3, log10(5), 100), logspace(-3, 0, 80)};

numfactors = {
0.001:0.001:0.01,
0.0001:0.0002:0.002,
0.01:0.01:0.15,
0.005:0.005:0.05
};

numfactors = {
logspace(-3, log10(5), 40),
logspace(-3, log10(5), 40),
logspace(-3, log10(5), 40),
logspace(-3, 0, 30);
};
initial_LearningRate=MaxLearningRate;

if numel(NumBatches)==0

  prefixpath = sprintf('%s/%s/CORA_%s_%s_N%d_Sliding%d_mqeMode%d_MAXLR_%g_SAM_%d/', directory, netType, experimentName, trainMode, numRepeats, slidingWindow, mqeMode, MaxLearningRate, NumSamples);
  system(sprintf('mkdir -p %s', prefixpath));

  DemoCompetitiveCORA(cora, prefixpath, numRepeats, numfactors, pre_load, makeResults, makeResults, 'online', netType, trainMode, mqeMode, slidingWindow, initial_LearningRate, MaxLearningRate, NumSamples);

else

  prefixpath = sprintf('%s/%s/CORA_%s_BATCH%d_%s_N%d_Sliding%d_mqeMode%d_MAXLR_%g_SAM_%d/', directory, netType, experimentName, NumBatches, trainMode, numRepeats, slidingWindow, mqeMode, MaxLearningRate, NumSamples);
  system(sprintf('mkdir -p %s', prefixpath));

  DemoCompetitiveCORABatches(cora, prefixpath, numRepeats, numfactors, pre_load, makeResults, makeResults, 'online', netType, trainMode, mqeMode, slidingWindow, initial_LearningRate, MaxLearningRate, NumSamples, NumBatches);

end
