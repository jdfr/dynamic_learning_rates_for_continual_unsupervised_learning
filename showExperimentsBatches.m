function showExperimentsBatches(varargin);

showExperiments1(varargin{:});
%showExperiments2(varargin{:});

function showExperiments2(pathS, pathM, measure, histSiz);

if false
showExperimentsBatches('experiments_10batches_100exps/competitive/CORA_exploration_BATCH10_singleLR_N100_Sliding1_mqeMode0_MAXLR_1_SAM_4/experiments.mat', 'experiments_10batches_100exps/competitive/CORA_exploration_BATCH10_multiLR_N100_Sliding1_mqeMode0_MAXLR_1_SAM_4/experiments.mat', 1, 1, 1);
end

load(pathS);
experimentsS = experiments;
load(pathM);
experimentsM = experiments;
clear experiments;

labnames = {'meanQuantizationErrors', 'CalinskiHarabasz', 'DaviesBouldin', 'Silhouette', 'TopographicError', 'DunnIndex', 'accuracy'};

funcnames = {'linear', 'quadratic', 'inverse', 'constant'};

histSizes = [10, 50, 100];

colores1 = {'r', 'g', 'b', 'k'};
colores2 = {'c', 'm', 'y', '#77AC30'};
ls1 = {'-', '-', '-', '-'};
ls2 = {'-', '-', '-', ':'};
ns1 = {'y=ax (S)', 'y=ax^2 (S)', 'y=a/x (S)', 'y=a (S)'};
ns2 = {'y=ax (M)', 'y=ax^2 (M)', 'y=a/x (M)', 'y=a (M)'};
allnames = {ns1{1}, ns2{1}, ns1{2}, ns2{2}, ns1{3}, ns2{3}, ns1{4}, ns2{4}};

h = figure('Name', sprintf('net=%s, measure=%s, histSize=%g', experimentsS{1}.netType, labnames{measure}, histSizes(histSiz)));

tiledlayout(2,1,'TileSpacing','compact');

nexttile;
hold on;
for k=1:4
  valsS = getfield(experimentsS{k}, labnames{measure});
  valsM = getfield(experimentsM{k}, labnames{measure});
  plot(experimentsS{k}.factors, reshape(nanmean(valsS(histSiz,:,:), 3), numel(experimentsS{k}.factors), 1), 'Color', colores1{k}, 'LineStyle', ls1{k}, 'LineWidth',1);
  plot(experimentsM{k}.factors, reshape(nanmean(valsM(histSiz,:,:), 3), numel(experimentsM{k}.factors), 1), 'Color', colores2{k}, 'LineStyle', ls2{k}, 'LineWidth',1);
end
set(gca, 'XScale', 'log');
legend(allnames);
title('mean of all batches')

nexttile;
hold on;
for k=1:4
  valsS = getfield(experimentsS{k}, labnames{measure});
  valsM = getfield(experimentsM{k}, labnames{measure});
  plot(experimentsS{k}.factors, reshape(nanmean(valsS(histSiz,:,end,:), 4), numel(experimentsS{k}.factors), 1), 'Color', colores1{k}, 'LineStyle', ls1{k}, 'LineWidth',1);
  plot(experimentsM{k}.factors, reshape(nanmean(valsM(histSiz,:,end,:), 4), numel(experimentsM{k}.factors), 1), 'Color', colores2{k}, 'LineStyle', ls2{k}, 'LineWidth',1);
end
set(gca, 'XScale', 'log');
%legend(allnames);
title('mean of end batches')

function showExperiments1(experiments, measure, func, histSiz, fac);

labnames = {'meanQuantizationErrors', 'CalinskiHarabasz', 'DaviesBouldin', 'Silhouette', 'TopographicError', 'DunnIndex', 'accuracy'};

funcnames = {'linear', 'quadratic', 'inverse', 'constant'};

histSizes = [10, 50, 100];

if false
showExperimentsBatches(experiments, 1, 1, 1, 1);
end

h = figure('Name', sprintf('net=%s, measure=%s, function=%s, histSize=%g, factor=%g', experiments{1}.netType, labnames{measure}, funcnames{func}, histSizes(histSiz), experiments{func}.factors(fac)));

vals = getfield(experiments{func}, labnames{measure});

tiledlayout(3,3,'TileSpacing','compact');
for k=1:8;
  nexttile;
  plot(1:10, reshape(vals(histSiz,fac,:,k), 1,10));
  title(sprintf('experiment %d', k));
end;
nexttile;
plot(1:10, reshape(mean(vals(histSiz,fac,:,:), 4), 1,10));
title('mean of 100 experiments');
