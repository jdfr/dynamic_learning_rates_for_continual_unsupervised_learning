function showAllPerNAndmulti(config)

%showAllPerNAndmulti(1); showAllPerNAndmulti(2); showAllPerNAndmulti(3); showAllPerNAndmulti(4); 

%interesan: config==1 y 4

if config==0

prefix = 'experiments_long_ranges/';
postFixS = 'CORA_exploration_singleLR_N10_Sliding1_mqeMode0_MAXLR_1_SAM_10832/';
postFixM = 'CORA_exploration_multiLR_N10_Sliding1_mqeMode0_MAXLR_1_SAM_10832/';
typeReduction = 'all';
expression = 'each point==mean of 10 experiments. ';
names = {'MQE', 'Cal-Har', 'Dav-Bou', 'Sil', 'Dunn', 'TopErr'};
labnames = {'meanQuantizationErrors', 'CalinskiHarabasz', 'DaviesBouldin', 'Silhouette', 'DunnIndex', 'TopographicError'};
titlename = '';
args = {...
  {'competitive', ...
   [prefix 'competitive/' postFixS], ...
   [prefix 'competitive/' postFixM], ...
   [prefix 'singlePlot2/'],...
   typeReduction, expression, names, labnames, titlename},...
  {'som', ...
   [prefix 'som/' postFixS], ...
   [prefix 'som/' postFixM], ...
   [prefix 'singlePlot2/'], ...
   typeReduction, expression, names, labnames, titlename},...
  {'gng', ...
   [prefix 'gng/' postFixS], ...
   [prefix 'gng/' postFixM], ...
   [prefix 'singlePlot2/'],...
   typeReduction, expression, names, labnames, titlename},...
};

elseif config==1

prefix = 'experiments_10batches_100exps/';
postFixS = 'CORA_exploration_BATCH10_singleLR_N100_Sliding1_mqeMode0_MAXLR_1_SAM_4/';
postFixM = 'CORA_exploration_BATCH10_multiLR_N100_Sliding1_mqeMode0_MAXLR_1_SAM_4/';
typeReduction = 'all';
expression = 'each point==mean of 10 batches in 100 experiments. ';
expression = 'ev. all batches from 100 exp./conf., each one with itself. ';
names = {'accuracy', 'Calinski-Harabasz', 'Silhouette', 'Dunn''s index', 'MSE', 'Davies-Bouldin', 'Topographic Error'};
labnames = {'accuracy', 'CalinskiHarabasz', 'Silhouette', 'DunnIndex', 'meanQuantizationErrors', 'DaviesBouldin', 'TopographicError'};
%ylims = {[0.13 0.58], [1.04 10.8], [-0.15 0.062], [0.13 0.56], [40.65 56.2], [1.4 7.05], [0 0.9]};
ylims = {[0.3 0.58], [1.04 2.15], [-0.15 0.062], {[0.285 0.56], [0.285 0.39], [0.285 0.56]}, {[40 55.2], [40 46.45], [40 46.1]}, {[1.9 4.1], [2.8 4.1], [1.4 4.1]}, {[], [0.65 0.89], [0 0.165]}};
yticks = cell(size(ylims));
yticks = {0.3:0.05:0.55, [], [], [], {40:2:56, 40:2:46, 40:2:46}, {2:0.5:4, 2:0.5:4, 1.5:0.5:4}, []};
titlename = '_allBatches_evalByBatch';
args = {...
  {'competitive', ...
   [prefix 'competitive/' postFixS], ...
   [prefix 'competitive/' postFixM], ...
   [prefix 'singlePlot_new/'],...
   typeReduction, expression, names, labnames, titlename},...
  {'som', ...
   [prefix 'som/' postFixS], ...
   [prefix 'som/' postFixM], ...
   [prefix 'singlePlot_new/'], ...
   typeReduction, expression, names, labnames, titlename},...
  {'gng', ...
   [prefix 'gng/' postFixS], ...
   [prefix 'gng/' postFixM], ...
   [prefix 'singlePlot_new/'],...
   typeReduction, expression, names, labnames, titlename},...
};

elseif config==4

prefix = 'experiments_nobatches_100exps/';
postFixS = 'CORA_exploration_singleLR_N100_Sliding1_mqeMode0_MAXLR_1_SAM_10832/';
postFixM = 'CORA_exploration_multiLR_N100_Sliding1_mqeMode0_MAXLR_1_SAM_10832/';
typeReduction = 'all';
expression = 'each point==mean of 100 experiments without batches. ';
expression = 'ev. 100 experiments per configuration  without batches. ';
names = {'accuracy', 'Calinski-Harabasz', 'Silhouette', 'Dunn''s index', 'MSE', 'Davies-Bouldin', 'Topographic Error'};
labnames = {'accuracy', 'CalinskiHarabasz', 'Silhouette', 'DunnIndex', 'meanQuantizationErrors', 'DaviesBouldin', 'TopographicError'};
%ylims = {[0.13 0.58], [1.04 10.8], [-0.15 0.062], [0.13 0.56], [422 612],    [1.4 7.05], [0 0.9]};
ylims = {[0.3 0.58], [1.04 10.8], [-0.15 0.062], {[0.13 0.28], [0.13 0.18] [0.13 0.47]}, {[420 612], [420 484], [420 484]}, {[2.8 7.05], [4.5 7.05], [1.4 7.05]}, {[], [0.65 0.89], [0 0.165]}};
yticks = {0.3:0.05:0.55, [], [], [], {420:40:600, 420:20:480, 420:20:480}, {3:0.5:7, 4.5:0.5:7, 1.5:0.5:7}, []};
titlename = '_noBatches';
args = {...
  {'competitive', ...
   [prefix 'competitive/' postFixS], ...
   [prefix 'competitive/' postFixM], ...
   [prefix 'singlePlot_new/'],...
   typeReduction, expression, names, labnames, titlename},...
  {'som', ...
   [prefix 'som/' postFixS], ...
   [prefix 'som/' postFixM], ...
   [prefix 'singlePlot_new/'], ...
   typeReduction, expression, names, labnames, titlename},...
  {'gng', ...
   [prefix 'gng/' postFixS], ...
   [prefix 'gng/' postFixM], ...
   [prefix 'singlePlot_new/'],...
   typeReduction, expression, names, labnames, titlename},...
};

elseif config==2

prefix = 'experiments_10batches_100exps/';
postFixS = 'CORA_exploration_BATCH10_singleLR_N100_Sliding1_mqeMode0_MAXLR_1_SAM_4/';
postFixM = 'CORA_exploration_BATCH10_multiLR_N100_Sliding1_mqeMode0_MAXLR_1_SAM_4/';
typeReduction = 'last';
expression = 'each point==mean of last batches in 100 experiments, evByBatch. ';
expression = 'ev. last batches from 100 exp./conf., each one with itself. ';
names = {'accuracy', 'Calinski-Harabasz', 'Silhouette', 'Dunn''s index', 'MQE', 'Davies-Bouldin', 'Topographic Error'};
labnames = {'accuracy', 'CalinskiHarabasz', 'Silhouette', 'DunnIndex', 'meanQuantizationErrors', 'DaviesBouldin', 'TopographicError'};
ylims = {[0.13 0.58], [1.04 10.8], [-0.15 0.062], [0.13 0.56], [40.65 56.2], [1.4 7.05], [0 0.9]};
titlename = '_lastBatches_evalByBatch';
args = {...
  {'competitive', ...
   [prefix 'competitive/' postFixS], ...
   [prefix 'competitive/' postFixM], ...
   [prefix 'singlePlot/'],...
   typeReduction, expression, names, labnames, titlename},...
  {'som', ...
   [prefix 'som/' postFixS], ...
   [prefix 'som/' postFixM], ...
   [prefix 'singlePlot/'], ...
   typeReduction, expression, names, labnames, titlename},...
  {'gng', ...
   [prefix 'gng/' postFixS], ...
   [prefix 'gng/' postFixM], ...
   [prefix 'singlePlot/'],...
   typeReduction, expression, names, labnames, titlename},...
};

elseif config==3

prefix = 'experiments_10batches_100exps/';
postFixS = 'CORA_exploration_BATCH10_singleLR_N100_Sliding1_mqeMode0_MAXLR_1_SAM_4/';
postFixM = 'CORA_exploration_BATCH10_multiLR_N100_Sliding1_mqeMode0_MAXLR_1_SAM_4/';
typeReduction = 'all';
expression = 'each point==mean of last batches in 100 experiments, evByBatch. ';
expression = 'ev. last batches from 100 exp./conf., each one with whole dataset. ';
names = {'accuracy', 'Calinski-Harabasz', 'Silhouette', 'Dunn''s index', 'MQE', 'Davies-Bouldin', 'Topographic Error'};
labnames = {'alldataset_accuracy', 'alldataset_CalinskiHarabasz', 'alldataset_Silhouette', 'alldataset_DunnIndex', 'alldataset_meanQuantizationErrors', 'alldataset_DaviesBouldin', 'alldataset_TopographicError'};
ylims = {[0.13 0.58], [1.04 10.8], [-0.15 0.062], [0.13 0.56], [422 612],    [1.4 7.05], [0 0.9]};
titlename = '_lastBatches_evalByWholeDataset';
args = {...
  {'competitive', ...
   [prefix 'competitive/' postFixS], ...
   [prefix 'competitive/' postFixM], ...
   [prefix 'singlePlot/'],...
   typeReduction, expression, names, labnames, titlename},...
  {'som', ...
   [prefix 'som/' postFixS], ...
   [prefix 'som/' postFixM], ...
   [prefix 'singlePlot/'], ...
   typeReduction, expression, names, labnames, titlename},...
  {'gng', ...
   [prefix 'gng/' postFixS], ...
   [prefix 'gng/' postFixM], ...
   [prefix 'singlePlot/'],...
   typeReduction, expression, names, labnames, titlename},...
};

end

%config=1
%Minimos / maximos MQE: [41.0548934240704 55.2116124485802]
%Minimos / maximos Calinski-Harabasz: [1.04897420244243 2.12893435704484]
%Minimos / maximos Davies-Bouldin: [1.47191842527805 4.04718697463857]
%Minimos / maximos Silhouette: [-0.131272382193343 0.0518956618544265]
%Minimos / maximos Dunn's index: [0.285594587672206 0.551672251473435]
%Minimos / maximos accuracy: [0.303084515511822 0.558307844745114]
%Minimos / maximos Topographic Error: [6.27442941096077e-05 0.888911644116446]
%config=2
%Minimos / maximos MQE: [40.6993921595995 56.1517186557053]
%Minimos / maximos Calinski-Harabasz: [1.10930584612874 2.57429997098607]
%Minimos / maximos Davies-Bouldin: [1.89181136729791 5.66696643568219]
%Minimos / maximos Silhouette: [-0.13340589180118 0.0613131556473284]
%Minimos / maximos Dunn's index: [0.277675446736854 0.558827354715812]
%Minimos / maximos accuracy: [0.299574142408091 0.576638923055897]
%Minimos / maximos Topographic Error: [0 0.911360393603936]
%config=3
%Minimos / maximos MQE: [424.184023333966 607.718584925431]
%Minimos / maximos Calinski-Harabasz: [1.25653924101928 10.2396304476605]
%Minimos / maximos Davies-Bouldin: [1.44860435789861 6.80350895581616]
%Minimos / maximos Silhouette: [-0.144379529139522 0.0267906725160614]
%Minimos / maximos Dunn's index: [0.135109703774302 0.467867373769756]
%Minimos / maximos accuracy: [0.30306129985229 0.545011078286558]
%Minimos / maximos Topographic Error: [0 0.913209010339734]
%config=4
%Minimos / maximos MQE: [422.698625395249 611.011219770725]
%Minimos / maximos Calinski-Harabasz: [1.23749892954768 10.7371314558007]
%Minimos / maximos Davies-Bouldin: [1.66813320764681 7.03248727239529]
%Minimos / maximos Silhouette: [-0.144253607272207 0.0124433622668958]
%Minimos / maximos Dunn's index: [0.134640433234182 0.4685594953078]
%Minimos / maximos accuracy: [0.303312407680945 0.561776218611521]
%Minimos / maximos Topographic Error: [1.10782865583456e-05 0.879113737075333]
%
%batches
%ylims = {[40.65 56.2], [1.04 2.6], [1.4 7.05], [-0.15 0.062], [0.13 0.56], [0.13 0.565], [0 0.9]};
%nobatches
%ylims = {[422 612], [1.04 10.8], [1.4 7.05], [-0.15 0.062], [0.13 0.56], [0.13 0.565], [0 0.9]};

%fprintf('config=%d\n', config);
%ylims = repmat({[inf, -inf]}, 7, 1);
for k=1:numel(args)
  ylims = showAllPerNAndmultiInstance(ylims, yticks, args{k}{:});
end
%for k=1:numel(names)
%  fprintf('Minimos / maximos %s: %s\n', names{k}, mat2str(ylims{k}));
%end


function ylims = showAllPerNAndmultiInstance(ylims, yticks, netType, prefixpathSingle, prefixpathMulti, prefixpathSave, typeReduction, expression, names, labnames, titlename)

oneFigureOnePlot = false;
oneFigurePerNAndNet = ~oneFigureOnePlot;

if ~(strcmp(netType, 'competitive') || strcmp(netType, 'som') || strcmp(netType, 'gng'))
  error(sprintf('netType must be som, competitive or gng'));
end
isCompetitive = strcmp(netType, 'competitive');
isSOM         = strcmp(netType, 'som');
isGNG         = strcmp(netType, 'gng');
isAll         = strcmp(typeReduction, 'all');


  load([prefixpathSingle 'experiments.mat']);
  experimentsS = experiments;
  load([prefixpathMulti 'experiments.mat']);
  experimentsM = experiments;
  clear experiments;

  numRepeats = experimentsS{1}.numRepeats;
  exp_show = numel(experimentsS);
  for k=1:numel(names)
    mn(k) = inf;
    mx(k) = -inf;
    for a = 1:exp_show
      v1 = getfield(experimentsS{a}, labnames{k});
      v2 = getfield(experimentsM{a}, labnames{k});
      mn(k) = min(mn(k), min(min(v1, [], 'all'), min(v2, [], 'all')));
      mx(k) = max(mx(k), max(max(v1, [], 'all'), max(v2, [], 'all')));
    end
  end
  mn(:)=0;
  mx(:)=0;
  if isCompetitive
    m = 'competitive';
    idxNet = 1;
  elseif isSOM
    m = 'SOM';
    idxNet = 2;
  elseif isGNG
    m = 'GNG';
    idxNet = 3;
  end
  for a=1:exp_show
    xvalues{a} = experimentsS{a}.factors;
  end
  system(sprintf('mkdir -p %s', prefixpathSave));
  for z1=1:1%numel(experimentsS{1}.HistorySizes)
    if oneFigurePerNAndNet
      if isCompetitive
        opos = [0 0 (0.5+0.5*0.5*0.62) 0.98];
      else
        opos = [0 0 0.5 1];
      end
      h = figure('Name', sprintf('%sN=%g, %s', expression, experimentsS{1}.HistorySizes(z1), m),...
                 'units','normalized','outerposition',opos);
      %tiledlayout(2,4,'TileSpacing','tight', 'Padding', 'tight');%'compact');
      %if isCompetitive
      %  tiledlayout(1,6,'TileSpacing','tight', 'Padding', 'normal');%'compact');
      %else
      %  tiledlayout(1,7,'TileSpacing','tight', 'Padding', 'normal');%'compact');
      %end
      tiledlayout(4,2,'TileSpacing','tight', 'Padding', 'tight');%'compact');
    end
    for z2=1:numel(names)
      for a=1:exp_show
        values1{a} = getfield(experimentsS{a}, labnames{z2});
        values2{a} = getfield(experimentsM{a}, labnames{z2});
      end
      if ~isCompetitive || z2<7 %~strcmp(labnames{z2}, 'TopographicError')
        if oneFigureOnePlot
          h = figure('Name', sprintf('%sN=%g, %s, %s', expression, experimentsS{1}.HistorySizes(z1), m, names{z2}));
        end
        if oneFigurePerNAndNet
          nexttile(z2);
        end
        putLegend = oneFigureOnePlot || z2==1;
        putXLabel = oneFigureOnePlot;
        if iscell(ylims{z2})
          mn(z2) = ylims{z2}{idxNet}(1);
          mx(z2) = ylims{z2}{idxNet}(2);
        else
          mn(z2) = ylims{z2}(1);
          mx(z2) = ylims{z2}(2);
        end
        if iscell(yticks{z2})
          ytck   = yticks{z2}{idxNet};
        else
          ytck   = yticks{z2};
        end
        %alsoConstant = z1==1;
        [mno, mxo] = doPlotSingle(mn(z2), mx(z2), ytck, xvalues, values1, values2, names{z2}, z1, isAll, putLegend, putXLabel, isCompetitive);%, alsoConstant);
        %ylims{z2}(1) = min(ylims{z2}(1), mno);
        %ylims{z2}(2) = max(ylims{z2}(2), mxo);
        if oneFigureOnePlot
          savefig(h, [prefixpathSave sprintf('singleplot%s_N%d_%s_measure_%s.fig', titlename, experimentsS{1}.HistorySizes(z1), m, names{z2})], 'compact');
          close(h);
        end
      end
    end
    if oneFigurePerNAndNet
      %savefig(h, [prefixpathSave sprintf('multiplot%s_N%d_%s.fig', titlename, experimentsS{1}.HistorySizes(z1), m)], 'compact');
      set(findobj(gcf,'type','axes'),'FontSize',16);%,'LineWidth', 2);
      exportgraphics(h, [prefixpathSave sprintf('multiplot%s_N%d_%s.pdf', titlename, experimentsS{1}.HistorySizes(z1), m)],'BackgroundColor','none', 'Resolution',50);
      close(h);
      %if isSOM && z1==1; print(jaja); end
    end
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mno, mxo] = doPlotSingle(mn, mx, ytck, xvalues, values1, values2, label, HistoryIndex, isAll, putLegend, putXLabel, isCompetitive)%, alsoConstant)

mno = inf;
mxo = -inf;
alsoM = true;

colores1 = {'r', 'r', 'r', 'k'};
colores2 = {'b', 'b', 'b', 'g'};
ls1 = {'-', '-.', '--', '-'};
ls2 = {'-', '-.', '--', ':'};
lw = 0.75;%1.5;%1;
ns1 = {'$\alpha(t)=aC(t)$ (S)', '$\alpha(t)=a\left(C\left(t\right)\right)^{2}$ (S)', '$\alpha(t)=a\left(C\left(t\right)\right)^{-1}$ (S)', '$\alpha(t)=a$'};
ns2 = {'$\alpha(t)=aC(t)$ (M)', '$\alpha(t)=a\left(C\left(t\right)\right)^{2}$ (M)', '$\alpha(t)=a\left(C\left(t\right)\right)^{-1}$ (M)', '$\alpha(t)=a$'};
%ns1 = {'y=ax (S)', 'y=ax^2 (S)', 'y=a/x (S)', 'y=a    '};
%ns2 = {'y=ax (M)', 'y=ax^2 (M)', 'y=a/x (M)', 'y=a (M)'};
if alsoM
  allnames = {ns1{1}, ns2{1}, ns1{2}, ns2{2}, ns1{3}, ns2{3}, ns1{4}};%, ns2{4}};
else
  allnames = {'$\alpha(t)=aC(t)$', '$\alpha(t)=a\left(C\left(t\right)\right)^{2}$', '$\alpha(t)=a\left(C\left(t\right)\right)^{-1}$', '$\alpha(t)=a$'};
end
hold on;
for k=1:numel(values1)
  if k==4
    ActualHistoryIndex = 1;
  else
    ActualHistoryIndex = HistoryIndex;
  end
  if isAll
    meanValues1 = reshape(nanmean(values1{k}(ActualHistoryIndex,:,:), 3), numel(xvalues{k}), 1);
    meanValues2 = reshape(nanmean(values2{k}(ActualHistoryIndex,:,:), 3), numel(xvalues{k}), 1);
  else
    meanValues1 = reshape(nanmean(values1{k}(ActualHistoryIndex,:,end,:), 4), numel(xvalues{k}), 1);
    meanValues2 = reshape(nanmean(values2{k}(ActualHistoryIndex,:,end,:), 4), numel(xvalues{k}), 1);
  end
  mno = min(mno, min(meanValues1, [], 'all'));
  mno = min(mno, min(meanValues2, [], 'all'));
  mxo = max(mxo, max(meanValues1, [], 'all'));
  mxo = max(mxo, max(meanValues2, [], 'all'));
  %if k<4 || alsoConstant
    plot(xvalues{k}, meanValues1, 'Color', colores1{k}, 'LineStyle', ls1{k}, 'LineWidth',lw);
  %end
  if alsoM && k<4
    plot(xvalues{k}, meanValues2, 'Color', colores2{k}, 'LineStyle', ls2{k}, 'LineWidth',lw);
  end
end
if mn<mx
  ylim([mn mx]);
end
if ~isempty(ytck)
  set(gca, 'YTick', ytck);
end
if putXLabel
  xlabel('a');
end
ylabel(label);
if putLegend
  l = legend(allnames, 'Location', 'layout', 'Interpreter', 'latex');
  if isCompetitive
    l.Layout.Tile = 'east';
  else
    l.Layout.Tile = 8;
  end
end
set(gca, 'XScale', 'log');
%set(gca, 'YScale', 'log');
xticks(10.^[-3 -2 -1 0]);


