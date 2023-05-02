function writeExperimentResults(config);

%writeExperimentResults(1); writeExperimentResults(2); writeExperimentResults(3);

%config=0;

if config==0;
  prefix = 'experiments_long_ranges/';
  pathsSaveMatrix = {[prefix 'singlePlot2/10repeats_per_config_summary_by_measure.txt']};
  pathsSaveFunc   = {[prefix 'singlePlot2/10repeats_per_config_summary_by_function.txt']};
  types           = {'all'};
  expressions     = {'across these experiments'};

  postFixS = 'CORA_exploration_singleLR_N10_Sliding1_mqeMode0_MAXLR_1_SAM_10832/experiments.mat';
  postFixM = 'CORA_exploration_multiLR_N10_Sliding1_mqeMode0_MAXLR_1_SAM_10832/experiments.mat';

  pathSs = {...
    [prefix 'competitive/' postFixS]; ...
    [prefix 'som/' postFixS]; ...
    [prefix 'gng/' postFixS]};
  pathMs = {...
    [prefix 'competitive/' postFixM]; ...
    [prefix 'som/' postFixM]; ...
    [prefix 'gng/' postFixM]};

  labnames    = {'meanQuantizationErrors', 'CalinskiHarabasz', 'DaviesBouldin', 'Silhouette', 'TopographicError', 'DunnIndex', 'accuracy'};
  initbest    = [inf,                      -inf,                inf,             -inf,          inf,               -inf,        -inf];
  labnames    = {'meanQuantizationErrors', 'CalinskiHarabasz', 'DaviesBouldin', 'Silhouette', 'TopographicError', 'DunnIndex'};
  initbest    = [inf,                      -inf,                inf,             -inf,          inf,               -inf,     ];

  doit = true;

elseif config==1

  prefix = 'experiments_10batches_100exps/';
  pathsSaveMatrix  = {[prefix 'singlePlot/summary_100repeats_evaluate_all_batches_each_one_with_itself_results%s%s_by_measure.txt'];...
                      [prefix 'singlePlot/summary_100repeats_evaluate_last_batches_each_one_with_itself_results%s%s_by_measure.txt']};
  pathsSaveFunc    = {[prefix 'singlePlot/summary_100repeats_evaluate_all_batches_each_one_with_itself_results%s%s_by_function.txt'];...
                      [prefix 'singlePlot/summary_100repeats_evaluate_last_batches_each_one_with_itself_results%s%s_by_function.txt']};
  pathsSummary     = {[prefix 'singlePlot/summary_all_batches_bybatch.txt'];...
                      [prefix 'singlePlot/summary_last_batches_bybatch.txt']};
  types           = {'all';...
                     'last'};
  expressions     = {'across all batches';...
                     'for last batches'};

  postFixS = 'CORA_exploration_BATCH10_singleLR_N100_Sliding1_mqeMode0_MAXLR_1_SAM_4/experiments.mat';
  postFixM = 'CORA_exploration_BATCH10_multiLR_N100_Sliding1_mqeMode0_MAXLR_1_SAM_4/experiments.mat';

  pathSs = {...
    [prefix 'competitive/' postFixS]; ...
    [prefix 'som/' postFixS]; ...
    [prefix 'gng/' postFixS]};
  pathMs = {...
    [prefix 'competitive/' postFixM]; ...
    [prefix 'som/' postFixM]; ...
    [prefix 'gng/' postFixM]};

  labnames    = {'accuracy', 'CalinskiHarabasz', 'Silhouette', 'DunnIndex', 'meanQuantizationErrors', 'DaviesBouldin', 'TopographicError'};
  initbest    = [-inf,       -inf,                -inf,         -inf,       inf,                       inf,             inf];

  doit = false;
  indexes = {1, 2, 3, 4, 5, 6, 7, [1, 2, 3, 4, 5, 6, 7]};
  for k=1:numel(indexes)
    if numel(indexes{k})==1
      pathsSaveMatrix2 = {sprintf(pathsSaveMatrix{1}, '_', labnames{indexes{k}}), sprintf(pathsSaveMatrix{2}, '_', labnames{indexes{k}})};
      pathsSaveFunc2   = {sprintf(pathsSaveFunc{1},   '_', labnames{indexes{k}}), sprintf(pathsSaveFunc{2},   '_', labnames{indexes{k}})};
      %doNet(pathSs, pathMs, pathsSaveMatrix2, pathsSaveFunc2, types, expressions, labnames(indexes{k}), initbest(indexes{k}), pathsSummary, false, true);
      %doNet(pathSs, pathMs, pathsSaveMatrix2, pathsSaveFunc2, types, expressions, labnames(indexes{k}), initbest(indexes{k}), pathsSummary, true, false);
      doNet(pathSs, pathMs, pathsSaveMatrix2, pathsSaveFunc2, types, expressions, labnames(indexes{k}), initbest(indexes{k}), pathsSummary, true, true);
    else
      pathsSaveMatrix2 = {sprintf(pathsSaveMatrix{1}, '', ''), sprintf(pathsSaveMatrix{2}, '', '')};
      pathsSaveFunc2   = {sprintf(pathsSaveFunc{1},   '', ''), sprintf(pathsSaveFunc{2},   '', '')};
      %doNet(pathSs, pathMs, pathsSaveMatrix2, pathsSaveFunc2, types, expressions, labnames(indexes{k}), initbest(indexes{k}), pathsSummary, false, true);
      %doNet(pathSs, pathMs, pathsSaveMatrix2, pathsSaveFunc2, types, expressions, labnames(indexes{k}), initbest(indexes{k}), pathsSummary, true, false);
      doNet(pathSs, pathMs, pathsSaveMatrix2, pathsSaveFunc2, types, expressions, labnames(indexes{k}), initbest(indexes{k}), pathsSummary, true, true);
    end
  end

elseif config==2

  prefix = 'experiments_10batches_100exps/';
  pathsSaveMatrix  = {[prefix 'singlePlot/summary_100repeats_evaluate_last_batches_with_whole_dataset_results%s%s_by_measure.txt']};
  pathsSaveFunc    = {[prefix 'singlePlot/summary_100repeats_evaluate_last_batches_with_whole_dataset_results%s%s_by_function.txt']};
  pathsSummary     = {[prefix 'singlePlot/summary_last_batches_whole.txt']};
  types           = {'all'};
  expressions     = {'across all batches (for all samples in dataset)'};

  postFixS = 'CORA_exploration_BATCH10_singleLR_N100_Sliding1_mqeMode0_MAXLR_1_SAM_4/experiments.mat';
  postFixM = 'CORA_exploration_BATCH10_multiLR_N100_Sliding1_mqeMode0_MAXLR_1_SAM_4/experiments.mat';

  pathSs = {...
    [prefix 'competitive/' postFixS]; ...
    [prefix 'som/' postFixS]; ...
    [prefix 'gng/' postFixS]};
  pathMs = {...
    [prefix 'competitive/' postFixM]; ...
    [prefix 'som/' postFixM]; ...
    [prefix 'gng/' postFixM]};

  labnames    = {'alldataset_accuracy', 'alldataset_CalinskiHarabasz', 'alldataset_Silhouette', 'alldataset_DunnIndex', 'alldataset_meanQuantizationErrors', 'alldataset_DaviesBouldin', 'alldataset_TopographicError'};
  initbest    = [-inf,       -inf,                -inf,         -inf,       inf,                       inf,             inf];

  doit = false;
  indexes = {1, 2, 3, 4, 5, 6, 7, [1, 2, 3, 4, 5, 6, 7]};
  for k=1:numel(indexes)
    if numel(indexes{k})==1
      pathsSaveMatrix2 = {sprintf(pathsSaveMatrix{1}, '_', labnames{indexes{k}}(12:end))};
      pathsSaveFunc2   = {sprintf(pathsSaveFunc{1},   '_', labnames{indexes{k}}(12:end))};
      %doNet(pathSs, pathMs, pathsSaveMatrix2, pathsSaveFunc2, types, expressions, labnames(indexes{k}), initbest(indexes{k}), pathsSummary, false, true);
      %doNet(pathSs, pathMs, pathsSaveMatrix2, pathsSaveFunc2, types, expressions, labnames(indexes{k}), initbest(indexes{k}), pathsSummary, true, false);
      doNet(pathSs, pathMs, pathsSaveMatrix2, pathsSaveFunc2, types, expressions, labnames(indexes{k}), initbest(indexes{k}), pathsSummary, true, true);
    else
      pathsSaveMatrix2 = {sprintf(pathsSaveMatrix{1}, '', '')};
      pathsSaveFunc2   = {sprintf(pathsSaveFunc{1},   '', '')};
      %doNet(pathSs, pathMs, pathsSaveMatrix2, pathsSaveFunc2, types, expressions, labnames(indexes{k}), initbest(indexes{k}), pathsSummary, false, true);
      %doNet(pathSs, pathMs, pathsSaveMatrix2, pathsSaveFunc2, types, expressions, labnames(indexes{k}), initbest(indexes{k}), pathsSummary, true, false);
      doNet(pathSs, pathMs, pathsSaveMatrix2, pathsSaveFunc2, types, expressions, labnames(indexes{k}), initbest(indexes{k}), pathsSummary, true, true);
    end
  end

elseif config==3

  prefix = 'experiments_nobatches_100exps/';
  pathsSaveMatrix = {[prefix 'singlePlot/summary_100repeats_without_batches_results%s%s_by_measure_again.txt']};
  pathsSaveFunc   = {[prefix 'singlePlot/summary_100repeats_without_batches_results%s%s_by_function_again.txt']};
  pathsSummary    = {[prefix 'singlePlot/summary.txt']};
  types           = {'all'};
  expressions     = {'across these experiments'};

  postFixS = 'CORA_exploration_singleLR_N100_Sliding1_mqeMode0_MAXLR_1_SAM_10832/experiments.mat';
  postFixM = 'CORA_exploration_multiLR_N100_Sliding1_mqeMode0_MAXLR_1_SAM_10832/experiments.mat';

  pathSs = {...
    [prefix 'competitive/' postFixS]; ...
    [prefix 'som/' postFixS]; ...
    [prefix 'gng/' postFixS]};
  pathMs = {...
    [prefix 'competitive/' postFixM]; ...
    [prefix 'som/' postFixM]; ...
    [prefix 'gng/' postFixM]};

  labnames    = {'accuracy', 'CalinskiHarabasz', 'Silhouette', 'DunnIndex', 'meanQuantizationErrors', 'DaviesBouldin', 'TopographicError'};
  initbest    = [-inf,       -inf,                -inf,         -inf,       inf,                       inf,             inf];

  doit = false;
  indexes = {1, 2, 3, 4, 5, 6, 7, [1, 2, 3, 4, 5, 6, 7]};
  for k=1:numel(indexes)
    if numel(indexes{k})==1
      pathsSaveMatrix2 = {sprintf(pathsSaveMatrix{1}, '_', labnames{k})};
      pathsSaveFunc2   = {sprintf(pathsSaveFunc{1},   '_', labnames{k})};
      %doNet(pathSs, pathMs, pathsSaveMatrix2, pathsSaveFunc2, types, expressions, labnames(indexes{k}), initbest(indexes{k}), pathsSummary, false, true);
      %doNet(pathSs, pathMs, pathsSaveMatrix2, pathsSaveFunc2, types, expressions, labnames(indexes{k}), initbest(indexes{k}), pathsSummary, true, false);
      doNet(pathSs, pathMs, pathsSaveMatrix2, pathsSaveFunc2, types, expressions, labnames(indexes{k}), initbest(indexes{k}), pathsSummary, true, true);
    else
      pathsSaveMatrix2 = {sprintf(pathsSaveMatrix{1}, '', '')};
      pathsSaveFunc2   = {sprintf(pathsSaveFunc{1},   '', '')};
      %doNet(pathSs, pathMs, pathsSaveMatrix2, pathsSaveFunc2, types, expressions, labnames(indexes{k}), initbest(indexes{k}), pathsSummary, false, true);
      %doNet(pathSs, pathMs, pathsSaveMatrix2, pathsSaveFunc2, types, expressions, labnames(indexes{k}), initbest(indexes{k}), pathsSummary, true, false);
      doNet(pathSs, pathMs, pathsSaveMatrix2, pathsSaveFunc2, types, expressions, labnames(indexes{k}), initbest(indexes{k}), pathsSummary, true, true);
    end
  end

end

if doit
  doNet(pathSs, pathMs, pathsSaveMatrix, pathsSaveFunc, types, expressions, labnames, initbest);
end



function doNet(pathSs, pathMs, pathsSaveMatrix, pathsSaveFunc, types, expressions, labnames, initbest, pathsSummary, useS, useM);

besterfuncs = cell(size(initbest));
directions  = cell(size(initbest));
for k=1:numel(initbest)
  if initbest(k)==inf
    besterfuncs{k} = @min;
    directions{k}  = 'ascend';
  else
    besterfuncs{k} = @max;
    directions{k}  = 'descend';
  end
end

histSizes = [10, 50, 100];

ns1 = {'y=ax (S)  ', 'y=ax^2 (S)', 'y=a/x (S) ', 'y=a (S)   '};
ns2 = {'y=ax (M)  ', 'y=ax^2 (M)', 'y=a/x (M) ', 'y=a (M)   '};
if useS && useM
  all_funcs = {ns1{1}, ns2{1}, ns1{2}, ns2{2}, ns1{3}, ns2{3}, ns1{4}};%, ns2{4}};
elseif useS
  %ns1 = {'y=ax  ', 'y=ax^2', 'y=a/x ', 'y=a   '};
  ns1{end} = 'y=a       ';
  ns2{end} = 'y=a       ';
  all_funcs = ns1;
  for k=1:numel(pathsSaveMatrix)
    pathsSaveMatrix{k} = [pathsSaveMatrix{k}(1:end-4) '_just_single.txt'];
    pathsSaveFunc{k}   = [  pathsSaveFunc{k}(1:end-4) '_just_single.txt'];
    pathsSummary{k}    = [   pathsSummary{k}(1:end-4) '_just_single.txt'];
  end
elseif useM
  %ns1 = {'y=ax  ', 'y=ax^2', 'y=a/x ', 'y=a   '};
  ns1{end} = 'y=a       ';
  ns2{end} = 'y=a       ';
  all_funcs = ns2;
  for k=1:numel(pathsSaveMatrix)
    pathsSaveMatrix{k} = [pathsSaveMatrix{k}(1:end-4) '_just_multi.txt'];
    pathsSaveFunc{k}   = [  pathsSaveFunc{k}(1:end-4) '_just_multi.txt'];
    pathsSummary{k}    = [   pathsSummary{k}(1:end-4) '_just_multi.txt'];
  end
end

for p=1:numel(types)

  isAll = strcmp(types{p}, 'all');
  f = fopen(pathsSaveMatrix{p}, 'w');
  configs_won = cell(size(all_funcs));
  for j=1:numel(configs_won)
    configs_won{j} = cell(0,1);
  end
  if numel(labnames)==1
    fsum = fopen(pathsSummary{p}, 'a');
  end

  for h = 1:numel(histSizes)
    fprintf(f, '####################\nRESULTS FOR WINDOW SIZE %f\n####################\n\n', histSizes(h));

    for n = 1:numel(pathSs)
      load(pathSs{n});
      experimentsS = experiments;
      load(pathMs{n});
      experimentsM = experiments;
      clear experiments;
      fprintf(f, '  ####################\n  RESULTS FOR NET %s\n  ####################\n\n', experimentsS{1}.netType);
      for m = 1:numel(labnames)
        if strcmp(experimentsS{1}.netType, 'competitive') && strcmp(labnames{m}(max(1, end-(numel('TopographicError')-1)):end), 'TopographicError')
          continue;
        end
        fprintf(f, '    ####################\n    results for measure %s\n    ####################\n\n', labnames{m});
        
        fprintf(f, '      ####################\n      best means %s\n      ####################\n\n', expressions{p});
        mnmn = initbest(m);bestname='';
        vals = zeros(0,1);%zeros(size(all_funcs));
        names = cell(0,1);
        for k=1:4
          valsS = getfield(experimentsS{k}, labnames{m});
          valsM = getfield(experimentsM{k}, labnames{m});
          if useM && useS
            % in this case, we want to collate results from all S simulations, and all M simulations except y=a
            doWithS = useS && (k<4 || h==1);
            doWithM = useM && k<4;
          elseif (useM && ~useS) || (~useM && useS)
            % in this case, we want to collate results from either all S simulations, or all M simulations and S simulation for y=a  (to always compare to the same set of y=a simulations)
            doWithS = (useS && (k<4 || h==1)) || (useM && k==4 && h==1);
            doWithM = useM && k<4;
          end
          if doWithS
            if isAll
              means = reshape(nanmean(valsS(h,:,:), 3), numel(experimentsS{k}.factors), 1);
            else
              means = reshape(nanmean(valsS(h,:,end,:), 4), numel(experimentsS{k}.factors), 1);
            end
            [mn, idx] = besterfuncs{m}(means);
            vals(end+1,1)  = mn;
            names{end+1,1} = ns1{k};
            if besterfuncs{m}(mn, mnmn)==mn
              mnmn = mn;
              bestname = ns1{k};
            end
            fprintf(f, '        Best value for %s is %f at a=%f\n\n', ns1{k}, mn, experimentsS{k}.factors(idx));
          end
          if doWithM
            if isAll
              means = reshape(nanmean(valsM(h,:,:), 3), numel(experimentsS{k}.factors), 1);
            else
              means = reshape(nanmean(valsM(h,:,end,:), 4), numel(experimentsS{k}.factors), 1);
            end
            [mn, idx] = besterfuncs{m}(means);
            vals(end+1,1)  = mn;
            names{end+1,1} = ns2{k};
            if besterfuncs{m}(mn, mnmn)==mn
              mnmn = mn;
              bestname = ns2{k};
            end
            fprintf(f, '        Best value for %s is %f at a=%f\n\n', ns2{k}, mn, experimentsS{k}.factors(idx));
          end
        end
        fprintf(f, '          ####################\n          best function %s: %s\n          ####################\n\n', expressions{p}, bestname);
        for j=1:numel(all_funcs)
          if strcmp(all_funcs{j}, bestname)
            [ordered, idxord] = sort(vals, directions{m});
            percentages = ordered/mnmn*100;
            strs = cell(numel(ordered)-1,1);
            for jj=2:numel(ordered)
              %if ~strcmp(names{idxord(jj)}, ns1{end}) || h==1
                strs{jj-1} = sprintf('\n    %s was %f: %f%% of %s', names{idxord(jj)}, ordered(jj), percentages(jj), bestname);
              %end
            end
            str = [strs{:}];
            configs_won{j}{end+1,1} = sprintf('For window size %f, net %s and measure %s, %s won with %f. other values: %s\n\n', histSizes(h), experimentsS{1}.netType, labnames{m}, bestname, mnmn, str);
            break;
          end
        end
      end
    end
  end

  fclose(f);

  f = fopen(pathsSaveFunc{p}, 'w');
  if numel(labnames)==1
    fprintf(fsum, '\nFor % 25s: ', labnames{1});
    %fprintf(fsum, '\n', labnames{1});
  end

  for j=[1 3 5 2 4 6 7]%1:numel(all_funcs)
    fprintf(f, '####################\nRESULTS FOR FUNCTION %s: it won %d times!\n####################\n\n', all_funcs{j}, numel(configs_won{j}));
    if numel(labnames)==1
      %fprintf(fsum, '%s won %d times, ', all_funcs{j}, numel(configs_won{j}));
      fprintf(fsum, '%d ', numel(configs_won{j}));
    end
    for jj=1:numel(configs_won{j})
      fprintf(f, '  %s', configs_won{j}{jj});
    end
  end

  fclose(f);

end
