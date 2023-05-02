function writeExperimentResultsBatches;

prefix = 'experiments_10batches_100exps/';
pathsSaveMatrix  = {[prefix '100repeats_per_config_summary_all_batches_by_measure.txt'];...
                   [prefix '100repeats_per_config_summary_last_batches_by_measure.txt']};
pathsSaveFunc    = {[prefix '100repeats_per_config_summary_all_batches_by_function.txt'];...
                   [prefix '100repeats_per_config_summary_last_batches_by_function.txt']};
types           = {'all';...
                   'last'};
postFixS = 'CORA_exploration_BATCH10_singleLR_N100_Sliding1_mqeMode0_MAXLR_1_SAM_4/experiments.mat';
postFixM = 'CORA_exploration_BATCH10_multiLR_N100_Sliding1_mqeMode0_MAXLR_1_SAM_4/experiments.mat';


pathSs = {...
  [prefix 'competitive/' postFixS]; ...
  [prefix 'som/' postFixS]; ...
  };%[prefix 'gng/' postFixS]};
pathMs = {...
  [prefix 'competitive/' postFixM]; ...
  [prefix 'som/' postFixM]; ...
  };%[prefix 'gng/' postFixM]};

doNet(pathSs, pathMs, pathsSaveMatrix, pathsSaveFunc, types);


function doNet(pathSs, pathMs, pathsSaveMatrix, pathsSaveFunc, types);

labnames    = {'meanQuantizationErrors', 'CalinskiHarabasz', 'DaviesBouldin', 'Silhouette', 'TopographicError', 'DunnIndex', 'accuracy'};
%labnames    = {'all_datasets_meanQuantizationErrors', 'all_datasets_CalinskiHarabasz', 'all_datasets_DaviesBouldin', 'all_datasets_Silhouette', 'all_datasets_TopographicError', 'all_datasets_DunnIndex', 'all_datasets_accuracy'};
besterfuncs = {@min,                     @max,                @min,            @max,          @min,              @max,        @max};
initbest    = [inf,                      -inf,                inf,             -inf,          inf,               -inf,        -inf];
directions  = {'ascend',                 'descend',           'ascend',        'descend',     'ascend',          'descend',   'descend'};

funcnames = {'linear', 'quadratic', 'inverse', 'constant'};

histSizes = [10, 50, 100];

ns1 = {'y=ax (S)  ', 'y=ax^2 (S)', 'y=a/x (S) ', 'y=a (S)   '};
ns2 = {'y=ax (M)  ', 'y=ax^2 (M)', 'y=a/x (M) ', 'y=a (M)   '};

all_funcs   = {ns1{1}, ns2{1}, ns1{2}, ns2{2}, ns1{3}, ns2{3}, ns1{4}, ns2{4}};

for p=1:numel(types)

  isAll = strcmp(types{p}, 'all');
  f = fopen(pathsSaveMatrix{p}, 'w');
  configs_won = cell(size(all_funcs));
  for j=1:numel(configs_won)
    configs_won{j} = cell(0,1);
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
        if strcmp(experimentsS{1}.netType, 'competitive') && strcmp(labnames{m}, 'TopographicError')
          continue;
        end
        fprintf(f, '    ####################\n    results for measure %s\n    ####################\n\n', labnames{m});

        if isAll
          fprintf(f, '      ####################\n      best means across all batches\n      ####################\n\n');
        else
          fprintf(f, '      ####################\n      best means for last batches\n      ####################\n\n');
        end
        mnmn = initbest(m);bestname='';
        vals = zeros(size(all_funcs));
        idxfun = 1;
        for k=1:4
          valsS = getfield(experimentsS{k}, labnames{m});
          valsM = getfield(experimentsM{k}, labnames{m});
          if isAll
            means = reshape(nanmean(valsS(h,:,:), 3), numel(experimentsS{k}.factors), 1);
          else
            means = reshape(nanmean(valsS(h,:,end,:), 4), numel(experimentsS{k}.factors), 1);
          end
          [mn, idx] = besterfuncs{m}(means);
          vals(idxfun) = mn; idxfun=idxfun+1;
          if besterfuncs{m}(mn, mnmn)==mn
            mnmn = mn;
            bestname = ns1{k};
          end
          fprintf(f, '        Best value for %s is %f at a=%f\n\n', ns1{k}, mn, experimentsS{k}.factors(idx));
          if isAll
            means = reshape(nanmean(valsM(h,:,:), 3), numel(experimentsS{k}.factors), 1);
          else
            means = reshape(nanmean(valsM(h,:,end,:), 4), numel(experimentsS{k}.factors), 1);
          end
          [mn, idx] = besterfuncs{m}(means);
          vals(idxfun) = mn; idxfun=idxfun+1;
          if besterfuncs{m}(mn, mnmn)==mn
            mnmn = mn;
            bestname = ns2{k};
          end
          fprintf(f, '        Best value for %s is %f at a=%f\n\n', ns2{k}, mn, experimentsS{k}.factors(idx));
        end
        if isAll
          fprintf(f, '          ####################\n          best function across all batches: %s\n          ####################\n\n', bestname);
        else
          fprintf(f, '          ####################\n          best function for last batches: %s\n          ####################\n\n', bestname);
        end
        for j=1:numel(all_funcs)
          if strcmp(all_funcs{j}, bestname)
            [ordered, idxord] = sort(vals, directions{m});
            percentages = ordered/mnmn*100;
            strs = cell(numel(ordered)-1,1);
            for jj=2:numel(ordered)
              strs{jj-1} = sprintf('\n    %s was %f: %f%% of %s', all_funcs{idxord(jj)}, ordered(jj), percentages(jj), bestname);
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

  for j=1:numel(all_funcs)
    fprintf(f, '####################\nRESULTS FOR FUNCTION %s: it won %d times!\n####################\n\n', all_funcs{j}, numel(configs_won{j}));
    for jj=1:numel(configs_won{j})
      fprintf(f, '  %s', configs_won{j}{jj});
    end
  end

  fclose(f);

end
