function results = getBestCORAResults%(config)

config=3;

if config==3

  prefix = 'experiments_nobatches_100exps/';
  pathsSummary    = {[prefix 'singlePlot/bests.txt']};
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

  results = doNet(pathSs, pathMs, types, expressions, labnames, initbest, pathsSummary);
end




function results = doNet(pathSs, pathMs, types, expressions, labnames, initbest, pathsSummary);

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

ns1 = {'y=ax (S)  ', 'y=ax^2 (S)', 'y=a/x (S) ', 'y=a'};
ns2 = {'y=ax (M)  ', 'y=ax^2 (M)', 'y=a/x (M) ', 'y=a (M)   '};
all_funcs = {ns1{1}, ns2{1}, ns1{2}, ns2{2}, ns1{3}, ns2{3}, ns1{4}};%, ns2{4}};

netnames = {'competitive', 'som', 'gng'};

results = cell(size(types));

for p=1:numel(types)

  best_values = zeros(size(labnames));
  best_values(:) = initbest(:);
  best_config = cell(size(labnames));
  isAll = strcmp(types{p}, 'all');
  for h = 1:numel(histSizes)

    for n = 1:numel(pathSs)
      load(pathSs{n});
      experimentsS = experiments;
      load(pathMs{n});
      experimentsM = experiments;
      clear experiments;
      %fprintf(f, '  ####################\n  RESULTS FOR NET %s\n  ####################\n\n', experimentsS{1}.netType);
      for m = 1:numel(labnames)
        if strcmp(experimentsS{1}.netType, 'competitive') && strcmp(labnames{m}(max(1, end-(numel('TopographicError')-1)):end), 'TopographicError')
          continue;
        end
        %fprintf(f, '    ####################\n    results for measure %s\n    ####################\n\n', labnames{m});
        
        %fprintf(f, '      ####################\n      best means %s\n      ####################\n\n', expressions{p});
        for k=1:4
          valsS = getfield(experimentsS{k}, labnames{m});
          valsM = getfield(experimentsM{k}, labnames{m});
          %if useM && useS
            % in this case, we want to collate results from all S simulations, and all M simulations except y=a
            doWithS = (k<4 || h==1);
            doWithM =  k<4;
          %elseif (useM && ~useS) || (~useM && useS)
          %  % in this case, we want to collate results from either all S simulations, or all M simulations and S simulation for y=a  (to always compare to the same set of y=a simulations)
          %  doWithS = (useS && (k<4 || h==1)) || (useM && k==4 && h==1);
          %  doWithM = useM && k<4;
          %end
          if doWithS
            if isAll
              means = reshape(nanmean(valsS(h,:,:), 3), numel(experimentsS{k}.factors), 1);
            else
              means = reshape(nanmean(valsS(h,:,end,:), 4), numel(experimentsS{k}.factors), 1);
            end
            [mn, idx] = besterfuncs{m}(means);
            if besterfuncs{m}(mn, best_values(m))==mn
              best_values(m) = mn;
              best_config{m} = {labnames{m}, mn, netnames{n}, histSizes(h), 'single', ns1{k}, experimentsS{k}.factors(idx)};
            end
            %fprintf(f, '        Best value for %s is %f at a=%f\n\n', ns1{k}, mn, experimentsS{k}.factors(idx));
          end
          if doWithM
            if isAll
              means = reshape(nanmean(valsM(h,:,:), 3), numel(experimentsS{k}.factors), 1);
            else
              means = reshape(nanmean(valsM(h,:,end,:), 4), numel(experimentsS{k}.factors), 1);
            end
            [mn, idx] = besterfuncs{m}(means);
            if besterfuncs{m}(mn, best_values(m))==mn
              best_values(m) = mn;
              best_config{m} = {labnames{m}, mn, netnames{n}, histSizes(h), 'multi', ns2{k}, experimentsS{k}.factors(idx)};
            end
            %fprintf(f, '        Best value for %s is %f at a=%f\n\n', ns2{k}, mn, experimentsS{k}.factors(idx));
          end
        end
      end
    end
  end

  results{p} = best_config;
  fsum = fopen(pathsSummary{p}, 'w');
  for m = 1:numel(labnames)
    fprintf(fsum, '%s, %f, %s, %d, %s, %s, %f\n', best_config{m}{:});
  end
  fclose(fsum);

end


