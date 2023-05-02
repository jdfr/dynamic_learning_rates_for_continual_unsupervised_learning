function results = searchGoodEnough;

kmeans = [0.352 4.121 -0.0724 0.206 431.4 2.107 inf]';
kmedoids = [0.425 8.109 -0.032 0.136 440.2 4.959 inf]';
dbscan = [0.307 1.002 -0.225 0.097 438.9 2.996 inf]';
competitive = [0.478 7.516 0.0084 0.160 427.9 3.886 inf]';
som = [0.517 9.089 -0.030 0.162 425.4 5.664 0.716]';
gng = [0.471 7.513 -0.085 0.148 439.3 6.985 0.122]';

morebetter = [true, true, true, true, false, false, false]';

competitors = [kmeans kmedoids dbscan competitive som gng];

bestvalscompetitors = zeros(size(morebetter));
for k=1:numel(morebetter)
  if morebetter(k)
    bestvalscompetitors(k) = max(competitors(k,:));
  else
    bestvalscompetitors(k) = min(competitors(k,:));
  end
end

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

  results = doNet(pathSs, pathMs, types, expressions, labnames, initbest, pathsSummary, bestvalscompetitors, competitors);
end

end




function results = doNet(pathSs, pathMs, types, expressions, labnames, initbest, pathsSummary, bestvalscompetitors, competitors);

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

numfactors = 40;

for p=1:numel(types)

  %dimensions: history sizes (3), network types (3), update type (2), function type (7), numfactors (40), metrics (7)
  all_values = zeros(numel(histSizes), numel(pathSs), numel(ns1)*2-1, numfactors, numel(labnames));
  all_is_inverse = zeros(numel(histSizes), numel(pathSs), numel(ns1)*2-1, numfactors);
  all_is_inverse(:)=nan;
  all_values(:) = nan;
  all_stds = zeros(size(all_values));
  all_stds(:) = nan;
  %all_dimensions = {histSizes, pathSs, [ns1 ns2(1:end-1)], labnames};
  all_configs = cell(numel(histSizes), numel(pathSs), numel(ns1)*2-1, numfactors);
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
          % in this case, we want to collate results from all S simulations, and all M simulations except y=a
          doWithS = (k<4 || h==1);
          doWithM =  k<4;
          if doWithS
            if isAll
              means = reshape(nanmean(valsS(h,:,:), 3), numel(experimentsS{k}.factors), 1);
              stds  = reshape(nanstd(valsS(h,:,:),0,3), numel(experimentsS{k}.factors), 1);
            else
              means = reshape(nanmean(valsS(h,:,end,:), 4), numel(experimentsS{k}.factors), 1);
              stds  = reshape(nanstd(valsS(h,:,end,:),0,4), numel(experimentsS{k}.factors), 1);
            end
            all_values(h,n,k,1:numel(means),m) = means;
            all_stds(  h,n,k,1:numel(means),m) = stds;
            if m==1
              all_is_inverse(h,n,k,1:numel(means)) = k==3;
              for kk=1:numel(experimentsS{k}.factors)
                all_configs{h,n,k,kk} = {netnames{n}, histSizes(h), 'single', ns1{k}, experimentsS{k}.factors(kk)};
              end
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
            all_values(h,n,k+4,1:numel(means),m) = means;
            all_stds(  h,n,k+4,1:numel(means),m) = stds;
            if m==1
              all_is_inverse(h,n,k+4,1:numel(means)) = k==3;
              for kk=1:numel(experimentsS{k}.factors)
                all_configs{h,n,k+4,kk} = {netnames{n}, histSizes(h), 'multi', ns2{k}, experimentsS{k}.factors(kk)};
              end
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

  all_diffs = all_values-reshape(bestvalscompetitors, [1, 1, 1, 1, numel(labnames)]);
  all_diffs(:,:,:,:,5:7) = -all_diffs(:,:,:,:,5:7);
  num_positives = nansum(all_diffs>0, 5);
  max_positives = max(num_positives(:));
  num_max = sum(num_positives(:)==max_positives);
  best_stds   = zeros(numel(labnames), num_max);
  best_ranks  = zeros(numel(labnames), num_max);
  best_values = zeros(numel(labnames), num_max);
  best_diffs  = zeros(numel(labnames), num_max);
  best_is_inverse = zeros(num_max,1);
  bests = cell(num_max,5);
  bests2 = cell(num_max,1);
  bb=1;
  all_ranks = zeros(size(all_values));
  all_ranks(:) = nan;
  all_autoranks = zeros(size(all_values));
  all_autoranks(:) = nan;
  competitors_by_metric = cell(numel(labnames),1);
  for h=1:numel(labnames)
    competitors_by_metric{h} = competitors(h,~isinf(competitors(h,:)));
  end
  %competitors=[kmeans kmedoids dbscan competitive som gng];

  
  %funcnames = {ns1{1:4}, ns2{1:3}};
  funcnames = cell(7,1);
  meanvalues_byfunction = zeros(7, numel(labnames));
  for k=1:7
    funcnames{k} = all_configs{1,1,k,1}{4};
    for gg=1:numel(labnames)
      values = all_values(:, :, k, :, gg);
      values = values(:);
      values = values(~isnan(values));
      value = mean(values);
      meanvalues_byfunction(k,gg) = value;
    end
  end
  

  for h = 1:numel(histSizes)
    for n = 1:numel(pathSs)
      for k=1:7
        for kk=1:numfactors
          for gg=1:numel(labnames)
            if ~isnan(all_values(h, n, k, kk, gg))
              values = competitors_by_metric{gg};
              values(end+1) = all_values(h, n, k, kk, gg);
              [ordered, idxord] = sort(values, directions{gg});
              idx = min(find(ordered==all_values(h, n, k, kk, gg)));
              all_ranks(h, n, k, kk, gg) = idx;
            end
          end
          if all_diffs(h,n,k,kk,1)>0 && all_diffs(h,n,k,kk,5)>0  %num_positives(h,n,k,kk)==max_positives
            best_values(:,bb) =  all_values(h,n,k,kk,:);
            best_stds(  :,bb) =    all_stds(h,n,k,kk,:);
            best_ranks(:,bb)  =   all_ranks(h,n,k,kk,:);
            best_diffs(:,bb)  =   all_diffs(h,n,k,kk,:);
            bests(bb,:)       = all_configs{h,n,k,kk};
            bests2{bb}        = all_configs{h,n,k,kk};
            best_is_inverse(bb) = all_is_inverse(h,n,k,kk);
            %fprintf('Config: %s %d %s %s %d\n', bests{bb,:});
            bb=bb+1;
          end
        end
      end
    end
  end
  sorted_all_values = reshape(all_values, [], 7);
  sorted_all_values = sorted_all_values(~isnan(sorted_all_values(:,1)),:);
  for h=1:7
    sorted_all_values(:,h) = sort(sorted_all_values(:,h), directions{h});
  end
  for h = 1:numel(histSizes)
    for n = 1:numel(pathSs)
      for k=1:7
        for kk=1:numfactors
          for gg=1:numel(labnames)
            if ~isnan(all_values(h, n, k, kk, gg))
              idx = min(find(sorted_all_values(:,gg)==all_values(h, n, k, kk, gg)));
              all_autoranks(h, n, k, kk, gg) = idx;
            end
          end
        end
      end
    end
  end
  these_autoranks = reshape(all_autoranks,[],7);
  these_inverse   = all_is_inverse(:);
  sum_autoranks = nansum(these_autoranks, 2);
  sum_autoranks(sum_autoranks==0) = inf;
  first_quartile = quantile(sum_autoranks(~isinf(sum_autoranks)), 0.25);
  zz = (these_inverse(sum_autoranks<first_quartile));
  [numel(zz)/numel(sum_autoranks(~isinf(sum_autoranks))), sum(zz)/numel(zz)]

  if false
  %here, we compute configs in the pareto front for all configs being better than competitors for 4 metrics in absolute terms
  pareto_values = best_values';
  pareto_diffs  = best_diffs';
  pareto_ranks  = best_ranks';
  pareto_bests  = bests2;
  pareto_inverse = best_is_inverse;
  pareto_values(:,5:7) = -pareto_values(:,5:7);
  [membership, member_value]=find_pareto_frontier(pareto_values);%(:,[3,4,6]));
  pareto_values = best_values(:,membership)';
  pareto_diffs  = best_diffs(:,membership)';
  pareto_ranks  = best_ranks(:,membership)';
  pareto_bests  = bests(membership,:);
  pareto_bests2 = bests2(membership);
  pareto_inverse = best_is_inverse(membership);
  end

  if true
  %here, we compute configs in the pareto front for all configs being better than competitors with some for 4 metrics by comfortable margins
  pareto_values = reshape(all_values, [], 7);
  pareto_diffs  = reshape(all_diffs,[],7);
  pareto_ranks  = reshape(all_ranks,[],7);
  pareto_bests  = all_configs(:);
  pareto_inverse = all_is_inverse(:);
  
  fprintf('Total number of configs: %d\n', sum(~isnan(pareto_values(:,1))));
  estudio = cell(0,10);
  for n=0:4
    combs = nchoosek(1:7, n);
    for k=1:size(combs,1)
      pareto_values = reshape(all_values, [], 7);
      pareto_diffs  = reshape(all_diffs,[],7);
      pareto_ranks  = reshape(all_ranks,[],7);
      pareto_bests  = all_configs(:);
      pareto_inverse = all_is_inverse(:);
      comb = combs(k,:);
      if n==0
        best_for_metric = ~isnan(pareto_inverse);
        comb = [];
      else
        best_for_metric = pareto_ranks(:,comb(1))==1;
        for g=2:numel(comb)
          best_for_metric = best_for_metric & pareto_ranks(:,comb(g))==1;
        end
      end
      if sum(best_for_metric)>0
        num_best_all = sum(best_for_metric);
        num_best_inverse_all = sum(pareto_inverse(best_for_metric));
        %fprintf('For ALL            configurations with rank 1 in metric(s) [%s], inverse ratio is %03d/%03d=%.04f\n', num2str(comb), num_best_inverse_all, num_best_all, num_best_inverse_all/num_best_all);
        pareto_values = pareto_values(best_for_metric,:);
        pareto_diffs  =  pareto_diffs(best_for_metric,:);
        pareto_ranks  =  pareto_ranks(best_for_metric,:);
        pareto_bests  =  pareto_bests(best_for_metric);
        pareto_inverse = pareto_inverse(best_for_metric);
        
        sum_ranks = nansum(pareto_ranks,2);
        best_ranksum = min(sum_ranks);
        best_by_ranksum = sum_ranks==best_ranksum;
        num_best_ranksum = sum(best_by_ranksum);
        num_inv_best_ranksum = sum(best_by_ranksum & pareto_inverse);
        
        pareto_values(:,5:7) = -pareto_values(:,5:7);
        [membership, member_value]=find_pareto_frontier(pareto_values);
        pareto_values = pareto_values(membership,:);
        pareto_diffs  = pareto_diffs(membership,:);
        pareto_bests  = pareto_bests(membership);
        pareto_inverse = pareto_inverse(membership);
        pareto_ranks  = pareto_ranks(membership,:);
        %fprintf('For PARETO-OPTIMAL configurations with rank 1 in metric(s) [%s], inverse ratio is %03d/%03d=%.04f\n', num2str(comb), sum(pareto_inverse), numel(pareto_inverse), sum(pareto_inverse)/numel(pareto_inverse));
        estudio(end+1,:) = {comb, num_best_inverse_all, num_best_all, num_best_inverse_all/num_best_all, sum(pareto_inverse), numel(pareto_inverse), sum(pareto_inverse)/numel(pareto_inverse), num_inv_best_ranksum, num_best_ranksum, num_inv_best_ranksum/num_best_ranksum};
      end
    end
  end
  fprintf('Considering all configurations:\n');
  estudio

  for n_type=1:3
    for h_size=1:3
      for update_type=1:3
        if update_type==1
          funcs = 1:3;
          update = 'single';
        elseif update_type==2
          funcs = 5:7;
          update = 'multi';
        else
          funcs = 4;
          update = 'constant';
        end
        pareto_values_temp  = reshape(all_values(h_size, n_type, funcs, :,:),[],7);
        pareto_diffs_temp   = reshape( all_diffs(h_size, n_type, funcs, :,:),[],7);
        pareto_ranks_temp   = reshape( all_ranks(h_size, n_type, funcs, :,:),[],7);
        pareto_bests_temp   = reshape(all_configs(h_size, n_type, funcs, :), [],1);
        pareto_inverse_temp = reshape(all_is_inverse(h_size, n_type, funcs, :), [], 1);
        fprintf('\n\nTotal number of configs with network %s, history size %d, update type %s: %d\n', netnames{n_type}, histSizes(h_size), update, sum(~isnan(pareto_values_temp(:,1))));

  sorted_all_values = pareto_values_temp;
  sorted_all_values = sorted_all_values(~isnan(sorted_all_values(:,1)),:);
  for h=1:7
    sorted_all_values(:,h) = sort(sorted_all_values(:,h), directions{h});
  end
  pareto_autoranks_temp = zeros(size(pareto_ranks_temp));
  pareto_autoranks_temp(:) = nan;
  for h = 1:size(pareto_autoranks_temp,1)
    for gg=1:numel(labnames)
      if ~isnan(pareto_values_temp(h, gg))
        idx = min(find(sorted_all_values(:,gg)==pareto_values_temp(h, gg)));
        pareto_autoranks_temp(h, gg) = idx;
      end
    end
  end
  sum_autoranks = nansum(pareto_autoranks_temp, 2);
  sum_autoranks(sum_autoranks==0) = inf;
  first_quartile = quantile(sum_autoranks(~isinf(sum_autoranks)), 0.25);
  zz = (pareto_inverse_temp(sum_autoranks<first_quartile));
  fprintf('  Ratio of configurations with inverse function in first quartile: %g\n', sum(zz)/numel(zz));
  %[numel(zz)/numel(sum_autoranks(~isinf(sum_autoranks))), sum(zz)/numel(zz)]

        estudio = cell(0,10);
        for n=0:4
          combs = nchoosek(1:7, n);
          for k=1:size(combs,1)
            pareto_values = pareto_values_temp;
            pareto_diffs  = pareto_diffs_temp;
            pareto_ranks  = pareto_ranks_temp;
            pareto_bests  = pareto_bests_temp;
            pareto_inverse = pareto_inverse_temp;
            comb = combs(k,:);
            if n==0
              best_for_metric = ~isnan(pareto_inverse);
              comb = [];
            else
              best_for_metric = pareto_ranks(:,comb(1))==1;
              for g=2:numel(comb)
                best_for_metric = best_for_metric & pareto_ranks(:,comb(g))==1;
              end
            end
            if sum(best_for_metric)>0
              num_best_all = sum(best_for_metric);
              num_best_inverse_all = sum(pareto_inverse(best_for_metric));
              %fprintf('For ALL            configurations with rank 1 in metric(s) [%s], inverse ratio is %03d/%03d=%.04f\n', num2str(comb), num_best_inverse_all, num_best_all, num_best_inverse_all/num_best_all);
              pareto_values = pareto_values(best_for_metric,:);
              pareto_diffs  =  pareto_diffs(best_for_metric,:);
              pareto_ranks  =  pareto_ranks(best_for_metric,:);
              pareto_bests  =  pareto_bests(best_for_metric);
              pareto_inverse = pareto_inverse(best_for_metric);
        
              sum_ranks = nansum(pareto_ranks,2);
              best_ranksum = min(sum_ranks);
              best_by_ranksum = sum_ranks==best_ranksum;
              num_best_ranksum = sum(best_by_ranksum);
              num_inv_best_ranksum = sum(best_by_ranksum & pareto_inverse);
        
              pareto_values(:,5:7) = -pareto_values(:,5:7);
              [membership, member_value]=find_pareto_frontier(pareto_values);
              pareto_values = pareto_values(membership,:);
              pareto_diffs  = pareto_diffs(membership,:);
              pareto_bests  = pareto_bests(membership);
              pareto_inverse = pareto_inverse(membership);
              pareto_ranks  = pareto_ranks(membership,:);
              %fprintf('For PARETO-OPTIMAL configurations with rank 1 in metric(s) [%s], inverse ratio is %03d/%03d=%.04f\n', num2str(comb), sum(pareto_inverse), numel(pareto_inverse), sum(pareto_inverse)/numel(pareto_inverse));
              estudio(end+1,:) = {comb, num_best_inverse_all, num_best_all, num_best_inverse_all/num_best_all, sum(pareto_inverse), numel(pareto_inverse), sum(pareto_inverse)/numel(pareto_inverse), num_inv_best_ranksum, num_best_ranksum, num_inv_best_ranksum/num_best_ranksum};
            end
          end
        end
        %estudio


      end
    end
  end


  best_acc_mse  = pareto_ranks(:,1)==1 & pareto_ranks(:,5)==1;
  best_acc_mse  = pareto_ranks(:,4)==1;
  pareto_values = pareto_values(best_acc_mse,:);
  pareto_diffs  =  pareto_diffs(best_acc_mse,:);
  pareto_ranks  =  pareto_ranks(best_acc_mse,:);
  pareto_bests  =  pareto_bests(best_acc_mse);
  pareto_inverse = pareto_inverse(best_acc_mse);

  %cond1 = pareto_values(:,1)>0.52;
  %cond2 = pareto_values(:,2)>9.1;
  %cond3 = pareto_values(:,5)<425;
  %cond4 = pareto_values(:,7)<0.1;% | isnan(pareto_values(:,7));
  %condA = cond1 & cond2 & cond3 & cond4;
  %pareto_values = pareto_values(condA,:);
  %pareto_diffs  = pareto_diffs(condA,:);
  %pareto_bests  = pareto_bests(condA);
  %pareto_inverse = pareto_inverse(condA);

  %pareto_values(:,5:7) = -pareto_values(:,5:7);
  %[membership, member_value]=find_pareto_frontier(pareto_values);%(:,[3,4,6]));
  %pareto_values = pareto_values(membership,:);
  %pareto_diffs  = pareto_diffs(membership,:);
  %pareto_bests  = pareto_bests(membership);
  %pareto_inverse = pareto_inverse(membership);
  %pareto_ranks  = pareto_ranks(membership,:);
  end

  if false
  %here, we find number of times a config is better than competitors
  pareto_values = reshape(all_values, [], 7);
  pareto_diffs  = reshape(all_diffs,[],7);
  pareto_ranks  = reshape(all_ranks,[],7);
  pareto_bests  = all_configs(:);
  pareto_inverse = all_is_inverse(:);
  
  pareto_values = best_values';
  pareto_diffs  = best_diffs';
  pareto_ranks  = best_ranks';
  pareto_bests  = bests2;
  pareto_inverse = best_is_inverse;

  num_wins      = zeros(size(pareto_values,1),1);
  competitors   = competitors';
  for a=1:size(competitors,1)
    for b=1:size(competitors,2)
      if initbest(b)==-inf && ~isinf(competitors(a,b))
        num_wins  = num_wins + double(pareto_values(:,b)>competitors(a,b));
      else
        num_wins  = num_wins + double(pareto_values(:,b)<competitors(a,b));
      end
    end
  end
  maxiwin     = max(num_wins)
  best_values = pareto_values(num_wins==maxiwin,:)
  bests       = pareto_bests(num_wins==maxiwin,:);
  for k=1:numel(bests)
    zz=bests{k,:}
  end
  end

  pareto_sum_ranks = nansum(pareto_ranks, 2);
  isCompetitive = any(isnan(pareto_ranks), 2);
  pareto_sum_ranks(pareto_sum_ranks==0)=inf;
  idxs=find(min(pareto_sum_ranks)==pareto_sum_ranks);
  pareto_sum_ranks2 = nansum(pareto_ranks, 2);
  pareto_sum_ranks2(pareto_sum_ranks==0)=inf;
  pareto_sum_ranks2(isCompetitive)=inf;
  idxs2=find(min(pareto_sum_ranks2)==pareto_sum_ranks2);
  m1=min(pareto_sum_ranks);
  v1=pareto_values(idxs,:);
  i1=pareto_inverse(idxs);
  b1=pareto_bests(idxs,:);
  r1=pareto_ranks(idxs,:);
  m2=min(pareto_sum_ranks2);
  v2=pareto_values(idxs2,:);
  b2=pareto_bests(idxs2,:);
  i2=pareto_inverse(idxs2);
  r2=pareto_ranks(idxs2,:);


  fprintf(' Total: %d, inverse: %d\n', numel(best_is_inverse), nansum(best_is_inverse));
  fprintf(' Total: %d, inverse: %d\n', numel(pareto_inverse), nansum(pareto_inverse));
  results{p} = {all_values, all_configs, bestvalscompetitors, best_values, best_diffs, bests};
  %fsum = fopen(pathsSummary{p}, 'w');
  %for m = 1:numel(labnames)
  %  fprintf(fsum, '%s, %f, %s, %d, %s, %s, %f\n', best_config{m}{:});
  %end
  %fclose(fsum);

end

end


function [membership, member_value]=find_pareto_frontier(input)
out=[];
arenan = isnan(input(:,end));
input(arenan,end)=4;
data=unique(input,'rows');
for i = 1:size(data,1)
    
    c_data = repmat(data(i,:),size(data,1),1);
    t_data = data;
    t_data(i,:) = Inf(1,size(data,2));
    smaller_idx = c_data>=t_data;
    
    idx=sum(smaller_idx,2)==size(data,2);
    if ~nnz(idx)
        out(end+1,:)=data(i,:);
    end
end
%arenan = isnan(input(:,end));
%input_nnan = input;
%input_nnan(arenan,end) = 0;
%out_nnan = out;
%out_nnan(arenan,end) = 0;
%membership = ismember(input_nnan,out_nnan,'rows');
membership = ismember(input,out,'rows');
member_value = out;
end

