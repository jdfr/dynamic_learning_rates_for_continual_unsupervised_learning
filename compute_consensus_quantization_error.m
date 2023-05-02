function [err, varargout]=compute_consensus_quantization_error(mqeMode, prototypes, samples)

[nevermind,nump]    = size(prototypes);
quantization_errors = zeros(size(samples, 2),1);
%assigned_neurons    = zeros(size(samples, 2),1);
%running_totals      = zeros(nump, 1);
%num_samples         = zeros(nump, 1);
receptive_field_sizes = zeros(nump,1);
outputWinners         = nargout>2;
winners               = zeros(nump, 1);
for k=1:size(samples, 2)
    RepMySample               = repmat(samples(:,k),1,nump);        
    MyDistances               = sqrt(sum((RepMySample-prototypes(:,:)).^2,1));            			    
    [~,NdxWinner]             = min(MyDistances);
    winners(k)                = NdxWinner;
    dstvec                    = samples(:,k)-prototypes(:,NdxWinner);
    dst                       = sqrt(sum(dstvec.*dstvec));
    receptive_field_sizes(NdxWinner) = receptive_field_sizes(NdxWinner)+1;
    %running_totals(NdxWinner) = running_totals(NdxWinner) + dst;
    %num_samples(NdxWinner)    = num_samples(NdxWinner) + 1;
    quantization_errors(k) = dst;
    %assigned_neurons(k)    = NdxWinner;
end

err = sum(quantization_errors)/nump;

if nargout>1
  varargout{1} = receptive_field_sizes;
end
if nargout>2
  varargout{2} = winners;
end

if false
switch mqeMode
  % median of all errors (without binning by assigned neuron)
  case 0
    err = median(quantization_errors);
  % mean   of all errors (without binning by assigned neuron)
  case 1
    err = mean(quantization_errors);
  % bin errors by neuron, discard neurons with 0 errors, compute mean error by neuron, compute median of these means
  case 2
    nonzero = num_samples~=0;
    values  = running_totals(nonzero)./num_samples(nonzero);
    err     = median(values);
  % bin errors by neuron, discard neurons with 0 errors, compute mean error by neuron, compute mean   of these means
  case 3
    nonzero = num_samples~=0;
    values  = running_totals(nonzero)./num_samples(nonzero);
    err     = mean(values);
  % bin errors by neuron, discard neurons with 0 errors, add up errors by neuron, compute median of these summations
  case 4
    nonzero = num_samples~=0;
    values  = running_totals(nonzero);
    err     = median(values);
  % bin errors by neuron, discard neurons with 0 errors, add up errors by neuron, compute mean   of these summations
  case 5
    nonzero = num_samples~=0;
    values  = running_totals(nonzero);
    err     = mean(values);
  % bin errors by neuron, add up errors by neuron, compute median of these summations
  case 6
    values  = running_totals;
    err     = median(values);
  % bin errors by neuron, add up errors by neuron, compute mean   of these summations
  case 7
    values  = running_totals;
    err     = mean(values);
end
end









%[nevermind,nump]=size(prototypes);

%distances = pdist2(prototypes', samples');

%[mins, min_idxs] = min(distances, [], 1);

%total_err = 0;

%for i=1:nump
%  %samples_for_prototype = samples(:,min_idxs==i);
%  quantization_error_for_i = sum(distances(i,min_idxs==i));
%  total_err = total_err + quantization_error_for_i;
%end

%err = total_err/nump;

