function plotDistanceCorrelation(cora, Model)

isSOM = numel(size(Model.Prototypes))==3;
isGNG = numel(size(Model.Prototypes))==2 && isfield(Model, 'Connections');

if isSOM
  connections = zeros(Model.NumRowsMap*Model.NumColsMap);
  for NdxNeuro=1:size(connections,1)
      connections(NdxNeuro,:) = (Model.TopolDist{NdxNeuro} == 1);
  end
else
  connections = double(logical(Model.Connections));
end
hops = hopDistances(connections);

allVecs = double(cora.paperWords);
allEuclidean = squareform(pdist(allVecs));
allDists = zeros(numel(cora.paperClass), numel(cora.paperClass), 2);
allDists(:,:,1) = allEuclidean;

[errs, receptive_field_sizes, winners] = compute_consensus_quantization_error(nan, Model.Prototypes(:,:), cora.paperWords');
for a=1:numel(cora.paperClass)
  for b=a+1:numel(cora.paperClass)
    wa = winners(a);
    wb = winners(b);
    allDists(a,b,2) = hops(wa, wb);
  end
end

X = allDists(:,:,1);
Y = allDists(:,:,2);
X = X(:);
Y = Y(:);
notsame = Y~=0;
Y=Y(notsame);
X=X(notsame);

pos = [X Y];
poss = unique(pos, 'rows');

step = 1;

plot(poss(1:step:end,2), poss(1:step:end,1), '.');

xlabel('distancias topológicas (entre centroides ganadores de cada par de muestras)');
ylabel('distancias euclideas (entre cada par de muestras)');

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

