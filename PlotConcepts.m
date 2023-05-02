function plotConcepts(cora, Model)

%Graphics in the paper:
%    bestrun SOM nobatches multi linear:
%          best_run.Model{3}
%          experiments_nobatches_100exps/som/CORA_exploration_multiLR_N100_Sliding1_mqeMode0_MAXLR_1_SAM_10832/experiments.mat
%    bestrun GNG nobatches multi linear:
%          best_run.Model{2}
%          experiments_nobatches_100exps/gng/CORA_exploration_multiLR_N100_Sliding1_mqeMode0_MAXLR_1_SAM_10832/experiments.mat

%classIds = cora.ClassIds;
classIds = {'GA', 'ReL', 'T', 'RuL', 'CB', 'PM', 'NN'};
colors   = {'r', 'g', 'b', 'c', 'm', 'k', '#D95319'};

isSOM = numel(size(Model.Prototypes))==3;
isGNG = numel(size(Model.Prototypes))==2 && isfield(Model, 'Connections');

if isSOM
  [pos connections] = prepareSOMPlot(Model);
  shft = 0.1;
elseif isGNG
  [pos connections] = prepareGNGPlot(Model);
  shft = 0.1;
end

prototype_freqs = computePrototypeClasses(cora, Model.Prototypes);

%drawNetwork(Model.Prototypes, pos, connections, classIds, prototype_freqs, colors, shft, 1);
drawNetwork(Model.Prototypes, pos, connections, classIds, prototype_freqs, colors, shft, 2);
%drawNetwork(Model.Prototypes, pos, connections, classIds, prototype_freqs, colors, shft, 3);
%drawNetwork(Model.Prototypes, pos, connections, classIds, prototype_freqs, colors, shft, 4);


function [pos Connections] = prepareSOMPlot(Model)

posX = repmat((1:Model.NumRowsMap), Model.NumColsMap, 1);
posY = repmat((1:Model.NumColsMap)', 1, Model.NumRowsMap);
pos = [posX(:), posY(:)];

Connections = zeros(Model.NumRowsMap*Model.NumColsMap);
for NdxNeuro=1:size(Connections,1)
    Connections(NdxNeuro,:) = (Model.TopolDist{NdxNeuro} == 1);
end

function [pos Connections] = prepareGNGPlot(Model)

G = graph(Model.Connections);

layout = 'force';
%layout = 'subspace';
%layout = 'layered';
%layout = 'circle';

h = matlab.graphics.chart.primitive.GraphPlot('BasicGraph', MLGraph(G), 'Layout', layout, 'Iterations', 100000, 'UseGravity', true);
% Get its XData and YData properties:
X = h.XData;
Y = h.YData;

pos = [X(:), Y(:)];
Connections = Model.Connections;

function drawNetwork(prototypes, prototypesPos, connections, classNames, prototype_freqs, colors, shft, MAX_WORDS)

figure;

hold on;

prototypes = prototypes(:,:);

for NdxUnit=1:size(prototypes,2)
    if isfinite(prototypes(1,NdxUnit))
        NdxNeighbors=find(connections(NdxUnit,:));
        for NdxMyNeigh=1:numel(NdxNeighbors)
            line([prototypesPos(NdxUnit,1) prototypesPos(NdxNeighbors(NdxMyNeigh),1)],...
                 [prototypesPos(NdxUnit,2) prototypesPos(NdxNeighbors(NdxMyNeigh),2)], 'Color', '#EDB120');
        end
    end
end

for k=1:numel(prototype_freqs)
  [freqs, classes_by_freq] = sort(prototype_freqs{k}, 'descend');
  if isfinite(prototypes(1,k))
    for NdxWord=1:MAX_WORDS
      MyConceptFreq = 0.3;
      ShowWord(prototypesPos(k,:),NdxWord,MAX_WORDS,classNames{classes_by_freq(NdxWord)},shft,MyConceptFreq);
    end
  end
end

maxfreqs = zeros(size(prototype_freqs));
for k=1:numel(prototype_freqs)
  [~, maxfreqs(k)] = max(prototype_freqs{k});
end
for k=1:numel(colors)
  plot(prototypesPos(maxfreqs==k,1), prototypesPos(maxfreqs==k,2), '.', 'Color', colors{k}, 'MarkerSize', 10);
end
%plot(prototypesPos(:,1), prototypesPos(:,2), '.b', 'MarkerSize', 8);

axis off;
axis equal;

function ShowWord(Position,NdxWord,MAX_WORDS,MyWord,shft,MyConceptFreq)
% Show a word around the position of a node (4 words are assumed as maximum).
% Inputs:
%   Position = coordinates of the node position
%   NdxWord = index of the word (1-4)
%   MyWord = word to represent
%   MyConceptFreq = frequency of the concept

if MAX_WORDS==1
  switch NdxWord
      case 1
          text(Position(1),Position(2),MyWord,'Color',[0 0 0],'FontSize',max(10*MyConceptFreq, 8), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Interpreter', 'none');
  end
elseif MAX_WORDS==2
  switch NdxWord
      case 1
          text(Position(1),Position(2),MyWord,'Color',[0 0 0],'FontSize',max(10*MyConceptFreq, 8), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Interpreter', 'none');
      case 2
          text(Position(1),Position(2),MyWord,'Color',[0 0 0],'FontSize',max(10*MyConceptFreq,8),'HorizontalAlignment','center','VerticalAlignment','top', 'Interpreter', 'none');
  end
else
  switch NdxWord
      case 1
          text(Position(1),Position(2),MyWord,'Color',[0 0 0],'FontSize',max(10*MyConceptFreq, 8), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Interpreter', 'none');
      case 2
          text(Position(1)+shft,Position(2),MyWord,'Color',[0 0 0],'FontSize',max(10*MyConceptFreq,8), 'Interpreter', 'none');
      case 3
          text(Position(1),Position(2),MyWord,'Color',[0 0 0],'FontSize',max(10*MyConceptFreq,8),'HorizontalAlignment','center','VerticalAlignment','top', 'Interpreter', 'none');
      case 4
          text(Position(1)-shft,Position(2),MyWord,'Color',[0 0 0],'FontSize',max(10*MyConceptFreq,8),'HorizontalAlignment','right', 'Interpreter', 'none');
  end
end



function prototype_freqs = computePrototypeClasses(cora, prototypes)

numClasses = max(cora.paperClass);

prototypes = prototypes(:,:);

[errs, receptive_field_sizes, winners] = compute_consensus_quantization_error(nan, prototypes, cora.paperWords');

prototype_freqs = cell(size(prototypes,2),1);

for k=1:numel(prototype_freqs)
  classes_receptive_field = cora.paperClass(winners==k);
  freqs = zeros(numClasses,1);
  for x=1:numClasses
    freqs(x) = sum(classes_receptive_field==x);
  end
  prototype_freqs{k} = freqs;
end







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


