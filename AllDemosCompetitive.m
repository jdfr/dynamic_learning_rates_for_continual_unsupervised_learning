function AllDemosCompetitive(experimentName, imagestopresent, trainingType, netType, mqeMode, slidingWindow, MaxLearningRate, NumSamplesPerImage);


%AllDemosCompetitive('_exploration', '2', 'online', 'gng', 0, true, 1, 40000);
%system('mv images/gng images/gng50');
%AllDemosCompetitive('_exploration', '2', 'online', 'gng', 0, true, 1, 40000);
%system('mv images/gng images/gng25');
%AllDemosCompetitive('_exploration', '2', 'online', 'competitive', 0, true, 1, 40000);
%AllDemosCompetitive('_exploration', '2', 'online', 'som', 0, true, 1, 40000);
%
%AllDemosCompetitive('_exploration_80000', '2', 'online', 'som', 0, true, 1, 80000);

%AllDemosCompetitive('_exploration', '2', 'online', 'gng', 0, true, 1, 40000);
%system('mv images/gng images/gng50');
%AllDemosCompetitive('_exploration', '2', 'online', 'competitive', 0, true, 1, 40000);
%AllDemosCompetitive('_exploration', '2', 'online', 'som', 0, true, 1, 40000);
%AllDemosCompetitive('_exploration', '2', 'online', 'gng', 0, true, 1, 40000);
%system('mv images/gng images/gng25');



%AllDemosCompetitive('images', '2', 'online', 'competitive', 0, true, 1.0, 10000); AllDemosCompetitive('images', '2', 'online', 'som', 0, true, 1.0, 10000);


%AllDemosCompetitive('_probando', 'n', 'online', 'competitive', 0, true, 1.0, 10000);

%PROBAR TODAS LAS COMBINACIONES DE mqeMode Y slidingWindow, Y JUGAR CON initial_learningrate EN LOS QUE NO SON SLIDING WINDOW


n={'images/Anillo.bmp', 'images/Circle.bmp', 'images/CuadradoHueco.bmp', 'images/CuadradoNoHueco.bmp', 'images/Eme.bmp', 'images/Equis.bmp', 'images/Ese.bmp', 'images/Haltera.bmp', 'images/Irregular.bmp', 'images/Ka.bmp', 'images/Ocho.bmp', 'images/RelojArena.bmp', 'images/circulito1.png', 'images/circulito2.png'};

if     strcmp(imagestopresent, '1')
  names = cellfun(@(x){x}, n, 'UniformOutput', false);
elseif strcmp(imagestopresent, '2')
  idxs = [1 2;
          2 1;
          3 4;
          4 3;
          5 6;
          6 5;
          8 9;
          9 8;
          2 6;
          6 2;
          1 6;
          6 1;
          2 7;
          7 2;
          11 12;
          12 11;
          7 10;
          10 7;
          ];
  idxs = [1 6;
          3 4];
  idxs2 = [13 14;
          14 13];
  names = arrayfun(@(x){n{idxs(x,:)}}, 1:size(idxs,1), 'UniformOutput', false);
elseif strcmp(imagestopresent, 'n')
  names = {{n{1}, n{2}}};
end

if strcmp(trainingType, 'online')
  %parpool(8);
  %parpool(12);
  %parpool(6);

  initial_LearningRate = MaxLearningRate;
  otherargs = {trainingType, netType, mqeMode, slidingWindow, initial_LearningRate, MaxLearningRate, NumSamplesPerImage};

  dr = sprintf('images/%s/bateria%s_I%s_Sliding%d_mqeMode%d_MAXLR_%g_SAM_%d', netType, experimentName, imagestopresent, slidingWindow, mqeMode, MaxLearningRate, NumSamplesPerImage);
  system(sprintf('mkdir -p %s', dr));
  %system(sprintf('mv images/*.mat images/*.fig %s/', dr));
  dr = [dr '/'];

  parfor k=1:length(names);
  %for k=1:length(names);
    DemoCompetitive(dr, names{k}, false, false, otherargs{:});
  end

  for k=1:length(names);
    DemoCompetitive(dr, names{k}, true, true, otherargs{:});
  end


elseif strcmp(trainingType, 'traditional')

  for k=1:length(names);
    DemoCompetitive(names{k}, true, true, trainingType, netType, MaxLearningRate, NumSamplesPerImage);
  end
  system(sprintf('mkdir -p images/%s/baseline',                      netType));
  system(sprintf('mv images/*.mat images/*.fig images/%s/baseline/', netType));

end


