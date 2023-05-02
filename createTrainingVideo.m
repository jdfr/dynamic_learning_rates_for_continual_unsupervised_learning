function createTrainingVideo(videoname, fps, skipframes, matfilepath, experiment_index, FileNames, NumSamplesPerImage)


%createTrainingVideo('images/Circle.Anillo.mp4', 15, 10, 'images/competitive/bateria_exploration_I1_Sliding1_mqeMode0_MAXLR_1_SAM_40000/Circle.experiments.mat', 1, {'images/Circle.bmp', 'images/Anillo.bmp'}, 5000);

%createTrainingVideo('images/Circle.Anillo.mp4', 15, 10, 'images/competitive/bateria3_MAXLR_1_SAM_10000/Circle.Anillo.experiments.mat', 1, {'images/Circle.bmp', 'images/Anillo.bmp'}, 10000);

%createTrainingVideo('images/Anillo.Circle.mp4', 15, 10, 'images/competitive/bateria_dosimgs_I2_Sliding1_mqeMode0_MAXLR_1_SAM_10000/Anillo.Circle.experiments.mat', 1, {'images/Anillo.bmp', 'images/Circle.bmp'}, 10000);
%createTrainingVideo('images/Circle.Anillo.mp4', 15, 10, 'images/competitive/bateria_dosimgs_I2_Sliding1_mqeMode0_MAXLR_1_SAM_10000/Circle.Anillo.experiments.mat', 1, {'images/Circle.bmp', 'images/Anillo.bmp'}, 10000);
%createTrainingVideo('images/Ocho.RelojArena.mp4', 15, 10, 'images/competitive/bateria_dosimgs_I2_Sliding1_mqeMode0_MAXLR_1_SAM_10000/Ocho.RelojArena.experiments.mat', 1, {'images/Ocho.bmp', 'images/RelojArena.bmp'}, 10000);
%createTrainingVideo('images/RelojArena.Ocho.mp4', 15, 10, 'images/competitive/bateria_dosimgs_I2_Sliding1_mqeMode0_MAXLR_1_SAM_10000/RelojArena.Ocho.experiments.mat', 1, {'images/RelojArena.bmp', 'images/Ocho.bmp'}, 10000);
%createTrainingVideo('images/Eme.Equis.mp4', 15, 10, 'images/competitive/bateria_dosimgs_I2_Sliding1_mqeMode0_MAXLR_1_SAM_10000/Eme.Equis.experiments.mat', 1, {'images/Eme.bmp', 'images/Equis.bmp'}, 10000);
%createTrainingVideo('images/Equis.Eme.mp4', 15, 10, 'images/competitive/bateria_dosimgs_I2_Sliding1_mqeMode0_MAXLR_1_SAM_10000/Equis.Eme.experiments.mat', 1, {'images/Equis.bmp', 'images/Eme.bmp'}, 10000);

%createTrainingVideo('images/Anillo.Circle.mp4', 15, 10, 'images/competitive/bateria3_MAXLR_1_SAM_10000/Anillo.Circle.experiments.mat', 1, {'images/Anillo.bmp', 'images/Circle.bmp'}, 10000);
%createTrainingVideo('images/Circle.Anillo.mp4', 15, 10, 'images/competitive/bateria3_MAXLR_1_SAM_10000/Circle.Anillo.experiments.mat', 1, {'images/Circle.bmp', 'images/Anillo.bmp'}, 10000);
%createTrainingVideo('images/Ocho.RelojArena.mp4', 15, 10, 'images/competitive/bateria3_MAXLR_1_SAM_10000/Ocho.RelojArena.experiments.mat', 1, {'images/Ocho.bmp', 'images/RelojArena.bmp'}, 10000);
%createTrainingVideo('images/RelojArena.Ocho.mp4', 15, 10, 'images/competitive/bateria3_MAXLR_1_SAM_10000/RelojArena.Ocho.experiments.mat', 1, {'images/RelojArena.bmp', 'images/Ocho.bmp'}, 10000);
%createTrainingVideo('images/Eme.Equis.mp4', 15, 10, 'images/competitive/bateria3_MAXLR_1_SAM_10000/Eme.Equis.experiments.mat', 1, {'images/Eme.bmp', 'images/Equis.bmp'}, 10000);
%createTrainingVideo('images/Equis.Eme.mp4', 15, 10, 'images/competitive/bateria3_MAXLR_1_SAM_10000/Equis.Eme.experiments.mat', 1, {'images/Equis.bmp', 'images/Eme.bmp'}, 10000);


%createTrainingVideo('images/Anillo.Circle.mp4', 15, 10, 'images/competitive/bateria_probando_In_Sliding1_mqeMode0_MAXLR_1_SAM_10000/Anillo.Circle.experiments.mat', 1, {'images/Anillo.bmp', 'images/Circle.bmp'}, 10000);

%createTrainingVideo('images/circulito1.circulito2.mp4', 15, 10, 'images/circulito1.circulito2.experiments.mat', 1, {'images/circulito1.png', 'images/circulito2.png'}, 100);

%createTrainingVideo('images/Anillo.Circle.mp4', 15, 10, 'images/gng/bateria_exploration_N50_I2_Sliding1_mqeMode0_MAXLR_1_SAM_40000/Anillo.Circle.experiments.mat', 1, {'images/Anillo.bmp', 'images/Circle.bmp'}, 10000);

load(matfilepath); % experiments, best_run

[Samples, NumSteps, FusedFileNames] = generateSamplesMultiImgs([], FileNames, NumSamplesPerImage);

sample_limits = 1:NumSamplesPerImage:(NumSteps+NumSamplesPerImage*2);

a = experiment_index;

HistorySize = best_run.parameters(a, 1);
factr       = best_run.parameters(a, 2)
HistorySize=1;
factr=10;
%initial_LearningRate=experiments{a}.initial_LearningRate;
%slidingWindow=experiments{a}.slidingWindow;
%mqeMode=experiments{a}.mqeMode;
%netType=experiments{a}.netType;
initial_LearningRate=1.0;
slidingWindow=true;
mqeMode=0;
netType='gng';
trainMode = 'singleLR';

[Model, quantization_errors, consensus_quantization_errors, learning_rates]=TrainOnline(trainMode, cat(2, Samples{:}),netType, experiments{a}.params,NumSteps, HistorySize, experiments{a}.type, factr, mqeMode, slidingWindow, 0, initial_LearningRate, experiments{a}.MaxLearningRate, experiments{a}.MinLearningRate, {videoname, fps, skipframes, sample_limits}, Samples{1}, size(Samples{1}, 2), Samples);

