function createSnapshotFigure(matfilepath, experiment_index, FileNames, NumSamplesPerImage)

%createSnapshotFigure('images/competitive/bateria_logspace_I2_Sliding1_mqeMode0_MAXLR_1_SAM_10000/Equis.Eme.experiments.mat', 3, {'images/Equis.bmp', 'images/Eme.bmp'}, 10000);
%createSnapshotFigure('images/competitive/bateria_logspace_I2_Sliding1_mqeMode0_MAXLR_1_SAM_10000/Equis.Eme.experiments.mat', 3, {'images/Equis.bmp', 'images/Eme.bmp'}, 10000);
%createSnapshotFigure('images/competitive/bateria_logspace2_I2_Sliding1_mqeMode0_MAXLR_1_SAM_10000/Anillo.Equis.experiments.mat', 4, {'images/Anillo.bmp', 'images/Equis.bmp'}, 10000);
%createSnapshotFigure('images/Equis.Eme.experiments.mat', 3, {'images/Equis.bmp', 'images/Eme.bmp'}, 10000);

load(matfilepath); % experiments, best_run

[Samples, NumSteps, FusedFileNames] = generateSamplesMultiImgs([], FileNames, NumSamplesPerImage);

sample_limits = 1:NumSamplesPerImage:(NumSteps+NumSamplesPerImage*2);

a = experiment_index;

HistorySize = best_run.parameters(a, 1);
factr       = best_run.parameters(a, 2);

imgs = {imread(FileNames{1}), imread(FileNames{end})};
for x=1:numel(imgs)
  if ndims(imgs{x})==3
    imgs{x} = mean(imgs{x}, 3);
  end
  imgs{x} = logical(imgs{x});
end

%initial_LearningRate=experiments{a}.initial_LearningRate;
%slidingWindow=experiments{a}.slidingWindow;
%mqeMode=experiments{a}.mqeMode;
%netType=experiments{a}.netType;
initial_LearningRate=1.0;
slidingWindow=true;
mqeMode=0;
netType='competitive';
trainMode = 'singleLR';
HistorySize
experiments{a}.type
factr
[Model, quantization_errors, consensus_quantization_errors, learning_rates]=TrainOnline(trainMode,cat(2, Samples{:}),netType, experiments{a}.params,NumSteps, HistorySize, experiments{a}.type, factr, mqeMode, slidingWindow, 0, initial_LearningRate, experiments{a}.MaxLearningRate, experiments{a}.MinLearningRate, 'cadena de caracteres cualquiera', [], [], {});
Samples = cat(2, Samples{:});

ks=[1 1 NumSamplesPerImage 'k' 1; 10 1 NumSamplesPerImage 'k' 1; 5000 1 NumSamplesPerImage 'k' 1; 10000 1 NumSamplesPerImage 'w' 1; 10010 NumSamplesPerImage+1 2*NumSamplesPerImage 'k' 2; 15000 NumSamplesPerImage+1 2*NumSamplesPerImage 'k' 2; 20000 NumSamplesPerImage+1 2*NumSamplesPerImage 'k' 2];
h=figure;
tiledlayout(3, size(ks,1));
fontsize=24;
for i=1:size(ks,1)
  nexttile;
  %plot(Samples(1,ks(i,2):ks(i,3)), Samples(2,ks(i,2):ks(i,3)), ['*' char(ks(i,4))], Model.all_prototypes{ks(i,1)}(1,:),Model.all_prototypes{ks(i,1)}(2,:),'or','LineWidth',2,'MarkerFaceColor',[1 0 0],'MarkerSize',7);
  hold on;
  if ks(i,4)=='w'
    colormap([1 1 1; 1 1 1]);
  else
    colormap([0 0 0; 1 1 1]);
  end
  imagesc([0 1], [1 0], imgs{ks(i,5)});
  plot(Model.all_prototypes{ks(i,1)}(1,:),Model.all_prototypes{ks(i,1)}(2,:),'or','LineWidth',2,'MarkerFaceColor',[1 0 0],'MarkerSize',7);
  axis off;
  title(sprintf('$t=%d$', ks(i,1)),'Interpreter','latex', 'FontSize', fontsize);
end
aaa=nexttile([2,size(ks,1)]);
t=1:numel(Model.all_prototypes);
plot(t, quantization_errors, 'c', t, consensus_quantization_errors, 'b', t, learning_rates, 'k', 'LineWidth',2);
legend({'$\mathrm{qe}(t)$', '$C(t)$', '$\alpha(t)$'},'Interpreter','latex', 'FontSize', fontsize);
xl = xlabel('$t$','Interpreter','latex');
set(xl, 'FontSize', fontsize);
aaa.XAxis.FontSize=fontsize;
aaa.YAxis.FontSize=fontsize;
%savefig(h, [FusedFileNames 'snapshots.fig'], 'compact');



