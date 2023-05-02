%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEMO FOR VIDEO SEGMENTATION USING THE GHNF
% Select the video and color space desired
% ------------------------------------------
% E.J. Palomo - November 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
warning off
rng('default');

% Initial Parameters
NdxColorSpace = 5; % color space to test
NdxVideo = 7; % video to test
NumBackgroundFrames = 100;
Tau = 0.1; % GHNF hierarchical growing parameter
Mu = 1; % penalization parameter for some dataset features
NumEpochs = 2;
MaxNeurons = 50; % Maximum number of neurons in each graph
Features = 3; % 1=RGB, 2=Normalized RGB, and 3=Normalized RGB median filtered

% The following values of the parameters are those considered in the
% original GNG paper by Fritzke (1995)
Lambda = 100;
EpsilonB = 0.2;
EpsilonN = 0.006;
Alpha = 0.5;
AMax = 50;
D = 0.995;

%% COLOR SPACES
ColorSpaces = {'HSL','HSV','Lab','Luv','RGB','YCbCr'};
NumColorSpaces = length(ColorSpaces);

%% VIDEOS TO BE PROCESSED
% Escalator, OneShopOneWait1cor no tienen frames de background
% LightSwitch tiene un cambio de iluminación
% Video4 tiene una planta que se mueve
Videos = {'Campus','Curtain','Escalator','Fountain','LevelCrossing','OneShopOneWait1cor',...
    'Video2','Video4','WaterSurface','Lobby','LightSwitch'};
NumVideos = length(Videos);
PathVideos = '../../Videos/';
SubPathFrames = '/frames/f%07d.bmp';
SubPathFramesLightSwitch = '/frames/f%07d.jpg';
SubPathGT = '/GT/GT%07d.bmp';
DeltaFrames = [999,20999,1397,999,0,0,0,0,999,999,0];
NumFrames = [1438,2963,3417,523,500,1376,749,819,632,1545,2799];
GT = [20,20,20,20,163,11,444,316,20,20];
PrimerFrameGT = [1199,22755,1398,1152,275,0,253,273,1480,1152,666];

%% DATASET GENERATION
MyColorSpace = ColorSpaces{NdxColorSpace};
fprintf('\nCOLOR SPACE: %s\n',MyColorSpace);    
switch MyColorSpace,
    case {'RGB'}
        MyColorSpaceConversion = '';
    otherwise
        MyColorSpaceConversion = ['RGB->' MyColorSpace];
end
MyVideo = Videos{NdxVideo};
fprintf('\n\tVIDEO: %s\n',MyVideo);
PathMyVideo = [PathVideos MyVideo SubPathFrames];
if NdxVideo==11, % jpg instead of bmp extension
    PathMyVideo = [PathVideos MyVideo SubPathFramesLightSwitch];
end

% Store the last frames belonging to the background
NdxFrameInicio = PrimerFrameGT(NdxVideo) - NumBackgroundFrames;
if NdxFrameInicio <= DeltaFrames(NdxVideo),
    NdxFrameInicio = DeltaFrames(NdxVideo)+1;
end
fprintf('\t\tDataset frames from %d to %d (First Frame=%d)\n',NdxFrameInicio,NdxFrameInicio+NumBackgroundFrames-1,DeltaFrames(NdxVideo)+1);
MyFrame = imread(sprintf(PathMyVideo,NdxFrameInicio));
if ~isempty(MyColorSpaceConversion),
    MyFrame = colorspace(MyColorSpaceConversion,MyFrame);
end
FirstFrames = zeros(size(MyFrame,1),size(MyFrame,2),size(MyFrame,3),NumBackgroundFrames);
FirstFrames(:,:,:,1) = MyFrame;
for NdxFrame=2:NumBackgroundFrames,
    
    NdxFrameInicio = NdxFrameInicio + 1;    
    MyFrame = imread(sprintf(PathMyVideo,NdxFrameInicio));
    if ~isempty(MyColorSpaceConversion),
        MyFrame = colorspace(MyColorSpaceConversion,MyFrame);
    end
    FirstFrames(:,:,:,NdxFrame) = MyFrame;
end

% Obtain the dataset as the mean (or the median) of the last frames
% MyFrame = mean(FirstFrames,4);
MyFrame = median(FirstFrames,4);
% MyFrame = mode(FirstFrames,4);
MyFrame = NormalizarEspColor(double(MyFrame), MyColorSpace);

% % EZEQUIEL: HERE WE DIVIDE RED, GREEN AND BLUE BY THE SUM OF RED, GREEN AND
% % BLUE
% Data(1:3,:) = Data(1:3,:)./repmat(sum(Data(1:3,:),1),[3 1]);
% % % HERE WE AVOID ERRORS CAUSED BY NaNs COMING FROM DIVIDING BY ZERO
% Data(~isfinite(Data))=0;
%--------------------------------------------------------------------------
if Features>1,  
    % NORMALIZED RGB
    SumFrame = sum(MyFrame,3);
    SumFrame(SumFrame==0) = 1;
    MyFrameNorm = MyFrame;
    MyFrameNorm(:,:,1) = squeeze(MyFrame(:,:,1))./SumFrame;        
    MyFrameNorm(:,:,2) = squeeze(MyFrame(:,:,2))./SumFrame;
    MyFrameNorm(:,:,3) = squeeze(MyFrame(:,:,3))./SumFrame;
    if Features==3,
        % NORMALIZED RGB MEDIAN FILTERED
        MyFrameNorm(:,:,1) = medfilt2(MyFrameNorm(:,:,1),[5 5]);
        MyFrameNorm(:,:,2) = medfilt2(MyFrameNorm(:,:,2),[5 5]);
        MyFrameNorm(:,:,3) = medfilt2(MyFrameNorm(:,:,3),[5 5]);
    end
    MyFrame = MyFrameNorm;
end
%--------------------------------------------------------------------------
MyMatrixFrame = reshape(shiftdim(MyFrame,2),3,[]);
[i,j] = ind2sub([size(MyFrame,1) size(MyFrame,2)],[1:length(MyMatrixFrame)]);
Data = [MyMatrixFrame;Mu*NormalizeData(i);Mu*NormalizeData(j)];

%% TRAINING AND TEST        
NumSamples = size(Data,2);        
% NumSteps = NumEpochs*NumSamples;                        

tic;
[Model] = TrainGHNF(Data,NumEpochs,MaxNeurons,Tau,Lambda,EpsilonB,EpsilonN,Alpha,AMax,D,1);
CpuTime = toc;
fprintf('\t\tTraining finished\n');     

% Test GHNF
[Winners,Errors] = TestGHNF(GetCentroidsGHNF(Model),Model.Samples);

% Evaluate Performance
PathMyGT = [PathVideos MyVideo SubPathGT];
fprintf('\t\tTesting Frames...\n');
PerformanceMeasures = DemoTestGHNFSegmentation(Model,Winners,Errors,...
    MyColorSpace,MyVideo,NumFrames(NdxVideo),DeltaFrames(NdxVideo),PathMyVideo,PathMyGT,Mu,Features);

% Time is the training time + the mean detecting foreground time           
Results.Time = CpuTime + PerformanceMeasures(1,end);
MeanPerformanceMeasures = mean(PerformanceMeasures(:,1:end-1),2);
StdPerformanceMeasures = std(PerformanceMeasures(:,1:end-1),0,2);
Results.S.Mean = MeanPerformanceMeasures(1);
Results.S.Std = StdPerformanceMeasures(1);
Results.FP.Mean = MeanPerformanceMeasures(2);
Results.FP.Std = StdPerformanceMeasures(2);
Results.FN.Mean = MeanPerformanceMeasures(3);
Results.FN.Std = StdPerformanceMeasures(3);
Results.Precision.Mean = MeanPerformanceMeasures(4);
Results.Precision.Std = StdPerformanceMeasures(4);
Results.Recall.Mean = MeanPerformanceMeasures(5);
Results.Recall.Std = StdPerformanceMeasures(5);
Results.Accuracy.Mean = MeanPerformanceMeasures(6)
Results.Accuracy.Std = StdPerformanceMeasures(6)
Results.Fmeasure.Mean = MeanPerformanceMeasures(7);
Results.Fmeasure.Std = StdPerformanceMeasures(7);

