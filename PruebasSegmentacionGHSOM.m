% Script para realizar las pruebas de segmentación del GHSOM, basado en el
% script 'PruebasSegmentacionGHNF'
% NOTA: Ejecutar desde la carpeta 'Pruebas Segmentación'
% E.J. Palomo, Mayo 2020
clear all
close all
warning off
rng('default');

NumEpocas = 2;
% MaxNeurons = 20; % Maximum number of neurons in each graph
% Tau = 0.1; 
TauGHSOM = [0.01 0.001]; %% GHSOM
ErrorType='abs'; % GHSOM
Mu = 1; % penalization parameter for some dataset features
% Beta = (0.1:0.1:2); % number of standard deviations to be considered foreground
Features = 3; % only for the RGB color space: 1=RGB, 2=Normalized RGB, and 3=Normalized RGB median filtered
WindowSizeRGBMedian = 5; % window size for the normalized RGB median filtered
NomFichModelos = 'Modelos';
NomFichResultados = 'ResultadosSegmentacionGHSOM';
PathDatos = 'Datos/Mediana 100 últimos sin NaN RGB Normalizado Mediana 5x5/';

% % The following values of the parameters are those considered in the
% % original GNG paper by Fritzke (1995)
% Lambda = 100;
% EpsilonB = 0.2;
% EpsilonN = 0.006;
% Alpha = 0.5;
% AMax = 50;
% D = 0.995;

%% COLOR SPACES
ColorSpaces = {'HSL','HSV','Lab','Luv','RGB','YCbCr'};
NumColorSpaces = length(ColorSpaces);

%% VIDEOS TO BE PROCESSED
Videos = {'Campus','Curtain','Escalator','Fountain','LevelCrossing','OneShopOneWait1cor',...
    'Video2','Video4','WaterSurface','Lobby','LightSwitch'};
NumVideos = length(Videos);
PathVideos = '../../Videos/';
SubPathFrames = '/frames/f%07d.bmp';
SubPathFramesLightSwitch = '/frames/f%07d.jpg';
SubPathGT = '/GT/GT%07d.bmp';
DeltaFrames = [999,20999,1397,999,0,0,0,0,999,999,0];
NumFrames = [1438,2963,3417,523,500,1376,749,819,632,1545,2799];
NumGTFrames = [20,20,20,20,163,11,444,316,20,20,21];
% VideosValidos = [2,6,7,9];
VideosValidos = 1:11;

%% TRAINING AND TEST
if exist(['./' NomFichResultados '.mat'],'file'),
    load(['./' NomFichResultados '.mat'],'Results');
else
    Results = cell(NumColorSpaces,NumVideos);
end
for NdxColorSpace=1:NumColorSpaces,
% for NdxColorSpace=5:5,
    
    MyColorSpace = ColorSpaces{NdxColorSpace};
    fprintf('\nCOLOR SPACE: %s\n',MyColorSpace);
    
    if exist(['./' NomFichModelos '_' MyColorSpace '.mat'],'file'),
        load(['./' NomFichModelos '_' MyColorSpace '.mat'],'Modelos');
    else
        Modelos = cell(1,NumVideos);
    end
    
    for NdxVideo=1:NumVideos,
        
        if isempty(find(VideosValidos==NdxVideo, 1)),
            continue;
        end                              
        MyVideo = Videos{NdxVideo};
        fprintf('\n\tVIDEO: %s\n',MyVideo);                        
                
        if isempty(Modelos{NdxVideo}),                        
            % Load dataset corresponding to the color space and video sequence
            load(['./' PathDatos 'Dataset_' MyColorSpace '_' MyVideo '.mat'],'Data');
            fprintf('\n\t\t Dataset: %s\n',PathDatos);
            Muestras = Data;
            NumMuestras = size(Data,2);
%             NumPasos = NumEpocas*NumMuestras;
        
            % GHSOM Training
            fprintf('\n\t\tENTRENAMIENTO GHSOM\n');
            fprintf('------------------------------------\n');            
            % Compute the mean of the samples (initial weight vector)
            IniWeight = full(nansum(Muestras,2)/NumMuestras);

            % Compute the initial Bregman quantization errror (qe0)
            IniQE = QuantizationError(IniWeight,Muestras,ErrorType);
            tic;
            [Modelo] = TrainGHSOM(Muestras,(1:NumMuestras),NumEpocas,TauGHSOM(1),TauGHSOM(2),ErrorType,IniQE,IniQE,[],1);
            CpuTime = toc;
            Modelo.Samples = Muestras;
            if isempty(Modelo),
                continue;
            end
            Modelos{NdxVideo} = Modelo;        
            save([NomFichModelos '_' MyColorSpace '.mat'],'Modelos');
        else
            fprintf('\n\t\tModelo GHSOM Entrenado\n');        
            Modelo = Modelos{NdxVideo};
            CpuTime = 0;
        end
        
        if (isempty(Results{NdxColorSpace,NdxVideo})),
            % Test GHNF
            [Winners,Errors] = TestGHNF(GetCentroidsGHSOM(Modelo),Modelo.Samples);

            % Evaluate Performance        
            if NdxVideo==11, % jpg instead of bmp extension
                PathMyVideo = [PathVideos MyVideo SubPathFramesLightSwitch];
            else
                PathMyVideo = [PathVideos MyVideo SubPathFrames];
            end
            PathMyGT = [PathVideos MyVideo SubPathGT];           
            fprintf('\t\tTesting Frames...\n');
            PerformanceMeasures = TestGHSOMSegmentation(Modelo,Winners,Errors,...
                MyColorSpace,MyVideo,NumFrames(NdxVideo),NumGTFrames(NdxVideo),...
                DeltaFrames(NdxVideo),PathMyVideo,PathMyGT,Mu,Features,WindowSizeRGBMedian,GetNumberNeuronsGHSOM(Modelo)); 

            % Time is the training time + the mean detecting foreground time           
            Results{NdxColorSpace,NdxVideo}.Time = CpuTime + PerformanceMeasures(1,end);    
            MeanPerformanceMeasures = mean(PerformanceMeasures(:,1:end-1),2);
            StdPerformanceMeasures = std(PerformanceMeasures(:,1:end-1),0,2);
            Results{NdxColorSpace,NdxVideo}.S.Mean = MeanPerformanceMeasures(1);
            Results{NdxColorSpace,NdxVideo}.S.Std = StdPerformanceMeasures(1);
            Results{NdxColorSpace,NdxVideo}.FP.Mean = MeanPerformanceMeasures(2);
            Results{NdxColorSpace,NdxVideo}.FP.Std = StdPerformanceMeasures(2);
            Results{NdxColorSpace,NdxVideo}.FN.Mean = MeanPerformanceMeasures(3);
            Results{NdxColorSpace,NdxVideo}.FN.Std = StdPerformanceMeasures(3);
            Results{NdxColorSpace,NdxVideo}.Precision.Mean = MeanPerformanceMeasures(4);
            Results{NdxColorSpace,NdxVideo}.Precision.Std = StdPerformanceMeasures(4);
            Results{NdxColorSpace,NdxVideo}.Recall.Mean = MeanPerformanceMeasures(5);
            Results{NdxColorSpace,NdxVideo}.Recall.Std = StdPerformanceMeasures(5);
            Results{NdxColorSpace,NdxVideo}.Accuracy.Mean = MeanPerformanceMeasures(6);
            Results{NdxColorSpace,NdxVideo}.Accuracy.Std = StdPerformanceMeasures(6);
            Results{NdxColorSpace,NdxVideo}.Fmeasure.Mean = MeanPerformanceMeasures(7); 
            Results{NdxColorSpace,NdxVideo}.Fmeasure.Std = StdPerformanceMeasures(7);

            save([NomFichResultados '.mat'],'Results');   
        else
            fprintf('\t\tTesteo Realizado\n');
        end
    end
end
