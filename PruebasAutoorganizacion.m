% Pruebas de autoorganización para los modelos (GHNF, GHNG, GNF).
% Se ejecuta en la carpeta raíz ('Pruebas Autoorganización'), donde
% debe haber una carpeta 'Datos' y otra 'Resultados'. Dentro de
% 'Resultados' debe haber una subcarpeta por modelo con el nombre del
% modelo a entrenar. El modelo a entrenar debe estar en el path.
% Modificación del 26/05/2020 para incorporar el modelo GHSOM
clearvars

TipoModelo = 'GHSOM'; % GHNF, GHNG, GNF, GHSOM
NumMuestras = 10000;
NumEpocas = 2;
MaxNeurons = 20; % Maximum number of neurons in each graph
Tau = [0 0.1 0.2];
NumTau = length(Tau);
TauGHSOM = [0.001 0.01; 0.01 0.001; 0.001 0.001]; %% GHSOM
ErrorType='abs'; % GHSOM
Datasets = {'Ring','Circle','Hollowed Square','Non-hollowed Square','M',... % 2D distributions
    'X','S','Barbell','Barbell2','K','Eight','Sand Clock',... 
    'Ball','Big Ball','Ellipse','Hexagon','Octagon','Smooth Ball','Star',... % 2D gradient distributions
    'Thick S','Thin S','Three Balls','Two Shapes','Irregular X',... 
    'SwissRoll','SwissHole','CornerPlanes','PuncturedSphere',... % 3D distributions
   'TwinPeaks','3D Clusters','ToroidalHelix','Gaussian','UedaSpiral',...
   'UnitCircle','UnitSquare','UnitSegment','UnitSphere',... % 3D simple distributions
   'Trefoil Knot','Cinquefoil Knot','Septafoil Knot',... % 3D Torus Knots
   'Hypersphere 4D','Hypersphere 5D','Hypersphere 6D'};  % D-dimensional stuffed hyperspheres
NumDatasets = length(Datasets);
DistribucionDatasets = [12,24,33,37,40];
DatasetsEntrenamiento = [1,6,7,16,29,30,36,41,42,43];
TorusParams = [2 3; 5 2; 7 2];
HypersphereDims = [4,5,6];
PathDatos = '/Datos/';
PathResultados = ['/Resultados/' TipoModelo '/'];
NomFichModelos = ['ModelosAutoorganizacion' TipoModelo '.mat'];
NomFichEvaluaciones = ['EvaluacionesAutoorganizacion' TipoModelo '.mat'];
NomFichResultados = ['ResultadosAutoorganizacion' TipoModelo '.mat'];
NumRepeticiones = 10;
NumIntentos = 10; % #intentos para entrenar una red que devuelve un modelo vacío

% The following values of the parameters are those considered in the
% original GNG paper by Fritzke (1995)
Lambda = 100;
% Lambda = 40; % this parameter has been changed in order to be compared with other models in terms of number of neurons.
EpsilonB = 0.2;
EpsilonN = 0.006;
Alpha = 0.5;
AMax = 50;
D = 0.995;

% Cargamos las evaluaciones o las inicializamos
if exist(['.' PathResultados NomFichEvaluaciones],'file') == 0
    Evaluaciones = cell(NumDatasets,NumTau,NumRepeticiones);    
else
    load(['.' PathResultados NomFichEvaluaciones],'Evaluaciones');
end

% Cargamos los resultados o los inicializamos
if exist(['.' PathResultados NomFichResultados],'file') == 0
    Resultados = cell(NumDatasets,NumTau);    
else
    load(['.' PathResultados NomFichResultados],'Resultados');
end

% Cargamos los modelos o los inicializamos
if exist(['.' PathResultados NomFichModelos],'file') == 0
    Modelos = cell(NumDatasets,NumTau);    
else
    load(['.' PathResultados NomFichModelos],'Modelos');
end

% for NdxDataset=1:NumDatasets
for NdxDataset=DatasetsEntrenamiento
        
    NomFich = Datasets{NdxDataset};    
    if NdxDataset <= DistribucionDatasets(2)
        if NdxDataset <= DistribucionDatasets(1)
            Extension = '.bmp';
        else
            Extension = '.jpg';
        end    
        Muestras = GenerateSamplesImg(['.' PathDatos NomFich Extension],NumMuestras);       
    elseif NdxDataset <= DistribucionDatasets(3)
        Muestras = Generate3DSamples(NdxDataset-DistribucionDatasets(2),NumMuestras);
    elseif NdxDataset <= DistribucionDatasets(4)
        Muestras = GenerateSimple3D(Datasets{NdxDataset},0,NumMuestras);
    elseif NdxDataset <= DistribucionDatasets(5)
        p = TorusParams(NdxDataset-DistribucionDatasets(4),1);
        q = TorusParams(NdxDataset-DistribucionDatasets(4),2);
        Muestras = GenerateSamplesTorusKnot(p,q,NumMuestras);        
    else
        Dimension = HypersphereDims(NdxDataset-DistribucionDatasets(5));
        Muestras = insideSphereGen(Dimension,NumMuestras);
    end
    
    for NdxTau=1:NumTau           
        MiTau = Tau(NdxTau);                 
        MiTau1 = TauGHSOM(NdxTau,1);
        MiTau2 = TauGHSOM(NdxTau,2);
        % For GNF, we compute the number of neurons as the average between
        % the mean number of neurons of the GHNF and the mean number of
        % neurons of the GNF. GHNF and GHNG must be previously trained
        if strcmp(TipoModelo,'GNF')
            S1 = load('./Resultados/GHNF/ResultadosAutoorganizacionGHNF.mat');
            S2 = load('./Resultados/GHNG/ResultadosAutoorganizacionGHNG.mat');
            MiMaxNeurons = round((S1.Resultados{NdxDataset,NdxTau}.NumNeurons.Media +...
                S2.Resultados{NdxDataset,NdxTau}.NumNeurons.Media)/2);
        end
        if ~isempty(Modelos{NdxDataset,NdxTau}) 
            fprintf('\tTrainings already done!\n');  
        else
            % Si un modelo está vacío, hay que hacer todos los
            % entrenamientos (no nos podemos ahorra ninguno ya hecho),
            % ya que solo se guarda el mejor modelo y para ello hay que
            % comparar todos los NumRepeticiones modelos generados.
            MejorMSE = Inf;            
            for NdxRepeticion=1:NumRepeticiones
                fprintf('\nDISTRIBUTION: %s\n',NomFich);
                fprintf('\nTAU = %g\n',MiTau);  
                fprintf('\n%s TRAINING %d\n',TipoModelo,NdxRepeticion);
                fprintf('------------------------------------\n');            
                Modelo = [];
                intentos = 0;
                while isempty(Modelo) && intentos<NumIntentos
                    if strcmp(TipoModelo,'GHNF')
                        tic;
                        [Modelo] = TrainGHNF(Muestras,NumEpocas,MaxNeurons,MiTau,Lambda,EpsilonB,EpsilonN,Alpha,AMax,D,1);
                        CpuTime = toc;  
                    elseif strcmp(TipoModelo,'GHNG')
                        tic;
                        [Modelo] = TrainGHNG3(Muestras,NumEpocas,MaxNeurons,MiTau,Lambda,EpsilonB,EpsilonN,Alpha,AMax,D,1);
                        CpuTime = toc;  
                    elseif strcmp(TipoModelo,'GNF')
                        tic;
                        [Modelo] = TrainGNF(Muestras,MiMaxNeurons,Lambda,EpsilonB,EpsilonN,Alpha,AMax,D,NumMuestras*NumEpocas);
                        CpuTime = toc;
                    elseif strcmp(TipoModelo,'GHSOM')
                        % Compute the mean of the samples (initial weight vector)
                        IniWeight = full(nansum(Muestras,2)/NumMuestras);

                        % Compute the initial Bregman quantization errror (qe0)
                        IniQE = QuantizationError(IniWeight,Muestras,ErrorType);
                        tic;
                        [Modelo] = TrainGHSOM(Muestras,(1:NumMuestras),NumEpocas,MiTau1,MiTau2,ErrorType,IniQE,IniQE,[],1);
                        CpuTime = toc;
                        Modelo.Samples = Muestras;
                    end
                    intentos = intentos + 1;
                end

                if isempty(Modelo)
                    fprintf('\tImpossible training!\n');
                else                        
                    fprintf('\tTraining performed in %d attempts!\n',intentos);
                    
                    % Calculamos las medidas de rendimiento
                    fprintf('\tComputing performance measures %d...\n',NdxRepeticion);
                    if strcmp(TipoModelo,'GHNF')
                        Centroids = GetCentroidsGHNF(Modelo);
                        [Winners,Errors] = TestGHNF(Centroids,Modelo.Samples);
                    elseif strcmp(TipoModelo,'GHNG')
                        Centroids = GetCentroidsGHNG(Modelo);
                        [Winners,Errors] = TestGHNG(Centroids,Modelo.Samples);
                    elseif strcmp(TipoModelo,'GNF')
                        Centroids = Modelo.Means;
                        [Winners,Errors] = TestGNG(Modelo,Modelo.Samples);
                    elseif strcmp(TipoModelo,'GHSOM')
                        Centroids = GetCentroidsGHSOM(Modelo);
                        [Winners,Errors] = TestGHNF(Centroids,Muestras);
                    end
%                     [DBI,MSE,Dunn,MeanSilh]=EvaluateClustering(Centroids,Modelo.Samples)
                    Evaluaciones{NdxDataset,NdxTau,NdxRepeticion}.MSE = mean(Errors);
                    MaxI = max(max(Modelo.Samples)) - min(min(Modelo.Samples));
                    Evaluaciones{NdxDataset,NdxTau,NdxRepeticion}.PSNR = 10*log10(MaxI^2/Evaluaciones{NdxDataset,NdxTau,NdxRepeticion}.MSE);
                    Evaluaciones{NdxDataset,NdxTau,NdxRepeticion}.DB = db_index(Modelo.Samples',Winners,Centroids');
                    Evaluaciones{NdxDataset,NdxTau,NdxRepeticion}.CS = nanmean(silhouetteMEX(Modelo.Samples,Winners'-1,max(Winners)));
                    if strcmp(TipoModelo,'GHSOM')
                        [~,Evaluaciones{NdxDataset,NdxTau,NdxRepeticion}.NumNeurons] = GetNumberNeuronsGHSOM(Modelo);
                    else
                        [~,Evaluaciones{NdxDataset,NdxTau,NdxRepeticion}.NumNeurons] = GetNumberNeuronsGHNF(Modelo);
                    end
                    Evaluaciones{NdxDataset,NdxTau,NdxRepeticion}.Time = CpuTime;
                    save(['.' PathResultados NomFichEvaluaciones],'Evaluaciones');

                    % Nos vamos quedando con el mejor modelo (aquél con menor MSE)
                    if Evaluaciones{NdxDataset,NdxTau,NdxRepeticion}.MSE < MejorMSE
                        MejorMSE = Evaluaciones{NdxDataset,NdxTau,NdxRepeticion}.MSE;
                        Modelos{NdxDataset,NdxTau} = Modelo;
                    end                    
                    clc
                end
            end            
            save(['.' PathResultados NomFichModelos],'Modelos');          
        end              
        
        if isempty(Modelos{NdxDataset,NdxTau})
            fprintf('\tEmpty model!\n');
        else            
            % Calculamos la media y la desviación típica de cada medida
            fprintf('\tComputing means and standard deviations...\n');
            Resultados{NdxDataset,NdxTau} = ValidacionCruzada(Evaluaciones(NdxDataset,NdxTau,:));                
            
            if strcmp(TipoModelo,'GHSOM')
                NomFichFigura = [TipoModelo '_' NomFich '_' num2str(MiTau1) '_' num2str(MiTau2)];                 
            else
                NomFichFigura = [TipoModelo '_' NomFich '_' num2str(MiTau)]; 
            end
            if (exist(['.' PathResultados NomFichFigura '.fig'],'file') == 0) || (exist(['.' PathResultados NomFichFigura '.pdf'],'file') == 0)
                % Plot the resulting model
                fprintf('\tPlotting results...\n');
                if NdxDataset <= DistribucionDatasets(2)
                    Handle = figure;
                    Image = fliplr(rot90(double(rgb2gray(imread(['.' PathDatos NomFich Extension])))/255,2)); 
                    h = imshow(Image);
                    hold on
                    Ejes = axis(gca);
                    Ejes([1 3]) = round(Ejes([1 3]));
                    Ejes([2 4]) = floor(Ejes([2 4]))-1;
                    axis tight;
                    axis xy;
                    if strcmp(TipoModelo,'GHNF')
                        PlotGHNF(Modelos{NdxDataset,NdxTau},Ejes);
                    elseif strcmp(TipoModelo,'GHNG')
                        PlotGHNG(Modelos{NdxDataset,NdxTau},Ejes);
                    elseif strcmp(TipoModelo,'GNF')
                        PlotGNFScaled(Modelos{NdxDataset,NdxTau},Ejes);
                    else
                        PlotGHSOMScaled(Modelos{NdxDataset,NdxTau},Ejes);
                    end
                    axis off
                elseif NdxDataset <= DistribucionDatasets(5)
                    if strcmp(TipoModelo,'GHNF')
                        Handle = Plot3GHNF(Modelos{NdxDataset,NdxTau});
                    elseif strcmp(TipoModelo,'GHNG')
                       Handle = Plot3GHNG(Modelos{NdxDataset,NdxTau});
                    elseif strcmp(TipoModelo,'GNF')                        
                        PlotGNF3D(Modelos{NdxDataset,NdxTau});
                    else
                        Handle = Plot3GHSOM(Modelos{NdxDataset,NdxTau});
                    end
                    hold on
                    h = plot3(Muestras(1,:),Muestras(2,:),Muestras(3,:),'*y');     
                    uistack(h,'bottom');
                end                                     
                hgsave(gcf,['./' PathResultados NomFichFigura '.fig']);
                % Save the plot in a PDF file
                set(gcf,'PaperUnits','centimeters');
                set(gcf,'PaperOrientation','portrait');
                set(gcf,'PaperPositionMode','manual');
                set(gcf,'PaperSize',[12 10]);
                set(gcf,'PaperPosition',[0 0 12 10]);
                axis off;
                saveas(gcf,['.' PathResultados NomFichFigura '.pdf'],'pdf');       
            end
        end
    end    
    close all;
end
save(['.' PathResultados NomFichResultados],'Resultados');