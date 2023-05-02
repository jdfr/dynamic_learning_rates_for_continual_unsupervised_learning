% Pruebas Autoorganización GHNF con 2D
clear all
NumMuestras = 10000;
NumEpocas = 2;
MaxNeurons = 20; % Maximum number of neurons in each graph
Tau = [0 0.1 0.2 0.3 0.4 0.5];
Extension = '.jpg';

% The following values of the parameters are those considered in the
% original GNG paper by Fritzke (1995)
Lambda = 100;
EpsilonB = 0.2;
EpsilonN = 0.006;
Alpha = 0.5;
AMax = 50;
D = 0.995;

% Leemos los conjuntos de datos del directorio
d = dir(['*' Extension]);
TiemposEntrenamiento = zeros(length(d),length(Tau));

for dataset=1:length(d),
        
    NomFich = strrep(d(dataset).name,Extension,'');
    fprintf('\nIMAGEN: %s\n',NomFich);
    Muestras = GenerateSamplesImg(sprintf(['%s' Extension],NomFich),NumMuestras);
    
    for NdxTau=1:length(Tau),
        
        MiTau = Tau(NdxTau);        
        fprintf('\nTau = %g\n',MiTau);                                          
        NomFichSalida = ['GHNF_' NomFich '_' num2str(MiTau)];         

        if exist(['./' NomFichSalida '.mat'],'file') == 0,

            fprintf('\nENTRENAMIENTO GHNF\n');
            fprintf('------------------------------------\n');            
            tic;
            [Modelo] = TrainGHNF(Muestras,NumEpocas,MaxNeurons,MiTau,Lambda,EpsilonB,EpsilonN,Alpha,AMax,D,1);
            TiemposEntrenamiento(dataset,NdxTau) = toc;
            save([NomFichSalida '.mat'],'Modelo');      
            save('TiemposEntrenamiento.mat','TiemposEntrenamiento');
        else
            fprintf('\nEntrenamiento realizado. Guardando figura como PDF...\n');           
            load(['./' NomFichSalida '.mat'],'Modelo');
            
            % Dibujamos el modelo resultante
            Handle = figure;
            Image = fliplr(rot90(double(rgb2gray(imread(d(dataset).name)))/255,2)); 
            h = imshow(Image);
            hold on
            Ejes = axis(gca);
            Ejes([1 3]) = round(Ejes([1 3]));
            Ejes([2 4]) = floor(Ejes([2 4]))-1;
            axis tight;
            axis xy;
            PlotGHNF(Modelo,Ejes);
            axis off
            hgsave(gcf,[NomFichSalida '.fig']);
            fig = open([NomFichSalida '.fig']);
            set(gcf,'PaperUnits','centimeters');
            set(gcf,'PaperOrientation','portrait');
            set(gcf,'PaperPositionMode','manual');
            set(gcf,'PaperSize',[12 10]);
            set(gcf,'PaperPosition',[0 0 12 10]);
            axis off;
            saveas(gcf,[NomFichSalida '.pdf'],'pdf');
        end        
    end    
    close all;
end
