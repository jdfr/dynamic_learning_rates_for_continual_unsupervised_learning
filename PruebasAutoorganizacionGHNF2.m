% Pruebas Autoorganización GHNF con Distribuciones 3D
clear all
NumMuestras = 10000;
NumEpocas = 2;
MaxNeurons = 20; % Maximum number of neurons in each graph
Tau = [0 0.1 0.2];
NomFichResultados = 'ResultadosAutoorganizacionGHNF.mat';

% The following values of the parameters are those considered in the
% original GNG paper by Fritzke (1995)
Lambda = 100;
EpsilonB = 0.2;
EpsilonN = 0.006;
Alpha = 0.5;
AMax = 50;
D = 0.995;

Nombres = {'SwissRoll','SwissHole','CornerPlanes','PuncturedSphere',...
   'TwinPeaks','3D Clusters','ToroidalHelix','Gaussian','UedaSpiral'};

% Cargamos los resultados o los inicializamos
if exist(['./' NomFichResultados],'file') == 0
    Time = nan(length(Nombres),length(Tau));
    MSE = nan(length(Nombres),length(Tau));
    PSNR = nan(length(Nombres),length(Tau));
    DB = nan(length(Nombres),length(Tau));
    CS = nan(length(Nombres),length(Tau));
else
    load(['./' NomFichResultados],'MSE','PSNR','DB','CS','Time');
end

for dataset = 1:length(Nombres)
        
    NomFich = strrep(Nombres{dataset},'.bmp','');
    fprintf('\nIMAGEN: %s\n',NomFich);
    Muestras = Generate3DSamples(dataset,NumMuestras);
    
    for NdxTau=1:length(Tau)
        
        MiTau = Tau(NdxTau);        
        fprintf('\nTau = %g\n',MiTau);                                          
        NomFichModelo = ['./GHNF_' NomFich '_' num2str(MiTau)];         

        if exist([NomFichModelo '.mat'],'file') == 0

            fprintf('\nGHNF TRAINING\n');
            fprintf('------------------------------------\n');            
            tic;
            [Modelo] = TrainGHNF(Muestras,NumEpocas,MaxNeurons,MiTau,Lambda,EpsilonB,EpsilonN,Alpha,AMax,D,1);
            Time(dataset,NdxTau) = toc;
            save([NomFichModelo '.mat'],'Modelo');                     
        else
            fprintf('\tTraining already done\n');           
            load([NomFichModelo '.mat'],'Modelo');
        end
        
        if isempty(Modelo)
            fprintf('\tEmpty model!\n');
        else
            if isnan(MSE(dataset,NdxTau)) || isnan(PSNR(dataset,NdxTau)) ||...
                    isnan(DB(dataset,NdxTau)) || isnan(CS(dataset,NdxTau))
                % Calculamos las medidas de rendimiento
                fprintf('\tComputing performance measures...\n');
                Centroids = GetCentroidsGHNF(Modelo);
                [Winners,Errors] = TestGHNF(Centroids,Modelo.Samples);
                MSE(dataset,NdxTau) = mean(Errors);
                MaxI = max(max(Modelo.Samples)) - min(min(Modelo.Samples));
                PSNR(dataset,NdxTau) = 10*log10(MaxI^2/MSE(dataset,NdxTau));
                DB(dataset,NdxTau) = db_index(Modelo.Samples',Winners,Centroids');
                CS(dataset,NdxTau) = nanmean(silhouetteMEX(Modelo.Samples,Winners'-1,max(Winners)));    
            end
            
            if (exist(['./' NomFichModelo '.fig'],'file') == 0) || (exist(['./' NomFichModelo '.pdf'],'file') == 0)
                % Dibujamos el modelo resultante
                Handle = Plot3GHNF(Modelo);
                hold on
                h = plot3(Muestras(1,:),Muestras(2,:),Muestras(3,:),'*y');     
                uistack(h,'bottom');
        %         axis([0 1 0 1])
                hgsave(gcf,[NomFichModelo '.fig']);
                fig = open([NomFichModelo '.fig']);
                set(gcf,'PaperUnits','centimeters');
                set(gcf,'PaperOrientation','portrait');
                set(gcf,'PaperPositionMode','manual');
                set(gcf,'PaperSize',[12 10]);
                set(gcf,'PaperPosition',[0 0 12 10]);
                axis off;
                saveas(gcf,[NomFichModelo '.pdf'],'pdf'); 
            end
        end
    end    
    close all;
end
save(['./' NomFichResultados],'MSE','PSNR','DB','CS','Time');