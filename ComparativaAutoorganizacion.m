% Script para comparar los resultados de autoorganización de distintos
% modelos jerárquicos. Ejecutar en la carpeta 'Pruebas autoorganización'.
clc
clearvars

Modelos = {'GHNF','GHNG'};
NumModelos = length(Modelos);
Datasets = {'Ring','Circle','Hollowed Square','Non-hollowed Square','M letter',...
    'X letter','S letter','Barbell','Barbell 2','K letter','Eight number','Sand clock',... % 2D distributions
    'Ball','Big Ball','Ellipse','Hexagon','Octagon','Smooth Ball','Star',...
    'ThickS','ThinS','Three Balls','Two Shapes','X',... % 2D gradient distributions
    'SwissRoll','SwissHole','CornerPlanes','PuncturedSphere',...
   'TwinPeaks','3D Clusters','ToroidalHelix','Gaussian','UedaSpiral'}; % 3D distributions
NumDatasets = length(Datasets);
Medidas = {'MSE','PSNR','DB','CS','NumNeurons','Time'};
NumMedidas = length(Medidas);
Tau = [0 0.1 0.2];
NumTau = length(Tau);
PathResultados = '/Resultados/';
NomFichResultados = 'ResultadosAutoorganizacion';

ResultadosModelos = cell(1,NumModelos);
for NdxModelo=1:NumModelos
    MiModelo = Modelos{NdxModelo};
    MiPath = [PathResultados MiModelo '/'];
    load(['.' MiPath NomFichResultados MiModelo '.mat'],'Resultados');
    ResultadosModelos{NdxModelo} = Resultados;
end

for NdxTau=1:NumTau    
    MiTau = Tau(NdxTau);
    fprintf('\nTAU=%d\n',MiTau);
    ResultadosFormateados = cell(NumDatasets,NumMedidas);
    
    for NdxMedida=1:NumMedidas    
        MiMedida = Medidas{NdxMedida};
        fprintf('\nMedida: %s\n',MiMedida);        

        for NdxDataset=1:NumDatasets
            MiDataset = Datasets{NdxDataset};
            fprintf('\nDataset: %s\n',MiDataset);
            
            MaxPSNR = -Inf;
            NdxMejorModelo = 1;
            for NdxModelo=1:NumModelos
                Resultados = ResultadosModelos{NdxModelo};
                if isempty(Resultados{NdxDataset,NdxTau})
                    fprintf('Resultados no encontrados para el dataset %s y tau=%d\n',MiDataset,MiTau);
                else
                    if Resultados{NdxDataset,NdxTau}.PSNR.Media > MaxPSNR
                        MaxPSNR = Resultados{NdxDataset,NdxTau}.PSNR.Media;
                        NdxMejorModelo = NdxModelo;
                    end                    
                end     
            end
            if NdxMejorModelo==1
                fprintf('\tMejor Modelo: %s!!!!!!!!!!!!!!!!!!!!!!\n',Modelos{NdxMejorModelo});
            else
                fprintf('\tMejor Modelo: %s\n',Modelos{NdxMejorModelo});
            end
        end
    end
%     ResultadosFormateados = [{''} Medidas; [Datasets(DatasetsSeleccionados)' ResultadosFormateados]];
%     xlswrite(['ResultadosAutoorganización' Modelo '.xlsx'],ResultadosFormateados,num2str(MiTau));
end
