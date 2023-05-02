% Script para generar los resultados de autoorganización para distintos
% modelos (GHNF, GHNG, GNF y GHSOM). 
% Ejecutar desde la carpeta donde se encuentre el fichero de resultados.
clc
clearvars

Modelo = 'GHSOM';
Datasets = {'Ring','Circle','Hollowed Square','Non-hollowed Square','M',... % 2D distributions
    'X','S','Barbell','Barbell2','K','Eight','Sand Clock',... 
    'Ball','Big Ball','Ellipse','Hexagon','Octagon','Smooth Ball','Star',... % 2D gradient distributions
    'Thick S','Thin S','Three Balls','Two Shapes','Irregular X',... 
    'SwissRoll','SwissHole','CornerPlanes','PuncturedSphere',... % 3D distributions
   'TwinPeaks','3D Clusters','ToroidalHelix','Gaussian','UedaSpiral',...
   'UnitCircle','UnitSquare','UnitSegment','UnitSphere',... % 3D simple distributions
   'Trefoil Knot','Cinquefoil Knot','Septafoil Knot',... % 3D Torus Knots
   'Hypersphere 4D','Hypersphere 5D','Hypersphere 6D'};  % D-dimensional stuffed hyperspheres
DatasetsSeleccionados = [1,6,7,16,29,30,36,41,42,43];
NumDatasets = length(Datasets);
Medidas = {'MSE','PSNR','DB','CS','NumNeurons','Time'};
NumMedidas = length(Medidas);
Tau = [0 0.1 0.2];
NumTau = length(Tau);

load(['ResultadosAutoorganizacion' Modelo '.mat'],'Resultados');

for NdxTau=1:NumTau    
    MiTau = Tau(NdxTau);
    fprintf('\nTAU=%d\n',MiTau);
    ResultadosFormateados = cell(NumDatasets,NumMedidas);
    
    for NdxMedida=1:NumMedidas    
        MiMedida = Medidas{NdxMedida};
        fprintf('\nMedida: %s\n',MiMedida);

%         for NdxDataset=1:NumDatasets
%         ContDataset = 1;
        for NdxDataset=DatasetsSeleccionados
            MiDataset = Datasets{NdxDataset};
            fprintf('\nDataset: %s\n',MiDataset);
            
            if isempty(Resultados{NdxDataset,NdxTau})
                fprintf('Resultados no encontrados para el dataset %s y tau=%d\n',MiDataset,MiTau);
            else
                MiCampo = Resultados{NdxDataset,NdxTau}.(MiMedida); 
                ResultadosFormateados{NdxDataset,NdxMedida} = PrettyNumbers(MiCampo.Media,MiCampo.DesvTip,0);
%                 ResultadosFormateados{MiDataset,MiTau} = sprintf('%6.2f (%2.2f)',MiCampo.Mean,MiCampo.Std);                   
            end     
%             ContDataset = ContDataset + 1;
        end
    end
%     ResultadosFormateados = [{''} Medidas; [Datasets(DatasetsSeleccionados)' ResultadosFormateados]];
    ResultadosFormateados = [{''} Medidas; [Datasets' ResultadosFormateados]];
    xlswrite(['ResultadosAutoorganización' Modelo '.xlsx'],ResultadosFormateados,num2str(MiTau));
end
