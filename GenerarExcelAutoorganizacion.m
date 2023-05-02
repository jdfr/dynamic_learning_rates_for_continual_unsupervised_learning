% Script que genera el excel de resultados de las pruebas de 
% autoorganización para el modelo indicado
% Ejecutar desde el directorio 'Pruebas de Autoorganización'
% E.J. Palomo. July 2017
clear
Modelo = 'GHNF';
Distribuciones = {'2D','2D densidad','3D'};
NumDistribuciones = length(Distribuciones);
Datasets = {'Ring','Circle','Hollowed Square','Non-hollowed Square','M letter',...
    'X letter','S letter','Barbell','Barbell 2','K letter','Eight number','Sand clock',... % 2D distributions
    'Ball','Big Ball','Ellipse','Hexagon','Octagon','Smooth Ball','Star',...
    'ThickS','ThinS','Three Balls','Two Shapes','X',... % 2D gradient distributions
    'SwissRoll','SwissHole','CornerPlanes','PuncturedSphere',...
   'TwinPeaks','3D Clusters','ToroidalHelix','Gaussian','UedaSpiral'}; % 3D distributions
NumDatasets = length(Datasets);
NumDatasetsDistribucion = [12,12,9];
Medidas = {'MSE','PSNR','DB','CS','Time'};
NumMedidas = length(Medidas);
Tau = [0 0.1 0.2];
NumTau = length(Tau);

for NdxTau=1:NumTau    
    MiTau = Tau(NdxTau);
    fprintf('\nTAU=%d\n',MiTau);    
    Resultados = cell(NumDatasets,NumMedidas);
    DesplazamientoDataset = 0;
    
    for NdxDistribucion=1:NumDistribuciones
        MiDistribucion = Distribuciones{NdxDistribucion};
        MiPathResultados = ['./' MiDistribucion '/'];
        S = load([MiPathResultados 'ResultadosAutoorganizacion' Modelo '.mat']);
        for NdxMedida=1:NumMedidas   
            MiMedida = Medidas{NdxMedida};
            M = getfield(S,MiMedida);
            MiNumDatasets = NumDatasetsDistribucion(NdxDistribucion);
            for NdxDataset=1:MiNumDatasets                                             
                Resultados{NdxDataset+DesplazamientoDataset,NdxMedida} = sprintf('%6.2f',M(NdxDataset,NdxTau));
            end            
        end
        DesplazamientoDataset = DesplazamientoDataset + MiNumDatasets;
    end
    Resultados = [{''} Medidas; [Datasets' Resultados]];
    xlswrite(['ResultadosAutoorganización' Modelo '.xlsx'],Resultados,num2str(MiTau));
end
        

