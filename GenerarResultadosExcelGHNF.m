clc
clear

Videos = {'Campus','Curtain','Escalator','Fountain','LevelCrossing','OneShopOneWait1cor',...
    'Video2','Video4','WaterSurface','Lobby','LightSwitch'};
EspaciosColor = {'HSL','HSV','Lab','Luv','RGB','YCbCr'};
EspaciosColorOrdenados = {'RGB','Lab','Luv','HSV','HSL','YCbCr'};
VideosOrdenados = {'Video2','Curtain','WaterSurface','Lobby','OneShopOneWait1cor','LevelCrossing'};
Medidas = {'S','FP','FN','Precision','Recall','Accuracy','Fmeasure','Time'};
VideosValidos = [2,5,6,7,9,10];

load('ResultadosSegmentacionGHNG.mat','Results');

for NdxMedida=1:length(Medidas),
    
    MiMedida = Medidas{NdxMedida};
    fprintf('\nMedida: %s\n',MiMedida);
    ResultadosFormateados = cell(length(EspaciosColorOrdenados),length(VideosValidos));
    
    for NdxEspacioColor=1:length(EspaciosColor),

        MiEspacioColor = EspaciosColor{NdxEspacioColor};
        fprintf('\nEspacio Color: %s\n',MiEspacioColor);
        ContVideosValidos = 0;
        for NdxVideo=1:length(Videos),

            if ismember(NdxVideo,VideosValidos),
                MiVideo = Videos{NdxVideo};      
                fprintf('\nVídeo: %s\n',MiVideo);                
                ContVideosValidos = ContVideosValidos+1;
                if isempty(Results{NdxEspacioColor,NdxVideo}),
                    fprintf('Resultados no encontrados para el video ''%s''\n',MiVideo);
                else
                    MiCampo = getfield(Results{NdxEspacioColor,NdxVideo},MiMedida);
                    [~,MiNdxEspacioColor] = ismember(MiEspacioColor,EspaciosColorOrdenados);
                    [~,MiNdxVideo] = ismember(MiVideo,VideosOrdenados);
                    
                    if (strcmp(MiMedida,Medidas{end})==0),                
                        ResultadosFormateados{MiNdxEspacioColor,MiNdxVideo} = sprintf('%6.2f (%2.2f)',MiCampo.Mean,MiCampo.Std);   
                    else
                        ResultadosFormateados{MiNdxEspacioColor,MiNdxVideo} = sprintf('%6.3f',MiCampo);
                    end
                end
            end
        end
    end
%     ResultadosFormateados = [{''} Videos; [EspaciosColor' ResultadosFormateados]];
%     ResultadosFormateados = ResultadosFormateados{:,VideosValidos};
    ResultadosFormateados = [{''} EspaciosColorOrdenados; [VideosOrdenados' ResultadosFormateados']];
    xlswrite('ResultadosSegmentaciónGHNG.xlsx',ResultadosFormateados,MiMedida);
end
