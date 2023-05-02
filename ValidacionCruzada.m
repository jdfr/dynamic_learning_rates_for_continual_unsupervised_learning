function [EvaluacionValidada]=ValidacionCruzada(Evaluaciones)
% Evaluar el rendimiento del clustering usando validación cruzada
% Entrada:
%   Evaluaciones=Evaluaciones obtenidas para un número indeterminado de
%   medidas de rendimiento
% Salida:
%   EvaluacionValidada=Resumen de las evaluaciones de esas medidas de
%   rendimiento

if isempty(Evaluaciones)
    EvaluacionValidada = [];
else
    NumDatos = length(Evaluaciones);
    Medidas = fieldnames(Evaluaciones{1,1});
    NumMedidas = length(Medidas);
    
    for NdxMedida=1:NumMedidas    
        MiMedida = Medidas{NdxMedida};
        Vector = zeros(1,NumDatos);
        for ndx=1:NumDatos
            if ~isempty(Evaluaciones{ndx})
                Vector(ndx) = Evaluaciones{ndx}.(MiMedida);
            else
                Vector(ndx) = NaN;
            end
        end
        EvaluacionValidada.(MiMedida).Media = mean(Vector(isfinite(Vector)));
        EvaluacionValidada.(MiMedida).DesvTip = std(Vector(isfinite(Vector)));
    end           
end

