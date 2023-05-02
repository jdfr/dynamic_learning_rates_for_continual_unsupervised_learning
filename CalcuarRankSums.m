% Cómo calcular los rank sums que aparecen en el paper en las tablas de
% detección del primer plano

Modelos = {'WrenGA','GrimsonGMM','GMMV1','GMMV2','FASOM','MFBM','GNF','GHSOM','GHNG','GHNF'};

% A es la tabla de datos
% A = [0.55 0.62 0.43 0.65 0.72 0.73 0.83 0.73 0.69 0.82;
%     0.54 0.44 0.46 0.54 0.78 0.77 0.73 0.64 0.61 0.77;
%     0.77 0.55 0.68 0.53 0.83 0.91 0.83 0.89 0.89 0.75; 
%     0.21 0.4 0.12 0.44 0.41 0.22 0.36 0.46 0.46 0.53;
%     0.67 0.38 0.53 0.37 0.61 0.66 0.64 0.31 0.32 0.67;
%     0.66 0.76 0.53 0.68 0.78 0.85 0.54 0.79 0.76 0.77];

A = [0.53 0.55 0.45 0.57 0.60 0.84 0.90 0.84 0.81 0.90;
    0.68 0.57 0.58 0.68 0.87 0.87 0.84 0.78 0.76 0.87;
    0.87 0.67 0.80 0.67 0.91 0.95 0.91 0.94 0.94 0.86;
    0.29 0.54 0.18 0.58 0.54 0.28 0.36 0.60 0.60 0.66;
    0.80 0.54 0.67 0.53 0.76 0.79 0.78 0.48 0.48 0.80;
    0.79 0.86 0.68 0.80 0.88 0.92 0.54 0.88 0.86 0.87];

% Calculas el tiedrank de la matriz. Si la medida es cuanto más alto mejor,
% como es el caso del Accuracy y la F-Measure, la matriz se pone negativa
% (si no, se deja positiva). En esta caso trasponemos la matriz para que
% calcule los rangos por vídeo (depende de como tengas dispuestos los datos
% en la matriz.
R = tiedrank(-A');

% Sumas las filas y ya tienes el rank sum para cada método
RankSums = sum(R,2);

fprintf('\nRANK SUMS\n');
for i=1:length(Modelos)
    fprintf('%s: %g\n',Modelos{i},RankSums(i));
end