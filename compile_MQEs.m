function [MQEs1, MQEs2] = compile_MQEs;

n={'images/Anillo.bmp', 'images/Circle.bmp', 'images/CuadradoHueco.bmp', 'images/CuadradoNoHueco.bmp', 'images/Eme.bmp', 'images/Equis.bmp', 'images/Ese.bmp', 'images/Haltera.bmp', 'images/Irregular.bmp', 'images/Ka.bmp', 'images/Ocho.bmp', 'images/RelojArena.bmp'};

MQEs1 = cell(length(n)+1, 10);

MQEs1{1,1} = 'NOMBRE';
MQEs1{1,2} = 'BASELINE';
for k=1:length(n)
  %pth = ['images/competitive/baseline/' n{k}(8:end) '.baseline.mat'];
  pth = ['images/som/baseline/' n{k}(8:end) '.baseline.mat'];
  load(pth);
  MQEs1{k+1, 1} = n{k};
  MQEs1{k+1, 2} = mqe;
end

MQEs1{1,3} = 'MAX LR 0.4, LINEAR';
MQEs1{1,4} = 'MAX LR 0.4, QUADRATIC';
MQEs1{1,5} = 'MAX LR 0.4, CUBIC';
MQEs1{1,6} = 'MAX LR 0.4, EXPONENTIAL';
MQEs1{1,7} = 'MAX LR 1.0, LINEAR';
MQEs1{1,8} = 'MAX LR 1.0, QUADRATIC';
MQEs1{1,9} = 'MAX LR 1.0, CUBIC';
MQEs1{1,10} = 'MAX LR 1.0, EXPONENTIAL';
for k=1:length(n)
  %pth = ['images/competitive/bateria2_MAXLR_0.4/' n{k}(8:end) '.experiments.mat'];
  pth = ['images/som/bateria2_MAXLR_0.4_40000/' n{k}(8:end) 'experiments.mat'];
  load(pth);
  MQEs1{k+1, 3} = best_run.mqe(1);
  MQEs1{k+1, 4} = best_run.mqe(2);
  MQEs1{k+1, 5} = best_run.mqe(3);
  MQEs1{k+1, 6} = best_run.mqe(4);
  %pth = ['images/competitive/bateria2_MAXLR_1/' n{k}(8:end) '.experiments.mat'];
  pth = ['images/som/bateria2_MAXLR_1_40000/' n{k}(8:end) 'experiments.mat'];
  load(pth);
  MQEs1{k+1, 7} = best_run.mqe(1);
  MQEs1{k+1, 8} = best_run.mqe(2);
  MQEs1{k+1, 9} = best_run.mqe(3);
  MQEs1{k+1,10} = best_run.mqe(4);
end

names = {
{n{1}, n{2}},
{n{2}, n{1}},
{n{3}, n{4}},
{n{4}, n{3}},
{n{5}, n{6}},
{n{6}, n{5}},
{n{8}, n{9}},
{n{9}, n{8}},
{n{2}, n{6}},
{n{6}, n{2}},
{n{1}, n{6}},
{n{6}, n{1}},
{n{2}, n{7}},
{n{7}, n{2}},
{n{11}, n{12}},
{n{12}, n{11}},
};

MQEs2 = cell(length(names)+1, 9);
MQEs2{1, 1} = 'SECUENCIA';
MQEs2{1, 2} = 'MAX LR 0.4, LINEAR';
MQEs2{1, 3} = 'MAX LR 0.4, QUADRATIC';
MQEs2{1, 4} = 'MAX LR 0.4, CUBIC';
MQEs2{1, 5} = 'MAX LR 0.4, EXPONENTIAL';
MQEs2{1, 6} = 'MAX LR 1.0, LINEAR';
MQEs2{1, 7} = 'MAX LR 1.0, QUADRATIC';
MQEs2{1, 8} = 'MAX LR 1.0, CUBIC';
MQEs2{1, 9} = 'MAX LR 1.0, EXPONENTIAL';
for k=1:length(names)
  FusedFileNames = cell(length(names{k}), 1);
  for j=1:length(names{k})
    p1 = strfind(names{k}{j}, '/');
    p2 = strfind(names{k}{j}, '.');
    FusedFileNames{j} = [names{k}{j}(p1(end)+1:p2(end)-1) '.'];
  end
  FusedFileNames = [FusedFileNames{:}];
  MQEs2{k+1, 1} = FusedFileNames;
  %pth = ['images/competitive/bateria3_MAXLR_0.4_SAM_10000/' FusedFileNames 'experiments.mat'];
  pth = ['images/som/bateria3_MAXLR_0.4_SAM_20000/' FusedFileNames 'experiments.mat'];
  load(pth);
  MQEs2{k+1, 2} = best_run.mqe(1);
  MQEs2{k+1, 3} = best_run.mqe(2);
  MQEs2{k+1, 4} = best_run.mqe(3);
  MQEs2{k+1, 5} = best_run.mqe(4);
  %pth = ['images/competitive/bateria3_MAXLR_1_SAM_10000/' FusedFileNames 'experiments.mat'];
  pth = ['images/som/bateria3_MAXLR_1_SAM_20000/' FusedFileNames 'experiments.mat'];
  load(pth);
  MQEs2{k+1, 6} = best_run.mqe(1);
  MQEs2{k+1, 7} = best_run.mqe(2);
  MQEs2{k+1, 8} = best_run.mqe(3);
  MQEs2{k+1, 9} = best_run.mqe(4);
end

%save('images/competitive/resumen.mat', 'MQEs1', 'MQEs2');
save('images/som/resumen.mat', 'MQEs1', 'MQEs2');


