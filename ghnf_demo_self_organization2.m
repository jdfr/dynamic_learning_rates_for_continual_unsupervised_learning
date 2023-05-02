% GHNF demo for self-organization with 'X' letter shape

clear all
NumSamples=10000;
MaxNeurons = 20; % Maximum number of neurons in each graph
Tau = 0.1;

% The following values of the parameters are those considered in the
% original GNG paper by Fritzke (1995)
Lambda=100;
Epochs=2;
EpsilonB=0.2;
EpsilonN=0.006;
Alpha=0.5;
AMax=50;
D=0.995;

% Generate data ('X' letter shape)
Samples = GenerateSamplesImg('X_letter.bmp',NumSamples);

% GHNF Training
[Model] = TrainGHNF(Samples,Epochs,MaxNeurons,Tau,Lambda,EpsilonB,EpsilonN,Alpha,AMax,D,1);

% Plot the Model
MyAxis = [0 1 0 1];
Handle = PlotGHNF(Model,MyAxis);
h = plot(Samples(1,:),Samples(2,:),'*y');
uistack(h,'bottom');
axis([0 1 0 1])
hold off

