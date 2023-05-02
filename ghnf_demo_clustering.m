% GHNF demo for clustering of Obama's tweets

clear all
MaxNeurons = 20; % Maximum number of neurons in each graph
Tau = 0.01;

% The following values of the parameters are those considered in the
% original GNG paper by Fritzke (1995)
Lambda = 100;
Epochs = 2;
EpsilonB = 0.2;
EpsilonN = 0.006;
Alpha = 0.5;
AMax = 50;
D = 0.995;

% Load data
load('ObamaTweets.mat','Samples');
load('Dictionary.mat','Dictionary');
Samples = NormalizeData(Samples);

% GHNF Training
[Model] = TrainGHNF(Samples,Epochs,MaxNeurons,Tau,Lambda,EpsilonB,EpsilonN,Alpha,AMax,D,1);

% Plot most frequent words
Handle = PlotGHNFWords(Model,Dictionary);


