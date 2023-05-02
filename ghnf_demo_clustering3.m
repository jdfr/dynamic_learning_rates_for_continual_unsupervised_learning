% GHNF demo for clustering of Obama's tweets (WITH LSI)

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
% Samples = NormalizeData(Samples);

% Perform a LSI
% GlobalMean = mean(Samples')';
% [U, SamplesLSI] = tcslsi(Samples-repmat(GlobalMean,1,size(Samples,2)));
[U, SamplesLSI] = tcslsi(Samples);

% GHNF Training
[Model] = TrainGHNF(SamplesLSI,Epochs,MaxNeurons,Tau,Lambda,EpsilonB,EpsilonN,Alpha,AMax,D,1);

% Plot most frequent concepts
PlotGHNFConcepts(Model,U,Dictionary);

% Undo LSI
% Model.Means = U*Model.Means+repmat(GlobalMean,1,size(Model.Means,2));
Model.Means = U*Model.Means;

% Plot most frequent words
PlotGHNFWords(Model,Dictionary);


