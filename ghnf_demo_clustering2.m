% GHNF demo for clustering of Obama's tweets (WITH PCA)

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

% Perform a global PCA
Dimension = 15;
GlobalMean = mean(Samples')';
GlobalCov = cov(Samples');
[Uq, Lambdaq] = eigs(GlobalCov,Dimension,'LM');
Lambdaq = diag(Lambdaq);
Samples_zq = Uq'*(Samples-repmat(GlobalMean,1,size(Samples,2)));
Samples = Samples_zq;

% GHNF Training
[Model] = TrainGHNF(Samples,Epochs,MaxNeurons,Tau,Lambda,EpsilonB,EpsilonN,Alpha,AMax,D,1);

% Undo PCA
Model.Means = Uq*Model.Means+repmat(GlobalMean,1,size(Model.Means,2));

% Plot most frequent words
Handle = PlotGHNFWords(Model,Dictionary);


