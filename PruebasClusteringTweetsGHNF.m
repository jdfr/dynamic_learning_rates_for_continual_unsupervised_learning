% Pruebas Clustering de Tweets para el GHNF
clear all

% Parameters
Tau = 0.01;
NumRepetitions = 10;
NumEpochs = 2;
MaxNeurons = 20; % Maximum number of neurons in each graph
Politicians = {'Obama','Clinton','Sanders','Trump'};
MeasuresName = {'MSE','TE','DBI','Dunn','CS','CPUTime'};
HigherIsBetter = [0 0 0 1 1 0];
NumMeasures = length(MeasuresName);
PathData = './Data/';
PathResults = './Results/v31/';

% The following values of the parameters are those considered in the
% original GNG paper by Fritzke (1995)
Lambda = 100;
EpsilonB = 0.2;
EpsilonN = 0.006;
Alpha = 0.5;
AMax = 50;
D = 0.995;

% Load all data and word dictionary
load([PathData 'Dictionary.mat'],'Dictionary');
load([PathData 'AllTweets.mat'],'Samples');
Labels = Samples(end,:);
Samples = Samples(1:end-1,:);
[NumTerms,NumDocs] = size(Samples);

% Perform a tf-idf weighting
Samples = tfidf(Samples);

% Perform a LSI
[U, SamplesLSI] = tcslsi(Samples);

for NdxPolitician=1:length(Politicians),    
    % Load data
    MyPolitician = Politicians{NdxPolitician};
    fprintf('\nTweets de %s\n',MyPolitician);
    TrainingSamples = SamplesLSI(:,Labels==NdxPolitician);
    NomFich = ['GHNF_' MyPolitician];       

    if exist([PathResults NomFich 'Results.mat'],'file'),
        fprintf('\tResults already obtained for %s\n',MyPolitician);
        load([PathResults NomFich 'Results.mat'],'Model','Measures');
    else
        Models = cell(1,length(NumRepetitions));
        Measures = zeros(NumMeasures,NumRepetitions);    
        for NdxRepetition=1:NumRepetitions,
            % GHNF Training
            fprintf('\tEntrenamiento GHNF %d\n',NdxRepetition);
            tic;
            Models{NdxRepetition} = TrainGHNF(TrainingSamples,NumEpochs,MaxNeurons,Tau,Lambda,EpsilonB,EpsilonN,Alpha,AMax,D,1);
            Measures(6,NdxRepetition) = toc;                        
            while isempty(Models{NdxRepetition}),
                % GHNF Training
                fprintf('\tRepitiendo Entrenamiento GHNF %d\n',NdxRepetition);
                tic;
                Models{NdxRepetition} = TrainGHNF(TrainingSamples,NumEpochs,MaxNeurons,Tau,Lambda,EpsilonB,EpsilonN,Alpha,AMax,D,1);
                Measures(6,NdxRepetition) = toc;                
            end     
            
            % Compute Performance Measures
            Model = Models{NdxRepetition};
            Prototypes = GetCentroidsGHNG(Model);
            [Measures(3,NdxRepetition),Measures(1,NdxRepetition),Measures(4,NdxRepetition),Measures(5,NdxRepetition)] = EvaluateClustering(Prototypes,TrainingSamples);
            Measures(2,NdxRepetition) = TopographicErrorGHNG(Model,TrainingSamples);
            
        end

        % Compute the Best Trained Model
        Ranks = cell(1,NumMeasures);
        for NdxMeasure=1:NumMeasures,
            if HigherIsBetter(NdxMeasure),
                Ranks{NdxMeasure} = tiedrank(-Measures(NdxMeasure,:));
            else
                Ranks{NdxMeasure} = tiedrank(Measures(NdxMeasure,:));
            end
        end   
        [Minimum,BestCandidate] = min(Ranks{1}+Ranks{2}+Ranks{3}+Ranks{4}+Ranks{5}+Ranks{6});
        Model = Models{BestCandidate};
        save([PathResults NomFich 'Results.mat'],'Model','Measures');
    end
    
    % Plot most frequent concepts
    Handle = PlotGHNFConcepts(Model,U,Dictionary,[PathResults NomFich '_Concepts_']);       
    close all;
    
    % Plot most frequent words
    Model.Means = U*Model.Means;
    Handle = PlotGHNFWords(Model,Dictionary,[PathResults NomFich '_Words_']);        
    close all;    
end


