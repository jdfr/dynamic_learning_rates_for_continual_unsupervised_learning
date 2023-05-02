% Demo of Latent Semantic Indexing, see http://www1.se.cuhk.edu.hk/~seem5680/lecture/LSI-Eg.pdf
% Also http://c2.com/cgi/wiki?LatentSemanticIndexing
clear all

% Loading NumWords x NumDocuments matrix
load('./Pruebas clustering/Data/ObamaTweets.mat','Samples');

% Singular value decomposition
[U,S,V] = svd(Samples);

% Remove residual imaginary parts
U=real(U);
S=real(S);
V=real(V);

% Nearly perfect recovery (no dimensions dropped)
norm(Samples-U*S*V')

% Dimensionality reduction (to NumConcepts dimensions)
NumConcepts=10;
ReducedU=U(:,1:NumConcepts);
ReducedS=S(1:NumConcepts,1:NumConcepts);
ReducedV=V(:,1:NumConcepts);
norm(Samples-ReducedU*ReducedS*ReducedV')

% Reduced dimensionality documents are in V', the importance of the concepts
% are in S, and the words that define each concept are in U. This means
% that S(i,i) is the importance of the i-th concept in the dataset, U(i,j)
% is the contribution of the i-th word to the j-th concept (positive values
% of U(i,j) mean that the i-th word appears with the j-th concept more than
% average, while negative values of U(i,j) mean that the i-th word appears
% with the j-th concept less than average), and V(i,j) is the contribution
% of the j-th concept to the i-th document.

% In order to reduce the dimensionality of a not previously seen document, you do the same transformation as in the
% training samples:
% ReducedSample=inv(ReducedS)*ReducedU'
ReducedSamples=ReducedV';
%ReducedSamples2=inv(ReducedS)*ReducedU'*Samples;
ReducedSamples2=(ReducedS\ReducedU')*Samples;
norm(ReducedSamples-ReducedSamples2)
TestSample=rand(size(Samples,1),1);
ReducedTestSample=(ReducedS\ReducedU')*TestSample;

% Then you can go back to the original high dimensional space of words.
% This is an approximate reconstruction, which gets closer to the original
% data as we increase the number of concepts.
ReconstructedSamples=ReducedU*ReducedS*ReducedSamples;
norm(Samples-ReconstructedSamples)
ReconstructedTestSample=ReducedU*ReducedS*ReducedTestSample;
norm(TestSample-ReconstructedTestSample)


