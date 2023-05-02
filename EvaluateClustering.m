function [DBI,MSE,Dunn,MeanSilh]=EvaluateClustering(Prototypes,Samples)
% Evaluate a clustering with unsupervised measures: 
%   The standard Davies-Boulding index (lower is better)
%   The Mean Squared Error (lower is better)
%   A Dunn index (higher is better). We consider the maximum of the means
%   of the pairwise intracluster distances in the denominator, and the minimum
%   of the distances among cluster means in the numerator. This version of the 
%   Dunn index is advocated in:
%   Bezdek and Pal (1995). Cluster validation with generalized Dunn's
%   indices.
%   The choice of the numerator is
%   recommended because it is robust and computationally efficient, and the
%   choice of the denominator is recommended because it is stable in
%   datasets with outliers.
%   The mean silhouette coefficient (higher is better), in the interval [-1,1]

% Compute distances from samples to prototypes and remove empty clusters
Distances=distanceMEX(Samples,Prototypes);
[MinDistances,InitialAssignments]=min(Distances,[],2);
[NonEmptyClusters,Trash,Assignments]=unique(InitialAssignments);
Prototypes=Prototypes(:,NonEmptyClusters);
[Dimension,NumPrototypes]=size(Prototypes);

MSE=mean(MinDistances.^2);
MeanSilh=mean(silhouetteMEX(Samples,Assignments-1,NumPrototypes));

StdDev=zeros(1,NumPrototypes);
MeanPairwiseDistances=zeros(1,NumPrototypes);
ClusterMeans=zeros(Dimension,NumPrototypes);
for NdxProto=1:NumPrototypes
    MyCluster=Samples(:,Assignments==NdxProto);
    ClusterSize=size(MyCluster,2);
    Diffs=MyCluster-repmat(Prototypes(:,NdxProto),[1 size(MyCluster,2)]);
    StdDev(NdxProto)=sqrt(mean(Diffs(:).^2));
    PairwiseDistances=distanceMEX(MyCluster,MyCluster);
    MeanPairwiseDistances(NdxProto)=sum(PairwiseDistances(:))/...
        (ClusterSize*(ClusterSize-1));
    ClusterMeans(:,NdxProto)=mean(MyCluster,2);
end
MaxMeanPairwise=max(MeanPairwiseDistances);

DBI=0;
MinInterDist=1.0e200;
for NdxProto1=1:NumPrototypes    
    % ClusterDiv(NdxProto2) is R_{NdxProto1,NdxProto2}
    ClusterDiv=zeros(1,NumPrototypes);
    for NdxProto2=1:NumPrototypes
        if NdxProto2~=NdxProto1            
            ClusterDiv(NdxProto2)=(StdDev(NdxProto1)+StdDev(NdxProto2))/...
                norm(Prototypes(:,NdxProto1)-Prototypes(:,NdxProto2));
            InterMeanDist=norm(ClusterMeans(:,NdxProto1)-ClusterMeans(:,NdxProto2));
            if InterMeanDist<MinInterDist
                MinInterDist=InterMeanDist;
            end
        end
    end
    % Maximum is D_NdxProto1
    [Maximum Index]=max(ClusterDiv);
    DBI=DBI+Maximum;
end
DBI=DBI/NumPrototypes;

Dunn=MinInterDist/MaxMeanPairwise;
