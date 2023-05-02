function x = insideSphereGen(Dimension,NumSamples)
% Sample NumSamples points at uniform inside the D-dimensional unit hypersphere.

% Gaussian data
x = randn(NumSamples,Dimension);

% Normalization
x = (x.*repmat(rand(NumSamples,1).^(1/Dimension)./sqrt(sum(x.^2,2)),1,Dimension))';  
