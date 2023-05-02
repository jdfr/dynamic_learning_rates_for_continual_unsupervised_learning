function D = hopDistances(connections)

% Adjacency matrix

D = nan(size(connections));
B = connections;
k = 1;
while any(isnan(D(:)))

    % Check for new walks, and assign distance
    D(B>0 & isnan(D)) = k;

    % Iteration
    k = k+1;
    B = B*connections;
end

for k=1:size(D,1)
  D(k,k) = 0;
end
