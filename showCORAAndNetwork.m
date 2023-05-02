function showCORAAndNetwork(modes, varargin)

if strcmp(modes,'pca')
  pcaCORA(varargin{:});
elseif strcmp(modes,'cmds')
  cmdsCORA(varargin{:});
elseif strcmp(modes,'ncmds')
  ncmdsCORA(varargin{:});
end

function ncmdsCORA(cora, centroids, connections)

N=1;

centroids = centroids(:,:)';

allPoints = [centroids; cora.paperWords];

dsts = pdist(allPoints, 'euclidean');

% 'metricstress' fails with N=2,3
%[Y,stress] = mdscale(dsts,2,'criterion','metricstress');

% 'metricsstress' is really weird with N=2,3
%[Y,stress] = mdscale(dsts,2,'criterion','metricsstress');

%sammon
%[Y,stress] = mdscale(dsts,2,'criterion','sammon');

%strain
%[Y,stress] = mdscale(dsts,2,'criterion','strain');

%nonmetric stress
%[Y,stress] = mdscale(dsts,2,'criterion','stress');

%nonmetric sstress
[Y,stress] = mdscale(dsts,2,'criterion','sstress');

transformed_centroids = Y(1:size(centroids,1),:);
transformed_samples   = Y(size(centroids,1)+1:end,:);

figure;

plot(transformed_samples(1:N:end,1),transformed_samples(1:N:end,2),'.b');
%plot3(transformed_samples(1:N:end,1),transformed_samples(1:N:end,2),transformed_samples(1:N:end,3),'.b');

hold on;

plot(transformed_centroids(:,1),transformed_centroids(:,2),'.r');
%plot3(transformed_centroids(:,1),transformed_centroids(:,2),transformed_centroids(:,3),'.r');

[a, b] = find(connections);
for x=1:numel(a)
  for y=1:numel(b)
    line([transformed_centroids(a(x),1) transformed_centroids(b(y),1)],...
         [transformed_centroids(a(x),2) transformed_centroids(b(y),2)],...
    ...%     [transformed_centroids(a(x),3) transformed_centroids(b(y),3)],...
         'Color', 'k');
  end
end

axis equal;


function cmdsCORA(cora, centroids, connections)

N=1;
centroids = centroids(:,:)';

allPoints = [centroids; cora.paperWords];

dsts = pdist(allPoints, 'euclidean');

[Y,eigvals] = cmdscale(dsts);

transformed_centroids = Y(1:size(centroids,1),:);
transformed_samples   = Y(size(centroids,1)+1:end,:);

figure;

plot(transformed_samples(1:N:end,1),transformed_samples(1:N:end,2),'.b');
%plot3(transformed_samples(1:N:end,1),transformed_samples(1:N:end,2),transformed_samples(1:N:end,3),'.b');

hold on;

plot(transformed_centroids(:,1),transformed_centroids(:,2),'.r');
%plot3(transformed_centroids(:,1),transformed_centroids(:,2),transformed_centroids(:,3),'.r');

[a, b] = find(connections);
for x=1:numel(a)
  for y=1:numel(b)
    line([transformed_centroids(a(x),1) transformed_centroids(b(y),1)],...
         [transformed_centroids(a(x),2) transformed_centroids(b(y),2)],...
    ...%     [transformed_centroids(a(x),3) transformed_centroids(b(y),3)],...
         'Color', 'k');
  end
end

axis equal;


function pcaCORA(cora, centroids, connections)

N=1;
[coeffs,score,~,~,explanation,mu] = pca(double(cora.paperWords));

figure;

plot(score(1:N:end,1),score(1:N:end,2),'.b');
%plot3(score(1:N:end,1),score(1:N:end,2),score(1:N:end,3),'.b');
hold on;

centroids = centroids(:,:)';
centroids = ((centroids-mu)*coeffs)';
plot(centroids(:,1),centroids(:,2),'.r');
%plot3(centroids(:,1),centroids(:,2),centroids(:,3),'.r');
[a, b] = find(connections);
for x=1:numel(a)
  for y=1:numel(b)
    line([centroids(a(x),1) centroids(b(y),1)],...
         [centroids(a(x),2) centroids(b(y),2)],...
    ...%     [centroids(a(x),3) centroids(b(y),3)],...
         'Color', 'k');
  end
end

axis equal;

