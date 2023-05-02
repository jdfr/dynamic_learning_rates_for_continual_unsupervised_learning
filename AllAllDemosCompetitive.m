function AllAllDemosCompetitive(sws, mqeModes);

%sws = {true, false};
%mqeModes = 0:7;

for i=1:numel(sws)
  for j=1:numel(mqeModes)
    AllDemosCompetitive('_exploration', '1', 'online', 'competitive', mqeModes(j), sws{i}, 1, 40000);
  end
end

