function [samples, originalfull] = GenerateSamplesCORA(cora, NumSamples);

if isstruct(cora)
  originalfull = double(cora.paperWords');
else
  originalfull = cora;
end

original = sparse(originalfull);

matrices = cell(0,1);

while NumSamples>0
  %p = randperm(size(original, 2));
  %if NumSamples<numel(p)
  %  p = p(1:NumSamples);
  %end
  %matrices{end+1,1} = original(:,p);
  if NumSamples<size(original,2)
    matrices{end+1,1} = original(:,1:NumSamples);
    NumSamples = NumSamples-NumSamples;
  else
    matrices{end+1,1} = original;
    NumSamples = NumSamples-size(original,2);
  end
end

samples = [matrices{:}];


