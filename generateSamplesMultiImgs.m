function [Samples, NumSteps, FusedFileNames] = generateSamplesMultiImgs(prefixpath, FileNames, NumSamplesPerImage)

FusedFileNames = cell(length(FileNames), 1);
Samples = cell(length(FileNames),1);
for k=1:length(FileNames)
  Samples{k} = GenerateSamplesImg(FileNames{k},NumSamplesPerImage);
  p1 = strfind(FileNames{k}, '/');
  p2 = strfind(FileNames{k}, '.');
  FusedFileNames{k} = [FileNames{k}(p1(end)+1:p2(end)-1) '.'];
end
NumSteps = length(FileNames)*NumSamplesPerImage;
if numel(prefixpath)==0
  FusedFileNames = [FileNames{1}(1:p1(end)) FusedFileNames{:}];
else
  FusedFileNames = [prefixpath FusedFileNames{:}];
end

