function [batches, classes_per_batch, batches_expanded] = GenerateSamplesCORABatches(paperWordsTransposed, paperClasses, NumEpochs, NumBatches);

NumSamples = numel(paperClasses);

BaseNumSamples = floor(NumSamples/NumBatches);
remainder = NumSamples-BaseNumSamples*NumBatches;
rndp = randperm(NumBatches);
selected = rndp(1:remainder);
NumSamplesPerBatch = zeros(NumBatches,1)+BaseNumSamples;
NumSamplesPerBatch(selected) = BaseNumSamples+1;

idxs = randperm(size(paperWordsTransposed,2));
paperWordsTransposed = paperWordsTransposed(:, idxs);
paperClasses = paperClasses(idxs);

batches = mat2cell(paperWordsTransposed, size(paperWordsTransposed, 1), NumSamplesPerBatch);
classes_per_batch = mat2cell(paperClasses, NumSamplesPerBatch, 1);
batches_expanded = cell(size(batches));
for k=1:numel(batches)
  batches_expanded{k} = repmat(batches{k}, 1, NumEpochs);
  batches_expanded{k} = batches_expanded{k}(:, randperm(size(batches_expanded{k},2)));
end


