function cora = slurpCORA;

PAPER_TO_CLASS  = 'cora/paper_to_class.txt';
PAPER_TO_WORD   = 'cora/paper_to_word.txt';
CITED_TO_CITING = 'cora/paper_cited_to_citing.txt';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = fopen(PAPER_TO_CLASS, 'r'); 

paperClass = zeros(0,2, 'int32');
classIds = cell(0,1);
while ~feof(f)
  paper = fscanf(f, '%d', 1);
  if feof(f)
    break;
  end
  cstr  = fscanf(f, '%s', 1);
  idx = find(strcmp(classIds, cstr), 1);
  if numel(idx)==0
    classIds{end+1,1} = cstr;
    idx = numel(classIds);
  end
  paperClass(end+1,:) = [paper idx];
end

fclose(f);

paperIds = paperClass(:,1);
paperClass = int8(paperClass(:,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = fopen(PAPER_TO_WORD, 'r'); 

paperWordsExtended = fscanf(f, '%d word%d', [2 Inf])';

fclose(f);

numWords = max(paperWordsExtended(:,2));

paperWords = false(numel(paperIds),numWords);
for k=1:size(paperWordsExtended,1)
  paper = paperWordsExtended(k,1);
  idx   = find(paperIds==paper,1);
  %not actually necessary
  if numel(idx)==0
    paperIds(end+1) = paper;
    idx = numel(paperIds);
  end
  word  = paperWordsExtended(k,2);
  paperWords(idx, word) = true;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = fopen(CITED_TO_CITING, 'r'); 

citesExtended = fscanf(f, '%d %d', [2 Inf])';

fclose(f);

%not actually necessary
if true
  for k=1:size(citesExtended,1)
    paper1 = citesExtended(k,1);
    idx1   = find(paperIds==paper1,1);
    if numel(idx1)==0
      paperIds(end+1) = paper1;
      idx1 = numel(paperIds);
    end
    paper2 = citesExtended(k,2);
    idx2   = find(paperIds==paper2,1);
    if numel(idx2)==0
      paperIds(end+1) = paper2;
      idx2 = numel(paperIds);
    end
  end
end

%cites(i,j)==true IFF paper i cites paper j
cites = false(numel(paperIds), numel(paperIds));

for k=1:size(citesExtended,1)
  paperCited = citesExtended(k,1);
  idxCited   = find(paperIds==paperCited,1);
  paperCiting = citesExtended(k,2);
  idxCiting   = find(paperIds==paperCiting,1);
  cites(idxCiting, idxCited) = true;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cora = struct('paperIds', paperIds, 'classIds', {classIds}, 'paperClass', paperClass, 'paperWords', paperWords, 'cites', cites);

