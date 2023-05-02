function Samples = GetTweetDataset2(Filename,NumSamples,NumFeatures)

Samples = zeros(NumFeatures,NumSamples);
fid = fopen(Filename);
if fid == -1,
    fprintf('Error opening document %s',Filename);
else                    
    
    fprintf('\nExtracting tweets datasets...\n');
    % Read the file in a line
    tline = fgetl(fid);    
    if ~ischar(tline),   
        return;
    end    
          
    for NdxSample=1:NumSamples,
        fprintf('\tSample %d\n',NdxSample);
        for NdxFeature=1:NumFeatures,
            fprintf('\t\tFeature %d\n',NdxFeature);
            
            % Data is stored as [id, freq] where 'id' will be the feature and
            % 'freq' the number of ocurrences of each distinct word. New vectors
            % and the entire matrix are between brackets.  
            [token, tline] = strtok(tline,[' ' ',' '[' ']']);                                                  
            [token, tline] = strtok(tline,[' ' ',' '[' ']']);                                                          
            Samples(NdxFeature,NdxSample) = str2double(token);                
        end
    end               
    fclose(fid);
end