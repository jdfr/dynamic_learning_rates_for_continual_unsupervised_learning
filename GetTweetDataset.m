function Samples = GetTweetDataset(Filename)

Samples = [];
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
    
    NdxSample = 0;     
    NumBrackets = 0;
    NumBracketsOld = 0;
    isID = 1;
    while 1,

        % Read a word of the line
        [token, tline] = strtok(tline,[' ' ',']);                                                  
        if isempty(token),
            break;
        end        
        
        % Data is stored as [id, freq] where 'id' will be the feature and
        % 'freq' the number of ocurrences of each distinct word. New vectors
        % and the entire matrix are between brackets.        
        while strcmp(token(1),'['),
            NumBracketsOld = NumBrackets;
            NumBrackets = NumBrackets + 1;
            if (NumBracketsOld == 1) && (NumBrackets == 2), % new sample
                NdxSample = NdxSample + 1;     
                fprintf('\tSample %d\n',NdxSample);
                if NdxSample > 1,
                    fprintf('\t\tFeatures=%d\n',NdxFeature);
                end
            end
            token = token(2:end);
        end
                
        while strcmp(token(end),']'),
            NumBracketsOld = NumBrackets;
            NumBrackets = NumBrackets - 1;
            if NumBrackets == 0, % end of file
                break;
            end
            token = token(1:end-1);
        end
        
        if isID,
            NdxFeature = str2double(token)+1; 
            fprintf('\t\tFeature %d\n',NdxFeature);
            isID = 0;
        else
            Samples(NdxFeature,NdxSample) = str2double(token);                
            isID = 1;                  
        end
    end                   
    fclose(fid);
end