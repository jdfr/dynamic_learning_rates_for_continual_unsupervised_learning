function PerformanceMeasures = DemoTestGNGSegmentation(Model,TrainingWinners,TrainingErrors,ColorSpaceName,VideoName,NumFrames,DeltaFrames,VideoPath,GTPath,Mu,Features)

MinArea = 10;
PerformanceMeasures = zeros(7,0);
Cputime = 0;
WindowSize = 1;
Beta = 1; % number of standard deviations considered foreground

switch ColorSpaceName,
    case {'RGB'}
        MyColorSpaceConversion = '';
    otherwise
        MyColorSpaceConversion = ['RGB->' ColorSpaceName];
end

for NdxFrame=1:NumFrames,
    
    if (mod(NdxFrame,100)==0),
        fprintf('\t\t\tFrame %d\n',NdxFrame);
    end
    
    % Compute Test Samples
    if exist(sprintf(GTPath,DeltaFrames+NdxFrame),'file')
        
        MyFrameOrig = imread(sprintf(VideoPath,DeltaFrames+NdxFrame)); 
        MyFrame = MyFrameOrig;
        if ~isempty(MyColorSpaceConversion),
            MyFrame = colorspace(MyColorSpaceConversion,MyFrame);
            MyFrame(isnan(MyFrame)) = 0;
        end
        
        % Build AuxFrames to compute the median of the frame window
        MyNumTestFrames = min(NumFrames - NdxFrame,(WindowSize-1)/2); % por si hay menos frames futuros que (WindowSize-1)/2
        MyNumTestFrames = MyNumTestFrames + (WindowSize-1)/2 + 1;
        AuxFrames = zeros(size(MyFrame,1),size(MyFrame,2),size(MyFrame,3),MyNumTestFrames);
        AuxFrames(:,:,:,(WindowSize-1)/2 + 1) = MyFrame;
        
        % Store the frame window
        ContWindow = 1;
        for NdxWindow=NdxFrame-(WindowSize-1)/2:NdxFrame+(WindowSize-1)/2,
                           
            if NdxWindow~=NdxFrame && exist(sprintf(GTPath,DeltaFrames+NdxWindow),'file'),
                MyFrame = imread(sprintf(VideoPath,DeltaFrames+NdxWindow));
                if ~isempty(MyColorSpaceConversion),
                    MyFrame = colorspace(MyColorSpaceConversion,MyFrame);
                    MyFrame(isnan(MyFrame)) = 0;
                end
                AuxFrames(:,:,:,ContWindow) = MyFrame;   
            end
            ContWindow = ContWindow+1;
        end
        
        MyFrame = median(AuxFrames,4);        
%         MyFrame = mean(AuxFrames,4); 
        MyFrame = NormalizarEspColor(double(MyFrame), ColorSpaceName);
        %------------------------------------------------------------------
        if Features>1,
            % NORMALIZED RGB
            SumFrame = sum(MyFrame,3);
            SumFrame(SumFrame==0) = 1;
            MyFrameNorm = MyFrame;
            MyFrameNorm(:,:,1) = squeeze(MyFrame(:,:,1))./SumFrame;
            MyFrameNorm(:,:,2) = squeeze(MyFrame(:,:,2))./SumFrame;
            MyFrameNorm(:,:,3) = squeeze(MyFrame(:,:,3))./SumFrame;
            if Features==3,
                % NORMALIZED RGB MEDIAN FILTERED
                MyFrameNorm(:,:,1) = medfilt2(MyFrameNorm(:,:,1),[5 5]);
                MyFrameNorm(:,:,2) = medfilt2(MyFrameNorm(:,:,2),[5 5]);
                MyFrameNorm(:,:,3) = medfilt2(MyFrameNorm(:,:,3),[5 5]);
            end
            MyFrame = MyFrameNorm;
        end
        %------------------------------------------------------------------
        MyMatrixFrame = reshape(shiftdim(MyFrame,2),3,[]);            
        [i,j] = ind2sub([size(MyFrame,1) size(MyFrame,2)],[1:length(MyMatrixFrame)]);
        TestSamples = [MyMatrixFrame;Mu*NormalizeData(i);Mu*NormalizeData(j)];

        % DETECT THE FOREGROUND
        tic;        
        [TestWinners,TestErrors] = TestGHNF(GetCentroidsGHNF(Model),TestSamples);
        
        % Criterion 1
    %     ForegroundPixels = ((TrainingWinners == TestWinners) == 0);
        % Criterion 2 (more restrictive)    
        ErrorDiff = abs(TrainingErrors - TestErrors);
        MyThreshold = nanmean(ErrorDiff)+Beta*nanstd(ErrorDiff);
    %     BackgroundPixels = (ErrorDiff <= mean(ErrorDiff));
    %     BackgroundPixels = (ErrorDiff > (mean(ErrorDiff)+std(ErrorDiff)));
        ForegroundPixels = (ErrorDiff > MyThreshold);

        % Foreground pixels are those that fulfill both criteria
    %     ForegroundPixels = (ForegroundPixels - BackgroundPixels) == 1;
        Cputime = Cputime + toc;
        imMask1 = reshape(ForegroundPixels,size(MyFrame,1),size(MyFrame,2));       
        imMask1 = double(imMask1 >= 0.5);

        % Fill holes (size 1) and remove objects with minimum area
        imMaskPost1 = bwmorph(imMask1,'majority');
%         imMask = bwmorph(imMask,'close',5);
%         imMask = bwmorph(imMask,'open',5);
        imMaskPost1 = removeSpuriousObjects(imMaskPost1, MinArea); 
        
        % Measure performance
        GroundTruthFrame = double(imread(sprintf(GTPath,DeltaFrames+NdxFrame)));        
        PerformanceMeasures(:,end+1) = EvaluatePerformance(imMaskPost1>40/255,GroundTruthFrame);
        
        % Plot frames
        fig = figure;
        set (fig, 'Units', 'normalized', 'Position', [0,0,1,1]);
        
        subplot(2,2,1)
        imshow(MyFrameOrig);
        title('Original Frame');
        
        subplot(2,2,3)
        imshow(GroundTruthFrame);
        title('Ground Truth');
        
        subplot(2,2,2)
        imshow(imMask1);
        title('Criterion 1');
        
        subplot(2,2,4)
        imshow(imMaskPost1);
        title('Postprocessed');
                
        fprintf('\nPress any key to continue\n');
        pause
        close
    end
end
PerformanceMeasures(1,end+1) = Cputime/NumFrames;