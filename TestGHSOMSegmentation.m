function PerformanceMeasures = TestGHSOMSegmentation(Model,TrainingWinners,TrainingErrors,ColorSpaceName,VideoName,NumFrames,NumGTFrames,DeltaFrames,VideoPath,GTPath,Mu,Features,WindowSizeRGBMedian,MaxNeurons,Beta)

MinArea = 10;
PerformanceMeasures = zeros(7,0);
Cputime = 0;
contGT = 0;
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
        
        MyFrame = imread(sprintf(VideoPath,DeltaFrames+NdxFrame));
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
                MyFrameNorm(:,:,1) = medfilt2(MyFrameNorm(:,:,1),[WindowSizeRGBMedian WindowSizeRGBMedian]);
                MyFrameNorm(:,:,2) = medfilt2(MyFrameNorm(:,:,2),[WindowSizeRGBMedian WindowSizeRGBMedian]);
                MyFrameNorm(:,:,3) = medfilt2(MyFrameNorm(:,:,3),[WindowSizeRGBMedian WindowSizeRGBMedian]);
            end
            MyFrame = MyFrameNorm;
        end
        %------------------------------------------------------------------
        MyMatrixFrame = reshape(shiftdim(MyFrame,2),3,[]);            
        [i,j] = ind2sub([size(MyFrame,1) size(MyFrame,2)],[1:length(MyMatrixFrame)]);
        TestSamples = [MyMatrixFrame;Mu*NormalizeData(i);Mu*NormalizeData(j)];

        % DETECT THE FOREGROUND
        tic;
        % Test GHNF
        [TestWinners,TestErrors] = TestGHNF(GetCentroidsGHSOM(Model),TestSamples);
        
        % Criterion 1
    %     ForegroundPixels = ((TrainingWinners == TestWinners) == 0);
        % Criterion 2 (more restrictive)
%         ErrorDiff = abs(TrainingErrors - TestErrors);
%         MyThreshold = nanmean(ErrorDiff)+Beta*nanstd(ErrorDiff);
%         ForegroundPixels = (ErrorDiff > MyThreshold);
        
%         % EZEQUIEL'S THRESHOLD FOR FILTERING TRAINING ERRORS
%         % --------------------------------------------
%          TrainingErrorsFiltered=filter2(fspecial('gaussian'),...
%             reshape(TrainingErrors,size(MyFrame,1),size(MyFrame,2)));        
%         MyThreshold = TrainingErrorsFiltered(:)+Beta*nanstd(TrainingErrors);
%         ForegroundPixels = (TestErrors > MyThreshold);
        
        % EZEQUIEL'S THRESHOLD 2 FOR FILTERING TRAINING ERRORS
        % --------------------------------------------
        TrainingErrorsImg=reshape(TrainingErrors,size(MyFrame,1),size(MyFrame,2));
        TrainingErrorsMedian=colfilt(TrainingErrorsImg,[5 5],'sliding',@nanmedian);
        TrainingErrorsStd=colfilt(TrainingErrorsImg,[5 5],'sliding',@nanstd);
        
        MyThreshold = TrainingErrorsMedian(:)+Beta*TrainingErrorsStd(:);
        
        TestErrorsImg=reshape(TestErrors,size(MyFrame,1),size(MyFrame,2));
        TestErrorsMedian=colfilt(TestErrorsImg,[5 5],'sliding',@nanmedian);
        ForegroundPixels = (TestErrorsMedian(:) > MyThreshold);

%         % NORMAL THRESHOLD
%         %-------------------------------------------------------
%         ErrorDiff = abs(TrainingErrors - TestErrors);
%         MyThreshold = nanmean(ErrorDiff)+Beta*nanstd(ErrorDiff);
%         ForegroundPixels = (ErrorDiff > MyThreshold);
       
        % Foreground pixels are those that fulfill both criteria
    %     ForegroundPixels = (ForegroundPixels - BackgroundPixels) == 1;
        Cputime = Cputime + toc;
        imMask = reshape(ForegroundPixels,size(MyFrame,1),size(MyFrame,2));       
        imMask = double(imMask >= 0.5);

        % Fill holes (size 1) and remove objects with minimum area
        imMask = bwmorph(imMask,'majority');
%         imMask = bwmorph(imMask,'close',5);
%         imMask = bwmorph(imMask,'open',5);
        imMask = removeSpuriousObjects(imMask, MinArea);  
    
        % Measure performance
        contGT = contGT+1;
        if contGT == round(NumGTFrames/2),
            save([ColorSpaceName '_' VideoName '_' num2str(DeltaFrames+NdxFrame) '_' num2str(Beta) '_Beta.mat'],'imMask');
        end
        GroundTruthFrame = double(imread(sprintf(GTPath,DeltaFrames+NdxFrame)));
        %         PerformanceMeasures(:,end+1)=EvaluatePerformance(imMask>0.5,GroundTruthFrame);
        PerformanceMeasures(:,end+1) = EvaluatePerformance(imMask>40/255,GroundTruthFrame);                       
    end
end
PerformanceMeasures(1,end+1) = Cputime/NumFrames;