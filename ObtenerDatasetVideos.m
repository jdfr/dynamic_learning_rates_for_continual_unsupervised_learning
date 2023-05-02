% Obtiene un dataset a partir de N frames para cada vídeo y espacio de color
clear all;  

NumBackgroundFrames = 100;
ColorSpaces = {'HSL','HSV','Lab','Luv','RGB','YCbCr'};
% ColorSpaces = {'RGB'};
NumColorSpaces = length(ColorSpaces);
Mu = 1; % penalization parameter for some dataset features
Features = 3; % 1=RGB, 2=Normalized RGB, and 3=Normalized RGB median filtered
WindowSizeRGBMedian = 5; % window size for the normalized RGB median filtered

%% VIDEOS TO BE PROCESSED
Videos = {'Campus','Curtain','Escalator','Fountain','LevelCrossing','OneShopOneWait1cor',...
    'Video2','Video4','WaterSurface','Lobby','LightSwitch'};
NumVideos = length(Videos);
PathVideos = '../../Videos/';
SubPathFrames = '/frames/f%07d.bmp';
SubPathFramesLightSwitch = '/frames/f%07d.jpg';
SubPathGT = '/GT/GT%07d.bmp';
DeltaFrames = [999,20999,1397,999,0,0,0,0,999,999,0];
NumFrames = [1438,2963,3417,523,500,1376,749,819,632,1545,2799];
GT = [20,20,20,20,163,11,444,316,20,20,21];
PrimerFrameGT = [1199,22755,1398,1152,275,0,253,273,1480,1152,666];

Data = [];
for NdxColorSpace=1:NumColorSpaces,
       
    MyColorSpace = ColorSpaces{NdxColorSpace};
    fprintf('\nESPACIO COLOR: %s\n',MyColorSpace);
    switch MyColorSpace,
        case {'RGB'}
            MyColorSpaceConversion = '';
        otherwise
            MyColorSpaceConversion = ['RGB->' MyColorSpace];
    end 
    
    for NdxVideo=1:NumVideos,        

        MyVideo = Videos{NdxVideo};
        fprintf('\n\tVIDEO: %s\n',MyVideo);
        PathMyVideo = [PathVideos MyVideo SubPathFrames];
        if NdxVideo==11, % jpg instead of bmp extension
            PathMyVideo = [PathVideos MyVideo SubPathFramesLightSwitch];
        end
        
        % Store the first frames belonging to the background
        NdxFrameInicio = PrimerFrameGT(NdxVideo) - NumBackgroundFrames;
        if NdxFrameInicio <= DeltaFrames(NdxVideo),
            NdxFrameInicio = DeltaFrames(NdxVideo)+1;
        end
        fprintf('\t\tFrames desde %d hasta %d (Primer Frame=%d)\n',NdxFrameInicio,NdxFrameInicio+NumBackgroundFrames-1,DeltaFrames(NdxVideo)+1);
%         MyFrame = imread(sprintf(PathMyVideo,DeltaFrames(NdxVideo)+1));        
        MyFrame = imread(sprintf(PathMyVideo,NdxFrameInicio));
        if ~isempty(MyColorSpaceConversion),
            MyFrame = colorspace(MyColorSpaceConversion,MyFrame); 
            MyFrame(isnan(MyFrame)) = 0;
        end
        FirstFrames = zeros(size(MyFrame,1),size(MyFrame,2),size(MyFrame,3),NumBackgroundFrames);
        FirstFrames(:,:,:,1) = MyFrame;             
        for NdxFrame=2:NumBackgroundFrames,
                        
            NdxFrameInicio = NdxFrameInicio + 1;
%             MyFrame = imread(sprintf(PathMyVideo,DeltaFrames(NdxVideo)+NdxFrame));
            MyFrame = imread(sprintf(PathMyVideo,NdxFrameInicio));
            if ~isempty(MyColorSpaceConversion),
                MyFrame = colorspace(MyColorSpaceConversion,MyFrame);
                MyFrame(isnan(MyFrame)) = 0;
            end
            FirstFrames(:,:,:,NdxFrame) = MyFrame;
        end
            
        % Obtain the dataset as the median of the first frames
%         MyFrame = mean(FirstFrames,4);
        MyFrame = median(FirstFrames,4);
        MyFrame = NormalizarEspColor(double(MyFrame), MyColorSpace);
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
                fprintf('\t\tFeatures: Normalized RGB Median Filtered\n');               
                MyFrameNorm(:,:,1) = medfilt2(MyFrameNorm(:,:,1),[WindowSizeRGBMedian WindowSizeRGBMedian]);
                MyFrameNorm(:,:,2) = medfilt2(MyFrameNorm(:,:,2),[WindowSizeRGBMedian WindowSizeRGBMedian]);
                MyFrameNorm(:,:,3) = medfilt2(MyFrameNorm(:,:,3),[WindowSizeRGBMedian WindowSizeRGBMedian]);
            else
                fprintf('\t\tFeatures: Normalized RGB\n');
            end
            MyFrame = MyFrameNorm;
        else
            fprintf('\t\tFeatures: Raw RGB\n');
        end
        %------------------------------------------------------------------
        MyMatrixFrame = reshape(shiftdim(MyFrame,2),3,[]);            
        [i,j] = ind2sub([size(MyFrame,1) size(MyFrame,2)],[1:length(MyMatrixFrame)]);
        Data = [MyMatrixFrame;Mu*NormalizeData(i);Mu*NormalizeData(j)];           
        
%         % EZEQUIEL: HERE WE DIVIDE RED, GREEN AND BLUE BY THE SUM OF RED, GREEN AND
%         % BLUE
%         Data(:,1:3) = Data(:,1:3)./repmat(sum(Data(:,1:3),1),[size(Data,1) 1]);
%         % % HERE WE AVOID ERRORS CAUSED BY NaNs COMING FROM DIVIDING BY ZERO
%         Data(~isfinite(Data)) = 0;
%         Data(:,1:3) = double(medfilt2(uint16(32768*Data),[5 5]))/32768;

        save(['Dataset_' MyColorSpace '_' MyVideo '.mat'],'Data'); 
        Data = [];
    end      
end
        