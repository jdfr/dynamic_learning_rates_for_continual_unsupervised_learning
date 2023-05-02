function [Handle] = PlotGHNFWords(Model,Dictionary,FileName)
% Plot a GHNF where each prototype can be displayed as the most frequent 
% word from a dictionary that represents the neuron.
% Inputs:
%   Model = the trained GHNF model
%   Dictionary = dictionary of words present the documents
%   FileName = name of the file to save the image

Handle = PlotGHNFWordsRec(Model,Dictionary,1,1,FileName);

function [Handle] = PlotGHNFWordsRec(Model,Dictionary,Level,Unit,FileName)
% Plot a GHNF layer where each prototype can be displayed as the most frequent 
% word from a dictionary that represents the neuron.
% Inputs:
%   Model = the trained GHNF model
%   Dictionary = dictionary of words present the documents
%   Level = level of the hierarchy
%   Unit = parent unit index
%   FileName = name of the file to save the image

Handle = [];
if ~isempty(Model) && (Level <= 3),
    
    Handle = figure('Name',[FileName num2str(Level) '_' num2str(Unit)],'units','normalized','outerposition',[0 0 1 1]);    
    hold on

    Colors = distinguishable_colors(4);
    MyColor = Colors(Level,:);
    
    SpanningTree=biograph(Model.SpanningTree);
    dolayout(SpanningTree);

    for NdxUnit=1:Model.MaxUnits
        if isfinite(Model.Means(1,NdxUnit))      
            
            MyPos=SpanningTree.Nodes(NdxUnit).Position;
            
            % Draw edges
            NdxNeighbors=find(Model.SpanningTree(NdxUnit,:));                        
            for NdxMyNeigh=1:numel(NdxNeighbors)
                line([MyPos(1) SpanningTree.Nodes(NdxNeighbors(NdxMyNeigh)).Position(1)],...
                    [MyPos(2) SpanningTree.Nodes(NdxNeighbors(NdxMyNeigh)).Position(2)],'Color',MyColor,'LineWidth',0.5);
            end
            
            % Draw nodes                        
            WordOK = 0;
            while ~WordOK,
                [MyWordFreq,MyWordNdx] = max(Model.Means(:,NdxUnit));
                Model.Means(MyWordNdx,NdxUnit) = -Inf;
                MyWord = Dictionary{MyWordNdx};
                if isnan(str2double(MyWord)),                        
                    WordOK = 1;
                end
            end
            plot(MyPos(1),MyPos(2),'o','LineWidth',1,'MarkerFaceColor',MyColor,'MarkerSize',max(15*MyWordFreq,8));
            text(MyPos(1),MyPos(2),MyWord,'Color',[0 0 0],'FontSize',max(8*MyWordFreq,8),'HorizontalAlignment','center');            
        end
    end
    
    hold off
    axis off
    Figure2pdf([FileName num2str(Level) '_' num2str(Unit)],12,10);
    
    % Draw children in different figures
    for NdxUnit=1:length(Model.Child)
        PlotGHNFWordsRec(Model.Child{NdxUnit},Dictionary,Level+1,NdxUnit,FileName);          
    end    
end

