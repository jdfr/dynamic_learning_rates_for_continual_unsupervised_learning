function [Handle] = PlotGHNFConcepts(Model,U,Dictionary,FileName)
% Plot a GHNF where each prototype can be displayed as the most frequent 
% concept that represents the neuron. These concepts can then be related to
% terms from a dictionary. Training samples fed to the GHNF had to be
% previously reduced by using LSI.
% Inputs:
%   Model = the trained GHNF model
%   U = factor matrix from LSI
%   Dictionary = dictionary of words present the documents
%   FileName = name of the file to save the image

Handle = PlotGHNFConceptsRec(Model,U,Dictionary,1,1,FileName);
end

function [Handle] = PlotGHNFConceptsRec(Model,U,Dictionary,Level,Parents,FileName)
% Plot a GHNF layer where each prototype can be displayed as the most frequent 
% concept that represents the neuron. These concepts can then be related to
% terms from a dictionary. Training samples fed to the GHNF had to be
% previously reduced by using LSI.
% Inputs:
%   Model = the trained GHNF model
%   U = factor matrix from LSI
%   Dictionary = dictionary of words present the documents
%   Level = level of the hierarchy
%   Parents = parent unit indexes
%   FileName = name of the file to save the image

Handle = [];
if ~isempty(Model) && (Level <= 1)
        
    Handle = figure('Name',[FileName num2str(Level) '_' num2str(Parents)],'units','normalized','outerposition',[0 0 1 1]);    
    hold on
    MAX_WORDS = 4;
    
    Colors = distinguishable_colors(4);
    MyColor = Colors(Level,:);
    
    SpanningTree=biograph(Model.SpanningTree);
    SpanningTree.LayoutType = 'equilibrium';
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
            [MyConceptFreq,MyConceptNdx] = max(Model.Means(:,NdxUnit));                
            if isempty(Model.Child{NdxUnit}) || Level>2                
                plot(MyPos(1),MyPos(2),'o','LineWidth',1,'MarkerEdgeColor',MyColor,'MarkerFaceColor',MyColor,'MarkerSize',max(18*MyConceptFreq,12)); % node with layer color                
            else               
                plot(MyPos(1),MyPos(2),'o','LineWidth',1,'MarkerEdgeColor','yellow','MarkerFaceColor','yellow','MarkerSize',max(18*MyConceptFreq,12)); % yellow node                
            end
%             text(MyPos(1),MyPos(2),num2str(NdxUnit),'Color',[0 0 0],'FontSize',12,'HorizontalAlignment','center');  % unit number                  
            U2 = full(U);
            for NdxWord=1:MAX_WORDS
                WordOK = 0;
                while ~WordOK
                    [~,MyWordNdx] = max(U2(:,MyConceptNdx));
                    U2(MyWordNdx,MyConceptNdx) = -Inf;
                    MyWord = Dictionary{MyWordNdx};
                    if isnan(str2double(MyWord))   % si no es número                 
                        WordOK = 1;
                    end
                end
                ShowWord(MyPos,NdxWord,MyWord,MyConceptFreq);                                
            end            
        end
    end

    hold off
    axis off
    Figure2pdf([FileName num2str(Level) '_' num2str(Parents)],12,10);
    
%     Figure2pdf([FileName num2str(Level) '_' num2str(Parents) '_numerated'],12,10);
    
    % Draw children in different figures
    for NdxUnit=1:length(Model.Child)
        PlotGHNFConceptsRec(Model.Child{NdxUnit},U,Dictionary,Level+1,[Parents NdxUnit],FileName);          
    end    
end
end

function ShowWord(Position,NdxWord,MyWord,MyConceptFreq)
% Show a word around the position of a node (4 words are assumed as maximum).
% Inputs:
%   Position = coordinates of the node position
%   NdxWord = index of the word (1-4)
%   MyWord = word to represent
%   MyConceptFreq = frequency of the concept

switch NdxWord   
    case 1       
        text(Position(1),Position(2),MyWord,'Color',[0 0 0],'FontSize',max(10*MyConceptFreq,8),'HorizontalAlignment','center','VerticalAlignment','bottom');                    
    case 2
        text(Position(1)-3*MyConceptFreq,Position(2),MyWord,'Color',[0 0 0],'FontSize',max(10*MyConceptFreq,8),'HorizontalAlignment','right');                    
    case 3        
        text(Position(1)+3*MyConceptFreq,Position(2),MyWord,'Color',[0 0 0],'FontSize',max(10*MyConceptFreq,8));                            
    case 4        
        text(Position(1),Position(2),MyWord,'Color',[0 0 0],'FontSize',max(10*MyConceptFreq,8),'HorizontalAlignment','center','VerticalAlignment','top');                                                   
end
end
