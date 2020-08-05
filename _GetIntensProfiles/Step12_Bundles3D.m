%% The task of the program is to find points corresponding to MT bundles in 3D
clear all;
close all;
%% Parameters
% Value that is given to all background regions on the image where 
% single cell is cropped
CytAver = 0;       % 200;
% Radius for dilation of the cell mask
Rdil = 3;
% Path
PathCellPixels = '_OutputGI/output_CellsPixels.mat';
%--------------------------------------------------------------------------
load(PathCellPixels);  % 'AllCellsPixels': one line is one image
AllMTs = cell(1);
CellNb = 1;

StackFilesBase = dir('_InputImages/STACK*');    % List of the image stacks we are going to play with  
StackFilesBase = StackFilesBase(1:length(StackFilesBase));      % Take off '.' and '..'
for i_Stack = 1:length(StackFilesBase)    
    CellNbOnImage = 1;
    AllMTsOneIm = cell(0,0);
    NextFileBase = StackFilesBase(i_Stack).name;
    Stack = load(['_InputImages/' StackFilesBase(i_Stack).name]);
    Stack = Stack.Stack;
    %% Load max projection image    
    NextFileBase = StackFilesBase(i_Stack).name;
    k1 = findstr(NextFileBase, 'K_');
    k2 = findstr(NextFileBase, '.mat');
%     InitIm_File = ['_OutputGI/WoBg_AVG_' NextFileBase(k1 + 2:k2 - 1) '.mat'];
    InitIm_File = ['_InputImages/MAX_' NextFileBase(k1 + 2:k2 - 1) '.mat'];
    InitImage = load(InitIm_File); 
    InitImage = InitImage.MaxProj;      
    %% Filtering all slices of the stack 
    NbSlices = length(Stack(1, 1, :));       
    %% Loop on the cells detected on the current image
    for i_cell = 1:length(AllCellsPixels{i_Stack})    % Loop on all cells that we have on this image
        close all        
        MTs = cell(1);
        PixelList = AllCellsPixels{i_Stack}{i_cell};
        s_OneCell = [max(max(PixelList(:, 2))) + Rdil, max(max(PixelList(:, 1))) + Rdil];
        %% Creating a mask that contains white pixels at the location of the current cell
        CellMask = zeros(s_OneCell);
        CellMask(sub2ind(size(CellMask), PixelList(:,2), PixelList(:,1))) = 1;
        % Dilate cell mask a bit, to have some extra pixels around the cell
        CellMaskDil = imdilate(CellMask, strel('disk', Rdil, 0));  
        %% Cropping the stack to have only the area of one cell  
        StackCell = zeros(s_OneCell(1), s_OneCell(2), NbSlices);
        for i_Slice = 1:NbSlices    
            Tmp = Stack(:,:, i_Slice);            
            Tmp = imcrop(Tmp, [1, 1, s_OneCell(2) - 1, s_OneCell(1) - 1]);
            Tmp = Tmp .* CellMask; 
            StackCell(:,:, i_Slice) = Tmp;                     
        end   
        %% Visualise all slices of the stack
%         for i = 1:NbSlices
%             figure();  
%             imshow(StackCell(:,:, i), []);
%         end
        %% MaxProj image for one cell 
        OneCellMaxProj = zeros(s_OneCell);
        OneCellMaxProj = imcrop(InitImage, [1, 1, s_OneCell(2) - 1, s_OneCell(1) - 1]);
        OneCellMaxProj = OneCellMaxProj .* CellMaskDil;
        imshow(OneCellMaxProj, []);
        %% Projection of nearby z slices (3 slices for every z position)
        %% ??? Max proj instead of aver proj
        LocalProj = cell(NbSlices, 1);
        LocalProj{1} = max(StackCell(:,:, 1), StackCell(:,:, 2));
        for i = 2:NbSlices - 1
            LocalProj{i} = max(max(StackCell(:,:,i - 1), StackCell(:,:,i)), StackCell(:,:,i + 1));
        end    
        LocalProj{NbSlices} = max(StackCell(:,:,NbSlices - 1), StackCell(:,:,NbSlices));           
        %% Getting potential MT bundles tips to work with   
        [MTTipsI, MTTipsJ] = f_BundlesTips2D(OneCellMaxProj, CellMask);      
        %% Finding continuous MTs in 3D, from starting points defined previously        
        for i_MT = 1:length(MTTipsI)       % Loop on all MT bundles of one cell
            [OneMT3D, h_Res] = f_MTpart3D(StackCell, MTTipsI(i_MT), MTTipsJ(i_MT), OneCellMaxProj, LocalProj);                
            if isempty(MTs{1})
                MTs{1} = OneMT3D;
            else
                MTs{1} = [MTs{1}, OneMT3D];
            end 
            % Visualise one MT bundle found
            figure; imshow(InitImage, []);
            hold on; plot(OneMT3D{1}(:, 1), OneMT3D{1}(:, 2), '-ro', 'MarkerSize', 5);
        end  
        %% Visualise MT bundles found for one cell
        h = figure(); imshow(InitImage, []);
        for i_MT = 1:length(MTTipsI)
            hold on; plot(MTs{1}{i_MT}(:, 1), MTs{1}{i_MT}(:, 2), '-r.', 'MarkerSize', 1);        
        end                
        saveas(h, ['_OutputGI/MTsCell_' int2str(CellNb) '.fig']);
        CellNb = CellNb + 1;
        %% Add bundles found in this cell to the total stock
        AllMTsOneIm = [AllMTsOneIm; MTs];
%         AllMTsOneIm(i_Stack, CellNbOnImage) = MTs;  
%         AllCellsPixels(i_File) = {CellsPixels};  % Each line corresponds to one image  
        CellNbOnImage = CellNbOnImage + 1;
    end
    AllMTs(i_Stack, 1) = {AllMTsOneIm};
end
% AllMTs{*} corresponds to an image
% a = AllMTs{*}(*) corresponds to a cell on this image
% a{*} corresponds to the list of MTs in this cell
save('_OutputGI/_BundlesPoints.mat', 'AllMTs');
















%% BackUp
% %% MaxProj image for one cell (projection- only at the cell pixels positions)
%         x0 = PixelList(1, 1);
%         y0 = min(PixelList(:, 2));
%         OneCellMaxProj = StackCell(:,:,1);
%         s = size(OneCellMaxProj);    
%         for x = x0:s(2)
%             for y = y0:s(1)
%                 Max = 0;
%                 for i_Slice = 1:NbSlices   
%                     Max = max(Max, StackCell(y, x, i_Slice));
%                 end
%                 OneCellMaxProj(y, x) = Max;
%             end
%         end     
%         imshow(OneCellMaxProj, []);

% For local max projection:
%         x0 = PixelList(1, 1);
%         y0 = min(PixelList(:, 2));
%         OneCellMaxProj = StackCell(:,:,1);
%         s = size(OneCellMaxProj);    
%         for x = x0:s(2)
%             for y = y0:s(1)
%                 Max = 0;
%                 for i_Slice = 1:NbSlices   
%                     Max = max(Max, StackCell(y, x, i_Slice));
%                 end
%                 OneCellMaxProj(y, x) = Max;
%             end
%         end     

%             for i_Pt = 1:length(PixelList)
%                 x = PixelList(i_Pt, 1);
%                 y = PixelList(i_Pt, 2);
%                 StackCell(y, x, i_Slice) = Stack(y, x, i_Slice);
%             end  
%             find(StackCell(:,:, i_Slice) == 0) = 200;

% save('_OutputGI/MTTips2D.mat', 'MTTipsI', 'MTTipsJ');