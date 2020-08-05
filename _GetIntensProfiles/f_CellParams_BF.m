%% The task of the program is to find cell cortex border and determine cell
%% parameters (cell ends, cell length, angle, profile of cell width) 
%% for one image
% 'GoodCells' array contains: cell ends(1,2,3,4 (x1,y1,x2,y2)), 
% cell axis angle in degrees(5), cell width close to cell end(6), cell length(7)
function [GoodCells, CellsPixels] = f_CellParams_BF(InitImage, FluoIm, OutOfZMask)
%--------------------------------------------------------------------------
%!!!--!!! Maximal distance between centers of two detected areas when they
%-- are still considered as the same cell detected twice
DistDeDoubling = 45;
%!!!--!!! Threshold on average intensity along the contour to decide if it
% is the 3rd contour
Thres3Cont = 3;
%!!!--!!! Minimal area of a cell, in pixels
MinArea = 1800;  
MaxArea = 100000;
MinAreaCont = 200;       % Min area of a contour to be kept
%!!!--!!! Defining the size limits (in pixels) for objects to be recognized as being
%-- S.pombe cells
MinCellWidth = 15;  
MaxCellWidth = 70; 
MinCellLength = 45; 
MaxCellLength = 200; 
%!!!--!!! Figure number to accumulate detected contours
AllBinCellsFigNb = 100;
%--------------------------------------------------------------------------
%% Initialisations
Nb_Cell = 1;    
CellEnds = [];
CellWidthEnds = [];
CellsPixels = cell(0,0);
GoodCells = []; 
DoneFlag = [];
AllCellsBinary = zeros(size(InitImage));
%% Subtract strongly filtered image
AverFiltered = medfilt2(InitImage, [30 30], 'symmetric');
figure, imshow(AverFiltered, []);
AverFiltered = InitImage - AverFiltered; 
% AverFiltered(AverFiltered < 0) = 0;
figure, imshow(AverFiltered, []);
%% Filtering of the initial image to smooth the background random changes in intensity    
AverFiltered = medfilt2(AverFiltered, [5 5]);
h = fspecial('average', 3);       
AverFiltered = imfilter(AverFiltered, h, 'replicate');      
% figure; imshow(AverFiltered, []);
%% Determining the edges of the cells using 'edge'       
[CannyResult, a] = edge(AverFiltered, 'canny');
% figure;  imshow(CannyResult);         
%% Taking off very short contours
BigContours = bwareaopen(CannyResult, MinAreaCont); 
% figure; imshow(BigContours, []); 
Labels = bwlabel(BigContours);
StatsCont = regionprops(Labels, 'PixelList');   
for i = 1:length(StatsCont)            % Loop on all contours    
    close all;
%% Creating black image with just one white contour on it    
    OneCellInNature = zeros(size(InitImage));
    Pxs = StatsCont(i).PixelList; 
    OneCellInNature(sub2ind(size(OneCellInNature), Pxs(:,2), Pxs(:,1))) = 1;
    % Visualise
%     figure, imshow(AverFiltered, []);
%     hold on, plot(Pxs(:,1), Pxs(:,2), 'o');
%% Dilation + filling
    % Dilate to fill in tiny holes in the contour
    Packed = bwpack(OneCellInNature);
    Packed = imdilate(Packed, strel('disk', 2, 0), 'ispacked', size(OneCellInNature,1));
    OneCellInNature = bwunpack(Packed, size(OneCellInNature, 1));
    % Fill in the contour
    OneCellInNature = imfill(OneCellInNature, 'holes');
    figure, imshow(OneCellInNature, []);
%% Measuring parameters of the filled cells
    Labels = bwlabel(OneCellInNature);
    Stats = regionprops(Labels, 'Area');           
%% Discarding cells that have an area that is too small or too big (on Bin1 BF image)
%% (unfilled contours will also disappear)
    if (Stats.Area < MinArea) || (Stats.Area > MaxArea)
        continue
    end        
%% Adding next 'layer' to the segmented image (purely for visualisation purposes): 
% next cell on the black BkGd added to the ones accumulated previously     
    AllCellsBinary = AllCellsBinary + OneCellInNature;
%     figure, imshow(AllCellsBinary, []);  
%% Check if current contour is external one. If yes, not analyse it.
% (based on the fact that external contour is above white region around cells on BF image)  
    if f_IsThirdCont(OneCellInNature, AverFiltered, Thres3Cont)
        continue
    end
%% Erode cell borders to have a closer to reality position of the border
    Packed = bwpack(OneCellInNature);
    % Should be eroded at least pixel more than previously done dilation of the contour
    % to get rid of eventual 'tails' coming from the bulk cell mask
    Packed = imerode(Packed, strel('disk', 6, 0), 'ispacked', size(OneCellInNature,1));
    OneCellInNature = bwunpack(Packed, size(OneCellInNature, 1));
    if max(max(OneCellInNature)) == 0       % If after erosion the whole object is gone (it was a line)
        continue
    end
%% Binning2 of the binary image to make cell mask correspond 
%% to cell's position on fluorescence image 
    [OneCellInNature] = f_Bin2(OneCellInNature);
%% Checking if the current cell is not in the region with fluorescence out of Z-range    
    if ~isempty(find(OneCellInNature + OutOfZMask == 2))
        continue
    end      
%% Finding the four cell tips                        
    [OneCell_CellEnds, OneCell_CellWidthEnds, OneCell_CellsPixels, OneCell_GoodCells] = ... 
        f_4Tips(OneCellInNature, AllBinCellsFigNb, MinCellWidth, MaxCellWidth, MinCellLength, MaxCellLength);             
    % If no good cell was detected or after image erosion 
    % two cells were created instead of one, ignore this object
    if isempty(OneCell_GoodCells) || isempty(OneCell_CellEnds) || size(OneCell_GoodCells, 1) > 1  
        continue
    end
%% Accumulation of data for all the cells in the image
    CellEnds = [CellEnds; OneCell_CellEnds];
    CellWidthEnds = [CellWidthEnds; OneCell_CellWidthEnds];
    CellsPixels = [CellsPixels; OneCell_CellsPixels];
    Stats = regionprops(Labels, 'Area'); % Tis time it is for Bin2 (fluo image)
    % 'GoodCells' contains:
    % CellNb|Cell_Lengths|Cell_Width|AxisAngle|Cell_Center: x1|y1|Cell_Tips: x1|x2|y1|y2|Area 
    GoodCells = [GoodCells; Nb_Cell, OneCell_GoodCells(2:6), OneCell_CellEnds(1:4), Stats.Area];

    CellCenter = OneCell_GoodCells(5:6);
    Nb_Cell = Nb_Cell + 1;                                                
%% Solving the problem of having two or three detected areas overlapping 
%% (steming out from inner and outer edges of a cell) with slightly different widths             
    [LenWE, a] = size(CellWidthEnds);  
    if LenWE < 2
        continue
    end
    % Check if this cell was not already treated 
    % (if it was, then the third inner contour is looked at and we want to avoid keeping it)
    s_F = size(DoneFlag);
    if ~isempty(DoneFlag)
        FlagDist = sqrt((DoneFlag(1:s_F(1), 1) - CellCenter(1)).^2 + (DoneFlag(1:s_F(1), 2) - CellCenter(2)).^2);                
        if min(FlagDist) <= DistDeDoubling
            % Taking away the last cell's information from geometry arrays
            [LinG, a] = size(CellEnds);
            CellWidthEnds(LinG, :) = [];
            CellEnds(LinG, :) = []; 
            GoodCells(LinG, :) = []; 
            CellsPixels(LinG) = [];
            continue
        end        
    end
    % Finding the good contour to keep
    for i_UnDoubl = (LenWE - 1):-1:1  % back count as more chances to have an overlap with one of the 'closely previous' areas             
        % Calculation of the distance between current cell center and the
        % cell number i_DeDoubling center
        CentersDist = sqrt((CellWidthEnds(i_UnDoubl, 5) - CellCenter(1))^2 + (CellWidthEnds(i_UnDoubl, 6) - CellCenter(2))^2);                
        if CentersDist <= double(DistDeDoubling)    % if centers of the two cells are very close
            CellWidthOld = GoodCells(i_UnDoubl, 3);     
            CellWidth = OneCell_GoodCells(3);            
            if CellWidth > CellWidthOld                                                
                % Remove those lines also in geometry matrixes                        
                CellWidthEnds(i_UnDoubl, :) = [];
                CellEnds(i_UnDoubl, :) = [];
                GoodCells(i_UnDoubl, :) = [];
                % Taking away from the pixels set
                CellsPixels(i_UnDoubl) = [];
                % Adding the cell center in the list of cells treated
                % (to avoid keeping the third (most inner) outline)
                DoneFlag = [DoneFlag; CellCenter];               
                break
            else                                               
                % Taking away the last cell's information from
                % geometry arrays
                [LinG, a] = size(CellEnds);
                CellWidthEnds(LinG, :) = [];
                CellEnds(LinG, :) = []; 
                GoodCells(LinG, :) = []; 
                CellsPixels(LinG) = [];
                % Adding the cell center in the list of cells treated
                % (to avoid keeping the third (most inner) outline)
                DoneFlag = [DoneFlag; CellCenter];
                break
            end
        end
    end          
end              
%% Representing borders of the segmented image overlaid with fluorescence image 
% Creating cell mask
FinalImage = zeros(size(FluoIm));
% Filling with 1s the places where 'final good' cells are detected
for i_vis = 1:length(CellsPixels)       
    FinalImage(sub2ind(size(FinalImage), CellsPixels{i_vis}(:, 2), CellsPixels{i_vis}(:, 1))) = 1;     
end
figure, imshow(FinalImage, []);

BWoutline = bwperim(FinalImage);   % To find perimeter pixels in binary image
Segout = FluoIm;   
Segout(BWoutline) = max(max(Segout));     
figure, imshow(Segout, []);    
% Putting geometrical parameters on this overlaid image
Len_Res = length(CellsPixels);     
for i_vis = 1:Len_Res    
    % Visualisation of the cell length 
    line([CellEnds(i_vis, 1), CellEnds(i_vis, 3)], [CellEnds(i_vis, 2), CellEnds(i_vis, 4)]);            
    % Visualisation of the cell width            
    line([CellWidthEnds(i_vis, 1), CellWidthEnds(i_vis, 3)], [CellWidthEnds(i_vis, 2), CellWidthEnds(i_vis, 4)]);                    
    % Visualisation of the cell center        
    line(GoodCells(i_vis, 5), GoodCells(i_vis, 6), 'Color', [.8 0 0], 'Marker', 'o');  
end    