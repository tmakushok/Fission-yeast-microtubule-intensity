%% This version of the function detects cell tips until the point where
%% cell axes cross cell outline (going inwards along these axes)
function [CellEnds, CellWidthEnds, CellsPixels, GoodCells] = f_4Tips(BigRegions, ImageFigNb, MinCellWidth, MaxCellWidth, MinCellLength, MaxCellLength)
%--------------------------------------------------------------------------
%!!!--!!! The value added to the automatically detected half cell width to
%-- get the initial position to obtain maximum values for width and length
    HWPlus = 5;
%--------------------------------------------------------------------------
    CellEnds = [];
    CellWidthEnds = []; 
    CellsPixels = [];
    GoodCells = [];
    
    [m, n] = size(BigRegions);
    Labels = bwlabel(BigRegions); 
    Stats = regionprops(Labels, 'Centroid', 'MajorAxisLength', 'MinorAxisLength', 'Orientation', 'PixelList');               
    if length(Stats) ~= 1
        CellEnds = []; 
        CellWidthEnds = []; 
        CellsPixels = []; 
        GoodCells = []; 
        return
    end
%% Measuring cell length and checking its conformity to normal values
% Calculation of the coordinates of cell ends in format [x1, y1; x2, y2]                        
    x0 = Stats.Centroid(1);
    y0 = Stats.Centroid(2);
    Angle = Stats.Orientation;    
    % Taking the intensity profile through the two ellipse ends
    HalfLength = (Stats.MajorAxisLength / 2) + HWPlus;
    x1 = x0 + HalfLength * cosd(Angle);
    y1 = y0 - HalfLength * sind(Angle); 
    x2 = x0 - HalfLength * cosd(Angle);
    y2 = y0 + HalfLength * sind(Angle);
    [IPrX, IPrY, IPr] = improfile(BigRegions, [x1 x2], [y1 y2]);
    figure, plot(IPr);
    % Finding steps between black and white
    Diff = IPr(1:length(IPr)-1) - IPr(2:length(IPr));
    ind = find(Diff);
    % If the line crosses more that 2 borders, there is a problem
    % or if the two borders are not of opposite nature
    if (length(ind) ~= 2) || (Diff(ind(1)) * Diff(ind(2)) ~= -1)     
        CellEnds = []; 
        CellWidthEnds = []; 
        CellsPixels = []; 
        GoodCells = []; 
        return
    end
    % Finding positions of cell ends
    E1x = IPrX(ind(1));     
    E2x = IPrX(ind(2)); 
    E1y = IPrY(ind(1)); 
    E2y = IPrY(ind(2));         
    % Check if cell length falls into the limits
    CellLength = sqrt((E2x - E1x)^2 + (E2y - E1y)^2);        
    if (CellLength < MinCellLength) | (CellLength > MaxCellLength)                               
        CellEnds = []; 
        CellWidthEnds = []; 
        CellsPixels = []; 
        GoodCells = []; 
        return
    end    
%% Measuring cell width and checking its conformity to normal values           
    HalfWidth = (Stats.MinorAxisLength / 2) + HWPlus;
    x1 = x0 - HalfWidth * sind(Angle);
    y1 = y0 - HalfWidth * cosd(Angle); 
    x2 = x0 + HalfWidth * sind(Angle);
    y2 = y0 + HalfWidth * cosd(Angle);
    [IPrX, IPrY, IPr] = improfile(BigRegions, [x1 x2], [y1 y2]);
    figure, plot(IPr);
    % Finding steps between black and white
    Diff = IPr(1:length(IPr)-1) - IPr(2:length(IPr));
    ind = find(Diff);
    % If the line crosses more that 2 borders, there is a problem
    % or if the two borders are not of opposite nature
    if (length(ind) ~= 2) || (Diff(ind(1)) * Diff(ind(2)) ~= -1)     
        CellEnds = []; 
        CellWidthEnds = []; 
        CellsPixels = []; 
        GoodCells = []; 
        return
    end
    % Finding positions of cell "width ends"
    W1x = IPrX(ind(1));     
    W2x = IPrX(ind(2)); 
    W1y = IPrY(ind(1)); 
    W2y = IPrY(ind(2)); 
    % Check if cell width falls into the limits
    CellWidth = sqrt((W2x-W1x)^2 + (W2y-W1y)^2);        
    if (CellWidth < MinCellWidth) || (CellWidth > MaxCellWidth)                 
        CellEnds = []; 
        CellWidthEnds = []; 
        CellsPixels = []; 
        GoodCells = []; 
        return
    end          
%% Preparing the output results
    CellCenter = [(E1x + E2x)/2, (E1y + E2y)/2];
    CellEnds = [E1x, E1y, E2x, E2y, CellCenter];
    CellWidthEnds = [W1x, W1y, W2x, W2y, CellCenter];
    CellsPixels = Stats.PixelList;
    GoodCells = [i, CellLength, CellWidth, Angle, CellCenter];
end
    
    
    
    
