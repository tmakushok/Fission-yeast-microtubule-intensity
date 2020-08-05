%% Taking off the cells that are not well segmented (with mouse clicks)
[CellOff_x, CellOff_y] = ginput(MaxBadCells);  % User mouse-clicks-in coordinates, 'Enter' to continue    
for i_Off = 1:length(CellOff_x)     % Loop on the cells clicked        
    i_Pix = 1;
    while i_Pix < (length(CellsPixels) + 1)            
        % Finding lines in CellsPixels{i_Pix} in which we have x =
        % clicked_x and y = clicked_y              
        LinCl = find(int16(CellsPixels{i_Pix}(:, 1)) == int16(CellOff_x(i_Off)) & int16(CellsPixels{i_Pix}(:, 2)) == int16(CellOff_y(i_Off)));                                    
        if ~isempty(LinCl) % if it is the cell number i_Pix that contain clicked coordinates                 
            CellEnds(i_Pix, :) = [];  
            CellWidthEnds(i_Pix, :) = [];
            GoodCells(i_Pix, :) = [];
            CellsPixels(i_Pix) = [];                
            break                          
        end
        i_Pix = i_Pix + 1;
    end
end  