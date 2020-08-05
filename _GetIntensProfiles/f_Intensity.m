%% The task of the function is to calculate total intensity and density of
%% intensity for all cells on the image
function [TotInt, DensInt, CellVolume, VolumeErr] = f_Intensity(FluoIm, Pixels, GoodCell)
TotInt = 0;         % For one cell
for i_px = 1:length(Pixels)      % Loop on the pixels of the cell 'i_c'
    xIn = Pixels(i_px, 1);
    yIn = Pixels(i_px, 2);                    
    TotInt = TotInt + double(FluoIm(yIn, xIn));                                        
end            
% What is written next is this: f_cellVolume(CellPixels, CellTips, Angle)
% VolumeErr is error on cell volume determination 
%   (normalised difference between final value and the two values for two cell halves)
[CellVolume, VolumeErr] = f_cellVolume(Pixels, [GoodCell(7), GoodCell(8);...
    GoodCell(9), GoodCell(10)], GoodCell(4));
DensInt = TotInt / CellVolume;      % Intensity of GFP per voxel inside the cell        
