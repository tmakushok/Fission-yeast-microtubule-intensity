%% The task of the function is to do binning 1/2 on a binary image
function [CropAreaNew] = f_BinHalf(CropArea)
    s = size(CropArea);
    CropAreaNew = zeros(s * 2);
    for i_Lin = 1:s(1)
        % Creating next line: line with 0 after each element
        NewLine = zeros(1, s(2)*2);
        NewLine(1:2:(s(2)*2)) = CropArea(i_Lin, :);                
        CropAreaNew(i_Lin * 2 - 1:i_Lin * 2, :) = [NewLine; NewLine];
    end
    imshow(CropAreaNew, [0 0.5]);
    % Average filtering to create continuous crop area instead of stripes
    h = fspecial('average', 3);        
    CropAreaNew = imfilter(CropAreaNew, h, 'replicate');     
    CropAreaNew(find(CropAreaNew > 0)) = 1;
    imshow(CropAreaNew, []);