%% The function determines if the contour is the most outer contour
%% by looking if a bit shrinked version of the contour lies above the white BF
%% region around current cell
% 'Binary' = black image with just one cell mask on it 
function [Res] = f_IsThirdCont(Binary, BF, BF_BkGd)
% Erode the picture a bit
Packed = bwpack(Binary);
% In the erosion parameter previous dilation is taken into account 
Packed = imerode(Packed, strel('disk', 6, 0), 'ispacked', size(Binary,1));
Binary = bwunpack(Packed, size(Binary, 1));
Binary = bwperim(Binary);
ind = find(Binary);
TestIntens = sum(BF(ind)) / length(ind);
if TestIntens > BF_BkGd
    Res = 1;
else
    Res = 0;
end