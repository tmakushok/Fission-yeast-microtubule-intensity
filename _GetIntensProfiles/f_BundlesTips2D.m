%% This function is designed to find potential MT bundles tips
function [MTTipsI, MTTipsJ] = f_BundlesTips2D(InitImage, CellMask)
%!!!--!!! Threshold for BkGd subtraction on LoG-filtered images
BkGdThres = 0.3;        %1.2;       
%!!!--!!! Min skeleton length to be considered as a bundle
MinSkeletLen = 4;
%--------------------------------------------------------------------------
BundlesLengths = [];        % To store length of each bundle with its number and number of the cell
TotalMTLengthAllCells = [];
TotalMTLengthAllCells_PerCellLengthUnit = [];
FigNb = 1;

InitFigNb = figure(), imshow(InitImage, []);
%% Using Laplacian of Gaussian for smoothing microtubules preserving their tips 
%-- Reveal MTs and eliminate BkGd noise
h = fspecial('log', 5, 2);  % 9 4
LapGaus = imfilter(InitImage, h, 'replicate'); 
% LapGaus = imfilter(ProcessedImage, h, 'replicate'); 
%% Invert colours on the image
LapGaus = -1 * LapGaus;
figure, imshow(LapGaus, []); 
%% Use non-dilated mask on the filtered image to get rid of the border
%% effect after LoG
LapGaus = LapGaus .* CellMask;
%% Binarisation using thresholding
BinaryIm = zeros(size(LapGaus));
BinaryIm(LapGaus > BkGdThres) = 1;   
figure, imshow(BinaryIm, []);  
%% Skeletonization
Skelet = bwmorph(BinaryIm, 'thin', Inf);
s = size(Skelet);
%% Keeping only long 'skelets'
Skelet = bwareaopen(Skelet, MinSkeletLen);  
% figure, imshow(Skelet);      
%% Finding all MT ends used later for 3D tracking
% Calculating an image where at a previously white on 'Skelet' image pixel 
% there is the sum of the pixel and of all surrounding pixels
WhiteSkelPos = find(Skelet);        % Linear indexes of all white pixels
s = size(Skelet);
ImMTEnds = zeros(s);
for i = 1:length(WhiteSkelPos)    
    [i_Lin, j_Lin] = ind2sub(s, WhiteSkelPos(i));
    ImMTEnds(i_Lin, j_Lin) = Skelet(i_Lin - 1, j_Lin - 1) + ...
        Skelet(i_Lin - 1, j_Lin) + Skelet(i_Lin - 1, j_Lin + 1) + ...
        Skelet(i_Lin, j_Lin - 1) + Skelet(i_Lin, j_Lin) + ...
        Skelet(i_Lin, j_Lin + 1) + Skelet(i_Lin + 1, j_Lin - 1) + ...
        Skelet(i_Lin + 1, j_Lin) + Skelet(i_Lin + 1, j_Lin + 1);
end
%% Visualisation of MT bundle tips found
[MTTipsI MTTipsJ] = find(ImMTEnds == 2);

figure, imshow(InitImage, []), hold on;
for i = 1:length(MTTipsI)  
    line(MTTipsJ(i), MTTipsI(i), 'Color', [.8 0 0], 'Marker', 'o', 'MarkerSize', 3);       
end
hold off






