%% The task of the program is to define a region (with threshold)
%-- where inhomogeneous illumination image has acceptable levels.
%-- This image is already normalised to 1.
%-- Final area is going to be ellipse-shaped
%-- the rest of the image becomes 0
function [InImage, GoodArea] = cropImage(BkGdImage, Thres, InImage)
GoodArea = zeros(size(BkGdImage));
GoodArea(find(BkGdImage > Thres)) = 1;
if nargin > 2      
    InImage = InImage .* GoodArea;
%     figure, imshow(GoodArea,[]);
%     figure, imshow(InImage,[]);      
else             % if no input image is supplied, crop region is returned
    InImage = GoodArea;    
end
