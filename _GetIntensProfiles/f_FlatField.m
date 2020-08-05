%% The task of the function is to do flat-field correction of one image
function [InitImage] = f_FlatField(InitImage, CCDBkGdImage, BkGdImage, CropMask)
%% Subtraction of already prepared CCD camera problem image from the current image   
    InitImage = InitImage - CCDBkGdImage; 
    figure, imshow(InitImage, []); 
%% Cropping the image in a lemon shape (using the mask received as a parameter)
    InitImage = InitImage .* CropMask;   
    figure, imshow(InitImage, []); 
%% Finding the background (zeros obtained from cropping are ignored)
% (this step is done before filtering, not to get artefact pics on intensity histogram)
    BkGdValue = f_BkGd(InitImage);      
%% Filtering (a bit) of the image (is useful when pixel size is small)
%     h = fspecial('gaussian');   % The default value for hsize is [3 3]; the default value for sigma is 0.5. 
%     InitImage = imfilter(InitImage, h);    
%% Subtracting the 'natural' background           
    InitImage = InitImage - BkGdValue;
    figure, imshow(InitImage, []); 
%% Put negative values to zero    
    InitImage(InitImage < 0) = 0;   
%     figure, imshow(InitImage, []);      
%% Divide the subtracted image by the background image
    InitImage = InitImage ./ BkGdImage;                    
    figure, imshow(InitImage, []);     
