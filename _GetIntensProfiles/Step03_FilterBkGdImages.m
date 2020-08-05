%% The task of the program
%% is to filter strongly the two images needed for background correction
close all;
clear;
%--------------------------------------------------------------------------
%!!!--!!! Paths
ImageBkGdFile = '_InputImages/WithLazer.tif';
ImageCCDBkGdFile = '_InputImages/WOLazer.tif';
BkGdImage_Output = '_OutputGI/Filtered_WithLazer.mat';
CCDBkGdImage_Output = '_OutputGI/Filtered_WOLazer.mat';
%--------------------------------------------------------------------------
CCDBkGdImage = double(imread(ImageCCDBkGdFile));        % CCD camera shading problem image
% !!! This step should not be done if both correstion images and all
% fluorescent images are taken with the same binning
BkGdImage = f_Bin2_GrayValues(double(imread(ImageBkGdFile)));     % Artificial bin2 
% BkGdImage = double(imread(ImageBkGdFile));     % Artificial bin2 

figure, imshow(BkGdImage, []);
figure, imshow(CCDBkGdImage, [180, 193]);
%% Strong filtering of both correction images
% First a bit of median filtering to get rid of some noise
BkGdImage = medfilt2(BkGdImage, [3 3]);
CCDBkGdImage = medfilt2(CCDBkGdImage, [3 3]);
figure, imshow(BkGdImage, []);
figure, imshow(CCDBkGdImage, [180, 193]);
% Then a strong average filtering
h = fspecial('average', 15);      
BkGdImage = imfilter(BkGdImage, h, 'replicate');
CCDBkGdImage = imfilter(CCDBkGdImage, h, 'replicate');
figure, imshow(BkGdImage, []);
figure, imshow(CCDBkGdImage, [180, 193]);
%% Subtraction of CCD camera shading problem image (filtered) from the illumination correction image   
BkGdImage = BkGdImage - CCDBkGdImage;
imshow(BkGdImage, []);
%% Normalising intensity levels of imhomogeneous lazer illumination image to 1
BkGdImage = BkGdImage / max(max(BkGdImage));
%% Output 
save(BkGdImage_Output, 'BkGdImage');  
%% Cropping of CCD problem image in the same way as analysed images are going to be cropped
% CCDBkGdImage = cropImage(BkGdImage, CCDBkGdImage);  % 'BkGdImage' is used as template for cropping
%% Visualisation
figure, imshow(BkGdImage, []);
figure, imshow(CCDBkGdImage, [180, 193]);   % Between [] the range of the image is defined for proper visualisation
%% Output
save(CCDBkGdImage_Output, 'CCDBkGdImage');  

