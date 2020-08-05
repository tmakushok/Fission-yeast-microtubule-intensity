close all;
clear;
%% Preparing image files
% Step00_DeleteNonTiffFiles;
% Step01_StacksAverMaxProj;
Step02_BrFieldProj;
%% Flat-field correction
Step03_FilterBkGdImages;
% Make one lemon-shaped mask (for the region that is properly illuminated with the lazer) 
% instead of calling f_cropImage every time
ImageBkGdFile = '_OutputGI/Filtered_WithLazer.mat';
BkGdImage = load(ImageBkGdFile);       
BkGdImage = BkGdImage.BkGdImage;
CropMask = f_cropImage(BkGdImage, 0.5);
Step05_IlluminationCorrection_AVG(CropMask);
%% Cell detection
% If there is a correction image for bright field then do flat field correction 
% for bright field instead of subtracting strongly filtered image
Step07_Cells_Tips_Outlines_BF(CropMask);















