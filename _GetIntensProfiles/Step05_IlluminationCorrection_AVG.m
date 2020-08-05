%% The task of the program is to correct average projection images for 
%% CCD camera shading problem and for inhomogeneous lazer illumination
function [] = Step05_IlluminationCorrection_AVG(CropMask)
%--------------------------------------------------------------------------
ImageBaseName = '_InputImages/';
ImageBkGdFile = '_OutputGI/Filtered_WithLazer.mat';
ImageCCDBkGdFile = '_OutputGI/Filtered_WOLazer.mat';
%--------------------------------------------------------------------------
load(ImageCCDBkGdFile);         % CCD camera shading problem image
load(ImageBkGdFile);            % Inhomogeneous illumination correction image
AverageImFiles = dir([ImageBaseName 'AVG_*']);      % List of average projection images
for i_File = 1:length(AverageImFiles)               % Loop on the image files to analyse  
    close all;        
    FileName = AverageImFiles(i_File).name;  
    FilePath = [ImageBaseName FileName];        
    InitImage = load(FilePath); 
    InitImage = InitImage.AverProj; 
    figure, imshow(InitImage, []); 
    %% Correction of the current file
    InitImage = f_FlatField(InitImage, CCDBkGdImage, BkGdImage, CropMask);    
    %% Output
    ImWoBkGd_File = [ImageBaseName 'WoBg_' FileName]; 
    save(ImWoBkGd_File, 'InitImage');       
end
