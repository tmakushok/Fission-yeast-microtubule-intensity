%% The task of the program is to correct maximum projection images for 
%% CCD camera shading problem and for inhomogeneous lazer illumination
%-- using coresponding images preferably collected on the same day
function [] = Step06_IlluminationCorrection_MAX(CropMask)
%--------------------------------------------------------------------------
PathOutput = '_OutputGI/';                
ImageBaseName = '_InputImages/';
ImageBkGdFile = '_OutputGI/Filtered_WithLazer.mat';
ImageCCDBkGdFile = '_OutputGI/Filtered_WOLazer.mat';
%--------------------------------------------------------------------------
CCDBkGdImage = load(ImageCCDBkGdFile);        % CCD camera problem image
CCDBkGdImage = CCDBkGdImage.CCDBkGdImage;
BkGdImage = load(ImageBkGdFile);       
BkGdImage = BkGdImage.BkGdImage;

AverageImFiles = dir([ImageBaseName 'MAX_*']);
for i_File = 1:length(AverageImFiles)              % Loop on the image files to analyse         
    FileName = AverageImFiles(i_File).name;  
    FilePath = [ImageBaseName FileName];        
    InitImage = load(FilePath);
    InitImage = InitImage.MaxProj;
    %% Subtraction of already prepared CCD camera problem image from the current image   
    InitImage = InitImage - CCDBkGdImage; 
    %% Cropping the image in a lemon shape
    InitImage = InitImage .* CropMask;
%     InitImage = f_cropImage(BkGdImage, CropThres, InitImage);    
    %% Finding the background (zeros obtained from cropping are ignored)
    % (before filtering, not to get artefact pics on intensity histogram)
    BkGdValue = f_naturalImageBkGd(InitImage);      
    %% Filtering (a bit) of the image (not necessary for bin 2)
    %% Actually, to get MT intensities out, better not to do filtering,
    %% otherwise values are lowered (because MTs are thin and intens. is smaller all around them)
%     h = fspecial('gaussian');   % The default value for hsize is [3 3]; the default value for sigma is 0.5. 
%     InitImage = imfilter(InitImage, h);    
    %% Subtracting the 'natural' background           
    InitImage = InitImage - BkGdValue;
%     InitImage(find(InitImage < 0)) = 0;  
%     figure, imshow(InitImage, []);       
    %% Divide the subtracted image by the background image
    InitImage = InitImage ./ BkGdImage;                    
%     figure, imshow(InitImage, []);     
%     pause(0.3);
    %% Output
    ImWoBkGd_File = [PathOutput 'WoBg_' FileName]; 
    save(ImWoBkGd_File, 'InitImage');       
end
