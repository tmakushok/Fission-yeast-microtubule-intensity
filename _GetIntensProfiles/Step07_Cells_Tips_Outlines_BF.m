%% The task of the program is to find cell borders from bright field
%% projection images and to obtain geometrical parameters of the cells
function [] = Step07_Cells_Tips_Outlines_BF(CropArea)
% Minimal area of a cell, in pixels
MinArea = 1000;       
% Defining the size limits for objects to be recognized as being
% S.pombe cells, in pixels
MinCellWidth = 30;  
MaxCellWidth = 100; 
MinCellLength = 90; 
MaxCellLength = 350;    
% Basic name for images where we'll get cells from 
ListFileName = '_InputImages/BF_AVG_*.mat';
% Paths for flat-field correction for finding out of Z-range cells
ImageBkGdFile = '_OutputGI/Filtered_WithLazer.mat';
ImageCCDBkGdFile = '_OutputGI/Filtered_WOLazer.mat';
% Paths for outputs
PathAllCellsPixels = '_OutputGI/output_CellsPixels.mat';
PathAllGoodCells = '_OutputGI/output_GoodCellsParams.mat';
PathTMP = '_OutputGI/__TMP_Pixels_GoodCells';
PathOutput = '_OutputGI/output_CellsCenterEndsIntens.txt'; 
%% 
ImFiles = dir(ListFileName);  % Obtaining the list of files of BF projections
load(ImageCCDBkGdFile);         % CCD camera shading problem image
load(ImageBkGdFile);            % Inhomogeneous illumination correction image
AllCellsPixels = cell(length(ImFiles), 1);  % Will be used for results storage
AllGoodCells = cell(length(ImFiles), 1);    % Will be used for results storage
%% Preparing the output file
fid = fopen(PathOutput, 'w');
fprintf(fid, '%s\n', 'File|CellNb|Cell_Length|Cell_Width|AxisAngle|Cell_Center: x1|y1|Cell_Tips: x1|x2|y1|y2|Intensity_Density|Total_Intensity|Volume(in px3)|Errors on volume');
% fclose(fid); 
%% Preparing crop mask for bright field images
% Do binning 1/2 on 'CropArea' (Fluo is bin2, BF is bin1)
CropAreaNew = f_BinHalf(CropArea);    
imshow(CropAreaNew, []);
for i_File = 1:length(ImFiles)              % Loop on BF image files to analyse 
    i_File
    close all;
    %% Open BF image    
    File = ImFiles(i_File).name;
    BrFieldIm = load(['_InputImages/' File]);    
    BrFieldIm = BrFieldIm.Image;
    imshow(BrFieldIm, []);
    %% Crop bright field image in the same way as fluo images were     
    BrFieldIm = BrFieldIm .* CropAreaNew;
    imshow(BrFieldIm, []);
    %% Open fluo projection image
    FileFl = ['_InputImages/' regexprep(File, 'BF_AVG_', 'WoBg_AVG_')];                  
    FluoIm = load(FileFl);
    FluoIm = FluoIm.InitImage;
    %% Finding regions where cells are out of Z-range
    % Open the stack, extract bottom and top images for finding regions that have intensity out of the Z-range
    CorrStack = load(regexprep(FileFl, 'WoBg_AVG', 'STACK'));
    CorrStack = CorrStack.Stack;
    CorrImage1 = CorrStack(:,:,1); 
    CorrImage2 = CorrStack(:,:,length(CorrStack(1,1,:)));  
    % Flat-field correct the top and bottom images
    CorrImage1 = f_FlatField(CorrImage1, CCDBkGdImage, BkGdImage, CropArea); 
    CorrImage2 = f_FlatField(CorrImage2, CCDBkGdImage, BkGdImage, CropArea);      
    OutOfZMask = f_OutOfZRange(CorrImage1, CorrImage2);
    %% Obtaining cell parameters from transmission image       
    % Function for bright field images, 
    % second parameter is for visualisation of cell outlines on fluo image                                          
    % third parameter is the mask showing the regions that are out of Z-range
    % (to subtract cells coming out of Z-range)
    % 'GoodCells' contains:
    % CellNb|Cell_Lengths|Cell_Width|AxisAngle|Cell_Center:
    % x1|y1|Cell_Tips: x1|y1|x2|y2|Area   
    [GoodCells, CellsPixels] = f_CellParams_BF(BrFieldIm, FluoIm, OutOfZMask);  
    AllCellsPixels(i_File) = {CellsPixels};  % Each line corresponds to one image  
    AllGoodCells(i_File) = {GoodCells};     % Each line corresponds to one image  
    % Save the outlines_on_top_of_fluo image as control 
    saveas(gcf, ['_OutputGI/CellOutlines_Frame' int2str(i_File) '.fig'])    % 'gcf' returns the handle of the current figure    
    for i_GC = 1:length(CellsPixels)        % Loop on all cells
        %% Calculating total intensity and density of intensity 
        [TotInt, DensInt, CellVolume, VolumeErr] = f_Intensity(FluoIm, CellsPixels{i_GC}, GoodCells(i_GC, :));              
        %% Output of the coordinates of the width, center, both tips and cell
        %% intensity for the current cell        
        % File|CellNb|Cell_Length|Cell_Width|AxisAngle|Cell_Center:
        % x1|y1|Cell_Tips: x1|x2|y1|y2|Average_Intensity|Total_Intensity|Volume(in px3)|Errors on volume'        
        fprintf(fid, '%s\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n', FileFl, i_GC, GoodCells(i_GC, 2), GoodCells(i_GC, 3), GoodCells(i_GC, 4), ...
            GoodCells(i_GC, 5), GoodCells(i_GC, 6), GoodCells(i_GC, 7), GoodCells(i_GC, 9), GoodCells(i_GC, 8), GoodCells(i_GC, 10), DensInt, TotInt, CellVolume, VolumeErr(1), VolumeErr(2));                                                 
    end    
    % Temporary saving of the data, to have the possibility to stop the program at any moment
    save(PathTMP, 'AllCellsPixels', 'AllGoodCells');    
    figure, imshow(BrFieldIm, []);
end            
save(PathAllCellsPixels, 'AllCellsPixels');
save(PathAllGoodCells, 'AllGoodCells');





