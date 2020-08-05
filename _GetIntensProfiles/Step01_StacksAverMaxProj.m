%% The task of the program is to create average and maximum projections 
%% of stacks. Z-stacks of images obtained from image sequences are saved also.

%% !!! Working with 3D matrixes representing Z-stacks is very memory-consuming.
%% Here it is done because 3D matrices are easier to use while tracing
%% MT bundles in 3D rather than opening files corresponding to different
%% z-positions each time. But, of course, one can use separate z-slices files 
%% and get rid of 3D matrices (this will save some RAM-space 
%% and probably some time for saving and loading huge 3D matrices).
clear;
close all;
%% Parameters 
InputDir = '_InputImages/';
OutputDir = '_InputImages/';
% Characteristic part of the name of image folders (to distinguish them from other files)
% !!! If all folders are image folders, then one can separate them from
% simple files using 'isdir' function instead of defining specific part of
% the name indicating that it is an image folder
FoldPart = 'DB*'; 
% Characteristic part of the name of bright field images (to distinguish them
% from other files and not analyse them in this script)
BF_Part = '_BF';
%--------------------------------------------------------------------------
AllFolders = dir([InputDir FoldPart]);
for i_Fold = 1:length(AllFolders)           % Loop on all folders    
    ImFolder = [InputDir AllFolders(i_Fold).name];
    % If the folder contains only bright field images, don't consider it
    if ~isempty(strfind(ImFolder, BF_Part))
        continue
    end          
    AllFiles = dir([ImFolder '/' '*.tif']);   
%% Read images into a stack  
    NbSlices = length(AllFiles);   % Number of Z-slices to be in the stack
    % Reading the first slice to know the size of the images
    FileName = [ImFolder '/' AllFiles(1).name];   
    InitImage = double(imread(FileName)); 
    s = size(InitImage);
    Stack = zeros(s(1), s(2), NbSlices);    % Preallocation
    % Filling preallocated space with actual images
    for i_Slice = 1:NbSlices 
        FileName = [InputDir AllFolders(i_Fold).name '/' AllFiles(i_Slice).name];     
        Stack(:,:,i_Slice) = double(imread(FileName));
%         imshow(Stack(:,:,i_Slice), []);
%         pause(0.1);
    end        
%% Performing projections of the stack    
    MaxProj = Stack(:,:,1);      % Initialisation with the first z-slice
    SumProj = Stack(:,:,1); 
    for i = 2:NbSlices 
        MaxProj = max(MaxProj, Stack(:,:,i));
        SumProj = SumProj + Stack(:,:,i);
    end  
    AverProj = SumProj / NbSlices;
    imshow(MaxProj, []);
%     figure, imshow(AverProj, []);    
%     pause(0.01);
%% Saving images
    % Zeros have to be added at the beginning of the file number string to have files
    % sorted in the right way (not to have 'File10' coming before 'File3')
    SaveFileName = AllFolders(i_Fold).name;
    [UnderPos] = regexp(SaveFileName, '_');    
    FileNb = sprintf('%03.0f \n', str2double(SaveFileName(UnderPos(1)+1:length(SaveFileName))));             
    
    save([OutputDir 'STACK_' FileNb(1:3) '.mat'], 'Stack'); 
    save([OutputDir 'MAX_' FileNb(1:3) '.mat'], 'MaxProj'); 
    save([OutputDir 'AVG_' FileNb(1:3) '.mat'], 'AverProj');       
end

