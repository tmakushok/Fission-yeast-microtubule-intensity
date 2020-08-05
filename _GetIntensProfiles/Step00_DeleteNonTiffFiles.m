%% The task of the program is to delete all non-tiff files from images
%% folder
clear;
close all;
%--------------------------------------------------------------------------
% Where folders to clean are situated
ToClean = '_InputImages/';
%--------------------------------------------------------------------------
AllFolders = dir(ToClean);
% The loop has to start with 3 because first folders in the list are '.'
% and '..': current folder and parent folder
for i_Fold = 3:length(AllFolders)   % Loop on all folders
    AllFiles = dir([ToClean AllFolders(i_Fold).name]);  % List of files in the folder
    for i_File = 3:length(AllFiles)      % Loop on the files inside the folder  
        FileName = AllFiles(i_File).name;  
        IsTiff = findstr(FileName, '.tif');     % Empty string if not tiff (or if the file does not have the right extension)
        if isempty(IsTiff)
            delete([ToClean  AllFolders(i_Fold).name '/' FileName]);
        end    
    end
end


