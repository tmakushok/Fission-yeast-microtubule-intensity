%% The task of the program is to create BF image that can be used to detect
%% cell borders (this involves finding middle Z-slice 
%% (which corresponds to the most 'greyish' image of the stack: no clear white or black rings around cells)
%% and making projection of some slices taken according to the position of
%% the middle slice.

%% To find the middle Z-plane histograms of image intensity are analysed:
%% middle slice corresponds to the most bottom (narrow) histogram
%% in the region on both sides of the middle strong peak

%% There is a check whether the active z slice does not go out
%% of Z-range (on the global level, not single-cell level: 
%% if middle plane is not too close to the borders of Z-range).
clear;
close all;
%% Parameters
Directory = '_InputImages/';
DirImOut = '_InputImages/';       % To output projection images
DirOut = '_Output/';
Nbins = 40;    % Nb of bins for histograms    
ZPlus = 1;      % Number of Z-slices up and down from the middle slice for the projection
%% Determining the number of z-slices in the movie: 'SlicesTotal'
SliceNbOld = -1;
SlicesTotal = 0;

AllFolders = dir([Directory '*BF']);
for i_Fold = 1:length(AllFolders)   % Loop on all folders
%     close all;
    ImFolder = [Directory AllFolders(i_Fold).name];
    Files = dir([ImFolder '/*.tif']);
    SlicesTotal = length(Files);    
    %% Collecting whole_image intensity histograms for images of a stack
    NhistAll = zeros(SlicesTotal, Nbins);
    for i_Slice = 1:SlicesTotal        
        FileName = [ImFolder '/' Files(i_Slice).name]; 
        InitImage = double(imread(FileName));         
%         figure, imshow(InitImage, []);        
        ImInLine = reshape(InitImage, 1, size(InitImage, 1) * size(InitImage, 2));
        [Nhist, Xhist] = hist(ImInLine, Nbins);
%         figure, bar(Xhist, Nhist);        
        NhistAll(i_Slice,:) = Nhist;    % Each line corresponds to one image of the stack
        XhistAll(i_Slice,:) = Xhist;           
    end
    %% Analysing collected histograms to find the most bottom one (at outer regions)
    %% that corresponds to the BF image taken axactly at the middle of the cell
    MiddlePlane = f_FindMiddlePlane(XhistAll, NhistAll) + 3;    
    % Show the middle plane image
    FileName = [ImFolder '/' Files(MiddlePlane).name]; 
    InitImage = double(imread(FileName));    
    figure, imshow(InitImage, []);
%     % Ask the user to confirm by pressing Enter
%     input('Please press Enter :-)');
    %% Do projection close to the middle plane
    FileName = [ImFolder '/'  Files(MiddlePlane - 2).name];     
    Image1 = double(imread(FileName)); 
    FileName = [ImFolder '/'  Files(MiddlePlane - 1).name]; 
    Image2 = double(imread(FileName)); 
    FileName = [ImFolder '/'  Files(MiddlePlane - 3).name]; 
    Image3 = double(imread(FileName)); 
    Image = (Image1 + Image2 + Image3) / 3;
%     figure, imshow(Image, []);
    %% Save the image
    SaveFileName = AllFolders(i_Fold).name;
    [UnderPos] = regexp(SaveFileName, '_');    
    FileNb = sprintf('%03.0f \n', str2num(SaveFileName(UnderPos(1)+1:UnderPos(2)-1)));             
    save([DirImOut 'BF_AVG_' FileNb(1:3) '.mat'], 'Image'); 
end










