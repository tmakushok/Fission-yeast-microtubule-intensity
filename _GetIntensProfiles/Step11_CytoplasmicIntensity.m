%% The task of the program is to measure cytoplasmic intensity 
%% (without the signal from microtubules) of the cells
%% detected previously by analysing the distribution of intensities of
%% voxels belonging to the cells
clear;      % To erase variables that are in memory
close all;
%% Parameters
MicronsInPixel = 0.13;   % RS bin2            % 0.065;   RS bin1
% Ratio between Z-step in stacks (in microns) and resolution in XY- plane 
ZspacePx = 0.25 / MicronsInPixel;       % Distance between Z-planes(in X-Y pixel size units)
% File names
% PathInput = '_OutputGI/output_CellsCenterEndsIntens.txt';
PathAllGoodCells = '_OutputGI/output_GoodCellsParams.mat';
PathCellPixels = '_OutputGI/output_CellsPixels.mat';
PathMidPlanes = '_OutputGI/AllMiddlePlanes.mat';
PathOutputCytoplMax = '_OutputGI/output_CytoplMax_AverProj.txt';
PathCytoplInt = '_OutputGI/_CytoplIntens.mat';
%%
load(PathCellPixels);  % 'AllCellsPixels': one line is one image
load(PathAllGoodCells);  % 'AllGoodCells': one line is one image
load(PathMidPlanes); % 'AllMiddlePlanes': middle Z-positions used 
StackFilesBase = dir('_InputImages/WoBg_STACK*');    % List intensity corrected stacks we are going to play with  
CytInt = [];
%% Reading cell parameters from the input file 
for i_Stack = 1:length(StackFilesBase)            % Loop on 3D stacks
    %% Open the 3D stack containing detected cells  
    StackFile = StackFilesBase(i_Stack).name;
    Stack = load(['_InputImages/' StackFile]);
    Stack = Stack.Stack;
    % Load max projection image        
    k1 = findstr(StackFile, 'K_');
    k2 = findstr(StackFile, '.mat');
    InitIm_File = ['_InputImages/MAX_' StackFile(k1 + 2:k2 - 1) '.mat'];
    MaxProj = load(InitIm_File); 
    MaxProj = MaxProj.MaxProj;    
    % Extracting info for the cells detected on the current image
    GoodCells = AllGoodCells{i_Stack};
    Zmiddle = AllMiddlePlanes(i_Stack);
    %% Loop on the cells detected on the current image
    for i_c = 1:length(AllCellsPixels{i_Stack})    % Loop on all cells that we have on this image
        close all                
        IntensZAll = [];    % To store all voxel intensities
        PixelList = AllCellsPixels{i_Stack}{i_c};
%         s_OneCell = [max(max(PixelList(:, 2))) + Rdil, max(max(PixelList(:, 1))) + Rdil];
        % Creating a mask that contains white pixels at the location of the current cell
        CellMask = zeros(size(MaxProj));
        CellMask(sub2ind(size(CellMask), PixelList(:,2), PixelList(:,1))) = 1;   
        % Creating a matrix of cell outline coordinates
        [PerimPts(:,2), PerimPts(:,1)] = find(bwperim(CellMask)); % 'Outline' is [X,Y]    
        figure, plot(PerimPts(:, 1), PerimPts(:, 2), 'r.');
        %% Determining the 3D coordinates of all cell voxels
        % Retrieve the parameters of the cell number i_cell        
        Angle = GoodCells(i_c, 4);        
        CellTips = [GoodCells(i_c, 7), GoodCells(i_c, 8); GoodCells(i_c, 9), GoodCells(i_c, 10)]; % [x1, y1; x2, y2]                        
        %line(CellsEnds_x, CellsEnds_y);     % Visualisation of the position of the cell tips                  
        %% Going into reference system with bottom cell tip as 0 and rotating
        % Choosing the right tip
        if CellTips(1, 2) > CellTips(2, 2)
            TheTip = CellTips(1,:);
        else
            TheTip = CellTips(2,:);
        end
        % Transforming coordinates
        PerimPtsRot = zeros(size(PerimPts));
        % figure, hold on
        AngleRot = -Angle;
        % Putting 0 of the reference system in the bottom cell tip
        PerimPts(:, 1) = PerimPts(:, 1) - TheTip(1);
        PerimPts(:, 2) = PerimPts(:, 2) - TheTip(2);
        % Rotating the cell outline to have cell axis lying on the Ox
        PerimPtsRot(:, 1) = PerimPts(:, 1) * cosd(AngleRot) + PerimPts(:, 2) * sind(AngleRot); 
        PerimPtsRot(:, 2) = - PerimPts(:, 1) * sind(AngleRot) + PerimPts(:, 2) * cosd(AngleRot);
        figure, plot(PerimPtsRot(:, 1), PerimPtsRot(:, 2), 'r.');
        %% Take intensities from the middle plane
        IntensZ = Stack(:, :, Zmiddle) .* CellMask;
        figure, imshow(IntensZ, []);
        IntensZ = IntensZ(find(IntensZ));
        figure, hist(IntensZ, 30);
        IntensZAll = [IntensZAll; IntensZ];
        %% Loop on the Z-planes of the stack, to get intensities from each plane        
        % For planes other than the middle one, starting from the ones
        % closest to the middle plane
        for i_Z = Zmiddle - 1 : -1 : 1  % Half of the stack, as the cell is symmetric   
            dZ = (Zmiddle - i_Z) * ZspacePx; % Difference (in X-Y pixel size units) between Z of the middle plane and the current plane
            PerimShrRot(:, 1) = PerimPtsRot(:, 1);
            for i_pt = 1:size(PerimPtsRot, 1)
                Y = PerimPtsRot(i_pt, 2);
                Sin = dZ / Y;
                if abs(Sin) > 1
                    PerimShrRot(i_pt, 2) = 0;
                else
                    PerimShrRot(i_pt, 2) = Y * cosd(asind(Sin));
                end
            end
            figure, plot(PerimShrRot(:, 1), PerimShrRot(:, 2), 'r.');
            % Going back to real image coordinates
            % Rotating the shrinked cell outline 
            PerimShrinked(:, 1) = PerimShrRot(:, 1) * cosd(Angle) + PerimShrRot(:, 2) * sind(Angle); 
            PerimShrinked(:, 2) = - PerimShrRot(:, 1) * sind(Angle) + PerimShrRot(:, 2) * cosd(Angle);
            figure, plot(PerimShrinked(:, 1), PerimShrinked(:, 2), 'r.');
            % Changing reference system origin
            PerimShrinked(:, 1) = round(PerimShrinked(:, 1) + TheTip(1));
            PerimShrinked(:, 2) = round(PerimShrinked(:, 2) + TheTip(2));
            figure, plot(PerimShrinked(:, 1), PerimShrinked(:, 2), 'r.');
            %% Creation of the cell mask for the current Z-plane (shrinked in comparison to the middle plane)                        
            Mask = zeros(size(Stack(:, :, i_Z)));            
            Mask(sub2ind(size(Mask), PerimShrinked(:,2), PerimShrinked(:,1))) = 1;
            % Dilation of the contour
            Packed = bwpack(Mask);
            Packed = imdilate(Packed, strel('disk', 1, 0), 'ispacked', size(Mask,1));
            Mask = bwunpack(Packed, size(Mask, 1));
            % Filling in the outline to produce cell mask for the current
            % Z-position
            Mask = imfill(Mask, 'holes');
            figure, imshow(Mask);
            %% Taking intensities at current Z-plane (plane less than middle) at the cell mask coordinates
            IntensZ = Stack(:, :, i_Z) .* Mask;
            figure, imshow(IntensZ, []);
            IntensZ = IntensZ(find(IntensZ));
            figure, hist(IntensZ, 30);
            IntensZAll = [IntensZAll; IntensZ];
            %% Taking intensities at current Z-plane (plane more than middle) at the cell mask coordinates            
            IntensZ = Stack(:, :, Zmiddle + (Zmiddle - i_Z)) .* Mask;
            figure, imshow(IntensZ, []);
            IntensZ = IntensZ(find(IntensZ));
            figure, hist(IntensZ, 30);
            IntensZAll = [IntensZAll; IntensZ];
        end          
    end
end