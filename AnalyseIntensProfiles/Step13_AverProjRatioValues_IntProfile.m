%% Taking values from aver. proj. at positions defined by previous MT tracking.
%% Background is subtracted at that position
clear; 
%% Parameters
% Taking intensity profiles along lines that are this much longer
% than the cell width determined automaticly
MoreWidth = 5;
% Parameter: how much the smoothing
% spline curve for initial data will be moved up: all points above are MTs
CutOffPlus = 10;
% Parameter defining the 'regidity' of the two fittings: 
% 1) to extract MTs; 2) to obtain values for cytoplasmic background
SmSplineParam1 = 0.001;       
SmSplineParam2 = 0.05;   
% File names
PathInputProfiles = '_OutputAI/output_Profiles_AverProj_AllCells_WoCenter.txt';  
PathLines = '_OutputAI/output_ProfileLines_WoCenter.txt';   
PathInputCytoplasmMax = '_OutputAI/output_CytoplMax_AverProj.txt';
PathMTs = '../__GetIntensProfiles/_OutputGI/UnDoubledMTs.mat';
PathCellsInfo = '../__GetIntensProfiles/_OutputGI/output_CellsCenterEndsIntens.txt';
ImagePathBase = '../__GetIntensProfiles/';
%%
FigNb = 1;
RatioMTsBkGd = [];
MaxAlongMTs3D = [];
MTsMaxWoBg_WholeCell = [];
Result = [];
%% Reading MT tracking results from a file
load(PathMTs);
%% Reading max cytoplasmic intensities from a file
fid = fopen(PathInputCytoplasmMax, 'r');
CytMaxima = textscan(fid, '%f');
fclose(fid);
%% Reading cell properties information from a file
fid = fopen(PathCellsInfo, 'r');
In_Cells = textscan(fid, '%s%f%f%f%f%f%f%f%f%f%f%f%f', 'headerLines', 1); 
CellNb = length(In_Cells{1});
Result_AverIntensOfMTmax = cell(CellNb, 1);  
Result_RatioOfMTmax = cell(CellNb, 1);  
%% Creating list of coordinates (x,y) corresponding to MTs tracked    
MTCoordList = f_MTsInCellCoordsList(UnDoubledMTs);    
for i_cell = 1:CellNb    % Loop on all the cells
    % Open the image where the cell is taken from
    a = In_Cells{1}(i_cell);
    FileName = [ImagePathBase a{1}];
    InitImage = load(FileName); 
    InitImage = InitImage.InitImage;      
    for i_pt = 1:length(MTCoordList{i_cell}(:,1))  % Loop on all MT points detected for that cell
        close all;  
        % Coordinates of MT point
        MT_Point = [MTCoordList{i_cell}(i_pt,1); MTCoordList{i_cell}(i_pt,2)];
        % Get info about cell number 'i_cell'
        CellLength = In_Cells{3}(i_cell);
        CellWidth = In_Cells{4}(i_cell) + MoreWidth;
        CellAxisAngle = In_Cells{5}(i_cell);
        CellCenter = [In_Cells{6}(i_cell), In_Cells{7}(i_cell)];
        CellsEnds_x = [In_Cells{8}(i_cell), In_Cells{9}(i_cell)];
        CellsEnds_y = [In_Cells{10}(i_cell), In_Cells{11}(i_cell)];    
        %% Take intens. profile across the cell, going through the point 'i_pt'
        % We are in the normal XY axis, but with the center in the
        % center of the cell        
        % Finding closest point on cell axis
        [X_CellAxis, Y_CellAxis, a] = improfile(InitImage, CellsEnds_x, CellsEnds_y);
        Dist = sqrt((X_CellAxis - MT_Point(1)) .^ 2 + (Y_CellAxis - MT_Point(2)) .^ 2);
        [a, i_min] = min(Dist);
        PlaneCenter1 = [X_CellAxis(i_min); Y_CellAxis(i_min)];        
        figure; imshow(InitImage, []);
%         line(CellsEnds_x, CellsEnds_y);       
        line(X_CellAxis, Y_CellAxis); 
%         line(MT_Point(1), MT_Point(2), 'Marker', 'o'); 
        line(PlaneCenter1(1), PlaneCenter1(2), 'Marker', 'o');                 
        Line1_x = [PlaneCenter1(1) - (CellWidth/2)*sind(CellAxisAngle), PlaneCenter1(1) + (CellWidth/2)*sind(CellAxisAngle)];
        Line1_y = [PlaneCenter1(2) + (CellWidth/2)*cosd(CellAxisAngle), PlaneCenter1(2) - (CellWidth/2)*cosd(CellAxisAngle)];
%         % Passage back to normal axes
%         Line1_x = Line1_x + CellCenter(1);
%         Line1_y = -Line1_y + CellCenter(2);               
        figure; imshow(InitImage, []);
        line(Line1_x, Line1_y);        % To check the position of the line                  
        line(PlaneCenter1(1), PlaneCenter1(2), 'Marker', 'o');        % To check the position of the plane center                 
        DistIntens = improfile(InitImage, Line1_x, Line1_y, CellWidth, 'bicubic'); 
        DistIntens = [(1:length(DistIntens))', DistIntens];
        %% First round of fitting (smoothing spline) of the cytoplasmic signal:
        %% to take off the points corresponding to MTs
    %     CutOffCurve = csaps(DistIntens(:,1), DistIntens(:,2), SmSplineParam1, DistIntens(:,1)) + CutOffPlus;    
        FitOptions = fitoptions('Method', 'SmoothingSpline', 'SmoothingParam', SmSplineParam1);
        FitType = fittype('smoothingspline');
        cfun = fit(DistIntens(:,1), DistIntens(:,2), FitType, FitOptions);
        CutOffCurve = cfun(DistIntens(:,1)) + CutOffPlus;
        figure;  %(FigNb); FigNb = FigNb + 1;
        hold on
        plot(DistIntens(:,1), DistIntens(:,2), 's', 'MarkerSize', 4);
        plot(DistIntens(:,1), CutOffCurve, '-r');
        title('Intensity profile with cut-off curve');
        xlabel('Position, in pixels');
        ylabel('Intensity');   
        i = 1;
        j = 1;
        Max = length(DistIntens(:,1));
        DistIntensWoMT = DistIntens;
        while i < Max + 1                   % Loop through the points of the intensity profile considered
            if (DistIntensWoMT(i,2) > CutOffCurve(j))            
                DistIntensWoMT(i,:) = [];           % Construction of an array with dots under the fit (taking MTs signal away)           
                i = i - 1;  
                Max = Max - 1;  
            end
            i = i + 1;  
            j = j + 1;         
        end
        %% Second round of fitting (smoothing spline) of the cytoplasmic signal:     
    %     [BkGd] = csaps(DistIntensWoMT(:,1), DistIntensWoMT(:,2), SmSplineParam2, DistIntens(:,1));
        FitOptions = fitoptions('Method', 'SmoothingSpline', 'SmoothingParam', SmSplineParam2);
        FitType = fittype('smoothingspline');
        cfun = fit(DistIntensWoMT(:,1), DistIntensWoMT(:,2), FitType, FitOptions);
        BkGd = cfun(DistIntens(:,1));

        figure;  %(FigNb); FigNb = FigNb + 1;  
        hold on;
        plot(1:length(DistIntens(:,1)), DistIntens(:,2), 's', 'MarkerSize', 4);    
        plot(1:length(DistIntens(:,1)), BkGd, '-r');
        hold off;
         %plot(DistIntens(:,1), BkGd, DistIntens(:,1), DistIntens(:,2), 's', 'MarkerSize', 4);         
        %% Taking off the background from the values obtained from
        %% average projection images at the positions for MT maxima determined before
        % Take intensity value not from MT maximum point, but from a square 3*3
        % around the point           
        [a, Pos_MaxBkGd] = max(BkGd);
        MaxWoBkGd = (DistIntens(Pos_MaxBkGd, 2) - BkGd(Pos_MaxBkGd));
        Result_AverIntensOfMTmax{CellNb} = [Result_AverIntensOfMTmax{CellNb}; MaxWoBkGd];
        % The same for ratios
        BasicBkGdMax = CytMaxima{1,1}(CellNb);
        RatioMTsBkGd = MaxWoBkGd / BasicBkGdMax;
        Result_RatioOfMTmax{CellNb} = [Result_RatioOfMTmax{CellNb}; RatioMTsBkGd];         
    end
end
save('_OutputAI/_ResAllIntensAllCells.mat', 'Result_AverIntensOfMTmax');  
save('_OutputAI/_ResAllRatiosAllCells.mat', 'Result_RatioOfMTmax');