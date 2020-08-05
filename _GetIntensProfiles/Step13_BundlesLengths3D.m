%% The task of the program is to select from all potential bundles paths the
%% ones that are the most likely to be ones 
%% Then the paths are smoothed and bundle lengths are extracted
clear;
close all;
%% Parameters
% Minimal number of points to have a bundle
MinPts = 4;
% Parameter to say 'bundle' gets curved too much to be a single
% bundle
dP_MT_Kick = 0.17;      %5;       % 2.1;      
PtsMTBending = 7;    %7   % number of points for linear fit for looking at how smooth is a MT
% Maximal distance at which bundles ends / beginnings are
% considered close
CloseDist = 10;  
% Maximal distance from the other bundle_to_delete end (opposite the close one)
% to the bundle that is going to stay
CloseDist_OtherEnd = 17;
% Maximal distance to take away overlapping part of bundles
CloseDist_Overlap = 6;
% Parameter for smoothing spline to smooth MT bundles
p_SmSpline = 0.2;      
% The amount of micrometers in one pixel of the camera and settings used
PixelTo_um = 0.13;   % RS bin2            % 0.065;   RS bin1
Px = 7.7; % 7.7 comes from RS having x-y resolution 0.065um for 100x and slice spacing being 0.5um
% Ratio between Z-step in stacks (in microns) and resolution in XY- plane 
RatioZ_ResXY = 0.25 / PixelTo_um;
%% 
UnDoubledMTs = {};
AllMTsSmoothed = {};
MTsLength = {};
MTsLength_um = {};
CellNb = 0;
load('_OutputGI/_BundlesPoints.mat');
% AllMTs{*} corresponds to an image
% a = AllMTs{*}(*) corresponds to a cell on this image
% a{*} corresponds to the list of MTs in this cell
Len_AM = length(AllMTs);
for i_Im = 1:Len_AM         % Loop on all images analysed
    MTs_WholeImage = AllMTs{i_Im}; 
    for i_Cell = 1:length(MTs_WholeImage)
        pause(0.5);
        close all;             
        MTs = MTs_WholeImage{i_Cell};
        if isempty(MTs)     % No bundles were found in that cell
            continue
        end
        CellNb = CellNb + 1;                
        % Take off [] at the end of the array containing potential bundles
        i_Off = length(MTs); 
        while 1
            if isempty(MTs{i_Off})
                i_Off = i_Off - 1;
                MTs = MTs(1:i_Off);
            else
                break
            end
        end
        %% Visualisation in 3 D
        h = figure(); plot3(MTs{1,1}(:, 2), MTs{1,1}(:, 1), MTs{1,1}(:, 3), '-o', 'MarkerSize', 3);
        grid on; xlabel('X'); ylabel('Y'); zlabel('Z'); 
  
        for i = 2:length(MTs) 
            hold all;  
            plot3(MTs{1,i}(:, 2), MTs{1,i}(:, 1), MTs{1,i}(:, 3), '-o', 'MarkerSize', 3);
        end          
        %% Saving initial 3D bundles representation  
        saveas(h, ['_OutputGI/MTs3DRaw_' int2str(CellNb) '.fig']);
        %% Taking off the bundle arrays that are too short to be a bundle
        i = 1; 
        while i < length(MTs)  
            if length(MTs{1,i}(:, 1)) < MinPts
                % Replace the too short bundle with the bundle from the end of the list
                MTs{1, i} = MTs{1, length(MTs)};
                MTs = MTs(1:length(MTs) - 1);
            end 
            % Check if after replacement length of the MT is big enough
            % (to avoid situations where short MTs are copied to this position)
            if length(MTs{1,i}(:, 1)) < MinPts
                i = i - 1;
            end
            i = i + 1;
        end
        %% Visualisation in 3 D
%         figure; plot3(MTs{1,1}(:, 2), MTs{1,1}(:, 1), MTs{1,1}(:, 3), '-o', 'MarkerSize', 3);
%         grid on; xlabel('X'); ylabel('Y'); zlabel('Z'); 
% 
%         for i = 2:length(MTs)  
%             hold all; 
%             plot3(MTs{1,i}(:, 2), MTs{1,i}(:, 1), MTs{1,i}(:, 3), '-o', 'MarkerSize', 3);
%         end
        %% Cutting in two bundle arrays that are too bent in Oxy to be a single bundle
        for i_MT = 1:length(MTs)                                         
            if length(MTs{1, i_MT}(:, 1)) < PtsMTBending + 1
                continue 
            end            
%             figure; plot3(MTs{1, i_MT}(:, 2), MTs{1, i_MT}(:, 1), MTs{1, i_MT}(:, 3), '-o', 'MarkerSize', 3); grid on;
%             hold on;                 
            x = MTs{1, i_MT}(1:PtsMTBending, 2);
            y = MTs{1, i_MT}(1:PtsMTBending, 1);                                          
            if (max(x) - min(x)) > (max(y) - min(y))
                p = polyfit(x, y, 1);  
                hold on
%                 plot3(x, p(1) * x + p(2), 8 * ones(PtsMTBending), '-rs', 'MarkerSize', 5);  
            else
                p = polyfit(y, x, 1);
                hold on
%                 plot3(p(1) * y + p(2), y, 8 * ones(PtsMTBending), '-gs', 'MarkerSize', 5);
            end              
            OldP = p(1);                                                                                                            
            for i = 2:length(MTs{1, i_MT}(:, 1)) - PtsMTBending + 1
                x = MTs{1, i_MT}(i:i + PtsMTBending - 1, 2);
                y = MTs{1, i_MT}(i:i + PtsMTBending - 1, 1);                    
                if (max(x) - min(x)) > (max(y) - min(y))
                    p = polyfit(x, y, 1);  
                    hold on
%                     plot3(x, p(1) * x + p(2), 8 * ones(PtsMTBending), '-rs', 'MarkerSize', 5);  
                else
                    p = polyfit(y, x, 1);
                    hold on
%                     plot3(p(1) * y + p(2), y, 8 * ones(PtsMTBending), '-gs', 'MarkerSize', 5);
                end                                
                if abs(p(1) - OldP) > dP_MT_Kick
                    % Copying last part of the bundle (after the strong curvature) to the end of the 'MTs' array
                    MTLenCurr = length(MTs);
                    MTs{1, MTLenCurr + 1} = MTs{1, i_MT}(i + floor(PtsMTBending/2) + 2:length(MTs{1, i_MT}(:, 1)), :);
                    MTs{1, i_MT}(i + floor(PtsMTBending/2) + 2:length(MTs{1, i_MT}(:, 1)), :) = [];
                    break
                end
                OldP = p(1);
            end 
        end            
        %% Take off empty bits (starting from the end of 'MTs')
        if length(MTs{1, length(MTs)}(:,1)) == 1    
            MTs = MTs(1:length(MTs) - 1);
        end
        i = length(MTs);
        while i > 0 
            if length(MTs{1,i}(:,1)) == 1
                MTs{1, i} = MTs{1, length(MTs)};
                MTs = MTs(1:length(MTs) - 1);
            end
            i = i - 1;
        end
        %% Visualisation in 3 D
        figure; plot3(MTs{1,1}(:, 2), MTs{1,1}(:, 1), MTs{1,1}(:, 3), '-o', 'MarkerSize', 3);
        grid on; xlabel('X'); ylabel('Y'); zlabel('Z'); 

        for i = 2:length(MTs) 
            hold all; 
            plot3(MTs{1,i}(:, 2), MTs{1,i}(:, 1), MTs{1,i}(:, 3), '-o', 'MarkerSize', 3);
        end
        %% Constructing an array containing distances between the beginings of
        %% potential bundles: one column corresponds to distances from this start to other starts
        DistStarts = zeros(length(MTs), length(MTs));
        for i = 1:length(MTs)
            for j = 1:length(MTs)  
                DistStarts(j, i) = sqrt((MTs{1,i}(1, 2) - MTs{1,j}(1, 2)).^2 + (MTs{1,i}(1, 1) - MTs{1,j}(1, 1)).^2 + ((MTs{1,i}(1, 3) - MTs{1,j}(1, 3)) * Px) .^2);                       
            end 
        end
        %% Taking off MT arrays that start at close locations
        MTsOld = MTs;       % to be able to compare with bundles that are already [0 0 0]
        for i_Col = 2:length(MTs)
            for i_Lin = 1:i_Col - 1
                if DistStarts(i_Lin, i_Col) < CloseDist           
                    % Take away shorter one
                    if length(MTsOld{1, i_Col}(:, 1)) > length(MTsOld{1, i_Lin}(:, 1))                
                        MTs{1, i_Lin} = [0 0 0];
                    else
                        MTs{1, i_Col} = [0 0 0];
                    end
                end
            end
        end
        %% Take off empty bits (starting from the end of 'MTs')
        if length(MTs{1, length(MTs)}(:,1)) == 1    
            MTs = MTs(1:length(MTs) - 1);
        end
        i = length(MTs);
        while i > 0 
            if length(MTs{1,i}(:,1)) == 1
                MTs{1, i} = MTs{1, length(MTs)};
                MTs = MTs(1:length(MTs) - 1);
            end
            i = i - 1;
        end
        %% Visualisation in 3 D
        figure; plot3(MTs{1,1}(:, 2), MTs{1,1}(:, 1), MTs{1,1}(:, 3), '-o', 'MarkerSize', 3);
        grid on; xlabel('X'); ylabel('Y'); zlabel('Z'); 

        for i = 2:length(MTs) 
            hold all; 
            plot3(MTs{1,i}(:, 2), MTs{1,i}(:, 1), MTs{1,i}(:, 3), '-o', 'MarkerSize', 3);
        end 
        %% Constructing an array containing distances between the ends of
        %% potential bundles: one column corresponds to distances from this start to other starts
        L = [];   
        for i = 1:length(MTs) 
            L(i, 1) = length(MTs{1,i}(:, 2));        % Number of points in each 'bundle'
        end        
        DistEnds = zeros(length(MTs), length(MTs));
        for i = 1:length(MTs)
            for j = 1:length(MTs)  
                DistEnds(j, i) = sqrt((MTs{1,i}(L(i, 1), 2) - MTs{1,j}(L(j, 1), 2)).^2 + (MTs{1,i}(L(i, 1), 1) - MTs{1,j}(L(j, 1), 1)).^2 + ((MTs{1,i}(L(i, 1), 3) - MTs{1,j}(L(j, 1), 3)) * Px) .^2);                       
            end 
        end
        %% Taking off MT arrays that end at close locations
        MTsOld = MTs;       % to be able to compare with bundles that are already [0 0 0]
        for i_Col = 2:length(MTs)
            for i_Lin = 1:i_Col - 1
                if DistEnds(i_Lin, i_Col) < CloseDist                                     
                    if length(MTsOld{1, i_Col}(:, 1)) > length(MTsOld{1, i_Lin}(:, 1))   % Finding shorter one  
                        NumOff = i_Lin;
                        NumLeft = i_Col;
                        Point = MTsOld{1, i_Lin}(1, :);     % Bundle start point, opposite the one that is close to the other MT
                        HalfMT_Dist = sqrt((Point(2) - MTs{1, i_Col}(:, 2)).^2 + (Point(1) - MTs{1, i_Col}(:, 1)).^2 + ((Point(3) - MTs{1, i_Col}(:, 3)) * 7.7) .^2);
                    else
                        NumOff = i_Col;
                        NumLeft = i_Lin;                        
                        Point = MTsOld{1, i_Col}(1, :);     % Bundle start point, opposite the one that is close to the other MT                          
                        HalfMT_Dist = sqrt((Point(2) - MTs{1, i_Lin}(:, 2)).^2 + (Point(1) - MTs{1, i_Lin}(:, 1)).^2 + ((Point(3) - MTs{1, i_Lin}(:, 3)) * 7.7).^2);
                    end
                    % Looking if the beginning of the shorter MT lies close to any region
                    % of the long one
                    if min(min(HalfMT_Dist)) < CloseDist_OtherEnd                                     
                        MTs{1, NumOff} = [0 0 0];    
                    else      % take off part of MT that goes closely along other MT 
%                             % start taking off from the end we know is
%                             % close in the two MTs                        
                        if (length(MTs{1, i_Col}) == 3) | (length(MTs{1, i_Lin}) == 3)      % then it is probably [0 0 0]
                            continue
                        end                        
                        x_PartOff_Col = length(MTs{1, i_Col}) - length(MTs{1, NumOff}) + 1 : length(MTs{1, i_Col});
                        x_PartOff_Lin = length(MTs{1, i_Lin}) - length(MTs{1, NumOff}) + 1 : length(MTs{1, i_Lin});
                        
                        Dist_CloseParts = sqrt((MTs{1, i_Col}(x_PartOff_Col, 2) - ...
                            MTs{1, i_Lin}(x_PartOff_Lin, 2)).^2 + (MTs{1, i_Col}(x_PartOff_Col, 1) - ...
                            MTs{1, i_Lin}(x_PartOff_Lin, 1)).^2 + ((MTs{1, i_Col}(x_PartOff_Col, 3) - ...
                            MTs{1, i_Lin}(x_PartOff_Lin, 3)) * 7.7).^2);
                        % Take off the shorter bundle the part that overlaps with the long one
                        MTs{1, NumOff}(find(Dist_CloseParts < CloseDist_Overlap), :) = [];                        
                    end                    
                end
            end
        end
        %% Take off empty bits (starting from the end of 'MTs')
        if length(MTs{1, length(MTs)}(:,1)) == 1    
            MTs = MTs(1:length(MTs) - 1);
        end
        i = length(MTs);
        while i > 0 
            if length(MTs{1,i}(:,1)) == 1
                MTs{1, i} = MTs{1, length(MTs)};
                MTs = MTs(1:length(MTs) - 1);
            end
            i = i - 1;
        end
        %% Visualisation in 3 D
        figure; plot3(MTs{1,1}(:, 2), MTs{1,1}(:, 1), MTs{1,1}(:, 3), '-o', 'MarkerSize', 3);
        grid on; xlabel('X'); ylabel('Y'); zlabel('Z'); 
        for i = 2:length(MTs) 
            hold all; 
            plot3(MTs{1,i}(:, 2), MTs{1,i}(:, 1), MTs{1,i}(:, 3), '-o', 'MarkerSize', 3);
        end 
        %% Constructing an array containing distances between the beginings and ends of
        %% potential bundles: one column corresponds to distances from this start to other starts
        DistStEnds = zeros(length(MTs), length(MTs));
        L = [];   
        for i = 1:length(MTs) 
            L(i, 1) = length(MTs{1,i}(:, 2));        % Number of points in each 'bundle'
        end 
        for i = 1:length(MTs)       % i- current bundle
            for j = 1:length(MTs)   % j- bundle it is compared to
                DistStEnds(j, i) = sqrt((MTs{1,i}(1, 2) - MTs{1,j}(L(j, 1), 2)).^2 + ...
                    (MTs{1,i}(1, 1) - MTs{1,j}(L(j, 1), 1)).^2 + ((MTs{1,i}(1, 3) - MTs{1,j}(L(j, 1), 3)) * RatioZ_ResXY) .^2);       
                % 7.7 comes from RS having x-y resolution 0.065um for 100x and slice spacing being 0.5um
            end
        end
        %% Taking off MT arrays that start and end at close locations
        DistMT_End = [];
        MTsOld = MTs;       % to be able to compare with bundles that are already [0 0 0]
        for i_Col = 1:length(MTs)
            for i_Lin = 1:length(MTs)
                if DistStEnds(i_Lin, i_Col) < CloseDist           
                    % Finding shorter one
                    L_Col = length(MTsOld{1, i_Col}(:, 1));
                    L_Lin = length(MTsOld{1, i_Lin}(:, 1));
                    if L_Col > L_Lin
                        NumOff = i_Lin;
                        NumLeft = i_Col;
                        Point = MTsOld{1, i_Lin}(1, :);     % The other end point, opposite the one that is close to the other MT
                        HalfMT_Dist = sqrt((Point(2) - MTs{1, i_Col}(:, 2)).^2 + (Point(1) - MTs{1, i_Col}(:, 1)).^2 + ((Point(3) - MTs{1, i_Col}(:, 3)) * 7.7).^2);
                    else
                        NumOff = i_Col;
                        NumLeft = i_Lin;                        
                        Point = MTsOld{1, i_Col}(L_Col, :);     % The other end point, opposite the one that is close to the other MT                          
                        HalfMT_Dist = sqrt((Point(2) - MTs{1, i_Lin}(:, 2)).^2 + (Point(1) - MTs{1, i_Lin}(:, 1)).^2 + ((Point(3) - MTs{1, i_Lin}(:, 3)) * 7.7).^2);
                    end
                    % Looking if the other end of the shorter MT lies close to any region
                    % of the long one
                    if min(min(HalfMT_Dist)) < CloseDist_OtherEnd                                     
                        MTs{1, NumOff} = [0 0 0];
%                     else                            
                    end
                end
            end
        end
        %% Take off empty bits
        if length(MTs{1, length(MTs)}(:,1)) == 1    
            MTs = MTs(1:length(MTs) - 1);
        end
        i = length(MTs);
        while i > 0 
            if length(MTs{1,i}(:,1)) == 1
                MTs{1, i} = MTs{1, length(MTs)};
                MTs = MTs(1:length(MTs) - 1); 
            end
            i = i - 1;
        end 
        % Check if the there is still some bundles left after the cleaning steps
        if isempty(MTs)
            AllMTsSmoothed = [AllMTsSmoothed; {0}];     % Written in the order of cells, not images with cells inside any more
            MTsLength = [MTsLength; {0}];
            MTsLength_um = [MTsLength_um; {0}];
            continue
        end
        %% Visualisation in 3 D
        h = figure; plot3(MTs{1,1}(:, 2), MTs{1,1}(:, 1), MTs{1,1}(:, 3), '-o', 'MarkerSize', 3);
        grid on; xlabel('X'); ylabel('Y'); zlabel('Z'); 

        for i = 2:length(MTs) 
            hold all; 
            plot3(MTs{1,i}(:, 2), MTs{1,i}(:, 1), MTs{1,i}(:, 3), '-o', 'MarkerSize', 3);
        end
        UnDoubledMTs = [UnDoubledMTs; {MTs}];     % Written in the order of cells, not images with cells inside any more
        %% Saving 3D bundles representation after single-making procedure 
        saveas(h, ['_OutputGI/MTs3D_ForLength_' int2str(CellNb) '.fig']);        
        %% Smoothing final version of bundles (only their z- coordinate) with smoothing spline        
        MTsSmoothed = MTs;
        MTsLength_Elem = [];
        for i = 1:length(MTs)
            ToSmooth = MTs{1, i}(:, 3);                                    
            x = 1:length(ToSmooth);            
            FitOptions = fitoptions('Method', 'SmoothingSpline', 'SmoothingParam', p_SmSpline);
            FitType = fittype('smoothingspline');
            cfun = fit(x', ToSmooth, FitType, FitOptions);            
            figure();
            plot(x, ToSmooth, '-bo'); 
            hold on;
            plot(x, cfun(x), '-r*');                 
            MTsSmoothed{1, i}(:, 3) = cfun(x);
            %% Measuring the length of the smoothed bundle
            NbPts = length(MTsSmoothed{1, i}(:, 2));
            MTsLength_Elem(i) = sum(sum(sqrt((MTsSmoothed{1, i}(1:NbPts-1, 2) - MTsSmoothed{1, i}(2:NbPts, 2)).^2 + ...
                (MTsSmoothed{1, i}(1:NbPts-1, 1) - MTsSmoothed{1, i}(2:NbPts, 1)).^2 + ...
                ((MTsSmoothed{1,i}(1:NbPts-1, 3) - MTsSmoothed{1, i}(2:NbPts, 3)) * RatioZ_ResXY) .^2)));                   
        end  
        AllMTsSmoothed = [AllMTsSmoothed; {MTsSmoothed}];     % Written in the order of cells, not images with cells inside any more
        MTsLength = [MTsLength; {MTsLength_Elem}];     
        MTsLength_um = [MTsLength_um; {MTsLength_Elem * PixelTo_um}]; 
    end
end 
save('_OutputGI/UnDoubledMTs.mat', 'UnDoubledMTs');
save('_OutputGI/SmoothedMTs.mat', 'AllMTsSmoothed');
save('_OutputGI/_MTsLengthsInPx.mat', 'MTsLength');
save('_OutputGI/_MTsLengths_um.mat', 'MTsLength_um');



