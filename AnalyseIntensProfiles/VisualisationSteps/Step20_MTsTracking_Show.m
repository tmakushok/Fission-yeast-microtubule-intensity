clear;
%--------------------------------------------------------------------------
%!!!--!!! Minimal number of points that constitute a MT
MinMTLength = 8;    %6;       
%!!!--!!! Parts of (presumably) MTs at this distance or smaller 
%-- will be connected into one MT
% 4.25 is bit more than max distance for 3px separation        
% 2.83 is a bit more than max distance for 2px separation
SmallDistance = 2.83; 
%!!!--!!! 'Rigidity' parameter for smoothing spline (for intensities along each MT)
%p_SingleMTFit = 0.3;
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
p_AllMTFit = 0.05; 
p_ResFit = 0.05;  % 0.05 for ratios
%!!!--!!! Coefficient by which the number of bins is devided
NbBinsCoeff = 30;    %1; for ratios
%!!!--!!! Putting that many points to 0 in the each cell values histogram 
%-- (to assure that fitting gets down a bit at 0)
%HistZeroPoints = 6;
%!!!--!!! Basic name of the file where 'Connected' array is stored
% For MT max values use 'ConnectedAver'
ValueBaseName = '_OutputImagesAndVariables\ConnectedAver';
ResOutput = '_OutputImagesAndVariables/_Result_MTmax';
% For ratio MT/cytoplasm values use 'ConnectedAverRatio'
%ValueBaseName = '_OutputImagesAndVariables/ConnectedAverRatio';
%--------------------------------------------------------------------------
Result = [];
MTs = [];                   % For the line 'clear MTs' to work
TotalMTLengthAllCells = [];
TotalMTLength = 0;
FileNb = 0;                 % The number of the first file - 1, for example, 0 for 'Image1.tif'
FigNb = 1;
while 1              % Loop on the 'Connected' matrixes to analyse  
    close all;    
    FileNb = FileNb + 1;
    str = sprintf('%0.5g', FileNb);     % Converts number to string
    FileName = strcat(ValueBaseName, str, '.mat');   
    % To check if the file with this name does exist
    fid = fopen(FileName, 'r');     
    if fid == -1        
        break;
    end            
    load(FileName);   
    Connected = ConnectedAver;
    clear ConnectedAver;
     
    [Lin, Col] = size(Connected);
%% Adding of a square of zeros all around the matrix (to be able to analyse border elements)
    Connected = [zeros(1, Col + 2); zeros(Lin, 1), Connected, zeros(Lin, 1); zeros(1, Col + 2)];
    Lin = Lin + 1;
    Col = Col + 1;   
%     figure(FigNb); FigNb = FigNb + 1;
%     Handle = waterfall(Connected);         % MaxAlongMTs3D(:,3), MaxAlongMTs3D(:,1), Z_3D);
%     set(Handle, 'Marker', 's', 'MarkerSize', 4)        
%     %grid on
%     title('Connected');        
%     xlabel('Position of MT max along the cell');
%     ylabel('Position of MT max across the cell');   
%     zlabel('MT maximum value');       
%     rotate3d on
%     rotate3d off
    
%% Filling the cell array 'MTs' with data
    MTs = [];      % Taking away previous 'MTs' from the memory, not to have weird values
    Analysed = zeros(Lin, Col);    
    NewMT_Flag = 1;             % Shows whether the next point will belong to the next MT
    MT_Number = 0;              % Number of the current element in cell array
    for i_Track = 2:(Lin - 1)
        for j_Track = 2:(Col - 1)
            Value = Connected(i_Track, j_Track);
            if (Value ~= 0) & (Analysed(i_Track, j_Track) == 0)     % and not analysed already
                if NewMT_Flag == 1      % Time to start new MT
                    NewMT_Flag = 0;
                    MT_Number = MT_Number + 1;
                    i_InsideMT = 1;             % Number of the current element to add inside current MT
                    
                    MTs{MT_Number}(i_InsideMT, :) = [i_Track, j_Track, Value];
                    Analysed(i_Track, j_Track) = 1;
                    i_InsideMT = i_InsideMT + 1;
                    % Searching for the next non-zero connected element
                    % !!! For the moment, branching is not taken into account
                    while 1
                        if (Connected(i_Track - 1, j_Track - 1) ~= 0) & (Analysed(i_Track - 1, j_Track - 1) == 0)
                            i_Track = i_Track - 1;
                            j_Track = j_Track - 1;
                            MTs{MT_Number}(i_InsideMT, :) = [i_Track, j_Track, Connected(i_Track, j_Track)];
                            Analysed(i_Track, j_Track) = 1;
                            i_InsideMT = i_InsideMT + 1;
                            continue
                        end
                        if (Connected(i_Track - 1, j_Track) ~= 0) & (Analysed(i_Track - 1, j_Track) == 0)
                            i_Track = i_Track - 1;                            
                            MTs{MT_Number}(i_InsideMT, :) = [i_Track, j_Track, Connected(i_Track, j_Track)];
                            Analysed(i_Track, j_Track) = 1;
                            i_InsideMT = i_InsideMT + 1;
                            continue
                        end                    
                        if (Connected(i_Track - 1, j_Track + 1) ~= 0) & (Analysed(i_Track - 1, j_Track + 1) == 0)
                            i_Track = i_Track - 1;
                            j_Track = j_Track + 1;
                            MTs{MT_Number}(i_InsideMT, :) = [i_Track, j_Track, Connected(i_Track, j_Track)];
                            Analysed(i_Track, j_Track) = 1;
                            i_InsideMT = i_InsideMT + 1;
                            continue                        
                        end
                        if (Connected(i_Track, j_Track + 1) ~= 0) & (Analysed(i_Track, j_Track + 1) == 0)                           
                            j_Track = j_Track + 1;
                            MTs{MT_Number}(i_InsideMT, :) = [i_Track, j_Track, Connected(i_Track, j_Track)];
                            Analysed(i_Track, j_Track) = 1;
                            i_InsideMT = i_InsideMT + 1;
                            continue                     
                        end
                        if (Connected(i_Track + 1, j_Track + 1) ~= 0) & (Analysed(i_Track + 1, j_Track + 1) == 0)
                            i_Track = i_Track + 1;
                            j_Track = j_Track + 1;
                            MTs{MT_Number}(i_InsideMT, :) = [i_Track, j_Track, Connected(i_Track, j_Track)];
                            Analysed(i_Track, j_Track) = 1;
                            i_InsideMT = i_InsideMT + 1;
                            continue                     
                        end
                        if (Connected(i_Track + 1, j_Track) ~= 0) & (Analysed(i_Track + 1, j_Track) == 0)
                            i_Track = i_Track + 1;                            
                            MTs{MT_Number}(i_InsideMT, :) = [i_Track, j_Track, Connected(i_Track, j_Track)];
                            Analysed(i_Track, j_Track) = 1;
                            i_InsideMT = i_InsideMT + 1;
                            continue                     
                        end
                        if (Connected(i_Track + 1, j_Track - 1) ~= 0) & (Analysed(i_Track + 1, j_Track - 1) == 0)
                            i_Track = i_Track + 1;
                            j_Track = j_Track - 1;
                            MTs{MT_Number}(i_InsideMT, :) = [i_Track, j_Track, Connected(i_Track, j_Track)];
                            Analysed(i_Track, j_Track) = 1;
                            i_InsideMT = i_InsideMT + 1;
                            continue                      
                        end
                        if (Connected(i_Track, j_Track - 1) ~= 0) & (Analysed(i_Track, j_Track - 1) == 0)                            
                            j_Track = j_Track - 1;
                            MTs{MT_Number}(i_InsideMT, :) = [i_Track, j_Track, Connected(i_Track, j_Track)];
                            Analysed(i_Track, j_Track) = 1;
                            i_InsideMT = i_InsideMT + 1;
                            continue                    
                        end    
                        NewMT_Flag = 1;                        
                        break
                    end                         
                end 
            end 
        end         
    end    

%----- 3D visualisation of 'short versions' of MT 
%     for MT_Number = 1:length(MTs)
%         if length(MTs{MT_Number}(:, 1)) >= MinMTLength
%             OneMT3D = zeros(Lin, Col);
%             for i = 1:length(MTs{MT_Number}(:, 1))
%                 OneMT3D(MTs{MT_Number}(i, 1), MTs{MT_Number}(i, 2)) = MTs{MT_Number}(i, 3); 
% %                 if MTs{MT_Number}(i, 3) ~= -1
% %                     AllMTsTogether = [AllMTsTogether; MTs{MT_Number}(i, 3)];
% %                 end
%             end

%             figure(FigNb); FigNb = FigNb + 1;
%             Handle = waterfall(OneMT3D);        
%             set(Handle, 'Marker', 's', 'MarkerSize', 4)        
%             %grid on
%             title('One MT tracked');        
%             xlabel('Position of MT max along the cell');
%             ylabel('Position of MT max across the cell');   
%             zlabel('MT maximum value');       
%             rotate3d on
%             rotate3d off
%         end
%     end
%% Measuring distances between all starts and ends of all tracked MTs.
%-- If distance is less than the parameter 'SmallDistance', regions 
%-- are connected to form a single MT that is longer.
%-- For the case where MT end is close to another MT's beginning  
    if (length(MTs) == 0)
        Result = [Result; 0]; 
        FigNb = FigNb + 1;
        continue
    end
    for i_LongerMT = 1:length(MTs)
        for j_LongerMT = 1:length(MTs)
%             if (length(MTs{i_LongerMT}) == 0) | (length(MTs{j_LongerMT}) == 0)
%                 continue;
%             end 
            if (MTs{j_LongerMT}(1,3) ~= -1) & (MTs{i_LongerMT}(1,3) ~= -1) & (j_LongerMT ~= i_LongerMT)         % Without this check it will measure distance between beginning and end of the same MT
                [LastJ, a] = size(MTs{j_LongerMT}(:,1));
                [LastI, a] = size(MTs{i_LongerMT}(:,1));
                DistFirstLastEnd = sqrt( (MTs{i_LongerMT}(1, 1) - MTs{j_LongerMT}(LastJ, 1))^2 + (MTs{i_LongerMT}(1, 2) - MTs{j_LongerMT}(LastJ, 2))^2 );
                DistFirstEndLast = sqrt( (MTs{i_LongerMT}(LastI, 1) - MTs{j_LongerMT}(1, 1))^2 + (MTs{i_LongerMT}(LastI, 2) - MTs{j_LongerMT}(1, 2))^2 );

                if DistFirstLastEnd < SmallDistance
                    MTs{i_LongerMT} = [MTs{j_LongerMT}; MTs{i_LongerMT}];                      
                    % Now we look at all possible previously located elements ends (going through all elements)
                    % when previous element found or if not found, we delete the element number 'j_LongerMT' that was last added in front                                       
                    i_LongRepeat = i_LongerMT;
                    j_LongRepOld = j_LongerMT;
                    j_LongRepeat = 1;
                    while j_LongRepeat <= length(MTs)       % With control of the index, will be repeated until there is no elements left to complement the chain 
                        if (MTs{j_LongerMT}(1,3) ~= -1) & (MTs{i_LongerMT}(1,3) ~= -1) & (j_LongRepeat ~= i_LongRepeat) & (j_LongRepeat ~= j_LongRepOld)         % Without this check it will measure distance between beginning and end of the same piece (big chain and last piece added)                               
                            [LastJRep, a] = size(MTs{j_LongRepeat}(:,1));
                            DistFirstLastEnd = sqrt( (MTs{i_LongRepeat}(1, 1) - MTs{j_LongRepeat}(LastJRep, 1))^2 + (MTs{i_LongRepeat}(1, 2) - MTs{j_LongRepeat}(LastJRep, 2))^2 );                            

                            if DistFirstLastEnd < SmallDistance % Accumulation from last to first continued
                                MTs{i_LongRepeat} = [MTs{j_LongRepeat}; MTs{i_LongRepeat}];                                  
                                MTs{j_LongRepOld} = [1 1 -1];      % This gives artificial pics in (1,1). With max MT value ~= -1 we can get rid of those            
                                j_LongRepOld = j_LongRepeat;
                                j_LongRepeat = 0;       % With the new piece added we should compare the new starting point of the assembled chain with all of the end points of pieces available                                                                                            
                            end                                
                        end
                        j_LongRepeat = j_LongRepeat + 1;
                    end
                    MTs{j_LongRepOld} = [1 1 -1];                                                                               
                end                     

                if DistFirstEndLast < SmallDistance % Accumulation from first to last
                    MTs{i_LongerMT} = [MTs{i_LongerMT}; MTs{j_LongerMT}];                       
                    % Now we look at all possible next elements starts (going through all elements)
                    % when next element found or if not found, we delete the element number 'j_LongerMT'                                        
                    i_LongRepeat = i_LongerMT;
                    j_LongRepOld = j_LongerMT;
                    j_LongRepeat = 1;
                    while j_LongRepeat <= length(MTs)       % With control of the index, will be repeated until there is no elements left to complement the chain 
                        if (MTs{j_LongerMT}(1,3) ~= -1) & (MTs{i_LongerMT}(1,3) ~= -1) & (j_LongRepeat ~= i_LongRepeat) & (j_LongRepeat ~= j_LongRepOld)         % Without this check it will measure distance between beginning and end of the same piece (big chain and last piece added)                               
                            [LastIRep, a] = size(MTs{i_LongRepeat}(:,1));
                            DistFirstEndLast = sqrt( (MTs{i_LongRepeat}(LastIRep, 1) - MTs{j_LongRepeat}(1, 1))^2 + (MTs{i_LongRepeat}(LastIRep, 2) - MTs{j_LongRepeat}(1, 2))^2 );

                            if DistFirstEndLast < SmallDistance % Accumulation from first to last continued
                                MTs{i_LongRepeat} = [MTs{i_LongRepeat}; MTs{j_LongRepeat}];                                                                                                              
                                MTs{j_LongRepOld} = [1 1 -1];      % This gives artificial pics in (1,1). With max MT value ~= -1 we can get rid of those            
                                j_LongRepOld = j_LongRepeat;
                                j_LongRepeat = 0;       % With the new piece added we should compare the new end of the assembled chain with all of the pieces available                                                                                            
                            end                                
                        end
                        j_LongRepeat = j_LongRepeat + 1;
                    end
                    MTs{j_LongRepOld} = [1 1 -1];                      
                end                     
            end
        end
    end 
%% Taking away the MTs in the 'MTs' cell array that were replaced with [1 1 -1]
%% and of those whose length is smaller than 'MinMTLength'
%% and of those that are strictly on the same horisontal line close to the border
%% (there are chances that it is an artifact of MT detection)
    i_Long = 1;
    for i = 1:length(MTs)   
        if (MTs{i}(1,3) ~= -1) & (length(MTs{i}(:,1)) >= MinMTLength) & (min(MTs{i}(:,1)) ~= simple_max(MTs{i}(:,1)))   
            MTs{i_Long} = MTs{i};       % Putting 'good' MTs more to the beginning
            i_Long = i_Long + 1;
        end        
    end
    for i = i_Long:length(MTs)               % Taking away [1 1 -1]s left at the end of 'MTs' cell array
        MTs{i} = [];               
    end            
    AllMTsTogether = [];        % An array containing all MT max positions and values
    for MT_Number = 1:length(MTs)
%% 3D visualisation of 'long versions' of MT (some close regions are connected)
        if length(MTs{MT_Number}) == 0
            break
        end
        % Array containing the values for all MTs inside a cell
        AllMTsTogether = [AllMTsTogether; MTs{MT_Number}(:, 3)];
        % Value representing total length of all MT inside one cell
        for i_InsideMT = 2:length(MTs{1, MT_Number})
            x1 = MTs{1, MT_Number}(i_InsideMT - 1, 1);
            y1 = MTs{1, MT_Number}(i_InsideMT - 1, 2);
            x2 = MTs{1, MT_Number}(i_InsideMT, 1);
            y2 = MTs{1, MT_Number}(i_InsideMT, 2);            
            TotalMTLength = TotalMTLength + sqrt((x2 - x1)^2 + (y2 - y1)^2);
        end
    end
    % Array containing values of total length of all cell's MTs for all cells
    TotalMTLengthAllCells = [TotalMTLengthAllCells, TotalMTLength];
    TotalMTLength = 0;
%% Analysing MT max values along all MTs in order to get a single value per cell
%% Fitting of the histogram of all values for a cell with a gaussian
    if (length(AllMTsTogether) == 0)
        continue
    end
    x = 1:length(AllMTsTogether);
%     figure(FigNb); FigNb = FigNb + 1;
%     %hold on
%     %plot(AllMTsTogether(:,3), 'o', 'MarkerSize', 3);
%     plot(AllMTsTogether, 'o', 'MarkerSize', 3);   
    
%     %p = polyfit(x', AllMTsTogether(:,3), 0);
%     p = polyfit(x', AllMTsTogether, 0);
%     ResultFit = p;
%     plot(x, ResultFit, 'o', 'MarkerSize', 2);%     
%     title('Horisontal Fit');        
%     xlabel('n');
%     ylabel('MT max values');   
%     hold off    

    NbBins = ceil(simple_max(AllMTsTogether)) / NbBinsCoeff;
    [HistNumbers, HistXout] = hist(AllMTsTogether, NbBins);
    
%% Artificial way of fighting with the 'first tail' of the fitting goin up      
    %HistNumbers = [0 0 0 0 0 0 HistNumbers 0 0 0 0 0 0];
    HistNumbers = [HistNumbers];
    x = -5:(length(HistNumbers) - 6);   
    Spl = csaps(x, HistNumbers, p_AllMTFit, x);
    
    figure(FigNb); FigNb = FigNb + 1;
    hold on;    
%!!!!!!!!!!!!!
    %plot(x, HistNumbers, 'o', 'MarkerSize', 2);
    bar(x, HistNumbers);
    plot(x, Spl, '-r', 'LineWidth', 2);
    hold off;       
    %pause(1);    
    
    [a, ResultFit] = simple_max(Spl);     % In HistXout are stored the centers of histogram bins
%     ResultFit = x(ResultFit);
%     ResultFit = HistXout(ResultFit);           
%% Representation of 'Results' matrix after each round instead of doing it in the end because the function 'load' exits with the error of being unable to read the next file which is inexistant    
    Result = [Result; ResultFit];       % Value obtained after several fittings for all MT analysed together is added    
%     figure(FigNb); FigNb = FigNb + 1;
%     plot(Result, 'bs', 'MarkerSize', 3);
%     grid on
%     title('MT max values obtained after several rounds of line fitting');        
%     xlabel('n');
%     ylabel('MT max values');    
end 

save(ResOutput, 'Result');
%% Visualisation of the results
figure(FigNb); FigNb = FigNb + 1;
plot(Result, 'bs', 'MarkerSize', 3);
grid on
xlabel('n');
saveas(FigNb-1, '_OutputImagesAndVariables/_ResMTmax');

NbBins = ceil(simple_max(Result) - min(Result))/20;
figure(FigNb); %FigNb = FigNb + 1;
[NbsHist, Xout] = hist(Result, NbBins);
bar(Xout, NbsHist);
ylabel('N');
saveas(FigNb-1, '_OutputImagesAndVariables/_HistRes_MTmax_Preliminary');

% figure(FigNb); FigNb = FigNb + 1;
% plot(NbsHist, 'bs', 'MarkerSize', 3);
% hold on
% x_Fit = 1:0.2:length(NbsHist);
% Spl = csaps(1:length(NbsHist), NbsHist, p_ResFit, x_Fit);
% plot(x_Fit, Spl, '-r');
% hold off
% grid on
% ylabel('N');
% 
% [NbsHist, Xout] = hist(TotalMTLengthAllCells, NbBins);
% bar(Xout, NbsHist);
% ylabel('N');
