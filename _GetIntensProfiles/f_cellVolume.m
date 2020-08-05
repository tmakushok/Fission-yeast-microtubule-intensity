function [CellVolumeAvg, Err] = f_cellVolume(CellPixels, CellTips, Angle)
%% Visualisation of the cell mask
FinalImage = [];
for i_OnePix = 1:length(CellPixels)     % Better to use 'size', but here it doesn't matter
    FinalImage(CellPixels(i_OnePix, 2), CellPixels(i_OnePix, 1)) = 1;
end
% Visualisation of cell axis
% figure, imshow(FinalImage);
% hold on;
% line([CellTips(1, 1), CellTips(2, 1)], [CellTips(1, 2), CellTips(2, 2)]); 
PerimIm = bwperim(FinalImage);
% figure, imshow(PerimIm);
PerimPts = [];
[PerimPts(:,2), PerimPts(:,1)] = find(PerimIm);     % PerimPts: [X, Y]
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
Angle = -Angle;
for i = 1:length(PerimPts)
    PerimPts(i, 1) = PerimPts(i, 1) - TheTip(1);
    PerimPts(i, 2) = PerimPts(i, 2) - TheTip(2);
%% Rotating the cell (the tip will stay unchanged)
    PerimPtsRot(i, 1) = PerimPts(i, 1) * cosd(Angle) + PerimPts(i, 2) * sind(Angle); 
    PerimPtsRot(i, 2) = - PerimPts(i, 1) * sind(Angle) + PerimPts(i, 2) * cosd(Angle); 
    % Visualisation
%     plot(PerimPtsRot(i, 1), PerimPtsRot(i, 2), 'ro');
end
% hold off
%% Extracting half outline (using the sign of converted Y)
ind1 = find(PerimPtsRot(:, 2) >= 0);
ind2 = find(PerimPtsRot(:, 2) < 0);
HalfOutline1 = PerimPtsRot(ind1, :); 
HalfOutline2 = PerimPtsRot(ind2, :); 
% Visualisation
% figure, hold on
% for i = 1:length(HalfOutline1)    
%     plot(HalfOutline1(i, 1), HalfOutline1(i, 2), 'ro');
% end    
% hold off    
%% Calculating cell volume using disk method (two values: one for each cell half)
CellVolume(1,1) = f_VolumeIntegral(HalfOutline1);
CellVolume(2,1) = f_VolumeIntegral(HalfOutline2);
CellVolumeAvg = (CellVolume(1,1) + CellVolume(2,1)) / 2;
Err(1,1) = (CellVolumeAvg - CellVolume(1,1)) / CellVolumeAvg;
Err(2,1) = (CellVolumeAvg - CellVolume(2,1)) / CellVolumeAvg;
 





























