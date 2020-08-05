function [CellVolume] = f_VolumeIntegral(CoordsXY)
CellVolume = 0;
dx = 1;
for x = floor(min(CoordsXY(:, 1))) + 1 : ceil(max(CoordsXY(:, 1)))
    ind = find((CoordsXY(:, 1) >= x - dx) .* (CoordsXY(:, 1) < x));
    if isempty(ind)
        ind = indOld;
    end
    % If more than one pt is in the X-interval, take the average value of the function
    Fx = (max(CoordsXY(ind, 2)) + min(CoordsXY(ind, 2))) / 2;     
    CellVolume = CellVolume + pi * Fx^2 * dx;  
    indOld = ind;
end