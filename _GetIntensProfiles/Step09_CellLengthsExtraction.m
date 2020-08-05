%% Lengths of all cells. This information is obtained from 'output_CellsCenterEndsIntens.txt'
clear;
% File names
PathInput = '_OutputGI/output_CellsCenterEndsIntens.txt';
PathCellLengths = '_OutputGI/_CellLengths.mat';
%% Reading all necessary data from the input file 
fid = fopen(PathInput, 'r');       
Input = textscan(fid, '%s%f%f%f%f%f%f%f%f%f%f%f%f', 'headerLines', 1);    
fclose(fid);
CellLengths = Input{3};
save(PathCellLengths, 'CellLengths');


