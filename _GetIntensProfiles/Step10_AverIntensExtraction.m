%% The aim of the program is to create a .mat file with array with lengths
%% of all cells. This information is obtained from 'output_CellsCenterEndsIntens.txt'
clear;
% File names
PathInput = '_OutputGI/output_CellsCenterEndsIntens.txt';
PathAverIntens = '_OutputGI/_AverIntens.mat';
%% Reading all necessary data from the input file 
fid = fopen(PathInput, 'r');       
Input = textscan(fid, '%s%f%f%f%f%f%f%f%f%f%f%f%f', 'headerLines', 1);    
fclose(fid);
AverIntens = Input{12};
TotIntens = Input{13}; 
save(PathAverIntens, 'AverIntens');


