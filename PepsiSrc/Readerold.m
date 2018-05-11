%This script reads the netcdf files.
%Author: Renato Frasson
%February 13, 2018
clear all

%pathtoncfiles='./Pepsi1/';
pathtoncfiles='./Pepsi2/';
Files=dir([pathtoncfiles '*.nc']);
numbfiles=size(Files,1);
for cf=1:numbfiles
    ReadRiver;
end
    