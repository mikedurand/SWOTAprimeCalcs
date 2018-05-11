%This script reads the netcdf files, filters to include only good reaches
%and (if the switches are equal to 1) plot reaches and nodes.
%Author: Renato Frasson
% v0: 2/13/18 RPF
% v1: 2/23/18 MTD. Adapted for Mac. Miscellaneous edits. 

clear all

pathtoncfiles='./Pepsi2/';
Files=dir([pathtoncfiles '*.nc']);
numbfiles=size(Files,1);
plotreaches=2; %0: no plots. 1: long profiles only. 2: time series only. 3: both.
plotnodes=0; %0: no plots. 1: long profiles .

for cf=1:numbfiles,
    
    ReadRiver;
    
    FilterRiver;
    
    if plotreaches==1 || plotreaches == 3
        PlotReachLongProfiles;
    end
    
    if plotreaches==2 || plotreaches == 3,
        PlotReachTimeseries;
    end
    
    if plotnodes==1 || plotreaches == 3
        PlotNodeLongProfiles;
    end
end



    