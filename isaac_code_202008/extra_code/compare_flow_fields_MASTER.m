%% Compare Flow Fields from All Datasets
% Isaac B Breinyn 2020
%
% This is a script that reads in all the data from comparing flow fields
% of various datasets and combines them in a single workspace. All data for
% this script should already exist. If some data is missing from a
% particular dataset, run the compare_flow_fields_____.m script that
% corresponds to that dataset. The only manual input in this script is time
% offset between datasets (calculated by qualitatively extracting fold
% depth using ImageJ)
%
%% Clean up Matlab and define Datasets

close all; clear all; clc ;

% We now have to collect all of the relevant information and data necessary
% to creating a single workspace. We will do this by organizing the
% following information into a structure with fields: data name, data set
% location, binned velocities, relative velocities, the full flow data, the
% phases, and the normed dot products. All of these are stored in the
% datasets' file directories.

% First, the base name of each dataset (will pop up in naming conventions)
data.dset{1} = '1644_folding' ;
data.dset{2} = '2038_e5' ;
data.dset{3} = 'midfold_e4' ;

% Next, the directories of each dataset
data.dir{1} = 'F:\Streichan\data\202003111630_mef2gal4klarUASCAAXmChHiFP_wo_great\1644_folding_30s_2um_la8_4p5zoom_63x\' ;
data.dir{2} = 'F:\Streichan\data\202003111830_mef2gal4klarUASCAAXmChHiFP_wo_63x_e4e5\2038_e5\' ;
data.dir{3} = 'F:\Streichan\data\202003111830_mef2gal4klarUASCAAXmChHiFP_wo_63x_e4e5\midfold_e4\' ;

% Just to save the information, let's store the (t,x,y,z) resolution
data.res{1} = [30 .0802 .0802 1.9994] ;
data.res{2} = [120 .0401 .0401 1.4996] ;
data.res{3} = [60 .0833 .0833 .9999] ;

for ii = 1:length(data.dset)
    try
        data.binvesl{ii} = importdata([data.dir{ii} 'binvels.m']) ; % The binned velocities
        data.relvels{ii} = importdata([data.dir{ii} 'relvels.m']) ; % The relative velocities
        data.fullflow{ii} = importdata([data.dir{ii} 'fullflow.m']) ; % The full flow data (as outputted by PIVlab)
        data.phases{ii} = importdata([data.dir{ii} 'phases.m']) ; % The phases
        data.surfaceArray{ii} = importdata([data.dir{ii} 'surfaceArray.m']) ;% The surface array (used to find time offsets)
    catch
        disp(['Data missing from dataset ' data.dset{ii} ', please rerun compare_flow_fields_' data.dset{ii} '.m now']) ;
    end
end

%% Use surface array to find time offset
% The surface array contains the surface in the FOV as a function of time.
% By comparing this surface, we can determine a reasonable time offset
% between datasets.

close all

for ii= 1:length(data.dset)
    sArray = data.surfaceArray{ii} ; % load the surface for this dset
    sArray = abs(sArray); % take the sqrt of the square of all values
    % now grab the max value for each TP
    for tt = 1:size(sArray, 1)
        slice = squeeze(sArray(tt,:,:)) ;
        slice = abs(medfilt2(slice, [2 2])) ;
        msArray(tt) = max(slice(:)) ;
        if max(slice(:)) > 1000
            print('LIMIT')
        end
        msArray = msArray.*data.res{ii}(4) ; % multiply the max height by the resolution in z
    end
    plot(msArray) ;
    hold on
end
legend

%% Combine data into working arrays
% Now that the data from each dataset has been consolidated into a single
% workspace, we can combine datasets into working arrays that can then be
% plotted and/or manipulated. In order to do so, we have to define a time
% offset between the datasets. This will be done by qualitatively looking
% at the fold depth in an image processing tool such as ImageJ. Will add
% more notes after doing so.

% Use surface array to plot surface as a functin of time. Use that as
% alignment between datasets.