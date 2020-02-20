function [map, flatmap] = buildLookupMap(lookupMetaFn, save_map) 
%BUILDLOOKUPMAP Principal function for lookupTable Class object construction
%
% Parameters
% ----------
% lookupMetaFn : str
%   file to find metadata for datasets: /mnt/data/analysis/lookupMeta.txt
% save_map : bool (optional, default=false)
%   whether to save the map object to disk
%
% Returns
% -------
% map : containers.Map object
%   The map from markers to datasets
% flatmap : #datasets x 1 cell array
%   an array of structs, each giving details on a dataset
%   
% Example Usage 
% -------------
% 
%
% NPMitchell 2020

%% Build a lookuptable for finding embryos at a given time of a given datatype
map = containers.Map() ;

%% Read all the labeltype folders from txt file
fid = fopen(lookupMetaFn, 'r');
dat = textscan(fid, '%s%d%d%s', 'HeaderLines',1) ;
fclose(fid) ;
markers = dat{1} ;
channels = dat{2} ;
nch = dat{3} ;
dirs = dat{4} ;
%% Obtain all markers
allMarkers = unique(markers) ;

%% Examine each channel in turn 
for mii = 1:length(allMarkers)
    marker = allMarkers{mii} ;

    first = true ;
    for qq = 1:length(markers) 
        % Is this dataset a marker match?
        if strcmp(markers{qq}, marker)
            if first
                substruct.folders = {dirs{qq}} ;
                substruct.channels = [channels(qq)] ;
                substruct.nchannels = [nch(qq)] ;
                % Load low-RAM data about the dataset here, like volume,
                % sa, centerline, etc        
                fn = fullfile(dirs{qq}, 'surfacearea_volume_stab.mat') ;
                if exist(fn, 'file')
                    load(fn, 'aas', 'vvs')
                    substruct.area = {aas} ;
                    substruct.volume = {vvs} ;
                else
                    disp(['Warning: could not load areas, volumes for ' dirs{qq}])
                    substruct.area = {} ;
                    substruct.volume = {} ;
                end
                first = false ;
            else
                substruct.folders{length(substruct.folders) + 1} = dirs{qq} ;
                substruct.channels(length(substruct.channels) + 1) = channels(qq) ;
                substruct.nchannels(length(substruct.nchannels) + 1) = nch(qq) ;
                
                % Load low-RAM data about the dataset here, like volume,
                % sa, centerline, etc
                fn = fullfile(dirs{qq}, 'surfacearea_volume_stab.mat') ;
                if exist(fn, 'file')
                    load(fn, 'aas', 'vvs')
                    substruct.area{length(substruct.area)+1} = aas ;
                    substruct.volume{length(substruct.volume)+1} = vvs ;
                else
                    disp(['Warning: could not load areas, volumes for ' dirs{qq}])
                    substruct.area{length(substruct.area)+1} = [] ;
                    substruct.volume{length(substruct.volume)+1} = [] ;
                end
            end
        end
    end
    % Add this marker info to the map
    map(marker) = substruct ;
end

% Flatten the map into a #datasets x 1 cell array
if nargout > 1
    % build flattened map
    flatmap = {} ;
    keys = map.keys ;
    for kk = 1:length(keys)
        key = keys{kk} ;
        tmp = map(key) ;
        tmp.marker = key ;
        flatmap{length(flatmap) + 1} = tmp ;
    end
end

if nargin > 1
    if save_map
        % Save the lookup Map
        save_map(fullfile(mutantDir, 'lookuptable_containersMap.mat'), 'map')
    end
end
disp('done building map(s)')


