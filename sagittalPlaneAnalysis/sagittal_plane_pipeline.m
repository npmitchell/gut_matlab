%% Options
clear all
overwrite = false ;
preview = true ;


%% Paths
addpath('/mnt/data/code/gut_matlab/sagittalPlaneAnalysis/')
addpath('/mnt/data/code/gut_matlab/addpath_recurse/')
addpath_recurse('/mnt/data/code/gptoolbox/external/toolbox_fast_marching/')
addpath_recurse('/mnt/data/code/gptoolbox/mex/')
addpath_recurse('/mnt/data/code/gut_matlab/distanceTransform/')
addpath_recurse('/mnt/data/code/gut_matlab/curve_functions/')
addpath_recurse('/mnt/data/code/gut_matlab/tiff_handling')
addpath_recurse('/mnt/data/code/gut_matlab/data_handling')
addpath_recurse('/mnt/data/code/imsane_for_git/imsane/generalfunctions/')
addpath_recurse('/mnt/data/code/gut_matlab/plotting/')
addpath_recurse('/mnt/data/code/gut_matlab/geometry/')

%% Settings
rootDir = '/mnt/data/confocal_data/gut/sagittal_plane/';

% GFP
% labelDir = fullfile(rootDir, '48YGAL4UASCAAXmChH2AGFP') ;
% dataDir = fullfile(labelDir, '202009142006_sagittal_48YGAL4klarUASCAAXmChH2AGFP_40x_halooil27immersol_left');
% dataDir = fullfile(labelDir, '202009151613_sagittal_48YGAL4klarUASCAAXmChH2AGFP_40x_halooil27immersol_right');
% dataDir = fullfile(labelDir, '20200913_sagittalVentral_48YGAL4UASCAAXmChH2AGFP_40x_water') ;

% HiFP
labelDir = fullfile(rootDir, '48YGAL4UASCAAXmChHiFP') ;
dataDir = fullfile(labelDir, '202009162117_sagittal_ventral_48YGal4klarUASCAAXmChHiFP_40x_haloOilImmersol_120spf') ;
% dataDir = fullfile(labelDir, '202009182053_48YGAL4klarUASCAAXmChHiFP_1p1x_40x_hcoilImmersol_120spf') ;

cd(dataDir)
settingsFn = fullfile(dataDir, 'masterSettings.mat') ;
if exist(settingsFn, 'file') && ~overwrite
    disp('Loading settings...')
    load(settingsFn, 'settings')
else
    name = '2117_sagittal_Z%d_T%03d' ;
    timePoints = 0:40 ;
    zplanes = 0:1 ;
    channels = [1,2]; 
    
    % 202009142006_sagittal_48YGAL4klarUASCAAXmChH2AGFP_40x_halooil27immersol_left
    % resolution = 0.284090909 ;  % um per pixel
    
    % 202009151613_sagittal_48YGAL4klarUASCAAXmChH2AGFP_40x_halooil27immersol_right
    % resolution = 0.25896480938 ;
    
    % 202009162117_sagittal_ventral_48YGal4klarUASCAAXmChHiFP_40x_haloOilImmersol_120spf
    % resolution = 0.284090909 ;  % um per pixel 
    
    % 202009182053_48YGAL4klarUASCAAXmChHiFP_1p1x_40x_hcoilImmersol_120spf
    % resolution = 0.257705962 ;  % um per pixel -- 1.1x 40x
    
    norm1 = 120 ; % intensity normalization for channel 1
    norm2 = 250 ; % intensity normalization for channel 1
    dt = 2 ;
    minsz_hole = 350 ;
    settings = struct('dataDir', dataDir, ...
        'name', name, ...
        'zplanes', zplanes, ...
        'timePoints', timePoints, ...
        'npts_skel', 10000, ...
        'resolution', resolution, ...
        'norms', [norm1, norm2], ...
        'channels', channels, ...
        'spaceUnits', '$\mu$m', ...
        'timeUnits', 'min', ...
        'dt', dt, ...
        'minsz_hole', minsz_hole) ;
    save(settingsFn, 'settings') ;
end

%% Unpacking
SE = SagittalExperiment(settings) ;
SE.plotting.preview = preview ;

%% Segment out the membrane
for zz = 0 % SE.zplanes
    
    % frame = dir(fullfile(zDir, [name(1:4) '*.tif'])) ;
    % probs = dir(fullfile(zDir, [name(1:4) '*Probabilities.h5'])) ;
    
    %% Segmentation // pathMasks // thicknessMasks
    for tidx = 1:length(SE.timePoints)
        tp = SE.timePoints(tidx) ;        
        disp(['Segmentation & masks: t = ' num2str(tp)])
        SE.setZPlane(zz) ;
        SE.setTime(tp) ; 
        
        im = SE.getIm() ;
        imL = SE.getImL() ;
        
        % Segment out the membrane/gut cross section
        SE.getBW() ;
    
        % Landmarks for curve extraction
        SE.getLandmarks() ;
    end
    
    %% Check that number of landmarks is correct
    for tidx = 1:length(SE.timePoints)
        tp = SE.timePoints(tidx) ;        
        disp(['Segmentation & masks: t = ' num2str(tp)])
        SE.setZPlane(zz) ;
        SE.setTime(tp) ; 
        lm = SE.getLandmarks() ;
        resave = false ;
        for qq = 1:length(lm)
            if length(find(lm{qq}.id == '/')) == length(lm{qq}.v)
                disp('good')
            else
                disp('bad')
                disp(lm{qq}.id)
                imL = SE.getImL() ;
                imshow(imL); hold on;
                scatter(lm{qq}.v(:, 1), lm{qq}.v(:, 2), 20, 'filled')
                correct = input('Correct landmarks: ', 's') ;
                lm{qq}.id = correct ;
            end
        end
        landmarks = lm ;
        if resave
            save(SE.filename.landmarks, 'landmarks')
        end
    end
    close all
    
    %% 
    for tidx = length(SE.timePoints)
        tp = SE.timePoints(tidx) ;        
        disp(['Segmentation & masks: t = ' num2str(tp)])
        SE.reset()
        SE.setZPlane(zz) ;
        SE.setTime(tp) ; 
        
        % Path masks for curve extraction
        SE.getPathMasks(false, true) ;
        
        % Thickness masks for thickness extraction via DT sampling
        SE.getThicknessMasks(true, 'fix') ;
    end

    %% Now extract curves & measure thickness
    clearvars skel_ss skels pathlength tp tidx im DTs landmarks skelFull
    zDir = fullfile(dataDir, ['z' num2str(zz)]) ; 
    zOutDir = [zDir '_results'] ;
    for tidx = 38:length(SE.timePoints)
        tp = SE.timePoints(tidx) ;
        disp(['Extracting curves: t=', num2str(tp)])
        SE.setZPlane(zz) ; 
        SE.setTime(tp) ;
        SE.getPathMasks() ;
        SE.getSkel() ;
        SE.getThicknessMasks(false, 'prev') ;
        SE.getThickness() ;
    end
    
    %% Fix pathMasks
    for tidx = 36:36 %length(SE.timePoints)
        tp = SE.timePoints(tidx) ;
        disp(['Extracting curves: t=', num2str(tp)])
        SE.reset()
        SE.setZPlane(zz) ; 
        SE.setTime(tp) ;
        SE.getPathMasks(true, false, 'fix') ;
    end
    
    %% Fix thicknessMasks
    for tidx = 36:36 %length(SE.timePoints)
        tp = SE.timePoints(tidx) ;
        disp(['Extracting curves: t=', num2str(tp)])
        SE.reset()
        SE.setZPlane(zz) ; 
        SE.setTime(tp) ;
        SE.getImL() ;
        SE.getIm() ;
        SE.getThicknessMasks(true, 'prev') ;
    end
    
    %% Check all landmark names
    for tidx = 40  % :length(SE.timePoints)
        tp = SE.timePoints(tidx) ;
        disp(['thickness: t = ', num2str(tp)])
        SE.reset()
        SE.setZPlane(zz) ; 
        SE.setTime(tp) ;
        landmarks = SE.getLandmarks() ;
        if isempty(landmarks{1}.id) || true 
            lmkfn = fullfile(SE.zDir.landmarks, ...
                sprintf('landmarks_%04d.mat', SE.currentTime)) ;
            load(lmkfn, 'landmarks') ;
    
            imL = SE.getImL() ;
            imshow(imL) ; hold on;
            for qq = 1:length(landmarks{1})
                plot(landmarks{1}.v(:, 1), landmarks{1}.v(:, 2), '.-')
            end
            
            disp('Old landmark Names: ')
            disp(landmarks{1}.id)
            landmarks{1}.id = input('New landmark names:', 's') ;
            error('here')
            save(lmkfn, 'landmarks')
        end
    end
    
    %% Measure thicknesses
    for tidx = 1:length(SE.timePoints)
        tp = SE.timePoints(tidx) ;
        disp(['thickness: t = ', num2str(tp)])
        SE.setZPlane(zz) ; 
        SE.setTime(tp) ;
        SE.getThickness() ;
    end

    %% Redo thickness unc
    % for tidx = 1:length(SE.timePoints)
    %     tp = SE.timePoints(tidx) ;
    %     disp(['thickness: t = ', num2str(tp)])
    %     SE.setZPlane(zz) ; 
    %     SE.setTime(tp) ;
    %     SE.getThickness() ;
    %     avg_thickness = SE.thickness.avg_thickness ;
    %     thickness_eval = SE.thickness.thickness_eval ;
    %     thickness_medial = SE.thickness.thickness_medial ;
    %     unc_thickness = SE.thickness.unc_thickness ;
    %     for qq = 1:length(avg_thickness)
    %         unc_thickness_um{qq} = unc_thickness{qq} * SE.resolution ;
    %     end
    %     xo = SE.thickness.xo ;
    %     yo = SE.thickness.yo ;
    %     unc_thickness = SE.thickness.unc_thickness ;
    %     avg_thickness_um = SE.thickness.avg_thickness_um ; 
    %     landmarkIds = SE.thickness.landmarkIds ;
    %     landmarkNames = SE.thickness.landmarkNames ;
    %     thicknessFn = fullfile(SE.zDir.thickness, ...
    %         sprintf('thickness_%04d.mat', SE.currentTime)) ;
    %     save(thicknessFn, 'thickness_eval', 'thickness_medial', ...
    %         'avg_thickness', 'unc_thickness', ...
    %         'avg_thickness_um', 'unc_thickness_um', ...
    %         'xo', 'yo', 'landmarkIds', 'landmarkNames')
    % end

    %% Plot intensity kymograph versus pathlength
    for tidx = 1:length(SE.timePoints)
        tp = SE.timePoints(tidx) ;
        disp(['skel intensity, t=', num2str(tp)])
        SE.reset() ;
        SE.setZPlane(zz) ; 
        SE.setTime(tp) ;
        
        % Load landmarks
        SE.getSkelIntensity(true) ;
    end
    
    %% Draw kymographs
    SE.getKymographs(false, true)
    
    %% Draw thickness kymographs
    SE.getThicknessKymographs(true)
        
    %% Draw curvature kymographs
    SE.reset() ;
    SE.setZPlane(zz) ; 
    SE.getSkelCurvature()
    
    %% Select lobe regions to include in measurement
    for tidx = 1:length(SE.timePoints) 
        tp = SE.timePoints(tidx) ;
        SE.setZPlane(zz) ; 
        SE.setTime(tp) ; 
        disp(['lobeMasks: t = ', num2str(tp)])
        SE.getLobeFilters() ;
    end
    
    %% 
    SE.getLobeThicknesses() ;
    
    %% Select fold regions to include in measurement
    SE.identifyFolds() ;
end
