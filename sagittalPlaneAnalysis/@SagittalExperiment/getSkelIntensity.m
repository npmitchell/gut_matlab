function skelIntensity = getSkelIntensity(SE, overwrite)


if nargin < 2
    overwrite = false ;
end

% Check that all fields are present in SE.thickness
allFieldsPresent = true ;
try
    SE.thickness.skels ;
    SE.thickness.DTs ;
    SE.thickness.skelFull ;
    SE.thickness.skel_ss ;
    SE.thickness.pathlength ;
    if ~isempty(SE.thickness.skels) || ~isempty(SE.thickness.DTs) ...
            || ~isempty(SE.thickness.skelFull) ...
            || ~isempty(SE.thickness.skel_ss) ...
            || SE.thickness.pathlength            
        allFieldsPresent = false ;
    end
catch
    allFieldsPresent = false ;
end

if allFieldsPresent && ~overwrite
    if nargout > 0
        skelIntensity = SE.skel ;
    end
else
    % Load or compute
    skelIntensityFn = fullfile(SE.zDir.skelIntensity, ...
        sprintf('skelIntensity_%04d.mat', SE.currentTime)) ;
    if exist(skelIntensityFn, 'file') && ~overwrite
        SE.skelIntensity = load(skelIntensityFn) ;
        if nargout > 0
            skelIntensity = SE.skelIntensity ;
        end
    else
        % Original image -- unnormalized
        fname = fullfile(SE.zDir.zIm, ...
            [sprintf(SE.name, SE.zplane, SE.currentTime) '.tif']) ;
        im = loadtiff(fname) ;
        
        % Split into channels
        ch1 = squeeze(im(:, :, 1)) ;
        ch2 = squeeze(im(:, :, 2)) ;

        % Clear prev vars / preallocate
        s1 = cell(0) ;
        s2 = cell(0) ;
        pathlength_um = cell(0) ;
        
        % Load skeleton
        skel = SE.getSkel() ;
        skel_ss = skel.skel_ss ;
        
        % Load thickness for landmarkIds and landmarkNames
        thickness = SE.getThickness() ;
        landmarkIds = thickness.landmarkIds ;
        landmarkNames = thickness.landmarkNames ;
        
        % Sift through landmarkNames and toss out utility names. Keep only
        % real landmarks
        for qq = 1:length(landmarkNames)
            keep = find(~contains(landmarkNames{qq}, 'u')) ;
            landmarkIds{qq} = landmarkIds{qq}(keep) ;
            landmarkNames{qq} = landmarkNames{qq}(keep) ;
        end
        
        for qq = 1:length(skel_ss)
            
            % Extract intensity at skel
            imIdx = sub2ind(size(ch1), round(skel_ss{qq}(:, 2)), ...
                round(skel_ss{qq}(:, 1))) ;
            s1{qq} = ch1(imIdx) ;
            s2{qq} = ch2(imIdx) ;

            % Save the skeleton intensities
            pathlength_um{qq} = skel.pathlength{qq} * SE.resolution ;
        end
        
        % Save output
        save(skelIntensityFn, 's1', 's2', 'pathlength_um', ...
            'landmarkIds', 'landmarkNames')
    end
end