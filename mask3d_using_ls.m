% MASK3D_USING_LS crops out a layer around a given level set surface, 
% +/- some specified thickness, from a 3d volume and plots it. 
% Note that an alternative approach is to use the level sets, dilate and
% erode them, then multiply that by the data. 
%
% Run this from projectDir (ie where the data lives)

% Parameters
rdilate = 4 ;
rerode = 3 ;
ssfactor = 4;

% Define directories
dataDir    =  cd; 
projectDir = cd ;

% Start by creating an experiment object, optionally pass on the project
% directory (otherwise it will ask), and change into the directory of the
% data.  This serves as a front-end for data loading, detection, fitting
% etc.
xp = project.Experiment(projectDir, dataDir);

fn = 'Time_%06d_c1_stab';
fileMeta                    = struct();
fileMeta.dataDir            = dataDir;
fileMeta.filenameFormat     = [fn, '.tif'];
fileMeta.nChannels          = 1;
fileMeta.timePoints         = 110:190;
fileMeta.stackResolution    = [.2619 .2619 .2619];
fileMeta.swapZT             = 1;
% first_tp is also required, which sets the tp to do individually.
first_tp = 1 ;
expMeta                     = struct();
expMeta.channelsUsed        = 1;
expMeta.channelColor        = 1;
expMeta.description         = 'Apical membrane in Drosophila gut';
expMeta.dynamicSurface      = 1;
expMeta.jitterCorrection    = 0;  % 1: Correct for sample translation
expMeta.fitTime             = fileMeta.timePoints(first_tp);
expMeta.detectorType        = 'surfaceDetection.integralDetector';
expMeta.fitterType          = 'surfaceFitting.meshWrapper';
% Now set the meta data in the experiment.
xp.setFileMeta(fileMeta);
xp.setExpMeta(expMeta);
xp.initNew();
clear fileMeta expMeta

%% Define msls directory
foreGroundChannel = 1;
ofn_ls = 'msls_apical_stab_%06d' ;
lambda1 = 1 ;
lambda2 = 1 ;   
smoothing = 0.1 ;
nu = 0.0 ;
pre_nu = -5 ;
pre_smoothing = 1 ;
post_nu = 2;
post_smoothing = 4 ;

% Name the output mesh directory ------------------------------------------
msls_exten = ['_prnu' strrep(strrep(num2str(pre_nu, '%d'), '.', 'p'), '-', 'n')];
msls_exten = [msls_exten '_prs' strrep(num2str(pre_smoothing, '%d'), '.', 'p') ];
msls_exten = [msls_exten '_nu' strrep(num2str(nu, '%0.2f'), '.', 'p') ];
msls_exten = [msls_exten '_s' strrep(num2str(smoothing, '%0.2f'), '.', 'p') ];
msls_exten = [msls_exten '_pn' num2str(post_nu, '%d') '_ps',...
    num2str(post_smoothing)];
msls_exten = [msls_exten '_l' num2str(lambda1) '_l' num2str(lambda2) ];
if projectDir(end) ~= '/'
    projectDir = [projectDir '/'];
end
mslsDir = [projectDir 'msls_output'];
mslsDir = [mslsDir msls_exten '/'] ;

% Make output directories for figures
exten = sprintf('_d%03d_e%03d', rdilate, rerode) ;
outfigdir = fullfile(mslsDir, ['mask3dls' exten]) ;
if ~exist(outfigdir, 'dir')
    mkdir(outfigdir)
end
outmip1dir = fullfile(outfigdir, ['mip1_raw' exten]) ;
outmip2dir = fullfile(outfigdir, ['mip2_raw' exten]) ;
outmip3dir = fullfile(outfigdir, ['mip3_raw' exten]) ;
tmp = {outmip1dir, outmip2dir, outmip3dir} ;
for i = 1:length(tmp)
    if ~exist(tmp{i}, 'dir')
        mkdir(tmp{i})
    end
end

% Load rotation, translation, and limits in APDV coords
meshdir = mslsDir ;
fns = dir(fullfile(meshdir, 'mesh_apical_stab_0*.ply')) ;
rotname = fullfile(meshdir, 'rotation_APDV') ;
transname = fullfile(meshdir, 'translation_APDV') ;
xyzlimname = fullfile(meshdir, 'xyzlim_APDV') ;

%% Load transformations
% Load the rotation matrix
rot = importdata([rotname '.txt']) ;
% Load the translation to put anterior to origin
trans = importdata([transname '.txt']) ;
% Load plotting limits
xyzlim = importdata([xyzlimname '.txt']) ;
xmin = (xyzlim(1) - rdilate * ssfactor) * resolution ;
ymin = (xyzlim(2) - rdilate * ssfactor) * resolution ;
zmin = (xyzlim(3) - rdilate * ssfactor) * resolution ;
xmax = (xyzlim(4) - rdilate * ssfactor) * resolution ;
ymax = (xyzlim(5) - rdilate * ssfactor) * resolution ;
zmax = (xyzlim(6) - rdilate * ssfactor) * resolution ;

%% Create structural element for dilation/erosion
se_d = strel('sphere', rdilate) ;
se_e = strel('sphere', rerode) ;

%% go through each level set. Dilate, erode and subtract
for tt = xp.fileMeta.timePoints
    name = ['msls_shell' exten] ;
    
    % Load levelset
    disp(['Loading LS for t = ' num2str(tt)])
    lsfn = fullfile(mslsDir, [sprintf(ofn_ls, tt) '.h5']) ;
    ls = h5read(lsfn, '/implicit_levelset') ;
    
    % Expand the ls to full resolution
    expand = ones(ssfactor, ssfactor, ssfactor);
    ls = superkron(double(ls), expand) ;
    ls = permute(ls, [2 1 3]) ;
    
    lsd = imdilate(ls, se_d) ;
    lse = imdilate(ls, se_e) ;
    mask = lsd - lse ;
    
    % Apply to the data
    xp.loadTime(tt);
    xp.rescaleStackToUnitAspect();
    tmp = xp.stack.image.apply() ;
    data = tmp{1} ;
    
    % Pad mask with zeros to match size of data
    if all(size(data) - size(mask) < 0)
        [xl, yl, zl] = size(data) ;
        padmask = mask(1:xl, 1:yl, 1:zl) ;
    else
        error('Unexpected size of mask array: smaller than data ')
        % padmask = padarray(mask, size(data) - size(mask), 0, 'post') ;
    end
    masked_data = uint16(padmask) .* data ;
    mdata16 = uint16(double(masked_data) / double(max(masked_data(:))) * (2^16 - 1)) ;
    
    % Make MIPs
    disp(['Making NMIPs for t = ' num2str(tt)])
    [m1, inds] = max(mdata16, [], 1) ;
    mip1 = squeeze() ;
    mip2 = squeeze(max(mdata16, [], 2)) ;
    mip3 = squeeze(max(mdata16, [], 3)) ;
    
    % Save mips as images
    close all
    imwrite(mip1, fullfile(outmip1dir, [name '.png'])) ;
    imwrite(mip2, fullfile(outmip2dir, [name '.png'])) ;
    imwrite(mip3, fullfile(outmip3dir, [name '.png'])) ;
end


%% Now rotate each datavolume, interpolate onto new grid, and save
for tt = xp.fileMeta.timePoints
    % rotate in 3d to align
    [mx, my, mz] = meshgrid(1:xl, 1:yl, 1:zl) ;
    mxv = reshape(mx, [xl*yl*zl, 1]) ;
    myv = reshape(my, [xl*yl*zl, 1]) ;
    mzv = reshape(mz, [xl*yl*zl, 1]) ;
    mxyzr = rot .* [mxv, myv, mzv] ;  
    
    % define query points
    % [Xq, Yq, Zq] = meshgrid(xmin:res:xmax, ymin:res:ymax, zmin:res:zmax) ;
    % interp3(X, Y, Z, masked_data, Xq, Yq, Zq)
    % bad = good
end

