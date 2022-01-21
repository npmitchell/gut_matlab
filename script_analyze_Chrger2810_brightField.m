%% Bright field image sequences of UAS-Chrger2810:
%  ==> Analyze their motion upon opto activation
rootdir = '/mnt/data/confocal_data/gut/mef2GAL4klarUASChrger2810/' ;
rootdir = [rootdir '202201071205_mef2G4klarUASChrger2810_3um_45spf_1p3x20x_12pc561'] ;

imseqDir = fullfile(rootdir, 'e01_0p5spf_excite12pc_1400');
imseqDir = fullfile(rootdir, 'e10/e10_0p5spf') ;


% Simp[le time series



%% Stacks of positions taken sequentially
addpath('/mnt/data/code/gut_matlab/addpath_recurse')
addpath_recurse('/mnt/data/code/gut_matlab')


rootdir = '/mnt/data/confocal_data/gut/mef2GAL4klarUASChrger2810/' ;
rootdir = [rootdir '202201071205_mef2G4klarUASChrger2810_3um_45spf_1p3x20x_12pc561'] ;
strobeDir = fullfile(rootdir, 'e04_strobed_10framesIn8sec_20sps') ;
imseqDir = fullfile(strobeDir, 'image_sequence') ;
imstkDir = fullfile(strobeDir, 'image_resequence') ;
imhisteqDir = fullfile(strobeDir, 'image_histequalize') ;
imoutDir = fullfile(strobeDir, 'image_overlays') ;
ensureDir(imstkDir) 
ensureDir(imhisteqDir) 
ensureDir(imoutDir) 

% Parameters
dt1 = 8/10 ; % seconds between one position snap and the next position's snap
dt2 = 20   ; % seconds between frames in a give position sequence
nPositions = 10 ;
histequalize = true ;
pix2um = 1./ 1.7590 ; % micron per pixel

% Compile list of all filenames as Cell of Cells
for pos = 1:nPositions
    fnpos = dir(fullfile(imseqDir, sprintf('e04_Position%03d_*.png', pos))) ;
    for qq = 1:length(fnpos)
        fnPosList{qq} = fullfile(fnpos(qq).folder, fnpos(qq).name) ;
        fnPosListNames{qq} = fnpos(qq).name ;
    end
    fns{pos} = fnPosList ;
    fnNames{pos} = fnPosListNames ;
    
end

%% PIVLAB piv
% Rename/save
dmyk = 1 ;
for jj = 1:ntps
    disp(['Resaving images for position stack, tp=' num2str(jj) ' / ' num2str(ntps)])
    % Get PIV from one position's snap to the next position's snap
    for ii = 1:nPositions
        outfn = fullfile(imstkDir, [sprintf('%06d_', dmyk) fnNames{ii}{jj}]) ;
        if ~exist(outfn, 'file')
            disp('creating file')
            im1 = imread(fns{ii}{jj}) ;
            imwrite(im1, outfn) ;
        end
        
        
        outfn = fullfile(imhisteqDir, [sprintf('%06d_', dmyk) fnNames{ii}{jj}]) ;
        if ~exist(outfn, 'file') || true
            disp('creating histeq file')
            im1 = imread(fns{ii}{jj}) ;
            im1 = histeq(im1) ;
            timestamp = sec2hms(dt1 * (ii - 1) + dt2*(jj-1), '%02d:%02d:%02.f') ;
            im1 = insertText(im1, ...
                [size(im1,2)-100, 5],timestamp,...
                'textColor', 'black', 'FontSize', 18, 'boxOpacity', 0) ;
            % scalebar
            % scbarLeft = size(im1, 2)*0.75 ;
            % scbarTop = size(im1, 1) * 0.75 ;
            % im1(scbarTop:scbarTop+8, scbarLeft:scbarLeft+50/pix2um, :) =0 
            imwrite(im1, outfn) ;
        end
        dmyk = dmyk + 1 ;
    end
end


%% load PIVLab result
piv = load(fullfile(strobeDir, 'piv_results.mat')) ;

%% Plot results from PIVlab
rmsspeed = zeros(ntps*nPositions,1) ;
nexttime = 0 ;
timestamps = rmsspeed ;
overwrite = false ;
dmyk = 1 ;
for jj = 1:ntps
    disp([num2str(jj) ': Resaving overlay velocity images for position stack, tp=' num2str(jj) ' / ' num2str(ntps)])
    % Get PIV from one position's snap to the next position's snap
    for ii = 1:nPositions
        
        im1 = imread(fns{ii}{jj}) ;
        if histequalize
           im1 = histeq(im1) ;
        end

        % associated piv
        if ii < nPositions
            vx = piv.u_filtered{dmyk} * pix2um / dt1 ;
            vy = piv.v_filtered{dmyk} * pix2um / dt1 ;
            timenow = nexttime ;
            nexttime = nexttime + dt1 ;
        else
            vx = piv.u_filtered{dmyk} * pix2um / (dt2-(nPositions-1)*dt1) ;
            vy = piv.v_filtered{dmyk} * pix2um / (dt2-(nPositions-1)*dt1) ;
            timenow = nexttime ;
            nexttime = nexttime + (dt2-(nPositions-1)*dt1) ;
        end
        speed = sqrt(vx.^2 + vy.^2) ;
        
        % Save imageoverlay
        overlay_outfn = fullfile(imoutDir, ...
            [sprintf('%06d_', dmyk) fnNames{ii}{jj}]) ;
        if ~exist(overlay_outfn, 'file') || overwrite

            XX = piv.x{dmyk} ;
            YY = piv.y{dmyk} ;
            xyfstruct = struct() ;
            tX = XX' ; tY = YY' ;
            xyfstruct.x = XX ;
            xyfstruct.y = YY ;
            vscale = 2 ;
            opts = struct() ;
            opts.outfn = overlay_outfn ;
            opts.qscale = 10 ;
            opts.qcolor = 'y' ;
            opts.label = '$|v|$ [$\mu$m/s]' ;
            opts.title = num2str(sec2hms(timenow, '%02d:%02d:%02.f')) ;
            vectorFieldHeatPhaseOnImage(im1, ...
                xyfstruct, vx, vy, vscale, opts) ;

            % % Check result as histogram
            % prodt = [vy(:), vx(:)] .* ...
            %     ([cos(pi/4)* ones(size(vx(:))), sin(pi/4)* ones(size(vx(:)))]) ;
            % prodt = prodt / (-2*sqrt(2)) ;
            % dotprodt = prodt(:, 1), + prodt(:, 2) ;
            % dotprodt = dotprodt(abs(vx) > 4e-1 | abs(vy) > 4e-1) ;
            % if length(dotprodt) < floor(104*48/64)
            %     dotprodt(floor(104*48/64)+1) = 0 ;
            %     dotprodt = dotprodt(1:floor(104*48/64)) ;
            % end
            % histogram(dotprodt, 100)
            % xlabel('$(v \cdot v_{input}) / |v_{input}^2|$', 'interpreter', 'latex')
            % ylabel('occurrence', 'interpreter', 'latex')
            % title('PIVlab velocity estimation')
            % saveas(gcf, fullfile(imoutDir, 'histogram_dotv_pivLab.png'))

        end
        
        timestamps(dmyk) = timenow ;
        rmsspeed(dmyk) = sqrt (mean (speed(:) .^2) ) ;
        
        dmyk = dmyk + 1 ;
    end
end


clf 
plot(timestamps(1:300), rmsspeed(1:300)) ;
xlabel('time [sec]')
ylabel('rms velocity')
title('ChRger 28-10 strobed illumination')
saveas(gcf, 


%% PIV using SJS code GetPIV and xcorrfft2
% edge length for PIV
eL = 12 ;
% Load first image to get imagesize
im1 = imread(fns{1}{1}) ;
step = eL*0.5 ;
[XX,YY] = meshgrid(eL/2:step:size(im1,1)-eL/2,eL/2:step:size(im1,2)-eL/2); 
uX = unique(XX) ;
uY = unique(YY) ;
ntps = length(fns{1}) ;  % number of timepoints per position

dmyk = 1 ;
VX = zeros(ntps * nPositions, size(XX, 1), size(XX, 2)) ;
VY = VX ;
for jj = 1:ntps
    disp(['Computing PIV for position stack, tp=' num2str(jj) ' / ' num2str(ntps)])
    % Get PIV from one position's snap to the next position's snap
    for ii = 1:nPositions
        im1 = imread(fns{ii}{jj}) ;
        if ii < nPositions
            im2 = imread(fns{ii+1}{jj}) ;
        else
            im2 = imread(fns{1}{jj+1}) ;
        end
       if histequalize
           im1 = histeq(im1) ;
           im2 = histeq(im2) ;
       end
        [vx,vy] = GetPIV(im1,im2,XX,YY,eL);
        VX(dmyk, :, :) = vx ;
        VY(dmyk, :, :) = vy ;        
        speed = sqrt(vx.^2 + vy.^2) ;
        
        % Save imageoverlay
        xyfstruct = struct() ;
        xyfstruct.x = YY' ;
        xyfstruct.y = XX' ;
        vscale = 2*sqrt(2) ;
        opts = struct() ;
        opts.outfn = fullfile(imoutDir, fnNames{ii}{jj}) ;
        vectorFieldHeatPhaseOnImage(ones(size(im1)), ...
            xyfstruct, vx', vy', vscale, opts)
        
        % % Check result as histogram
        % prodt = [vy(:), vx(:)] .* ...
        %     ([cos(pi/4)* ones(size(vx(:))), sin(pi/4)* ones(size(vx(:)))]) ;
        % prodt = prodt / (-2*sqrt(2)) ;
        % dotprodt = prodt(:, 1), + prodt(:, 2) ;
        % dotprodt = dotprodt(abs(vx) > 1e-5 | abs(vy) > 1e-5) ;
        % if length(dotprodt) < floor(104*48/36)
        %     dotprodt(floor(104*48/36)+1) = 0 ;
        %     dotprodt = dotprodt(1:floor(104*48/36)) ;
        % end
        % histogram(dotprodt, 100)
        % xlabel('$(v \cdot v_{input}) / |v_{input}^2|$', 'interpreter', 'latex')
        % ylabel('occurrence', 'interpreter', 'latex')
        % saveas(gcf, fullfile(imoutDir, 'histogram_dotv_SJSpiv.png'))
        
        % Save imageoverlay -- scalarVectorFieldsOnImage
        % subf = 5 ;
        % opts = struct() ;
        % opts.qsubsample = 1 ;
        % opts.qscale = 10 ;
        % opts.sscale = 5 ;
        % opts.outfn = fullfile(imoutDir, fnNames{ii}{jj}) ;
        % subsample the quiver
        % vxsub = imresize(vx, 1/subf);
        % vysub = imresize(vy, 1/subf) ;
        % Xsub = imresize(XX, 1/subf) ;
        % Ysub = imresize(YY, 1/subf) ;
        % scalarVectorFieldsOnImage(fliplr(im1')', uY, uX, speed', ...
        %     Ysub, Xsub, vysub(:), vxsub(:), opts) ;
        
        dmyk = dmyk + 1 ;
    end
end

