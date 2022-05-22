function [ssrM, minddssr_cc, minname_cc, minerror_cc, ...
    minweights_cc, ssr_minimum_cc ] = ...
    computeSSRMatrixMeshes(mca, clist, lastTP, rlist, options)
% [ssrM, minddssr_cc, minname_cc, minerror_cc, ...
%    minweights_cc, ssr_minimum_cc ] = ...
%    computeSSRMatrixMeshes(mca, clist, lastTP, rlist, refExptID, options)
% 
% Parameters
% ----------
% mca : cell of structs
%   mesh cell array, built as:
%       for jj = 1:ndatasets
%           meshes = dir(fullfile(alignedMeshDir, 'mesh_*_APDV_um.ply')) ;
%           for kk = 1:length(meshes)
%               mca{jj, kk} = meshes(kk);
%           end
%       end
% lastTP : int
% rlist : reference dataset index list
% options : struct with fields
%   ssrccDir : char path to output directory
%   overwrite : bool (default=false)
%       overwrite results on disk
%   cExptID : string
%       comparison dataset name, such as 'caax201902072000'
%   refExptID : string
%       reference dataset name, such as 'caax201902072000'
%
% Returns 
% -------
% ssrM : nTps_cc x length(rlist)
%   matrix of sum of squared Residuals
% minddssr_cc : nTps_cc x 1 float
%   minimum ssr index for each timepoint in comarison timeline cc
% minname_cc : nTps_cc x 1 cell;
%   reference timepoints that match each (cc, ii)
% minerror_cc : nTps_cc x 1 float
%   uncertainty in minimum, to use later in plotting
% minweights_cc : nTps_cc x 1 float
%   1/dix^2, for weighted statistics/fitting
% ssr_minimum_cc : nTps_cc x 1 float
%   ssr at minimum timestamp in comparison timeline cc
% 
% How to handle output
% --------------------
% minddssr(cc, ii) = rlist(ix);
% % Place reference timepoint that matches this (cc, ii)
% % tuple into minname
% minname{cc, ii} = mca{1, rlist(ix)}.name;
% minerror(cc, ii) = dix ; % store error to use later in plotting
% minweights(cc, ii) = 1/dix^2 ;
% ssr_minimum(cc,ii) = min(ssr) ;
% 
%
% NPMitchell 2020-2022

% Default options and unpack options struct
overwrite = false ;
rsubsampling = 1 ;
ssample_factor = 40 ;
ssfactorMedianThres = 8 ;
fitrange= 20 ;
t0ref = rlist(1) ;
if nargin < 5 || ~isfield(options, 'ssrccDir')
    error('must supply options with fields ssrccDir, cc (id for comparison in mca), and refID (id for reference in mca)')
else
    ssrccDir = options.ssrccDir ;
    cc = options.cc ;
    refID = options.refID ;
    cExptID = options.cExptID ;
    refExptID = options.refExptID ;
end

if isfield(options, 'rsubsampling')
    rsubsampling = options.rsubsampling ;
end
if isfield(options, 'ssample_factor')
    ssample_factor = options.ssample_factor ;
end
if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
end
if isfield(options, 'overwrite_figures')
    overwrite_figures = options.overwrite_figures ;
else
    overwrite_figures = overwrite ;
end
if isfield(options, 'ssfactorMedianThres')
    ssfactorMedianThres = options.ssfactorMedianThres ;
end
if isfield(options, 'fitrange')
    fitrange = options.fitrange ;
end
if isfield(options, 't0ref')
    t0ref = options.t0ref ;
end


% Preallocation
ssrM = [] ;
cid2do = clist(clist < (lastTP + 1)) ;
nTps_cc = length(cid2do) ;
minddssr_cc = zeros(nTps_cc, 1);
% Place reference timepoint that matches this (cc, ii)
% tuple into minname
minname_cc = cell(nTps_cc, 1);
minerror_cc = zeros(nTps_cc, 1) ; % store error to use later in plotting
minweights_cc = zeros(nTps_cc, 1) ;
ssr_minimum_cc = zeros(nTps_cc, 1) ;

% Consider each dataset TP (t_c) and match to reference mesh
for ii = cid2do % index of dataset thats being parsed against reference

    tic 

    % Consider the mesh only if mca{cc,ii} is populated with
    % struct that has field 'name'.
    if isfield(mca{cc, ii}, 'name') == 1

        % Check if SSR(cc, ii, rlist) has already been saved
        ssrii_fn = fullfile(ssrccDir, ...
            sprintf('ssr_ref%s_tp%04d_rsub%03d_ptsub%03d.mat',...
            refExptID, ii, rsubsampling, ssample_factor)) ;

        % Assume that we must redo calculation unless proven
        % that otherwise if file exists and overwrite == false
        % ---------------------------------------------------
        redoii = true ;
        if exist(ssrii_fn, 'file') && ~overwrite
            disp(['Loading ssr from: ' ssrii_fn])
            tmp = load(ssrii_fn, 'rlist', 'ssr', 'numpts1_cc', 'numpts2_cc') ;
            % Check that the rlist is indeed the same
            if length(tmp.rlist) == length(rlist)
                if all(tmp.rlist == rlist)
                    % Calculation exists for the same rlist for this ii
                    redoii = false ;
                    % assign ssr for this cc, ii
                    ssr = tmp.ssr ;
                    try
                        numpts1(cc,:) = tmp.numpts1_cc ;
                        numpts2(cc,:) = tmp.numpts2_cc ;
                    catch
                        disp('no numpts on disk!')
                    end
                end
            end
        end

        % (Re)Compute the SSR(cc, ii, rr) if we could not load ssr
        if redoii
            % prepare the dataset pointcloud
            cCloud = pcread(fullfile(mca{cc, ii}.folder, mca{cc, ii}.name)); % the mesh in question
            cxyz = cCloud.Location(1:ssample_factor:end,:) ;

            % preallocate ssr
            ssr = zeros(length(rlist), 1) ;
            tforms = zeros(length(rlist), 4, 4) ;

            % Consider each reference mesh to find match
            % Use dummy index qq to allow subsampled rlist
            for qq = 1:length(rlist)
                rr = rlist(qq) ;
                % why would ntimepoints_ref be necessary here?
                % if rr < ntimepoints_ref + 1
                if mod(qq, 50) == 0 || qq == 1
                    disp(['dataset_c=',sprintf('%d',cc), ...
                        ' timepoint_i=',sprintf('%d/%d',ii,lastTP), ...
                        ' against refTP_r=',sprintf('%d',rr), ...
                        ' of reference dataset=', sprintf('%s', refExptID)])
                end

                % Load the CAAX mesh used for comparison
                refCloud = pcread(fullfile(mca{refID, rr}.folder, mca{refID, rr}.name));
                refxyz = refCloud.Location;
                refxyz = refxyz(1:ssample_factor:end,:);
                tform = pcregistericp(cCloud,refCloud,'Extrapolate',true); % extracting the 4x4 transform from cmesh to refmesh
                cxyzpad = horzcat(cxyz, (zeros(length(cxyz),1)+1)); % pad the locations with a ones vector
                cxyztrans = cxyzpad*tform.T;    % apply the transform
                cxyztrans = cxyztrans(:,1:3);   % remove padding

                numpts1(cc,rr) = length(refxyz) ;
                numpts2(cc,rr) = length(cxyztrans) ;

                % find sum of squared residuals using two different ref orders
                % and finding geometric mean to penalize
                % dissimilarities in both directions of time
                indx1 = pointMatch(refxyz, cxyztrans);
                indx2 = pointMatch(cxyztrans, refxyz);

                % Normalize by the number of points considered in each
                % mesh separately to get true mean distance between
                % correspondence points
                ssr1 = sum(sum((cxyztrans(indx1,:)-refxyz).^2,2)) / length(indx1) ;
                ssr2 = sum(sum((refxyz(indx2,:)-cxyztrans).^2,2)) / length(indx2) ;
                % Geometric mean of the two SSRs
                ssr(qq) = sqrt(ssr1*ssr2);
                tforms(qq, :, :) = tform.T ;
                % end

            end
            % save matching curve for this cc tp == ii
            disp(['saving tforms: ' ssrii_fn])
            numpts1_cc = numpts1(cc, :) ;
            numpts2_cc = numpts2(cc, :) ;
            save(ssrii_fn, 'rlist', 'ssr', 'tforms', ...
                'numpts1_cc', 'numpts2_cc')
        end

        ssrM(ii, :) = ssr ;                

        % -------------------------------------------------
        % 1D FITTING TO SSR CURVE AS FUNCTION OF REFTIME
        % -------------------------------------------------
        % Apply median filter if we have a fine enough sampling
        if rsubsampling < ssfactorMedianThres
            ssr = movmedian(ssr, 11);
        end
        [minssr_value,ix] = min(ssr);

        %find polyfit of min to find error of match
        if ix + fitrange > length(ssr)
            endmax = length(ssr) ;
        else
            endmax = ix + fitrange ;
        end
        if ix > fitrange && (ix + fitrange) <= length(rlist)
            xpoly = rlist(ix-fitrange:endmax) ;
            ypoly = ssr(ix-fitrange:endmax) ;
        elseif ix < fitrange
            xpoly = rlist(1:2*fitrange);
            ypoly = ssr(1:2*fitrange);
        elseif ix + fitrange > length(rlist)
            xpoly = rlist(length(rlist)-2*fitrange:end) ;
            ypoly = ssr(length(rlist)-2*fitrange:end) ;
        elseif ix == fitrange && (ix + fitrange) <= length(rlist)
            xpoly = rlist(ix-fitrange+1:endmax) ;
            ypoly = ssr(ix-fitrange+1:endmax) ;
        else
            error('could not assign xpoly & ypoly')
        end

        % Fit to parabola
        [f, S] = polyfix(xpoly', ypoly, 2, rlist(ix), min(ssr), rlist(ix), 0) ;
        [poly, delta] = polyval(f, xpoly, S) ;
        Rinv = inv(S.R) ;
        covm = (Rinv*Rinv')*S.normr^2/S.df ;
        vara = covm(1,1) ;
        varb = covm(2,2) ;
        varc = covm(3,3) ;
        dix = ( ix* 0.5 ) * sqrt(vara/f(1)^2 + varb/f(2)^2) ;
        disp(['Uncertainty = ' num2str(dix)])
        
        % store minimum index as minddssr
        minddssr_cc (ii) = rlist(ix);
        % Place reference timepoint that matches this (cc, ii)
        % tuple into minname
        minname_cc{ii} = mca{1, rlist(ix)}.name;
        minerror_cc(ii) = dix ; % store error to use later in plotting
        minweights_cc(ii) = 1/dix^2 ;
        ssr_minimum_cc(ii) = min(ssr) ;
        assert(abs(polyval(f, rlist(ix)) - min(ssr))<1e-7)

        % make figure invisible for speedup
        figDirOut = fullfile(ssrccDir, ['figs_' sprintf('rsub%03d',rsubsampling)]) ;
        if ~exist(figDirOut, 'dir')
            mkdir(figDirOut)
        end
        figoutfn = fullfile(figDirOut, ['ICPPlot2_Dataset_', ...
            sprintf('%s',cExptID),'_', ...
            sprintf('%s',refExptID),'_', ...
            sprintf('%03d',ii),'.png']) ;    
        if redoii || ~exist(figoutfn, 'file') || overwrite_figures
            fig = figure('visible', 'off') ;
            plot(rlist - t0ref, ssr, '.');
            hold on
            plot(xpoly - t0ref, ypoly, '-') ;
            plot(rlist(ix) - t0ref, minssr_value, 'o')
            plot(xpoly - t0ref,poly, 'color', 'magenta', 'linewidth', 1) ;
            ylim([min(ssr)-10 max(ssr)+10])
            title(['dataset ', cExptID, ' vs ', refExptID,...
                ': $t=$', sprintf('%03d', ii), '$\rightarrow$', ...
                sprintf('%03.1f', rlist(ix) - t0ref), '$\pm$', ...
                sprintf('%0.1f', dix), ' min'], ...
                'Interpreter', 'Latex')
            xlabel('Reference Timepoint [min]', 'Interpreter', 'Latex');
            labstr = 'Mismatch, $\langle \Sigma |\mathrm{min}(\vec{x}_c - \vec{x}_{\mathrm{ref}})|^2 \rangle$ [$\mu$m]$^2$' ;
            ylabel(labstr, 'Interpreter', 'Latex');
            disp(['Saving ' figoutfn]);
            saveas(fig, figoutfn);
            close all
        end
    end
    toc
end