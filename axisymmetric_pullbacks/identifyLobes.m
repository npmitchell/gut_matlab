function [folds, ssfold, ssfold_frac, ssmax, rmax, fold_onset] = ...
    identifyLobes(timePoints, sphiBase, guess123, max_wander,...
    visualize, method, first_tp_allowed)
%idenitfyLobes(timePoints, sphiBase, guess1, guess2, guess3, visualize) Find folds in meshes 
%   Load each spcutMesh, find the local minima in radius, and mark these as 
%   fold locations. Track the location of those local minima over time, and
%   mark fold locations before the appearance of a fold as the location
%   where it appears. The pathlength is computed from the mean centerline's
%   pathlength ('mcline'), not ringpath_ss, which uses the surface
%   pathlength.
%
% Parameters
% ----------
% timePoints : N x 1 float/int array
%   The timestamps of each file to load (1 per minute is assumed)
% sphiBase : string
% guess123 : array of 3 floats, each between 0 and 1
%   Guesses for fractional position along U of each fold
% max_wander : float
%   maximum distance a fold location can wander in a given timepoint once 
%   identified, in units of pathlength (ss)
% visualize : bool
% method : ('ringpath' or 'avgpts')
%   which method to use to find local minima, 'avgpts' recommended
% first_tp_allowed : int, default=-1
%   First timepoint in which we allow detected local minumum to be 
%   considered a true, physical fold
%
% Returns
% -------
% fold_onset : #folds x 1 float
%   timestamps (not indices) of fold onset
% folds : #timepoints x #folds int
%   indices of nU sampling of folds
% ssmax : #timepoints x 1 float
%   maximum length of the centerline for each timepoint
% ssfold : #timepoints x #folds float
%   positional pathlength along centerline of folds
% rmax : #timepoints x 1 float
%   maximum radius of any part of the surface at each timepoint
%         
%
%
% NPMitchell 2019

first1 = true ;
first2 = true ;
first3 = true ;
% pathlength preallocations
ssfold = zeros(length(timePoints) + 1, 3) ;
ssfold_frac = zeros(length(timePoints) + 1, 3) ;
ssmax = zeros(length(timePoints), 1) ;
    
% Convert method to a boolean
if strcmp(method, 'avgpts')
    method = true ;
elseif strcmp(method, 'ringpath')
    method = false ;
end

% Store maximum radius if desired as output
if nargout > 4
    rmax = zeros(length(timePoints), 1) ;
end

% Fold onset timepoint indices
k1 = 0 ;
k2 = 0 ;
k3 = 0 ;

% Consider each timepoint
for kk = 1:length(timePoints)
    % Translate to which timestamp
    t = timePoints(kk) ;
    load(sprintf(sphiBase, t), 'spcutMesh') ;
    
    % Extract centerline pathlength from DVhoop mean curve
    if method
        ss = spcutMesh.avgpts_ss ;
    else
        ss = spcutMesh.ringpath_ss ;
    end
    
    maxss = max(ss) ;
    ssmax(kk, :) = maxss ;
    
    if kk == 1
        % Initialize the fold positions to be approximately .3, .5, .8
        folds = ones(length(timePoints) + 1, 3) ;
        folds = [guess123(1) * length(spcutMesh.phi0s), ...
                        guess123(2) * length(spcutMesh.phi0s), ...
                        guess123(3) * length(spcutMesh.phi0s)] .* folds;
    end
    
    % Find the minima in the radius. First make radius 1d
    rad = mean(spcutMesh.radii_from_mean_uniform_rs, 2) ;
    minidx = islocalmin(rad) ;
    if any(minidx) && (t > (first_tp_allowed - 1))
        minidx = find(minidx) ;
        % Identify each minimum which fold it may be
        dd = zeros(length(minidx), 1) ;
        which_fold = zeros(length(minidx), 1) ;
        for jj = 1:length(minidx)
            [dd(jj), which_fold(jj)] = min(abs(minidx(jj) - folds(kk, :))) ;
        end
        
        % If there is an increase of two folds, name them separately
        % if (length(minidx) - nnz(folds)) > 1 && nnz(folds) < 2
        %    % There are now at least two more folds so two are assigned to 
        %    % the same fold index: assigned to nnz(folds) + 1.
        %    % If there are more than two more new assignments, this is no
        %    % good
        % ALTERNATIVE:
        % We instead say 0.3 is fold 1, 0.5 is fold 2, 0.8 is fold 3
        
        % Mark if this is the first appearance of a fold
        if first1
            if any(which_fold == 1)
                first1 = false ;
                k1 = kk ;
            end
        end
        if first2
            if any(which_fold == 2)
                first2 = false ;
                k2 = kk ;
            end
        end
        if first3 
            if any(which_fold == 3)
                first3 = false ;
                k3 = kk ;
            end
        end
        
        % Sort indices to folds
        if length(minidx) > 3
            disp('there are more minima than folds. Sorting...')
            % Find closest fold candidate for each fold (1,2,3)
            for pp = 1:3
                if length(find(which_fold == pp)) > 1
                    disp(['more than one fold #' num2str(pp) ' found. Choosing closest....'])
                    [~, closest_idx] = min(dd(which_fold == 1)) ;
                    options = find(which_fold ==pp) ;
                    choice1 = options(closest_idx) ;
                    to_remove = setdiff(options, choice1) ;
                    to_keep = setdiff(1:length(minidx), to_remove) ;
                    % Remove the further matches
                    minidx = minidx(to_keep) ;
                    which_fold = which_fold(to_keep) ;
                    dd = dd(to_keep) ;
                end
            end
            % 
            % % fold 2
            % if length(find(which_fold == 2)) > 1
            %     disp('more than one fold #2 found. Choosing....')
            %     [~, which2] = min(dd(which_fold == 2)) ;
            %     options = find(which_fold == 2) ;
            %     choice2 = options(which2) ;
            %     to_remove = setdiff(options, choice2) ;
            %     to_keep = setdiff(1:length(minidx), to_remove) ;
            %     % Remove the further matches
            %     minidx = minidx(to_keep) ;
            %     which_fold = which_fold(to_keep) ;
            %     dd = dd(to_keep) ;
            % end
            % 
            % % fold 3
            % if length(find(which_fold == 3)) > 1
            %     disp('more than one fold #3 found. Choosing....')
            %     [~, which3] = min(dd(which_fold == 3)) ;
            %     options = find(which_fold == 3) ;
            %     choice3 = options(which3) ;
            %     to_remove = setdiff(options, choice3) ;
            %     to_keep = setdiff(1:length(minidx), to_remove) ;
            %     % Remove the further matches
            %     minidx = minidx(to_keep) ;
            %     which_fold = which_fold(to_keep) ;
            %     dd = dd(to_keep) ;
            % end
        end
        folds_kk = folds(kk, :) ;
        folds_kk(which_fold) = minidx ;
        
        % Ensure that no folds have wandered farther than max_wander unless
        % this is the first appearance of that fold
        firsts = [first1, first2, first3] ;
        onset_kks = [k1, k2, k3] ;
        for pp = 1:3
            if ~firsts(pp)
                % If the fold has already appeared, consider its motion
                if kk > onset_kks(pp)
                    % not first appearance, so check distance
                    if abs(ss(folds_kk(pp)) - ssfold(kk, pp)) > max_wander
                        % too far, mark as unchanged
                        folds_kk(pp) = folds(kk, pp) ;
                    end
                end
            end
        end
        
        % Append the updated fold positions
        folds(kk + 1, :) = folds_kk ;
        ssfold(kk + 1, :) = ss(folds_kk) ;
        ssfold_frac(kk + 1, :) = ss(folds_kk) / maxss ; 
    end
    
    % Store maximum radius if desired as output
    if nargout > 4
        rmax(kk) = max(spcutMesh.radii_from_mean_uniform_rs(:)) ;
    end

    % Plot the indices of the identified folds
    if visualize
        clf;
        plot(folds(:, 1)); hold on;
        plot(folds(:, 2)); hold on;
        plot(folds(:, 3)); hold on;
        title('Fold locations')
        xlabel('time [min]')
        pause(0.0001)
    end
end
% Truncate extra first entry which was the initial guess
folds = folds(2:end, :) ;
ssfold = ssfold(2:end, :) ;
ssfold_frac = ssfold_frac(2:end, :) ;

% Identify location prior to appearance of each fold as first location
folds(1:(k1-1), 1) = folds(k1, 1) ;
folds(1:(k2-1), 2) = folds(k2, 2) ;
folds(1:(k3-1), 3) = folds(k3, 3) ;
ssfold_frac(1:(k1-1), 1) = ssfold_frac(k1, 1) ;
ssfold_frac(1:(k2-1), 2) = ssfold_frac(k2, 2) ;
ssfold_frac(1:(k3-1), 3) = ssfold_frac(k3, 3) ;

% Go back and find the ss (pathlength) value for pre-fold indices
for qq = 1:max([k1, k2, k3])
    t = timePoints(qq) ;
    load(sprintf(sphiBase, t), 'spcutMesh') ;
    if method
        ss = spcutMesh.avgpts_ss ;
    else
        ss = spcutMesh.ringpath_ss ;
    end
    
    % Assign the correct ss value to each
    if qq < k1 
        assert(ss(folds(qq, 1)) > 0)
        ssfold(qq, 1) = ss(folds(qq, 1)) ;
        assert(ssfold(qq, 1) ~= 0)
    end
    if qq < k2 
        assert(ss(folds(qq, 2)) > 0)
        ssfold(qq, 2) = ss(folds(qq, 2)) ;
        assert(ssfold(qq, 2) ~= 0)
    end
    if qq < k3 
        assert(ss(folds(qq, 3)) > 0)
        ssfold(qq, 3) = ss(folds(qq, 3)) ;
        assert(ssfold(qq, 3) ~= 0)
    end
end

fold_onset = timePoints([k1, k2, k3]) ;

end

