function pathMasks = queryPathMasks(bw, landmarks, options)
%pathMasks = queryPathMasks(bw, options)
%
% Parameters
% ----------
% bw : NxM bool image
%   image on which to find landmarks
% options : struct with fields
%   reference_pathMasks : #previous pathMasks x 1 cell
%       reference masks for creating new masks
%   ref_priority : bool (default = true)
%       If true, use the reference mask if no polygon drawn. If no
%       reference mask, use previous polygon drawn.
%       If false, use the previous mask drawn if no polygon drawn. If no
%       previous polygon stored, use reference mask. 
%
% Returns
% -------
% pathMasks : (#lmkStrings)x1 cell of (#landmarks)x1 cell of NxM booleans
%   binary masks for path extraction between subsequent landmarks
%
% NPMitchell 2020

% Initialize 
pathMasks = cell(0); 
doneLmks = false ;
refs = {} ;
msz = 50 ;
titlestr = 'Identify path masks: polygon + <Enter>, or <Esc> for no mask' ;
ref_priority = true ;
default_required = false ;
forceSelection = false ;

% Unpack options
if nargin > 1
    if isfield(options, 'reference_pathMasks')
        refs = options.reference_pathMasks ;
    end
    if isfield(options, 'title')
        titlestr = options.title ;
    end
    if isfield(options, 'bgim')
        bgim = options.bgim;
    end
    if isfield(options, 'ref_priority')
        ref_priority = options.ref_priority;
    end
    if isfield(options, 'default_required')
        default_required = options.default_required;
    end
    if isfield(options, 'promptString')
        promptstr = options.promptString;
    elseif isfield(options, 'prompt')
        promptstr = options.prompt;
    else
        if default_required
            promptstr = 'Are pathMasks required? [n/Y]' ;
        else
            promptstr = 'Are pathMasks required? [N/y]' ;
        end
    end
    if isfield(options, 'forceSelection')
        forceSelection = options.forceSelection ;
    end
end

close all
set(gcf, 'visible', 'on')
    
%% Show image
% Convert to RGB if needed
if length(size(bw)) == 2 || size(bw, 3) == 1
    bw2 = zeros(size(bw, 1), size(bw, 2), 3) ;
    bw2(:, :, 1) = bw ;
    bw2(:, :, 2) = bw ;
    bw2(:, :, 3) = bw ;
    bw = bw2 ;
    clearvars bw2
end

% Show image for selection
if ~isempty(bgim)
    imshow(bgim);  hold on;
    h = imshow(bw) ;
    set(h, 'AlphaData', 0.2 * double(bw(:, :, 1)));
else
    imshow(bw); hold on ;
end
title(titlestr) 

% Show all landmarks
for qq = 1:length(landmarks)
    colors = lines(size(landmarks{qq}.v, 1)) ;
    scatter(landmarks{qq}.v(:, 1), landmarks{qq}.v(:, 2), msz, colors)
end

%% Start populating pathMasks
% lmId is index for string of continuously connected landmarks
set(gcf, 'visible', 'on')
if forceSelection
    skip_mask = false ;
else
    do_masks = input(promptstr, 's') ;
    if default_required
        skip_mask = strcmpi(do_masks, 'n') ;
    else
        % default is to skip the mask
        skip_mask = strcmpi(do_masks, 'n') || isempty(do_masks) ;
    end
end

if skip_mask
    for lmId = 1:length(landmarks)
        %% Iterate over landmark pairs (paths)
        pathMasks{lmId} = cell(0) ;
        cycle = landmarks{lmId}.v ;
        for qq = 1:length(cycle) 
            pathMasks{lmId}.mask{qq} = [] ;
            pathMasks{lmId}.xroi{qq} = [] ;
            pathMasks{lmId}.yroi{qq} = [] ;
        end
    end
else
    for lmId = 1:length(landmarks)
        %% Iterate over landmark pairs (paths)
        pathMasks{lmId} = cell(0) ;
        cycle = landmarks{lmId}.v ;
        for qq = 1:length(cycle) 

            % Show image for selection
            if ~isempty(bgim)
                figure(2)
                imshow(bgim)
                figure(1)
                clf
                imshow(bgim);  hold on;
                h = imshow(bw) ;
                set(h, 'AlphaData', 0.2 * double(bw(:, :, 1)));
            else
                imshow(bw)
            end

            % Identify pairs of landmarks for curves
            lm2 = mod(qq + 1, size(cycle, 1)) ;
            if lm2 == 0
                lm2 = size(cycle, 1);
            end
            pair = landmarks{lmId}.v([qq, lm2], :) ;
            % plot landmark pairs
            colors = lines(2) ;
            scatter(pair(:,1), pair(:, 2), msz, colors)
            refName = landmarks{lmId}.id ;

            colormap(colors)
            cb = colorbar() ;
            ylabel(cb, 'landmark ID')
            title(titlestr) 

            % Show reference mask if available
            plot_previous_mask = false ;
            % Does reference mask take priority over previous mask?
            if ref_priority
                % Plot reference mask if exists, otherwise plot previous if
                % it exists
                if length(refs) >= lmId
                    disp('refs available')
                    if length(refs{lmId}.xroi) >= qq
                        disp('roi available')
                        if ~isempty(refs{lmId}.xroi{qq})
                            disp('plotting roi')
                            h1 = plot(refs{lmId}.xroi{qq}, refs{lmId}.yroi{qq}, '.-') ;
                        elseif qq > 1 && ~isempty(pathMasks{lmId}.xroi{qq-1})
                            plot_previous_mask = true ;
                        end
                    elseif qq > 1 && ~isempty(pathMasks{lmId}.xroi{qq-1})
                        plot_previous_mask = true ;
                    end
                elseif qq > 1 && ~isempty(pathMasks{lmId}.xroi{qq-1})
                    plot_previous_mask = true ;
                end
            else
                % Plot previous mask if it exists, otherwise plot reference
                if qq > 1 && ~isempty(pathMasks{lmId}.xroi{qq-1})
                    plot_previous_mask = true ;
                elseif length(refs) >= lmId
                    disp('refs available')
                    if length(refs{lmId}.xroi) >= qq
                        disp('roi available')
                        if ~isempty(refs{lmId}.xroi{qq})
                            disp('plotting roi')
                            h1 = plot(refs{lmId}.xroi{qq}, refs{lmId}.yroi{qq}, '.-') ;
                        end
                    end
                end
            end
            
            % Plot the previous mask if no reference mask available
            if plot_previous_mask
                h1 = plot(pathMasks{lmId}.xroi{qq-1}, ...
                    pathMasks{lmId}.yroi{qq-1}, '.-') ;
            end

            %% Identify masks manually
            disp('Identify a closed polygon for masking this landmark pair, then <Enter>, or <Esc> for no mask')
            [bwout, xroi, yroi] = roipoly() ;
            if ref_priority
                if isempty(bwout) && length(refs) >= lmId && length(refs{lmId}.mask) >= qq
                    % if no landmarks identified, preserve reference marks as current
                    pathMasks{lmId}.mask{qq} = refs{lmId}.mask{qq} ;
                    pathMasks{lmId}.xroi{qq} = refs{lmId}.xroi{qq} ;
                    pathMasks{lmId}.yroi{qq} = refs{lmId}.yroi{qq} ;
                elseif isempty(bwout) && plot_previous_mask
                    % if no landmarks identified, preserve reference marks as
                    % previously cropped
                    pathMasks{lmId}.mask{qq} = pathMasks{lmId}.mask{qq-1} ;
                    pathMasks{lmId}.xroi{qq} = pathMasks{lmId}.xroi{qq-1} ;
                    pathMasks{lmId}.yroi{qq} = pathMasks{lmId}.yroi{qq-1} ;
                else
                    pathMasks{lmId}.mask{qq} = bwout ;
                    pathMasks{lmId}.xroi{qq} = xroi ;
                    pathMasks{lmId}.yroi{qq} = yroi ;
                end
            else
                if isempty(bwout) && plot_previous_mask
                    % if no landmarks identified, preserve reference marks as
                    % previously cropped
                    pathMasks{lmId}.mask{qq} = pathMasks{lmId}.mask{qq-1} ;
                    pathMasks{lmId}.xroi{qq} = pathMasks{lmId}.xroi{qq-1} ;
                    pathMasks{lmId}.yroi{qq} = pathMasks{lmId}.yroi{qq-1} ;
                elseif isempty(bwout) && length(refs) >= lmId
                    % if no landmarks identified, preserve reference marks as current
                    pathMasks{lmId}.mask{qq} = refs{lmId}.mask{qq} ;
                    pathMasks{lmId}.xroi{qq} = refs{lmId}.xroi{qq} ;
                    pathMasks{lmId}.yroi{qq} = refs{lmId}.yroi{qq} ;
                else
                    pathMasks{lmId}.mask{qq} = bwout ;
                    pathMasks{lmId}.xroi{qq} = xroi ;
                    pathMasks{lmId}.yroi{qq} = yroi ;
                end
            end
            
            % plot it briefly
            clf
            imshow(pathMasks{lmId}.mask{qq})
            pause(0.1)             
        end
    end
end

