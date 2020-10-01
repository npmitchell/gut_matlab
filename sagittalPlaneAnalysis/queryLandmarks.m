function landmarks = queryLandmarks(bw, options)
%landmarks = queryLandmarks(bw, options)
%
% Parameters
% ----------
% bw : NxM bool image
%   image on which to find landmarks
% options : struct with fields
%   reference_landmarks : #previous landmarks x 2 float array
%       locations of reference points to indicate when choosing new
%       landmarks
%
% Returns
% -------
% landmarks : #landmarks x 2 float array
%   location of landmarks in bw image
%
% NPMitchell 2020

% Initialize 
landmarks = cell(0); 
doneLmks = false ;
refs = [] ;
msz = 50 ;
titlestr = 'Identify landmarks' ;

% Unpack options
if nargin > 1
    if isfield(options, 'reference_landmarks')
        refs = options.reference_landmarks ;
    end
    if isfield(options, 'title')
        titlestr = options.title ;
    end
    if isfield(options, 'bgim')
        bgim = options.bgim;
    end
end

% Start populating landmarks
lmId = 1 ;
while ~doneLmks

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
    if ~isempty(bgim)
        imshow(bgim);  hold on;
        h = imshow(bw) ;
        set(h, 'AlphaData', 0.5 * double(bw(:, :, 1)));
    else
        imshow(bw)
    end
    if ~isempty(refs)
        if lmId <= length(refs) 
            colors = lines(size(refs{lmId}.v, 1)) ;
            scatter(refs{lmId}.v(:, 1), refs{lmId}.v(:, 2), msz, ...
                colors)
            refName = refs{lmId}.id ;
        else
            lastId = length(refs) ;
            colors = lines(size(refs{lastId}.v, 1)) ;
            scatter(refs{lastId}.v(:, 1), refs{lastId}.v(:, 2), msz, ...
                colors)
            refName = refs{length(refs)}.id ; 
        end
        colormap(colors)
        cb = colorbar() ;
        ylabel(cb, 'landmark ID')
    else
        refName = '' ;
    end
    title(titlestr) 

    %% Identify landmarks manually
    landmarks{lmId} = struct() ;
    disp('Identify a closed string of landmarks (counterclockwise), then <Enter>')
    [xi, yi] = getpts() ;
    if isempty(xi)
        % if no landmarks identified, preserve reference marks as current
        landmarks{lmId}.v = refs{lmId}.v ;
    else        
        landmarks{lmId}.v = [xi, yi] ;
    end
    landmarks{lmId}.v
    
    %% Display for classification
    colors = lines(size(landmarks{lmId}.v, 1)) ;
    hold on; 
    scatter(landmarks{lmId}.v(:, 1), ...
        landmarks{lmId}.v(:, 2), msz, colors, 'filled') 
    cb = colorbar() ;
    colormap(colors)
    ylabel(cb, 'landmark ID')
    title(titlestr) 
    
    %% Classify landmark
    % First grab reference names if given  
    classified = false ;
    msg = 'Name #LM+1 landmarks [a/v1/v2/v3/p/d3/d2/d1/a]' ;
    if ~isempty(refName)
        msg = [msg ' (default=' refName ')'] ;
    end
    landmarkName = input(msg, 's') ;
    % Keep previous names if none entered
    if isempty(landmarkName)
        landmarkName = refName ;
    end
    
    if any(contains(landmarkName, {'a', 'p', 'd1', 'd2', 'd3', 'v1', 'v2', 'v3', 'u'}))
        landmarks{lmId}.id = landmarkName ;
        classified = true ;
    else    
        if isempty(landmarkName)
            landmarks{lmId}.id = refName ;
        else
            error('Could not classify landmark')
        end
    end
    
    %% Decide on another or exit
    moreLmks = input('Done with landmarks? [Enter/N]', 's') ;
    if ~isempty(moreLmks)
        if strcmpi(moreLmks, 'n')
            doneLmks = false ;
        else
            doneLmks = true ;
        end
    else
        doneLmks = true ;
    end
    lmId = lmId + 1 ;
end