function [skel, DT, DD] = extractCurve(bw, startpt, endpt, options)
% [skel, DT, DD] = extractCurve(bw, startpt, endpt, options)
%   Extract a 1d curve from a black & white segmented 2D image using the
%   distance transform.
%
% Parameters
% ----------
% bw : NxM bool array
% startpt : 1x2 float array
% endpt : 1x2 float array
% options : struct with fields
%   exponent
%   eps
%   preview
%   weight
%   DD
%   DT
%   
% Returns
% -------
% skel : 
% DT : size(bw) int array
%   distance transform of bw
% DD : size(bw) float array
%   normalized distance transform exponentiated to the power of
%   options.exponent used for curve extraction via geodesic marching

%% Options
exponent = 1;
eps = 1e-7 ;
preview = false ;
weight = 0.1 ;
DD = [] ;
DT = [] ;

if nargin > 3
    if isfield(options, 'exponent')
        exponent = options.exponent ;
    end    
    if isfield(options, 'eps')
        eps = options.eps ;
    end    
    if isfield(options, 'weight')
        weight = options.weight ;
    end    
    if isfield(options, 'preview')
        preview = options.preview ;
    end    
    if isfield(options, 'preview')
        preview = options.preview ;
    end    
    if isfield(options, 'DT')
        DT = options.DT ;
    end    
    if isfield(options, 'DD')
        DD = options.DD ;
    end    
end
 
%% Check that segmentation is a single solid
% props = regionprops(bw, 'Area') ;
% ndilate = 0 ;
% while length(props.Volume) > 1
%     disp('dilating inside volume by two voxels')
%     [xb,yb,zb] = ndgrid(-4:4);
%     se = strel(sqrt(xb.^2 + yb.^2 + zb.^2) <=4);
%     inside = imdilate(inside, se) ;
%     props = regionprops(inside, 'Area') ;
%     ndilate = ndilate + 1;
%     if ndilate > 10
%         error('could not dilate array to connect components')
%     end
% end

%% use the distanceTransform from Yuriy Mishchenko
if isempty(DD)
    if isempty(DT)
        disp('Computing DT') ;
        outside = 1 - bw ;
        DT = bwdistsc(outside) ;
    end

    % DD = max(DD(:)) - DD ;
    DD = (DT + eps) ./ (max(DT(:)) + eps) ;
    % DD = 1 - DD ;
    DD = DD.^(exponent) ; 
    DD(logical(1 - bw)) = eps ;
    disp('> Computed normalized distance transform')
end

if preview
    % Preview DD
    close all ;
    disp('Previewing the distance transform')
    imagesc(DD')
    axis equal
    title('DT')
    colorbar
    pause(0.1)        
end

%% use Peyre's fast marcher
disp('> Computing centerline ');
tic
% From example (DD is W, with low values being avoided)
options.heuristic = weight * DD ;
% Convert here to the gridspacing of xx,yy,zz
% startpt_transposed = [startpt(2), startpt(1)]' ;
% endpt_transposed = [endpt(2), endpt(1)]' ;
[D2, S] = perform_fast_marching(DD, startpt', options);
path = compute_geodesic(D2, endpt');

% Show the intermediate result
disp('> Found skel via geodesic fast marching')  

% Convert skeleton's rows to columns and flip start/end
skel = fliplr(path)' ;

if preview
    clf
    imshow(bw'); hold on;
    plot(skel(:, 1), skel(:, 2), '.')
    pause(1)
end
