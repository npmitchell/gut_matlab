function [wr, wr_local, wr_nonlocal] = polarWrithe(xyz, ss)
%POLARWRITHE Compute the polar writhe, its local and nonlocal components
%   Compute as defined by "The writhe of open and closed curves, by 
%   Mitchell A Berger and Chris Prior, 2006.
%
% Parameters
% ----------
% ss : optional, pathlength parameterization
% xyz :
%
% Returns
% -------
% wr : 
% wr_local :
% wr_global :
%
% Example Usage
% -------------
% t = 0:0.01:1 ;
% xx = sin(2*pi*t) ;
% yy = cos(2*pi*t) ;
% zz = 2 * pi * (t - 0.5) ;
% xx = [xx xx(1:end) ];
% yy = [yy fliplr(yy(1:end)) ];
% zz = [zz fliplr(zz(1:end)) ];
% xyz = [xx ; yy; zz]';
% [wr, wr_local, wr_nonlocal] = polarWrithe(xyz) ;
% 
% NPMitchell 2019

if length(ss) < 1
    % get distance increment
    ds = vecnorm(diff(xyz), 2, 2) ;
    % get pathlength at each skeleton point
    ss = [0; cumsum(ds)] ;
end

% Divide the curve into segments
dz = gradient(xyz(:, 3)) ;
ds = gradient(ss) ;
dzds = dz ./ ds ;
% where does this change sign?
sgnz = sign(dzds) ;
turns = find(abs(diff(sgnz)) > 0) ;

if isempty(turns)
    disp('One single segment detected')
    % Compute the local writhe for each segment
    wr_local = local_writhe(ss, xyz) ;
else
    disp([ num2str(length(turns) + 1) ' segments detected'])
    % Compute the local writhe for each segment
    wr_local = local_writhe(ss, xyz) ;
    
    % Compute the nonlocal contribution of the polar writhe 
    % add contribution from zmin to zmax
    
    % Segments are defined 1:turns(1)-1, turns(1):turns(2)-1, etc
    % Compute sigmas, which tell which direction curve is trending
    sigma = zeros(length(turns) + 1, 1) ;
    for ii = 1:length(turns)
        if ii == 1
            segii = 1:turns(ii) - 1 ;
        else
            segii = turns(ii - 1):turns(ii) - 1 ;
        end
        [~, minID] = min(xyz(segii), 3) ;
        [~, maxID] = max(xyz(segii), 3) ;
        sigma(ii) = (minID < maxID) * 2 - 1 ;
    end
    
    % Now compute the twisting rate of the vector joining each segment to
    % every other
    for ii = 1:length(turns)
        if ii == 1
            segii = 1:turns(ii) - 1 ;
        else
            segii = turns(ii - 1):turns(ii) - 1 ;
        end
        [zmini, minID] = min(xyz(segii), 3) ;
        [zmaxi, maxID] = max(xyz(segii), 3) ;
        for jj = 1:length(turns)
            if ii ~= jj
                otherseg = turns(jj - 1):turns(jj) - 1 ;
                zminj = min(xyz(otherseg), 3) ;
                zmaxj = max(xyz(otherseg), 3) ;

                zmin = max(zmini, zminj) ;
                zmax = min(zmaxi, zmaxj) ;
                
                % draw rvec pointing from ii to jj curve along (zmin, zmax)
                % Interpolate segment i
                segi_intp =  0;
                % Interpolate segment j
                segj_intp =  0;
                
                % Compute the twisting of i around j
                dTheta_dz = 0;
                
                contrib = sigma(ii) * sigma(jj) * sum(dTheta_dz .* dz) ;
                wr_nonlocal = wr_nonlocal + contrib ;
            end
        end    
    end
    
end



end

