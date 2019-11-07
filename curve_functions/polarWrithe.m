function [wr, wr_local, wr_nonlocal, turns] = polarWrithe(xyz, ss)
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
% Handle the case where there is a stagnation point, ex sgnz=[..1, 0, -1..]
% by setting all negative signs to zero, so that now sgnz=[..1, 0, 0...]
sgnz(sgnz < 0) = 0 ;
turns = find(abs(diff(sgnz)) > 0) + 1;

if isempty(turns)
    disp('One single segment detected')
    % Compute the local writhe for each segment
    wr_local = local_writhe(ss, xyz) ;
    wr_nonlocal = [] ;
else
    disp([ num2str(length(turns) + 1) ' segments detected'])
    % Compute the local writhe for each segment
    wr_local = local_writhe(ss, xyz) ;
    wr_nonlocal = [] ;
    % Compute the nonlocal contribution of the polar writhe 
    % add contribution from zmin to zmax
    
    % Segments are defined 1:turns(1)-1, turns(1):turns(2)-1, etc
    % Compute sigmas, which tell which direction curve is trending
    sigma = zeros(length(turns) + 1, 1) ;
    for ii = 1:length(turns)+1
        if ii == 1
            segii = 1:turns(ii) - 1 ;
        elseif ii < length(turns) + 1
            segii = turns(ii - 1):turns(ii) - 1 ;
        else
            segii = turns(ii - 1):size(xyz, 1) ;
        end
        [~, minID] = min(xyz(segii, 3)) ;
        [~, maxID] = max(xyz(segii, 3)) ;
        sigma(ii) = (minID < maxID) * 2 - 1 ;
    end
    
    % Now compute the twisting rate of the vector joining each segment to
    % every other
    for ii = 1:length(turns) + 1
        disp(['considering segment ' num2str(ii)])
        
        % obtain segment indices of ii 
        if ii == 1
            segii = 1:turns(ii) - 1 ;
        elseif ii < length(turns) + 1
            segii = turns(ii - 1):turns(ii) - 1 ;
        else
            segii = turns(ii - 1):size(xyz, 1) ;
        end
        
        [zmini, minIDi] = min(xyz(segii, 3)) ;
        [zmaxi, maxIDi] = max(xyz(segii, 3)) ;
        for jj = 1:length(turns) + 1
            disp(['comparing segment ' num2str(ii) ' to segment ' num2str(jj)])
            % skip self-energy terms
            if ii ~= jj
                % obtain segment indices of jj
                if jj == 1
                    segjj = 1:turns(jj) - 1 ;
                elseif jj < length(turns) + 1
                    segjj = turns(jj - 1):turns(jj) - 1 ;
                else
                    segjj = turns(jj - 1):size(xyz, 1) ;
                end
                
                % disp('segjj = ')
                % disp(segjj)
                
                % find the z range over which to interpolate
                [zminj, minIDj] = min(xyz(segjj, 3)) ;
                [zmaxj, maxIDj] = max(xyz(segjj, 3)) ;

                zmin = max(zmini, zminj) ;
                zmax = min(zmaxi, zmaxj) ;
                
                % draw rvec pointing from ii to jj curve along (zmin, zmax)
                % Interpolate segment i
                dzz = (zmax - zmin) / 100 ;
                tt = zmin:dzz:zmax ;
                segi_x = interp1(xyz(segii, 3), xyz(segii, 1), tt)';
                segi_y = interp1(xyz(segii, 3), xyz(segii, 2), tt)';
                
                % Interpolate segment j
                segj_x = interp1(xyz(segjj, 3), xyz(segjj, 1), tt)';
                segj_y = interp1(xyz(segjj, 3), xyz(segjj, 2), tt)';
                
                % Compute the twisting of i around j
                rij = [segj_x - segi_x, segj_y - segi_y, zeros(size(segj_x))] ;
                rijprime = [gradient(rij(:, 1), dzz), gradient(rij(:, 2), dzz), gradient(rij(:, 3), dzz)] ;
                numprod = cross(rij, rijprime)  ;
                dTheta_dz = numprod(:, 3) ./ vecnorm(rij, 2, 2).^2;
                
                contrib = sigma(ii) * sigma(jj) * nansum(dTheta_dz .* dzz) ;
                wr_nonlocal = [wr_nonlocal, contrib] ;
            end
        end 
    end
    for ii = 1:length(wr_nonlocal)
        % normalize the nonlocal contribs by 2pi
        wr_nonlocal(ii) = wr_nonlocal(ii) / (2 * pi) ;
    end
end

disp('wr_nonlocal = ')
disp(wr_nonlocal)
wr = nansum(wr_local .* abs(dz)) + nansum(wr_nonlocal) ;

end

