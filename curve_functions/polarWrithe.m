function [wr, wr_local, wr_global] = polarWrithe(ss, xyz)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here


% Divide the curve into segments
dz = gradient(xyz(:, 3)) ;
ds = gradient(ss) ;
dzds = dz ./ ds ;
% where does this change sign?
sgnz = sign(dzds) ;
turns = find(abs(gradient(sgnz)) > 0) ;

if isempty(turns)
    disp('One segment detected')
    % Compute the local writhe for each segment
    wr_local = local_writhe(ss, xyz) ;
else
    disp([ num2str(length(turns) + 1) ' segments detected'])
    % Compute the local writhe for each segment
    wr_local = local_writhe(ss, xyz) ;
    
    % Compute the nonlocal contribution of the polar writhe 
end



end

