function [zmeans, zs] = binData2dGrid(uvz, uminmax, vminmax, nU, nV)
% Bin 3d data into 2d grid and take means
%
%
% NPMitchell 2020 adaptation of Walter Robinson from reference: 
% https://www.mathworks.com/matlabcentral/answers/
%       322113-binning-data-with-2d-coordinates

u_vals = linspace(uminmax(1), uminmax(2), nU) ;
v_vals = linspace(vminmax(1), vminmax(2), nV) ;

% Round uvz(:, [1,2]) onto u_vals, v_vals
urounded = uvz(:, 1) - uminmax(1)

% Now find indices of the rounded data to accumulate
[uu, ~, xidx] = unique(urounded);
[vv, ~, yidx] = unique(vrounded);

%count the number of points at each unique x/y combination
counts = accumarray([xidx(:), yidx(:)], 1);  
%average the z that fall into each unique x/y combination
sums = accumarray([xidx(:), yidx(:)], uvz(:,3).');
zmeans = sums / counts ;

if nargout > 1
    %create a list of the z that fall into each unique x/y combination
    zs = accumarray([xidx(:), yidx(:)], uvz(:,3).', [], @(V) {V}, {});
end
