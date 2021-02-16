function dilation = jacobian2dilation(jac)
% dilation = jacobian2dilation(jac)
% Given an MxN jacobian, return the dilation factor 
% todo: generalize to D-dim to D'-dim mappings.
% 
% Parameters
% ----------
% jac : 2 x 3 float array
%   jacobian matrix taking 3d vector to 2d plane according to some mapping
%
% Returns
% -------
% dilation : float
%
% Example Usage
% -------------
% >> jac = jacobian3Dto2DMesh(v2d, v3d, faces) ;
% >> jac{1}
% ans = 
%    -0.3584    0.0816    2.4829
%    -0.0001   -0.0064    0.0013
%
% >> dilation = zeros(size(faces, 1), 1) ;
% >> for f = 1:length(faces)
% >>     qg = jac{faces(f)} * jac{faces(f)}' ;
% >>     g_ab(f, :, :) =  qg ;
% >>     dilation(f) = sqrt(det(qg)) ;
% >> end
% 
%
% NPMitchell 2021

qg = jac * jac' ;
dilation = sqrt(det(qg)) ;