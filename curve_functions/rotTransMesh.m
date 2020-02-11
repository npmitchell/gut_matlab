function [ssdiff] = rotTransMesh(vars, curmesh, refmesh)
% Given two meshes, find Euler angles and translation vector from first
% mesh to second mesh using point matching
%
% Inputs
% ------
% mesh1&2 = N x 3 matrx of vertex coords
%
% Parameters
% ----------
% vars : a vector containing the variables to be minimized
%
% vars1 : eul1
%   The first euler angle for rotation
% vars 2 : eul2
%   The second euler angle for rotation
% vars 3 : eul3
%   The third euler angle for rotation
% vars 4 : xoff
%   The x-component of the translation vector
% vars 5 : yoff
%   The y-component of the translation vector
% vars 6 : zoff
%   The z-component of the translation vector
% 
% Returns 
% -------
% ssdiff : float
%   the sum of the squared differences between mesh1 and mesh2 as a
%   function of a rotation and translation matrix

eul1 = vars(1); % euler angle 1
eul2 = vars(2); % euler angle 2
eul3 = vars(3); % euler angle 3
xoff = vars(4); % x-component of translational vector
yoff = vars(5); % y-component of translational vector
zoff = vars(6); % z-component of tranlational vector

% Creating Proper Euler transformation matrix (X1.Z2.X3)

Eul = zeros(4,4) ;
Eul(1,1) = cos(eul2);
Eul(1,2) = -cos(eul3)*sin(eul2);
Eul(1,3) = sin(eul2)*sin(eul3);
Eul(2,1) = cos(eul1)*sin(eul2);
Eul(2,2) = cos(eul1)*cos(eul2)*cos(eul3)-sin(eul1)*sin(eul3);
Eul(2,3) = -cos(eul3)*sin(eul1)-cos(eul1)*cos(eul2)*sin(eul3);
Eul(3,1) = sin(eul1)*sin(eul2);
Eul(3,2) = cos(eul1)*sin(eul3)+cos(eul2)*cos(eul3)*sin(eul1);
Eul(3,3) = cos(eul1)*cos(eul3)-cos(eul2)*sin(eul1)*sin(eul3);
Eul(1,4) = xoff;
Eul(2,4) = yoff;
Eul(3,4) = zoff;
Eul(4,4) = 1;

tranrot = inv(Eul);

% Apply transforms

curmesh = curmesh * tranrot ;

end

