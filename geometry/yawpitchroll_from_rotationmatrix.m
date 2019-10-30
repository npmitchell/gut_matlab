function [yaw, pitch, roll] = yawpitchroll_from_rotationmatrix(rot)
% YAWPITCHROLL_FROM_ROTATIONMATRIX find angles from rotation matrix
% 
%
% Parameters
% ----------
% rot : 3 x 3 rotation matrix
% (r11 r12 r13 )
% (r21 r22 r23 )
% (r31 r32 r33 )
%
% Returns
% -------
% yaw : float, counterclockwise rotation about the z axis
% pitch : float, counterclockwise rotation about the y axis
% roll : float, counterclockwise rotation about the x axis
%

r11 = rot(1) ;
r21 = rot(2) ;
r31 = rot(3) ;
r32 = rot(6) ;
r33 = rot(9) ;

if r11 == 0
    error('r11 cannot be zero')
elseif r33 == 0
    error('r33 cannot be zero')
end

yaw = atan2(r11, r21) ;
pitch = atan2(-r31, sqrt(r32^2 + r33^2)) ;
roll = atan2(r32, r33) ;
    
return


% Note a similar procedure was already online:
% function [ alpha, beta, gamma ] = ypr2abg( Q )
% % this function is used to get rotation angles from euler rotation matrix
% %% Coded by
% % Mohamed Mohamed El-Sayed Atyya
% % mohamed.atyya94@eng-st.cu.edu.eg
% %% INPUTS:
% % Q  : ypr rotation matrix
% %% OUTPUTS:
% % alpha      : rotation about x
% % beta        : rotation about y
% % gamma   : rotation about z
% % ---------------------------------------------------------------------------------------------------------------------------------------------------------
% alpha=atand(Q(1,2)/Q(1,1));
% if Q(1,2) < 0 && Q(1,1) > 0 && alpha <0
%     alpha=alpha+360;
% elseif Q(1,1) < 0 && Q(1,2) > 0 && alpha <0
%     alpha=alpha+180;
% elseif Q(1,2) < 0 && Q(1,1) < 0
%     alpha=alpha+180;
% end
% beta=asind(-Q(1,3));
% gamma=atand(Q(2,3)/Q(3,3));
% if Q(2,3) < 0 && Q(3,3) > 0 && gamma < 0
%     gamma=gamma+360;
% elseif Q(3,3) < 0 && Q(2,3) > 0 && gamma < 0
%     gamma=gamma+180;
% elseif Q(3,3) < 0 && Q(2,3) < 0
%     gamma=gamma+180;
% end
% end