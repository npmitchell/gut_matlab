function c = continueButtonOnFigure(fig, msg)
% continueButtonOnFigure(fig, msg)
%
% Parameters
% ----------
% fig : figure instance
%   figure on which to create button
% msg : str 
%   message for title and display/command window
%
% Returns
% -------
% c : uicontrol instance
%
% NPMitchell 2022

if nargin < 1
    fig = gcf ; 
end
if nargin < 2
    msg = 'press Continue button to move on' ;
end

title(msg)
disp(msg)
c = uicontrol('String','Continue','Callback','uiresume(gcf)');
uiwait(fig)
