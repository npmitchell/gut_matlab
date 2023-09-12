function fillfig(ax)
% Set an axis to fill a figure fully

if nargin < 1
    ax = gca ;
end
drawnow;
InSet = get(ax, 'TightInset');
set(ax, 'Position', [InSet(1:2), 1-InSet(1)-InSet(3), 1-InSet(2)-InSet(4)])