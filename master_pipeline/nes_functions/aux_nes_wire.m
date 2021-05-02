function aux_nes_wire(QS, eIDx, VV, eL, eL0, clims, wirefn, titlestr)
% aux_nes_wire(QS, eIDx, VV, eL, eL0, clims, wirefn)
%
% NPMitchell 2021

if isempty(clims)
    cmin = min(eL ./ eL0) ;
    cmax = max(eL ./ eL0) ;
else
    cmin = clims(1) ;
    cmax = clims(2) ;
end
disp('Creating figure with wire frame bonds')
figure('visible', 'off')
colormap bwr
cmap = bwr ;
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 16 10]);  
sratio = eL ./ eL0 ;
cID = max(1, sum(sratio > linspace(cmin, cmax, length(cmap)), 2)) ;
ecolors = cmap(cID, :) ;
lsegs = [VV(eIDx(:,1), :), VV(eIDx(:, 2), :)] ;
plotColoredLinesegs(lsegs, ecolors) ;
c = colorbar ;
caxis([cmin, cmax])
c.Color = 'w' ;
c.Label.Interpreter = 'latex' ;

c.Label.String = '$\ell/\ell_0$' ;
% Figure properties
set(gca, 'color', 'k', 'xcol', 'w', 'ycol', 'w')
set(gcf, 'color', 'k')
title(titlestr, 'interpreter', 'latex', 'color', 'w'); 
axis equal
view(0,0)
axis off
% save figure
disp('Saving fig')
export_fig(wirefn, '-nocrop') ;

disp('done')
% set(gcf, 'visible', 'on')
% waitfor(gcf)    