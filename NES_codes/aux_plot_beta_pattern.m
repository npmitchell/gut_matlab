function aux_plot_beta_pattern(betaCut, Vc, eIDxCut, outfn)

cmin_beta = 0; cmax_beta = 1 ;
cmap = parula ;
colormap parula
cID = max(1, sum(abs(sin(betaCut)) > linspace(cmin_beta, cmax_beta, length(cmap)), 2)) ;
ecolors = cmap(cID, :) ;
figure('visible', 'off')
plotColoredLinesegs([Vc(eIDxCut(:,1), :), Vc(eIDxCut(:, 2), :)], ecolors, ...
    'linewidth', 10^3 / length(Vc)) ;
c = colorbar ;
caxis([cmin_beta, cmax_beta])
c.Color = 'w' ;
c.Label.Interpreter = 'latex' ;
c.Label.String = '$|\sin(\beta)|$' ;
% Figure properties
set(gca, 'color', 'k', 'xcol', 'w', 'ycol', 'w')
set(gcf, 'color', 'k')
title('axial angle $\beta$ definition', 'interpreter', 'latex', 'color', 'w'); 
axis equal
axis off
% save figure
export_fig(outfn, '-r300', '-nocrop') ;