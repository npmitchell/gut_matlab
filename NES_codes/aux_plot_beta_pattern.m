function aux_plot_beta_pattern(beta, V, eIDx, outfn)

cmin_beta = 0; cmax_beta = 1 ;
cmap = jet ;
cID = max(1, sum(abs(sin(beta)) > linspace(cmin_beta, cmax_beta, length(cmap)), 2)) ;
ecolors = cmap(cID, :) ;
figure('visible', 'off')
plotColoredLinesegs([V(eIDx(:,1), :), V(eIDx(:, 2), :)], ecolors, ...
    'linewidth', 10^3 / length(V)) ;
c = colorbar ;
colormap jet
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
drawnow
% save figure
export_fig(outfn, '-r300', '-nocrop') ;