function aux_plot_beta_pattern(beta, V, eIDx, outfn)

cmin_beta = 0; cmax_beta = 1 ;
cmap = cubehelix(128,1.23,2.98,1.35,1.77,[0.17,0.98],[0.96,0.51]) ; 
% cmap = cubehelix(128,1.23,1.36,1.35,1.77,[0.14,0.96],[0.96,0.27]) ; 
cID = max(1, sum(abs(sin(beta)) > linspace(cmin_beta, cmax_beta, length(cmap)), 2)) ;
ecolors = cmap(cID, :) ;
figure('visible', 'off')
plotColoredLinesegs([V(eIDx(:,1), :), V(eIDx(:, 2), :)], ecolors, ...
    'linewidth', 10^3 / length(V)) ;
c = colorbar ;
colormap(cmap)
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