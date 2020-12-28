function aux_plot_strain_pattern(V, eIDx, eL0, distribution,...
    strainstyle, strainStruct, cmin, cmax, outfn, Rad)


midx = 0.5 * (V(eIDx(:, 1), 1) + V(eIDx(:, 2), 1)) ;
midy = 0.5 * (V(eIDx(:, 1), 2) + V(eIDx(:, 2), 2)) ;
midz = 0.5 * (V(eIDx(:, 1), 3) + V(eIDx(:, 2), 3)) ;
midpt = [midx(:), midy(:), midz(:)] ;

beta = strainStruct.beta ;

switch distribution
    case 'gaussian'
        mag = exp(- midx.^2/ (2 * strainStruct.sigma^2)) ;
end

switch strainstyle
    case 'isopulse'
        scalefactor = 1 - mag ;  
    case 'hoopstrain'
        scalefactor = 1 - mag .* abs(sin(beta)) ;
    case 'halfhoop'
        % target strain is zero at midplane and above, varies to mag at z=-Rad
        scalefactor = 1 + mag .* abs(sin(beta)) .* (midy/Rad);
        scalefactor(scalefactor > 1) = 1 ;
end
eL = eL0 .* scalefactor ;
estrain = (eL - eL0) ./ eL0 ;
cmap = bwr ;
colormap bwr
cID = max(1, sum(estrain > linspace(cmin, cmax, length(cmap)), 2)) ;
ecolors = cmap(cID, :) ;
figure('visible', 'off')
plotColoredLinesegs([V(eIDx(:,1), :), V(eIDx(:, 2), :)], ecolors, ...
    'linewidth', 10^3 / length(V)) ;
drawnow
c = colorbar ;
colormap bwr
caxis([cmin, cmax])
c.Color = 'w' ;
c.Label.Interpreter = 'latex' ;
c.Label.String = '$\varepsilon$' ;
% Figure properties
set(gca, 'color', 'k', 'xcol', 'w', 'ycol', 'w')
set(gcf, 'color', 'k')
title('bond strain pattern', 'interpreter', 'latex', 'color', 'w'); 
axis equal
axis off
drawnow
% save figure
export_fig(outfn, '-r300', '-nocrop') ;




%%%%%%%