function aux_plot_strain_pattern(V, eIDx, eL0, distribution,...
    strainstyle, strainStruct, cmin, cmax, outfn)


midpts = 0.5 * (V(eIDx(:, 1), 1) + V(eIDx(:, 2), 1)) ;
% midy = 0.5 * (V(eIDx(:, 1), 2) + V(eIDx(:, 2), 2)) ;
beta = strainStruct.beta ;

switch distribution
    case 'gaussian'
        mag = exp(- midpts.^2/ (2 * strainStruct.sigma^2)) ;
end

switch strainstyle
    case 'isopulse'
        scalefactor = 1 - mag ;  
    case 'hoopstrain'
        disp('strain style is hoopstrain')
        scalefactor = 1 - mag .* abs(sin(beta)) ;
        % scatter(midpts + rand(size(midpts)), scalefactor, 50, beta, 'filled')
        % colorbar
        % waitfor(gcf)
    case 'halfhoop'
        scalefactor = 1 - mag .* abs(sin(beta)) ;
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
c = colorbar ;
caxis([cmin, cmax])
c.Color = 'w' ;
c.Label.Interpreter = 'latex' ;
c.Label.String = '$\epsilon$' ;
% Figure properties
set(gca, 'color', 'k', 'xcol', 'w', 'ycol', 'w')
set(gcf, 'color', 'k')
title('bond strain pattern', 'interpreter', 'latex', 'color', 'w'); 
axis equal
view(2);
axis off
% save figure
export_fig(outfn, '-r300', '-nocrop') ;
%%%%%%%