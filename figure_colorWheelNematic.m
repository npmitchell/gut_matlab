% figure_colorWheelNematic


% aux_paths_and_colors
xx = linspace(0, pi) ;
[mag, theta] = meshgrid(xx, xx) ;

close all
cbs{1} = phasebar('colormap', phasemap, ...
    'location', [0.1, 0.01, 0.8, 0.8], ...
    'style', 'nematicFillGradient') ;
axis off
axis equal
set(gcf, 'color', 'w')
        
print('/mnt/data/analysis/colorNematicWheel_white.png', '-dpng', '-r300')

set(gcf, 'color', 'k')
print('/mnt/data/analysis/colorNematicWheel_black.png', '-dpng', '-r300')
export_fig('/mnt/data/analysis/colorNematicWheel_black.png', '-dpng', '-r300')


