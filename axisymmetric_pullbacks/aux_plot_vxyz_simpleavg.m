function aux_plot_vxyz_simpleavg(im, vsm_ii, xx, yy, time_ii, vscale, ...
    pivSimpleAvgImXDir, pivSimpleAvgImYDir, pivSimpleAvgImZDir)
% AUX_PLOT_VXYZ_SIMPLEAVG(vsm_ii, gridsz, time_ii, ...
%    pivSimpleAvgImXDir, pivSimpleAvgImYDir, pivSimpleAvgImZDir)
%   Auxiliary function for Generate_Axisymmetric_Pullbacks_Orbifold.m, to
%   be used after smoothing of meshes and simple averaging of velocity
%   fiber
%
% vsm_ii
% gridsz : 
% time_ii : time(i)
% pivSimpleAvgImXDir : 
% pivSimpleAvgImYDir : 
% pivSimpleAvgImZDir : 
% 
% 
% NPMitchell 2020

alphaVal = 0.5 ;
gridsz = [length(xx), length(yy)] ;
vx = reshape(vsm_ii(:, 1), gridsz) ;
vy = reshape(vsm_ii(:, 2), gridsz) ;
vz = reshape(vsm_ii(:, 3), gridsz) ;

% x axis
close all
fig = figure('units', 'normalized', ...
    'outerposition', [0 0 1 1], 'visible', 'off') ;
scalarFieldOnImage(im, xx, yy, vx, alphaVal, vscale, 'v_x [\mum/min]')
ylim([size(im, 2) * 0.25, size(im, 2) * 0.75])
axis off
saveas(gcf, fullfile(pivSimpleAvgImXDir, [sprintf('%04d', time_ii) '.png']))

% y axis
close all
fig = figure('units', 'normalized', ...
    'outerposition', [0 0 1 1], 'visible', 'off') ;
scalarFieldOnImage(im, xx, yy, vy, alphaVal, vscale, 'v_y [\mum/min]')
ylim([size(im, 2) * 0.25, size(im, 2) * 0.75])
axis off
saveas(gcf, fullfile(pivSimpleAvgImYDir, [sprintf('%04d', time_ii) '.png']))

% z axis
close all
fig = figure('units', 'normalized', ...
    'outerposition', [0 0 1 1], 'visible', 'off') ;
scalarFieldOnImage(im, xx, yy, vz, alphaVal, vscale, 'v_z [\mum/min]')
ylim([size(im, 2) * 0.25, size(im, 2) * 0.75])
axis off
saveas(gcf, fullfile(pivSimpleAvgImZDir, [sprintf('%04d', time_ii) '.png']))