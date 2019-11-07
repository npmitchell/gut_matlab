% Demonstrate the polar Writhe definition on example curves
outdir = '~/Desktop/polarwrithe/' ;
if ~exist(outdir, 'dir')
    mkdir(outdir) ;
end
addpath('../plotting/')

%% First example: half twist, then twist back
t = 0:0.01:1 ;
xx = sin(2*pi*t) ;
yy = cos(2*pi*t) ;
zz = t ;
xx = [xx xx(1:end) ];
yy = [yy fliplr(yy(1:end)) ];
zz = [zz fliplr(zz(1:end)) ];
xyz = [xx ; yy; zz]';
[wr, wr_local, wr_nonlocal, turns] = polarWrithe(xyz, []) ;

subplot(2, 2, 1)
% scatter3(xx', yy', zz', 10, yy) 
cvals = 1:length(xx) ;
scatter3(xx', yy', zz', 10, cvals) 
xlabel('x')
ylabel('y')
zlabel('z')
title('Curve')

% Plot the segments colored
subplot(2,2,2)
seg1 = 1:turns(1) ;
seg2 = turns(1):length(xx) ;
plot3(xx(seg1), yy(seg1), zz(seg1), 's') ;
hold on
plot3(xx(seg2), yy(seg2), zz(seg2), 'o') ;
title('Segments')
% Plot the local writhe
subplot(2, 2, 3)
plot(zz, wr_local)
title('local writhe')
xlabel('position, z')
ylabel('local writhe')
% Plot the nonlocal writhe
subplot(2, 2, 4)
plot(1:length(wr_nonlocal), wr_nonlocal, 'o')
title('nonlocal writhe')
xlabel('segment index')
ylabel('nonlocal writhe')
axes( 'Position', [0, 0, 1, 1] ) ;
set( gca, 'Color', 'None', 'XColor', 'White', 'YColor', 'White' ) ;
titletext = ['Wr = ' num2str(wr)] ;
text( 0.5, 0.95, titletext, 'FontSize', 14', 'FontWeight', 'Bold', ...
      'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;
axis off
saveas(gcf, fullfile(outdir, 'twoseg_helix.png'))
close all

%% Half twist
t = 0:0.002 :1 ;
xx = sin(2*pi*t) ;
yy = cos(2*pi*t) ;
zz = t ;
xyz = [xx ; yy; zz]';
[wr, wr_local, wr_nonlocal, turns] = polarWrithe(xyz, []) ;

subplot(2, 2, 1)
% scatter3(xx', yy', zz', 10, yy) 
cvals = 1:length(xx) ;
scatter3(xx', yy', zz', 10, cvals) 
xlabel('x')
ylabel('y')
zlabel('z')
title('Curve')

% Plot the segments colored
subplot(2,2,2)
seg1 = 1:length(xyz) ;
plot3(xx(seg1), yy(seg1), zz(seg1), 's') ;
title('Segments')
% Plot the local writhe
subplot(2, 2, 3)
plot(zz, wr_local)
title('local writhe')
xlabel('position, z')
ylabel('local writhe')
% Plot the nonlocal writhe
subplot(2, 2, 4)
plot(1:length(wr_nonlocal), wr_nonlocal, 'o')
title('nonlocal writhe')
xlabel('segment index')
ylabel('nonlocal writhe')
axes( 'Position', [0, 0, 1, 1] ) ;
set( gca, 'Color', 'None', 'XColor', 'White', 'YColor', 'White' ) ;
titletext = ['Wr = ' num2str(wr)] ;
text( 0.5, 0.95, titletext, 'FontSize', 14', 'FontWeight', 'Bold', ...
      'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;
axis off
saveas(gcf, fullfile(outdir, 'oneseg_helix.png'))
close all


%% Half twist reverse
t = 0:0.01:1 ;
xx = sin(2*pi*t) ;
yy = cos(2*pi*t) ;
zz = t ;
xx = [xx ];
yy = [fliplr(yy) ];
zz = [fliplr(zz) ];
xyz = [xx ; yy; zz]';
[wr, wr_local, wr_nonlocal, turns] = polarWrithe(xyz, []) ;

subplot(2, 2, 1)
% scatter3(xx', yy', zz', 10, yy) 
cvals = 1:length(xx) ;
scatter3(xx', yy', zz', 10, cvals) 
xlabel('x')
ylabel('y')
zlabel('z')
title('Curve')

% Plot the segments colored
subplot(2,2,2)
seg1 = 1:length(xyz) ;
plot3(xx(seg1), yy(seg1), zz(seg1), 's') ;
title('Segments')
% Plot the local writhe
subplot(2, 2, 3)
plot(zz, wr_local)
title('local writhe')
xlabel('position, z')
ylabel('local writhe')
% Plot the nonlocal writhe
subplot(2, 2, 4)
plot(1:length(wr_nonlocal), wr_nonlocal, 'o')
title('nonlocal writhe')
xlabel('segment index')
ylabel('nonlocal writhe')
axes( 'Position', [0, 0, 1, 1] ) ;
set( gca, 'Color', 'None', 'XColor', 'White', 'YColor', 'White' ) ;
titletext = ['Wr = ' num2str(wr)] ;
text( 0.5, 0.95, titletext, 'FontSize', 14', 'FontWeight', 'Bold', ...
      'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;
axis off
saveas(gcf, fullfile(outdir, 'oneseg_helix_reverse.png'))
close all


%% half twist, then straight back
t = 0:0.01:1 ;
xx = sin(2*pi*t) ;
yy = cos(2*pi*t) ;
zz = t ;
xback = linspace(xx(end), xx(1), 10) ;
yback = linspace(yy(end), yy(1), 10) ;
zback = linspace(zz(end), zz(1), 10) ;
xx = [xx xback ];
yy = [yy yback ];
zz = [zz zback ];
xyz = [xx ; yy; zz]';
[wr, wr_local, wr_nonlocal, turns] = polarWrithe(xyz, []) ;

subplot(2, 2, 1)
% scatter3(xx', yy', zz', 10, yy) 
cvals = 1:length(xx) ;
scatter3(xx', yy', zz', 10, cvals) 
xlabel('x')
ylabel('y')
zlabel('z')
title('Curve')

% Plot the segments colored
subplot(2,2,2)
seg1 = 1:turns(1) ;
seg2 = turns(1):length(xx) ;
plot3(xx(seg1), yy(seg1), zz(seg1), 's') ;
hold on
plot3(xx(seg2), yy(seg2), zz(seg2), 'o') ;
title('Segments')
% Plot the local writhe
subplot(2, 2, 3)
plot(zz, wr_local)
title('local writhe')
xlabel('position, z')
ylabel('local writhe')
% Plot the nonlocal writhe
subplot(2, 2, 4)
plot(1:length(wr_nonlocal), wr_nonlocal, 'o')
title('nonlocal writhe')
xlabel('segment index')
ylabel('nonlocal writhe')
axes( 'Position', [0, 0, 1, 1] ) ;
set( gca, 'Color', 'None', 'XColor', 'White', 'YColor', 'White' ) ;
titletext = ['Wr = ' num2str(wr)] ;
text( 0.5, 0.95, titletext, 'FontSize', 14', 'FontWeight', 'Bold', ...
      'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;
axis off
saveas(gcf, fullfile(outdir, 'twoseg_helix.png'))
close all


%% traditional writhe:
ssx = ss_from_xyz(xyz) ;
[tangent, normal, binormal] = frenetSeretFrame(ssx, xyz(:, 1), xyz(:, 2), xyz(:, 3)) ;
ds = gradient(ssx) ;
wr = zeros(length(ssx), 1) ;
for jj=1:length(ssx)
    oind = setdiff(1:length(ssx), jj) ;
    rmr = xyz(jj, :) - xyz(oind, :) ;
    rmrmag = vecnorm(rmr')' ;
    txt = cross(tangent(jj,:) .* ones(length(oind), 3), tangent(oind, :));
    % Take row-wise inner product
    integrand = sum(sum(txt .* rmr, 2) ./ (rmrmag.^3 .* ones(size(txt))), 2) ;
    % Writhe per unit length is wr
    wr(jj) = sum(integrand) ;
end
wrtot = nansum(wr .* ds)

%% half twist that turns back
t = 0:0.01:1 ;
xx = sin(2*pi*t) ;
yy = cos(2*pi*t) ;
zz = 1.25 - (t - 0.5).^2 ;
xyz = [xx ; yy; zz]';
[wr, wr_local, wr_nonlocal, turns] = polarWrithe(xyz, []) ;

subplot(2, 2, 1)
% scatter3(xx', yy', zz', 10, yy) 
cvals = 1:length(xx) ;
scatter3(xx', yy', zz', 10, cvals) 
xlabel('x')
ylabel('y')
zlabel('z')
title('Curve')

% Plot the segments colored
subplot(2,2,2)
seg1 = 1:turns(1) ;
seg2 = turns(1):length(xx) ;
plot3(xx(seg1), yy(seg1), zz(seg1), 's') ;
hold on
plot3(xx(seg2), yy(seg2), zz(seg2), 'o') ;
title('Segments')
% Plot the local writhe
subplot(2, 2, 3)
plot(zz, wr_local)
title('local writhe')
xlabel('position, z')
ylabel('local writhe')
% Plot the nonlocal writhe
subplot(2, 2, 4)
plot(1:length(wr_nonlocal), wr_nonlocal, 'o')
title('nonlocal writhe')
xlabel('segment index')
ylabel('nonlocal writhe')
 axes( 'Position', [0, 0, 1, 1] ) ;
 set( gca, 'Color', 'None', 'XColor', 'White', 'YColor', 'White' ) ;
titletext = ['Wr = ' num2str(wr)] ;
text( 0.5, 0.95, titletext, 'FontSize', 14', 'FontWeight', 'Bold', ...
      'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;
saveas(gcf, fullfile(outdir, 'benthelix.png'))
close all

%% Torus wound curve
t = 0:0.01:1 ;
us = {t} ;
vs = {2 * t} ;
for ii = 1:length(us) 
    u = us{ii} ;
    v = vs{ii} ;
    yy = cos(2*pi * v) .* (2 + cos(2 * pi * u)) ;
    zz = sin(2*pi * v) .* (2 + cos(2 * pi * u)) ;
    xx = sin(2*pi * u) ;
    xyz = [xx ; yy; zz]';
    [wr, wr_local, wr_nonlocal, turns, segs, segpairs] = polarWrithe(xyz, []) ;

    colors = define_colors(length(segs)) ;
    markers = define_markers(length(segs)) ;
    
    subplot(2, 2, 1)
    % scatter3(xx', yy', zz', 10, yy) 
    cvals = 1:length(xx) ;
    scatter3(xx', yy', zz', 10, cvals) 
    xlabel('x')
    ylabel('y')
    zlabel('z')
    axis equal
    title('Curve')

    % Plot the segments colored
    subplot(2,2,2)
    for jj = 1:length(segs)
        plot3(xx(segs{jj}), yy(segs{jj}), zz(segs{jj}), '.', 'Color', colors(jj, :)) ;
        hold on
    end
    xlabel('x')
    ylabel('y')
    zlabel('z')
    axis equal
    title('Segments')
    % Plot the local writhe
    subplot(2, 2, 3)
    plot(zz, wr_local)
    title('local writhe')
    xlabel('position, z')
    ylabel('local writhe')
    % Plot the nonlocal writhe
    subplot(2, 2, 4)
    for jj = 1:length(wr_nonlocal)
        segments = segpairs{jj} ;
        segi = segments(1) ;
        segj = segments(2) ;
        plot(segi, wr_nonlocal(jj), markers{segj}, 'Color', colors(segj, :))
        hold on
    end
    title('nonlocal writhe')
    xlabel('segment pair index')
    ylabel('nonlocal writhe')
    axes( 'Position', [0, 0, 1, 1] ) ;
    set( gca, 'Color', 'None', 'XColor', 'White', 'YColor', 'White' ) ;
    titletext = ['Wr = ' num2str(wr)] ;
    text( 0.5, 0.95, titletext, 'FontSize', 14', 'FontWeight', 'Bold', ...
          'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;
    axis off
    saveas(gcf, fullfile(outdir, ['torusknot_' num2str(ii) '.png']))
    close all
end