% Demonstrate the polar Writhe definition on example curves
outdir = '~/Desktop/polarwrithe/' ;
mkdir(outdir) ;

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
suptitle(['Wr = ' num2str(wr)])
saveas(gcf, fullfile(outdir, 'twoseg_helix.png'))
close all

%% Half twist
t = 0:0.01:1 ;
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
suptitle(['Wr = ' num2str(wr)])
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
suptitle(['Wr = ' num2str(wr)])
saveas(gcf, fullfile(outdir, 'oneseg_helix_reverse.png'))
close all


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
suptitle(['Wr = ' num2str(wr)])
saveas(gcf, fullfile(outdir, 'twoseg_helix.png'))
close all
