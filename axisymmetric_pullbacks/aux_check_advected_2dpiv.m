% evaluate flow field at (x, Ysz - y). Add (u,v) to (x, Ysz - y).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the advected mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the advected tissue on top of the next frame
% interpolate velocities (x0, y0) onto mesh (mesh0x, mesh0y)
uuinterp = griddedInterpolant(x0', y0', uu') ; % in inverted Y
vvinterp = griddedInterpolant(x0', y0', vv') ; 
umeshpix = uuinterp(mesh0x, ysz - mesh0y) ; % interpolate in inverted Y
vmeshpix = vvinterp(mesh0x, ysz - mesh0y) ;
addx = umeshpix ;
addy = vmeshpix ;
mesh0adv_pix = [ mesh0x + addx, (ysz - mesh0y) + addy] ;

% Texture image options
Options.imSize = ceil( 1000 .* [ a_fixed 1 ] );  % row, col pix
Options.yLim = [0 ysz];
Options.xLim = [0 1000];
% original image RED
im0 = texture_patch_to_image(mesh0.f, [mesh0x, ysz - mesh0y], ...
    mesh0.f, mesh0.v(:, [ 2 1 3]), IV0, Options );
% im0 = imread(fullfile(fns(i).folder, fns(i).name)) ;
% advected image GREEN
patchIma = texture_patch_to_image(mesh0.f, mesh0adv_pix, ...
    mesh0.f, mesh0.v(:, [ 2 1 3]), IV0, Options );
% patchIma = adapthisteq(patchIma, 'NumTiles', [round(a_fixed * ntiles), round(2 * ntiles)]) ;
% Next timepoint BLUE
im1 = texture_patch_to_image(mesh1.f, [mesh1x, ysz - mesh1y], ...
    mesh1.f, mesh1.v(:, [2 1 3]), IV1, Options );
% Load next timepoint BLUE
% im1 = imread(fullfile(fns(i+1).folder, fns(i+1).name)) ;
% patchImRGB = cat(3, im0, uint8(patchIma * 255), im1) ;
patchImRGB = cat(3, im0, patchIma, im1) ;

% Plot the image and advected image, in inverted Y space but
% with increasing Y going down
close all; 
fig1 = figure(1) ; 
h = imshow( patchImRGB );
hold on;  
quiver(mesh0x, ysz - mesh0y, addx, addy, 0, 'color', yellow, 'linewidth', 3)
% plot(mesh0x, ysz - mesh0y, 'o')
% plot(mesh0adv_pix(:, 1), mesh0adv_pix(:, 2), 's')
triplot(mesh0.f, mesh0x, ysz - mesh0y, 'color', red, 'linewidth', 2)
triplot(mesh0.f, mesh0x + addx, ysz - mesh0y + addy, 'color', green, 'linewidth', 2)
% axis equal

% Check the displacement by toggling between frames
% fig = figure(2) ;
% pressed_enter = false ; 
% kk = 0 ;
% while ~pressed_enter
%     if kk == 0
%         imshow(flipud(im0))
%         %set( gca, 'YDir', 'Normal' );
%         titlestr = '0: <Enter> to exit, <-> to toggle' ;
%     else
%         imshow(flipud(im1))
%         %set( gca, 'YDir', 'Normal' );
%         titlestr = '1: <Enter> to exit, <-> to toggle' ;
%     end
%     title(titlestr)
%     was_a_key = waitforbuttonpress;
%     left = strcmp(get(fig, 'CurrentKey'), 'leftarrow');
%     rght = strcmp(get(fig, 'CurrentKey'), 'rightarrow') ;
%     if was_a_key && strcmp(get(fig, 'CurrentKey'), 'return')
%         pressed_enter = true ;
%     elseif was_a_key && left || rght 
%         kk = mod(kk + 1, 2) ;
%     end
% end

% Check differences
d0 = mat2gray(im1, [0, double(max(im1(:)))]) - mat2gray(im0, [0, double(max(im0(:)))]);
d0pos = 0*d0 ;
d0pos(d0 > 0) = d0(d0 > 0) ;
d0neg = 0*d0 ;
d0neg(d0 < 0) = abs(d0(d0 < 0)) ;
pd0 = cat(3, d0pos, 0*d0, d0neg) ; 

% diff between next and advected 
da = mat2gray(im1, [0, double(max(im1(:)))]) - patchIma ;
dapos = 0*da ;
dapos(da > 0) = da(da > 0) ;
daneg = 0*da ;
daneg(da < 0) = abs(da(da < 0)) ;
pda = cat(3, dapos, 0*da, daneg) ;
% plot both
fig1 = figure(1) ;
imshow(pd0) 
title(['\langle|t1 - t0|\rangle = ' num2str(100 * mean(pd0(:))) '%'])
waitfor(fig1)

fig2 = figure(2) ;
imshow(pda) 
title(['\langle|t1 - advected t0|\rangle = ' num2str(100 * mean(pda(:))) '%'])
waitfor(fig2)