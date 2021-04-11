% script_explore_FundFormNematic

% Create disk that transforms from spherical cap to saddle, 
% spins around, and back

outdir = '/mnt/data/analysis/demo_HopfDifferential_halfAngle/' ;
if ~exist(outdir, 'dir')
    mkdir(outdir)
end

RR = 0.5 ;
xy = vogelDisk(5000) * RR ;
faces = delaunay(xy) ;
Rphi = 1e10;
nn = 40 ;
% exponent = 0.3 ;
% r2do0 = [fliplr(1 ./ (logspace(-5, 1, nn)).^(exponent)), ...
%     -fliplr((logspace(-1, 5, nn)).^(exponent))] ;

exponent = 3 ;
r2do0 = [linspace(0.3^(1/exponent), 1, nn).^(exponent), 5, 10, -5, ...
    -fliplr(linspace(0.3^(1/exponent), 1, nn).^(exponent))] ;

extraN = 8*nn ;
r2do = [r2do0, r2do0(end) * ones(1,extraN)] ;
rot = [0 * r2do0, linspace(0, 4*pi, extraN) ] ;

id2do = 1:10:length(r2do) ;
id2do = [id2do, setdiff(1:length(r2do), id2do)] ;

for ii = id2do
    disp(['ii = ' num2str(ii)])
    Rz = r2do(ii) ;
    
    rot0 = rot(ii) ;
    Rot = [cos(rot0), sin(rot0); -sin(rot0), cos(rot0)];
    xy0 = (Rot * xy')' ;
    zz = catenoid(Rphi, Rz, xy0) ;
    zz = zz - mean(zz) ;
    
    % trisurf(faces, xy(:, 1), xy(:, 2), zz, ...
    %     'edgecolor', 'none')
    
    xyz = cat(2, xy, zz) ;
    [gg, bb] = constructFundamentalForms(faces, xyz, xy) ;
    QQ = hopfDifferential(bb) ;    

    mags = sqrt(abs(QQ)) ;
    thetas = -atan2(imag(QQ), real(QQ))*0.5 ;
    
    % Mask out boundary faces
    bnd = freeBoundary(triangulation(faces, xyz))  ;
    edgeFaces = find(any( ismember( faces, bnd(:, 1) ), 2 ));
    mags(edgeFaces) = 0 ;
    thetas(edgeFaces) = 0 ;
    
    % Create mesh
    m3d = struct() ;
    m3d.f = faces ;
    m3d.v = xyz ;
    m2d = m3d ;
    m2d.v(:, 3) = 0 ;
    
    bc2d = barycenter(m2d.v, m2d.f) ;
    nearBnd = sqrt(bc2d(:, 1).^2 + bc2d(:, 2).^2) > 0.95 * RR ;
    mags(nearBnd) = 0 ;
    thetas(nearBnd) = 0 ;
    
    %% Mean curvature
    DEC = DiscreteExteriorCalculus(m3d.f, m3d.v) ;
    m3d.vn = per_vertex_normals(m3d.v, m3d.f) ;
    H3d = sum(m3d.vn .* DEC.laplacian(m3d.v), 2) * 0.5 ;
    
    %% Plot as colored mesh with quiver
    clf
    opts.clims = {[0,0.5], [0, median(mags(:, 1))], [-2, 2]} ;
    opts.axisOff = true ;
    opts.labels = {'', 'Hopf differential, $\theta=-\frac{1}{2}\tan^{-1}(\Im Q/\Re Q)$', ''} ;
    opts.views = {[45, 45 ], [0, 90 ], [45,45]};
    opts.phasebarPosition = [0.33, 0.20, 0.1, 0.135] ;
    opts.xyzlims = 0.5 * [-1, 1;-1,1;-1,1] ;
    opts.cbarlabels = {['Hopf differential, $\sqrt{|Q|}$ [L$^{-1}$]'], ...
                       ['Hopf differential, $\sqrt{|Q|}$ [L$^{-1}$]'], ...
                       ['Mean curvature, $H$ [L$^{-1}$]']} ;
    opts.polarStyle = 'nematic';
    opts.visible = 'on'; 
    close all
    set(gcf, 'visible', 'off')
    opts.visible = 'off' ;
    axs = nFieldsOnSurface({m3d, m2d, m3d}, ...
        {{mags, thetas}, {mags, thetas}, H3d}, opts) ;
    hold on
    set(gcf, 'currentAxes', axs{1})
    bc3d = barycenter(m3d.v, m3d.f) ;
    Qf3d = pushVectorField2Dto3DMesh(mags .* [cos(thetas), sin(thetas)], ...
        m2d.v, m3d.v, m3d.f, 1:length(m3d.f)) ;
    Qf2d = mags .* [cos(thetas), sin(thetas), 0*thetas] ;
    
    idx = 1:50:length(bc3d) ;
    quiver3(bc3d(idx, 1), bc3d(idx, 2), bc3d(idx, 3), ...
        Qf3d(idx, 1), Qf3d(idx, 2), Qf3d(idx, 3), 1, 'w', 'ShowArrowHead', 'off') 
    quiver3(bc3d(idx, 1), bc3d(idx, 2), bc3d(idx, 3), ...
        -Qf3d(idx, 1), -Qf3d(idx, 2), -Qf3d(idx, 3), 1, 'w', 'ShowArrowHead', 'off') 
    % axis on
    % xlabel('x');
    % ylabel('y')
    % zlabel('z')
    
    %% Plot in 2d on right
    set(gcf, 'currentAxes', axs{2})
    quiver3(bc2d(idx, 1), bc2d(idx, 2), bc2d(idx, 3), ...
        Qf2d(idx, 1), Qf2d(idx, 2), Qf2d(idx, 3), 1, 'w', 'ShowArrowHead', 'off') 
    quiver3(bc2d(idx, 1), bc2d(idx, 2), bc2d(idx, 3), ...
        -Qf2d(idx, 1), -Qf2d(idx, 2), -Qf2d(idx, 3), 1, 'w', 'ShowArrowHead', 'off') 
    
    %% Title and output'-r150'
    if median(H3d) > -1e-19
        Rzsigned = abs(Rz) ;
    else
        Rzsigned = - abs(Rz) ;
    end
    sgtitle(['$R_1=$' sprintf('%0.3f',Rzsigned)], 'interpreter', 'latex')
    outfn = fullfile(outdir, sprintf('hopfDiff_%04d.png', ii)) ;
    set(gcf, 'color', 'w')
    export_fig(outfn, '-nocrop', '-r150')  % '-transparent', '-r150'

end