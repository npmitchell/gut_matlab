% for handGAL4 handGFP histGFP, define as
% options.lambda = 0 ;
% options.lambda_err = 0 ;
% options.lambda_mesh = 0 ;
QS.dir.metricKinematics.pathline.measurements = ...
    '/mnt/data/handGAL4klarHandGFPhistGFP/202105072030_1mWGFP/deconvolved_16bit/msls_output/gridCoords_nU0100_nV0100/metricKinematics_dtOffBy2/lambda0p000_lmesh0p000_lerr0p000_modes07w02/pathline_0001t0/measurements';
tSpan = 45 ;
div2 = true ;

%% Other wise use:
tSpan = 90 ;
div2 = false ;

%% Figure for eLife rebuttal

try
    apKymos = load(fullfile(sprintf(...
       QS.dir.metricKinematics.pathline.measurements, QS.t0set()),...
       'apKymographMetricKinematics.mat')) ;
catch
    error('apKymos not saved')
end
firstpass = true ;
tps = QS.t0set():QS.t0set()+tSpan ;
gdotM = cell(1,1) ;
divvM = cell(1,1) ;
H2vnM = cell(1,1) ;
for tt = tps
    divv = load(fullfile(sprintf(...
        QS.dir.metricKinematics.pathline.measurements, QS.t0set()),...
        sprintf('divv_pathline%04d_%06d.mat', QS.t0set(),tt))) ;
    H2vn = load(fullfile(sprintf(...
        QS.dir.metricKinematics.pathline.measurements, QS.t0set()),...
        sprintf('H2vn_pathline%04d_%06d.mat', QS.t0set(),tt))) ;
    gdot = load(fullfile(sprintf(...
        QS.dir.metricKinematics.pathline.measurements, QS.t0set()),...
        sprintf('gdot_pathline%04d_%06d.mat', QS.t0set(),tt))) ;
    
    if firstpass
        gdotAll = gdot.gdot ;
        divvAll = divv.divv ;
        H2vnAll = H2vn.H2vn ;
        firstpass = false ;
    else
        gdotAll = gdotAll + gdot.gdot ;
        divvAll = divvAll + divv.divv ;
        H2vnAll = H2vnAll + H2vn.H2vn ;
    end
    gdotM{tt} = gdot.gdot ;
    divvM{tt} = divv.divv ;
    H2vnM{tt} = H2vn.H2vn ;
    
end
gdotAll = gdotAll ./ length(tps) ;
divvAll = divvAll ./ length(tps) ;
H2vnAll = H2vnAll ./ length(tps) ;

[~,~,~,xyzlims] = QS.getXYZLims() ;
xyzlims = xyzlims + [-5, 5] ;
gdotAll(1:3, :) = 0 ;
divvAll(1:3, :) = 0 ;
H2vnAll(1:3, :) = 0 ;
gdotAll(end-3:end, :) = 0 ;
divvAll(end-3:end, :) = 0 ;
H2vnAll(end-3:end, :) = 0 ;
QS.setTime(QS.t0set())
mesh = QS.getCurrentSPCutMeshSmRS() ;

% for histRFP: lambda0p010_lmesh0p000_lerr0p010_modes07w01
% for handGAL4handhistGFP: lambda0p000_lmesh0p000_lerr0p000_modes07w02
if div2 
    gdotAll = gdotAll * 0.5 ;
    divvAll = divvAll * 0.5 ;
    H2vnAll = H2vnAll * 0.5 ;
end

%% Keep certain time window
gdot90 = apKymos.gdot_apM(QS.xp.tIdx(QS.t0set()):QS.xp.tIdx(QS.t0set()+tSpan),:) ;

if div2 
    gdot90 = gdot90 * 0.5 ;
end

close all
fig = figure('Position', [0 0 1200 900], 'Units', 'pixels') ;

    set(fig, 'Position', [0, 0, 400, 400])
imagesc(linspace(0, 1, QS.nU), 1:90, gdot90)
caxis(0.05 * [-1, 1])
colorbar


% COLORMAP
% colormap( brewermap(256, '*RdBu'))
colormap(bwr) 
cmapStr = '_bwr'

ylabel(['time [' QS.timeUnits ']'], 'Interpreter', 'latex') ;
xlabel('ap position in material frame [$s/L$]', 'Interpreter', 'latex') ;
axis square
axis tight
outdir = fullfile(sprintf(...
    QS.dir.metricKinematics.pathline.root, QS.t0set())) ;
outfn = fullfile(outdir, ['mean_gdot_apKymograph' cmapStr '.pdf']) ;
saveas(gcf, outfn) ;
outfn = fullfile(outdir, ['mean_gdot_apKymograph' cmapStr '.png']) ;
saveas(gcf, outfn) ;
axis off
FF = getframe() ;
outfn = fullfile(outdir, ['mean_gdot_apKymograph_cdata' cmapStr '.png']) ;
imwrite(FF.cdata, outfn)

%% Quantify asymmetry
ww = 10 ;
ventH2vn = mean(H2vnAll(:, [1:ww, QS.nV-ww:QS.nV-1]), 2) ;
dorsH2vn = mean(H2vnAll(:, 0.5*QS.nV-ww:QS.nV*0.5+ww), 2)  ;
leftH2vn = mean(H2vnAll(:, 0.25*QS.nV-ww:QS.nV*0.25+ww), 2)  ;
rightH2vn = mean(H2vnAll(:, 0.75*QS.nV-ww:QS.nV*0.75+ww), 2)  ;

close all
subtightplot(2, 1, 1)
plot(linspace(0,1,QS.nU), dorsH2vn) ; hold on;
plot(linspace(0,1,QS.nU), ventH2vn) ;
plot(linspace(0,1,QS.nU), dorsH2vn-ventH2vn) ;

subtightplot(2, 1, 2)
plot(linspace(0,1,QS.nU), leftH2vn) ; hold on;
plot(linspace(0,1,QS.nU), rightH2vn) ;
plot(linspace(0,1,QS.nU), leftH2vn-rightH2vn) ;

%ofn = fullfile('/mnt/data/analysis/tubular/gut_asymmetries', 'kinematics_histRFP_0p01_0p00_0p01.mat') ;
ofn = fullfile('/mnt/data/analysis/tubular/gut_asymmetries', 'kinematics_handGAL4handGFPhistGFP_0p00_0p00_0p00.mat') ;
save(ofn, 'H2vnAll', 'divvAll', 'gdotAll', 'dorsH2vn')

%%
fields = {H2vnAll, divvAll, gdotAll} ;
fnames = {'H2vn', 'divv', 'gdot' } ;
titles = {'$\langle2Hv_n\rangle$', ...
    '$\langle\nabla\cdot \mathbf{v}_\parallel\rangle$', ...
    '$\langle\frac{1}{2}\mathrm{Tr} \left[ \mathbf{g}^{-1} \dot{\mathbf{g}} \right]\rangle$', ...    
     } ;
resStr = '-r600' ;
cmax = 0.3 ;
for fieldID = 1:3
    field = fields{fieldID} ;
    fname = fnames{fieldID} ;
    close all
    fig = figure('Position', [0 0 1200 900], 'Units', 'pixels') ;
    set(fig, 'color', 'w')
    trisurf(triangulation(mesh.f, mesh.v), field, 'edgecolor', 'none');
    colorbar()
    colormap(bwr)
    view([0, 0])
    caxis(cmax * [-1, 1])
    grid off
    axis equal ;
    xlim(xyzlims(1, :))
    ylim(xyzlims(2, :))
    zlim(xyzlims(3, :))
    title(titles{fieldID}, 'interpreter', 'latex')
    outdir = fullfile(sprintf(...
        QS.dir.metricKinematics.pathline.root, QS.t0set()), ...
        strrep(sprintf('clim%0.2f', cmax), '.', 'p')) ;
    mkdir(outdir)
    outfn = fullfile(outdir, ['mean_' fname cmapStr '_lateralL.pdf']) ;
    disp(['Saving ' outfn])
    export_fig(fig, outfn, '-nocrop',resStr)
    view([0, -90])
    outfn = fullfile(outdir, ['mean_' fname cmapStr '_ventral.pdf']) ;
    disp(['Saving ' outfn])
    export_fig(fig,  outfn,'-nocrop', resStr)
    view([0, -180])
    export_fig(fig,  fullfile(outdir, ['mean_' fname cmapStr '_lateralR.pdf']),'-nocrop',  resStr)
    view([0, -270])
    export_fig(fig,  fullfile(outdir, ['mean_' fname cmapStr '_dorsal.pdf']), '-nocrop', resStr)
    
    set(fig, 'Position', [0, 0, 200, 100])
    export_fig(fig, fullfile(outdir, ['mean_' fname cmapStr '_scalebar.pdf']), resStr)
end
