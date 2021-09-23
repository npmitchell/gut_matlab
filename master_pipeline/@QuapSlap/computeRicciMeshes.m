function computeRicciMeshes(QS, options)
% computeRicciMeshes(QS, options) 
%   Compute all Ricci meshes (if flow converges) and plot aspect ratio for isothermal PB over time
%
% Parameters
% ----------
% options : optional struct with fields
%   resample : bool (default=true) ;
%   
%
% NPMitchell 2021

% Default options
maxIter = 200 ;
resample = true ;

if nargin < 2
    options = struct();
end
if isfield(options, 'maxIter')
    maxIter = options.maxIter ;
end
if isfield(options, 'resample')
    resample = options.resample ;
end

% Compute each Ricci flow mesh
opts = struct() ;
opts.maxIter = maxIter ;
opts.resample = resample ;

t0 = QS.t0set() ;
timePoints = QS.xp.fileMeta.timePoints ;
startID = 87 ;
tidx2do1 = startID:50:length(timePoints) ;
tidx2do1 = [tidx2do1, setdiff(startID:30:length(timePoints), tidx2do1)] ;
tidx2do2 = setdiff(startID:20:length(timePoints), tidx2do1) ;
tidx2do3 = setdiff(setdiff(startID:10:length(timePoints), ...
    tidx2do1), tidx2do2) ;
tidx2do123 = [tidx2do1, tidx2do2, tidx2do3] ;
tidx2do1234 = [tidx2do123, ...
    setdiff(startID:2:length(timePoints), tidx2do123)] ;
tidx2do = [tidx2do1234, ...
    setdiff(1:length(timePoints), tidx2do1234)] ;
tidx0 = QS.xp.tIdx(t0) ;
tidx2do = [tidx0, setdiff(tidx2do, [tidx0], 'stable')] ;
% check that there are no repeats
assert(length(tidx2do) == length(unique(tidx2do)))
assert(all(tidx2do == unique(tidx2do, 'stable')))

for tiix = 1:length(tidx2do)
    tidx = tidx2do(tiix) ;
    tp = timePoints(tidx) ;
    disp(['t = ', num2str(tp)])
    QS.setTime(tp)
    if resample
        QS.generateRicciMeshTimePoint(tp, opts) 
    else
        try
            QS.generateRicciMeshTimePoint(tp, opts) 
        catch
            disp('could not generate Ricci mesh -- self intersections?')
        end
    end
end

aratio_r = zeros(length(QS.xp.fileMeta.timePoints), 1) ;
aratio_a = zeros(length(QS.xp.fileMeta.timePoints), 1) ;
rads = zeros(length(QS.xp.fileMeta.timePoints), 1) ;
for tidx = 1:length(QS.xp.fileMeta.timePoints)
    tp = QS.xp.fileMeta.timePoints(tidx) ;
    disp(['t = ', num2str(tp)])
    QS.setTime(tp)
    try
        rmesh = QS.loadCurrentRicciMesh() ;
        aratio_r(tidx) = max(rmesh.rectangle.u(:,1)) / (2*pi) ;
        aratio_a(tidx) = max(vecnorm(rmesh.annulus.u, 2, 2)) / min(vecnorm(rmesh.annulus.u, 2, 2)) ;
    catch
        aratio_r(tidx) = NaN ;
        aratio_a(tidx) = NaN ;
    end
    [~, radius_cline] = QS.getRadii() ;
    rads(tidx) = mean(radius_cline) ;
end
Length_t = QS.measureLength() ;
lengs = Length_t.lengths ;

%% plot it
clf
t0 = QS.t0set() ;
tidx0 = QS.xp.tIdx(t0) ;
tps = (QS.xp.fileMeta.timePoints - t0) * QS.timeInterval ;
min2hr = strcmpi(QS.timeUnits, 'min') && (length(tps) > 90) ;
if min2hr
    tps = tps / 60 ;
end

plot(tps, aratio_r, '.-')
hold on
plot(tps, log10(aratio_a), '.-')
plot(tps, lengs / lengs(tidx0), '.-')
plot(tps, lengs ./ (2*pi*rads), '.-')
legend({'Ricci $L_\zeta/L_\phi$', ...
    'Ricci $\log_{10}(\textrm{max}\rho/\textrm{min}\rho)$', ...
    '$L_{\Re^3}/L_{\Re^3}^{0}$', ...
    '$L_{\Re^3}/(2\pi\langle r \rangle))$'}, ...
    'interpreter', 'latex', 'location', 'best')
title('Aspect ratios over time', 'interpreter', 'latex')
if min2hr
    xlabel('time, $t$ [hr]', 'interpreter', 'latex')
else
    xlabel(['time, $t$ [' QS.timeUnits ']'], 'interpreter', 'latex')
end
saveas(gcf, fullfile(QS.dir.ricci.data, 'aspect_ratios.pdf'))