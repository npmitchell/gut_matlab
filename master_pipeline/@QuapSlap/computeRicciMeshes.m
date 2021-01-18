function computeRicciMeshes(QS, options)
% computeRicciMeshes(QS, options) 
%   Compute all Ricci meshes (if flow converges) and plot aspect ratio for isothermal PB over time
%
% Parameters
% ----------
% options : optional struct with fields
%   
%
% NPMitchell 2021

% Default options
maxIter = 200 ;

if nargin < 2
    options = struct();
end
if isfield(options, 'maxIter')
    maxIter = options.maxIter ;
end

% Compute each Ricci flow mesh
opts = struct() ;
opts.maxIter = maxIter ;

for tp = QS.xp.fileMeta.timePoints
    disp(['t = ', num2str(tp)])
    QS.setTime(tp)
    try
        QS.generateRicciMeshTimePoint(tp, opts) 
    catch
        disp('could not generate Ricci mesh -- self intersections?')
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