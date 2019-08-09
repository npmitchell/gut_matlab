% Plot the gaussian curvature & metric tensor components
% NPM 2019

KDir = fullfile(soiDir, 'Kcurv/') ; 
if ~exist(KDir, 'dir')
    mkdir(KDir)
end
cmp = diverging_cmap(0:1e-3:1, 1, 2) ;

fig = figure('Visible', 'off') ;

for tp = 1:length(fileMeta.timePoints)
    if all(size(SOI.fields{2}(tp).patches) == 0)
        SOI.NCalcInducedMetric(fileMeta.timePoints(tp))
    end
    mm = SOI.fields{2}(tp).patches{1}.apply() ;
    EE = mm{1} ;
    FF = mm{2} ; 
    GG = mm{4} ;
    % Brioschi formula gives K = det()- det() / (EG - F^2)^2
    % see https://en.wikipedia.org/wiki/Gaussian_curvature
    disp('taking gradients')
    [Eu, Ev] = gradient(EE); [Fu, Fv] = gradient(FF); [Gu, Gv] = gradient(GG);
    [Euu, Evu] = gradient(Eu) ; 
    [Euv, Evv] = gradient(Ev) ;
    [Fuu, Fvu] = gradient(Fu) ;
    [Fuv, Fvv] = gradient(Fv) ;
    [Guu, Gvu] = gradient(Gu) ;
    [Guv, Gvv] = gradient(Gv) ;
    a11 = -0.5 * Evv + Fuv - 0.5 * Guu ;
    a12 = 0.5 * Eu ;
    a13 = Fu - 0.5 * Ev ;
    a21 = Fv - 0.5 * Gu ;
    a22 = EE ;
    a23 = FF ;
    a31 = 0.5 * Gv ;
    a32 = FF ;
    a33 = GG ;
    disp('defining first term in numerator')
    aa = reshape(arrayfun(@(x) det([a11(x), a12(x), a13(x); ...
                            a21(x), a22(x), a23(x); ...
                            a31(x), a32(x), a33(x)]), 1:numel(EE)), size(EE)) ;
    disp('defining second term in numerator')
    bb = reshape(arrayfun(@(x) det([ 0, 0.5 * Ev(x), 0.5 * Gu(x); ...
           0.5 * Ev(x), EE(x), FF(x); ...
           0.5 * Gu(x), FF(x), GG(x) ]), 1:numel(EE)), size(EE)) ; 
    disp('defining denominator')
    cc = reshape(arrayfun(@(x) (EE(x) .* GG(x) - FF(x) .* FF(x)) .^2, 1:numel(EE)), size(EE)) ;
    disp('computing curvature')
    KK = (aa - bb) ./ cc;

    % An alternate form is below for reference
    % KK = (ee .* gg - f .* f) ./ (EE .* GG - FF .* FF) ;

    % Filter the image
    sigma = 10 ;
    KK = imgaussfilt(KK, sigma) ;

    % Plot it
    % Define the colormap
    imagesc(KK)
    colormap(cmp)
    ax = gca ;
    ax.CLim = [-1e-3, 1e-3] ;
    % ax.XLim([1, 1000])
    % ax.YLim([1, 2000])
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    axis equal
    h = colorbar ;
    ylabel(h, 'K') ;
    name = sprintf('Kcurv_%05d.png', fileMeta.timePoints(tp));
    disp(['saving file: ', KDir, name]) 
    set(gcf, 'PaperPositionMode', 'auto');
    saveas(gcf, fullfile(KDir, name)) ;
    clf
end