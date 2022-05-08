% Figure drafting

foldIDfn3 = '/mnt/data/handGAL4klarHandGFPhistGFP/202105072030_1mWGFP/deconvolved_16bit/msls_output/gridCoords_nU0100_nV0100/lobes/fold_locations_sphi_nU0100_nV0100_avgpts.mat' ;
tmp = load(foldIDfn3) ;

datdir = '/mnt/data/analysis/tubular/gut_asymmetries/';
tmp1 = load(fullfile(datdir, 'kinematics_caax.mat')) ;
tmp2 = load(fullfile(datdir, 'kinematics_histRFP_0p01_0p00_0p01.mat')) ;
tmp3 = load(fullfile(datdir, 'kinematics_handGAL4handGFPhistGFP_0p00_0p00_0p00.mat')) ;
dats = {tmp1, tmp2, tmp3} ;
Ndat = 3 ;

%% nonlinear alignment
folds = [35,56,80;29,57,82;29,47,67];
sA = cell(Ndat, 1) ;
sB = cell(Ndat, 1) ;
sC = cell(Ndat, 1) ;
ss = cell(Ndat, 1) ;
Nsample = 50 ;

for qq = 1:Ndat
    % get mean position of first fold
    Apos = mean(double(folds(:, 1)-1)./(nU-1)) ;
    sA{qq} = linspace(0, Apos, folds(qq, 1)) ;
    
    % get mean position of second fold
    Bpos = mean(double(folds(:, 2)-1)./(nU-1)) ;
    sBtmp = linspace(Apos, Bpos, folds(qq, 2)-folds(qq, 1)+1) ;
    sB{qq} = sBtmp(2:end) ;
    
    % get mean position of third fold
    Cpos = mean(double(folds(:, 3)-1)./(nU-1)) ;
    sCtmp = linspace(Bpos, Cpos, folds(qq, 3)-folds(qq, 2)+1) ;
    sC{qq} = sCtmp(2:end) ;
    
    % get mean position of third fold
    Dpos = 1 ;
    sDtmp = linspace(Cpos, Dpos, nU-folds(qq, 3)+1) ;
    sD{qq} = sDtmp(2:end) ;
    
    ss{qq} = [sA{qq}, sB{qq}, sC{qq}, sD{qq}] ;
end

%% Obtain <= 1/4 circumference averages
ww = 10 ;
v2Hvn = cell(Ndat, 1) ;
d2Hvn = cell(Ndat, 1) ;
l2Hvn = cell(Ndat, 1) ;
r2Hvn = cell(Ndat, 1) ;
for qq = 1:Ndat
    dat = dats{qq} ;
    d2Hvn{qq} = mean(dat.H2vnAll(:, [1:ww, QS.nV-ww:QS.nV-1]), 2) ;
    v2Hvn{qq} = mean(dat.H2vnAll(:, 0.5*QS.nV-ww:QS.nV*0.5+ww), 2)  ;
    l2Hvn{qq} = mean(dat.H2vnAll(:, 0.25*QS.nV-ww:QS.nV*0.25+ww), 2)  ;
    r2Hvn{qq} = mean(dat.H2vnAll(:, 0.75*QS.nV-ww:QS.nV*0.75+ww), 2)  ;
    ddivv{qq} = mean(dat.divvAll(:, [1:ww, QS.nV-ww:QS.nV-1]), 2) ;
    vdivv{qq} = mean(dat.divvAll(:, 0.5*QS.nV-ww:QS.nV*0.5+ww), 2)  ;
    ldivv{qq} = mean(dat.divvAll(:, 0.25*QS.nV-ww:QS.nV*0.25+ww), 2)  ;
    rdivv{qq} = mean(dat.divvAll(:, 0.75*QS.nV-ww:QS.nV*0.75+ww), 2)  ;
end

%% Collate
sAll = ss{1} ;
dvAll2Hvn = d2Hvn{1} - v2Hvn{1} ;
lrAll2Hvn = l2Hvn{1} - r2Hvn{1} ;
dvAlldivv = ddivv{1} - vdivv{1} ;
lrAlldivv = ldivv{1} - rdivv{1} ;
for qq = 2:Ndat
    sAll = [sAll, ss{qq}];
    dvAll2Hvn = [dvAll2Hvn; d2Hvn{qq}-v2Hvn{qq}];
    lrAll2Hvn = [lrAll2Hvn; l2Hvn{qq}-r2Hvn{qq}];
    dvAlldivv = [dvAlldivv; ddivv{qq}-vdivv{qq}];
    lrAlldivv = [lrAlldivv; ldivv{qq}-rdivv{qq}];
end

NN = 90 ;
[midx, meanDV_2Hvn, stdDV_2Hvn, ~, steDV_2Hvn] = binDataMeanStd(sAll,dvAll2Hvn(:),linspace(0, 1, NN)) ;
[midx, meanLR_2Hvn, stdLR_2Hvn, ~, steLR_2Hvn] = binDataMeanStd(sAll,lrAll2Hvn(:),linspace(0, 1, NN)) ;

[midx, meanDV_divv, stdDV_divv, ~, steDV_divv] = binDataMeanStd(sAll,dvAlldivv(:),linspace(0, 1, NN)) ;
[midx, meanLR_divv, stdLR_divv, ~, steLR_divv] = binDataMeanStd(sAll,lrAlldivv(:),linspace(0, 1, NN)) ;

%% Save data -- raw and processed

save('DV_LR_asymmetries_raw.mat', 'ss', 'd2Hvn','v2Hvn','l2Hvn','r2Hvn',...
    'ddivv','vdivv','ldivv','rdivv',...
    'ww', 'folds', 'dats', ...
    'midx', 'meanDV_2Hvn', 'meanLR_2Hvn', 'stdDV_2Hvn', 'steDV_2Hvn', 'stdLR_2Hvn', 'steLR_2Hvn',...
    'midx', 'meanDV_divv', 'meanLR_divv', 'stdDV_divv', 'steDV_divv', 'stdLR_divv', 'steLR_divv')



%% Plot result -- raw
dfields = {d2Hvn, ddivv} ;
vfields = {v2Hvn, vdivv} ;
lfields = {l2Hvn, ldivv} ;
rfields = {r2Hvn, rdivv} ;
fnames = {'2Hvn', 'divv'} ;
titles = {'$2Hv_n$', '$\nabla \cdot \mathbf{v}_{\parallel}$'} ;

for pp = 1:2
    dd = dfields{pp};
    vv = vfields{pp};
    ll = lfields{pp};
    rr = rfields{pp};
    fname = fnames{pp} ;
    
    close all
    fig = gcf ;
    sh1 = subplot(2, 1, 1) ;
    % plot(linspace(0,1,QS.nU), dorsH2vn) ; hold on;
    % plot(linspace(0,1,QS.nU), ventH2vn) ;
    for qq = 1:Ndat
        dv = dd{qq} - vv{qq} ;
        plot(ss{qq}, dv) ;
        hold on;
    end
    xlabel('position, s/L')
    ylabel('DV asymmetry')
    ylims1 = ylim ;

    sh2 = subplot(2, 1, 2) ;
    for qq = 1:Ndat
        lr = ll{qq} - rr{qq} ;
        plot(ss{qq}, lr) ;
        hold on;
    end
    xlabel('position, s/L')
    ylabel('LR asymmetry')
    ylims2 = ylim ;
    ymax = max(abs([ylims1,  ylims2])) ;
    ylim(ymax * [-1,1])
    set(fig, 'currentaxes', sh1);
    ylim(ymax * [-1,1])
    
    sgtitle(titles{pp}, 'interpreter', 'latex') 

    saveas(gcf, fullfile(datdir, ...
        sprintf([fname '_DV_LR_asymmetries_raw_w%03d.pdf'], ww)))
end

%% Plot processed
colors = define_colors ;

DVmeanfields = {meanDV_2Hvn, meanDV_divv} ;
DVstdfields = {stdDV_2Hvn, stdDV_divv} ;
DVstefields = {steDV_2Hvn, steDV_divv} ;
LRmeanfields = {meanLR_2Hvn, meanLR_divv} ;
LRstdfields = {stdLR_2Hvn, stdLR_divv} ;
LRstefields = {steLR_2Hvn, steLR_divv} ;
fnames = {'2Hvn', 'divv'} ;
titles = {'$2Hv_n$', '$\nabla \cdot \mathbf{v}_{\parallel}$'} ;

for pp = 1:2
    mDV = DVmeanfields{pp};
    stdDV = DVstdfields{pp};
    steDV = DVstefields{pp};
    mLR = LRmeanfields{pp};
    stdLR = LRstdfields{pp};
    steLR = LRstefields{pp};
    fname = fnames{pp} ;

    close all
    fig = gcf ;
    sh1 = subplot(2, 1, 1) ;
    lineProps = {'-','color', colors(1, :)} ;
    shadedErrorBar(midx, mDV, stdDV, 'lineProps', lineProps);
    xlabel('position, s/L')
    ylabel('DV asymmetry')
    ylims1 = ylim ;

    sh2 = subplot(2, 1, 2) ;
    lineProps = {'-','color', colors(2, :)} ;
    shadedErrorBar(midx, mLR, stdLR, 'lineProps', lineProps);
    xlabel('position, s/L')
    ylabel('LR asymmetry')
    ylims2 = ylim ;
    ymax = max(abs([ylims1,  ylims2])) ;
    ylim(ymax * [-1,1])
    set(fig, 'currentaxes', sh1);
    ylim(ymax * [-1,1])

    sgtitle(titles{pp}, 'interpreter', 'latex') 
    saveas(gcf, fullfile(datdir, ...
        sprintf([fname '_DV_LR_asymmetries_binstats_w%03d.pdf'], ww)))

end