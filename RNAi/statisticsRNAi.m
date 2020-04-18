% Read embryo fate data from WT, mef2, and WASp

root = '/mnt/data/RNAi' ;
WTDir = fullfile(root, 'Control_mef2Gal4klar') ;
RokDir = fullfile(root, 'RokRNAi20_mef2GAL4klar_BF') ;
WASpDir = fullfile(root, 'WASpRNAi20_48YGAL4klar_BF') ; 
dirsWT = subdirs(WTDir) ;
dirsRok = subdirs(RokDir) ;
dirsWASp = subdirs(WASpDir) ;

wcell = buildExpCell(dirsWT) ;
rcell = buildExpCell(dirsRok) ;
scell = buildExpCell(dirsWASp) ;

% Find all WT in each stage, count by classification
buildFateHist(wcell) ;
title('Wild type')
saveas(gcf, fullfile(root, 'wildtype_fatehist.png'))
close all

buildFateHist(rcell) ;
title('Rok RNAi (VALIUM 20) x Mef2GAL4 klar')
saveas(gcf, fullfile(root, 'rok20_Mef2gal4klar_fatehist.png'))
close all

buildFateHist(scell) ;
title('WASp RNAi (VALIUM 20) x 48YGAL4 klar')
saveas(gcf, fullfile(root, 'rok20_48Ygal4klar_fatehist.png'))
close all


% Also plot as lumped
bins = [1, 10, 14, 17] ;
buildFateHist(wcell, bins) ;
title('Wild type')
saveas(gcf, fullfile(root, 'wildtype_fatehist_bunch.png'))
close all

buildFateHist(rcell, bins) ;
title('Rok RNAi (VALIUM 20) x Mef2GAL4 klar')
saveas(gcf, fullfile(root, 'rok20_Mef2gal4klar_fatehist_bunch.png'))
close all

buildFateHist(scell, bins) ;
title('WASp RNAi (VALIUM 20) x 48YGAL4 klar')
saveas(gcf, fullfile(root, 'rok20_48Ygal4klar_fatehist_bunch.png'))
close all


