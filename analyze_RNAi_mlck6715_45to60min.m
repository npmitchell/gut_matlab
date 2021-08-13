% analyze_RNAi_results
% Here, 
%   13 is bands forming, late 13 is bands formed
%   14a is closing
%   14b is mostly closed but not fully
%   15a is after closure
%   15b is after/during first fold

% 0: tub67;tub15
% 0.1: 48YGAL4 (no UAS)
% 0.2: 48YGAL4 klar x UAS-CAAX mCh
% 0.3: mef2Gal4 (no UAS)
% 0.4: UAS-SERCA (no GAL4)
% 1: tub67;tub15 x UAS-MLCK RNAi TRiP#4
% 2: mef2G4klar x UAS-SERCA.R751Q
x = 100 ;
scoreLowerLimit = 1 ;
scoreUpperLimit = 6 ;
outdir = '/mnt/data/analysis/mlckRNAi/';

%% 45min-60min 37C heatshock
% Date, EmbryoID, Genotype, Stage, Result
date = 1 ;
embryoID = 2 ; 
genotype = 3 ;
stage = 4;
score = 5; 
minStage = 13.15 ;
maxStage = 15.2 ;


%% 202107211653: 67;15 x MLCK RNAi #4 (TRiP 20) 37C for 45min
% 0 --> WT
% 1 --> disrupted but folds
% 2 --> incomplete Fold
% 3 --> missing one Fold
% 4 --> missing two folds
% 5 --> missing all folds
% 6 --> no gut closure
% 7 --> dead
stage_result_mlck = [...
    202107211653, 1, 1,  14, 2;
    202107211653, 2, 1,  14, 3;
    202107211653, 3, 1,  15.2, 1;
    202107211653, 4, 1,  11, 6;
    202107211653, 5, 1,  13, 3;
    202107211653, 6, 1,  12, 5;
     ... % EXCLUDE SINCE DID NOT GET FAR ALONG ENOUGH 12, 5; % <-- 7. haltClosure
    202107211653, 8, 1,  15.3, 0;
    202107211653, 9, 1,  15.2, 0;
    202107211653, 10, 1,  15.1, 4;
    202107211653, 11, 1,  14, 5;  % e11
    202107211653, 12, 1,  15.2, 0;
    202107211653, 13, 1,  14, 4;
    202107211653, 14, 1,  14, 6;
    202107211653, 15, 1,  13, 5;
    202107211653, 16, 1,  13, 7;
    202107211653, 17, 1,  15.1, 4;
    202107211653, 18, 1,  10, 7;
    202107211653, 19, 1,  14, 7;
    202107211653, 20, 1,  11, 7;
    202107211653, 21, 1,  13, 6;
    202107211653, 22, 1,  15.1, 4;
    202107211653, 23, 1,  11, 6;
    202107211653, 24, 1,  13, 4; % <--- 24. how many missing? 2 folds missing (posterior anterior)
    202107211653, 25, 1,  15.1, 2;
    202107211653, 26, 1,  15.1, 7;
    202107211653, 27, 1,  15.2, 1;
    202107211653, 28, 1,  15.1, 1;
    202107211653, 29, 1,  15.2, 0;
    202107211653, 30, 1,  15.1, 2; % <-- 30. check, slowly folds
    202107211653, 31,  1, 13, 6;
     ... % EXCLUDE NOT FAR ENOUGH ALONG 12, 5; % <--- check, rename 30-> 32
    202107211653, 33, 1,  12, 5; % <-- check 33
    202107211653, 34, 1,  15.1, 0;
    202107211653, 35, 1,  12, 7;
    202107211653, 36, 1,  12, 7;
    202107211653, 37, 1,  13, 7;
    202107211653, 38, 1,  15.1, 5;
    202107211653, 39, 1,  15.2, 4; %<-- check seems like one fold  39 (already had one fold!)
    202107211653, 40, 1,  15.2, 3; 
    202107211653, 41, 1,  15.1, 7; % 41
    202107211653, 42, 1,  15.2, 0;
    202107211653, 43, 1,  12, 6;
    202107211653, 44, 1,  15.1, 1;
    202107211653, 45, 1,  14, 4; % 45
    202107211653, 46, 1,  15.1, 2; 
    202107211653, 47, 1,  13, 5;
    202107211653, 48, 1,  15.3, 0;
... %  202107221231_6715UASMLCKRNAi4_37Csince1220_2p5x10x_0p3pc_4mpf_15um
202107221231, 1, 1,     15.1,     0;
202107221231, 2, 1,     15.1,     0;
202107221231, 3, 1,     13,     0   ;
202107221231, 4, 1,     13,     0 ;
202107221231, 5, 1,     16.1,     0 ;
202107221231, 6, 1,     15.1,     0;
202107221231, 7, 1,     10,        7 ;
202107221231, 8, 1,     15.1,    1 ;
202107221231, 9, 1,     13,        0 ;
202107221231, 10, 1,     15.2,    2 ;
202107221231, 11, 1,     15.1,     0;
202107221231, 12, 1,     14.1,     0;
202107221231, 13, 1,     15,         0;
202107221231, 14, 1,     13,         x;
202107221231, 15, 1,     14.2,     0;
202107221231, 16, 1,     13,         1;
202107221231, 17, 1,     15.1,     2;
202107221231, 18, 1,     13,         1;
202107221231, 19, 1,     15.1,     4;
202107221231, 20, 1,     14.2,    4 ;
202107221231, 21, 1,     15.1,     0;
202107221231, 22, 1,     12,         6;
202107221231, 23, 1,     14.2,     3;
202107221231, 24, 1,     13,         1;
202107221231, 25, 1,     13,          1;  %<-- might be 2, can't see middle fold but shape seems ok
202107221231, 26, 1,     14.2 ,      3;  % <-- might be 4, can't see posterior fold, but may be there. Middle fold missing
202107221231, 27, 1,     14,          1;  % <-- highly disrupted shape
202107221231, 28, 1,     15.2,      0;
202107221231, 29, 1,     13 ,         0 ;
202107221231, 30, 1,     15.1,         1;
202107221231, 31, 1,     13,           0;
202107221231, 32, 1,     15.1 ,        0 ;
202107221231, 33, 1,     15.1,     0 ;
202107221231, 34, 1,     13.2 ,      0 ;
202107221231, 35, 1,     15.1,      0;
202107221231, 36, 1,     13,         5 ;
202107221231, 37, 1,     15.1,      0;
202107221231, 38, 1,     15.1 ,    0 ;
202107221231, 39, 1,     14.1,         0;
202107221231, 40, 1,     14.1,        0;
202107221231, 41, 1,     13,         0;
202107221231, 42, 1,     15.2 ,  0 ;
202107221231, 43, 1,     15.2  , 0 ;
202107221231, 44, 1,     15.1 , 0 ;
202107221231, 45, 1,     14.2, 0 ;
202107221231, 46, 1,     16.1,  0;
    ];

%%  202107241336 CONTROL 55min --  date, genotype, embryoID, starting stage, result
stage_result_6715WT = [...
202107241336, 1, 0,  13, 0; ...
202107241336, 2, 0,  15.1,0 ; ...
202107241336, 3, 0,  13, 4; ...
202107241336, 4, 0,  14.2, 0; ...
202107241336, 5, 0,  15.1, 0; ...
202107241336, 6, 0,  14.1, 2; ...
202107241336, 7, 0,  10, 5 ; ...
202107241336, 8, 0,  15.1, 0; ...
202107241336, 9, 0,  11, 5; ...
202107241336, 10, 0,  14.2, 0 ; ...
202107241336, 11, 0,  16.1, 0; ...
202107241336, 12, 0,  13, 3; ...
202107241336, 13, 0,  x, x; ...
202107241336, 14, 0,  15.1, 0; ...
202107241336, 15, 0,  14.2, 3; ...
202107241336, 16, 0,  11, 0; ...
202107241336, 17, 0,  11, 0; ...
202107241336, 18, 0,  14.2, 3; ...
202107241336, 19, 0,  12 , 5; ... % <--- initially upside down within vitelline
202107241336, 20, 0,  15.1, 0; ...
202107241336, 21, 0,  15.2, 0; ...
202107241336, 22, 0,  15.2, 2; ... % <-- incomplete anterior forld
202107241336, 23, 0,  15.1, 0; ...
202107241336, 24, 0,  12,   3; ... % <-- irregular folds, as if gut is sideways but embryo is not
202107241336, 25, 0,  15.1, 0; ...
202107241336, 26, 0,  13,   5; ...
202107241336, 27, 0,  15.1, 0; ...
202107241336, 28, 0,  x ,   x; ...
202107241336, 29, 0,  14.1, 3; ...
202107241336, 30, 0,  16.1, 3; ...
202107241336, 31, 0,  15.1, 0; ...
202107241336, 32, 0,  15.1, 0; ...
202107241336, 33, 0,  15.1, 0 ; ...
202107241336, 34, 0,  12,   2;  ...% <-- folds out of order, anterior is very late as to be incomplete
202107241336, 35, 0,  15.1, 0 ; ...
202107241336, 36, 0,  15.1, 0; ...
202107241336, 37, 0,  15.1, 0 ; ...
202107241336, 38, 0,  15.1,1 ; ...
202107241336, 39, 0,  12, 0; ...
202107241336, 40, 0,  13, 0; ...
202107241336, 41, 0,  15.1, 0 ; ...
202107241336, 42, 0,  15.1, 0; ...
202107241336, 43, 0,  11, 3 ;...
202107241336, 44, 0,  10, 0 ;...
202107241336, 45, 0,  15.1, 0 ;...
202107241336, 46, 0,  13, 7 ;...
202107241336, 47, 0,  15.2, 0;...
];
    
results = [stage_result_mlck; stage_result_6715WT];
    

%% STATS
inBin = find(results(:, stage) > minStage & ...
    results(:, stage) < maxStage) ;
mutantIdx = intersect(inBin, find(results(:, genotype) == 1)) ;
controlIdx = intersect(inBin, find(results(:, genotype) < 1)) ;

% Clean out those that are x or 7
keepM = find(results(mutantIdx, score) < scoreUpperLimit) ;
keepC = find(results(controlIdx, score) < scoreUpperLimit) ;
mutantIdx = mutantIdx(keepM) ;
controlIdx = controlIdx(keepC) ;

% N mutant and N control
nM = length(mutantIdx) ;
nC = length(controlIdx) ;

misMutant = intersect(mutantIdx, find(results(:, score) > scoreLowerLimit)) ;
foldMutant = intersect(mutantIdx, find(results(:, score) <= scoreLowerLimit)) ;
misControl = intersect(controlIdx, find(results(:, score) > scoreLowerLimit)) ;
foldControl = intersect(controlIdx, find(results(:, score) <= scoreLowerLimit)) ;
fracBadFolds_m = length(misMutant) / nM ;
fracBadFolds_c = length(misControl) / nC ;

assert(nM == length(misMutant) + length(foldMutant))
assert(nC == length(misControl) + length(foldControl))

nfoldM = length(foldMutant) ;
nmisM = length(misMutant) ;
nfoldC = length(foldControl) ;
nmisC = length(misControl) ;

% Significance
xx = table([nfoldM; nmisM],[nfoldC;nmisC],...
    'VariableNames',{'RNAi','control'},'RowNames',{'folded','misfolded'}) ;

[h,pval,stats] = fishertest(xx) ;

success = [nfoldM / nM, nfoldC/nC] ;
bar(success)
hold on;
yerr0 = sqrt((success .* (1-success)) ./ [nM, nC]) ;
yneg = min(yerr0, success) ;
ypos = min(yerr0, 1-success) ;
errorbar([1,2], success, yneg, ypos , 'LineStyle','none')
set(gca,'xticklabel',{'UAS-MLCK RNAi', 'control'});
ylabel('probablity of forming three folds')
    
if ~exist(outdir, 'dir')
    mkdir(outdir)
end
saveas(gcf, fullfile(outdir, 'rnai_60min37C.pdf')) ;
    

%% superbar plot
hf = figure('Position', [100 100 400 400], 'units', 'centimeters');
clf;
Y = success;
E = cat(3, yneg, ypos) ;

Colors = [
    0.90    0.55    0.55
    0.62    0.76    0.84
    0.89    0.10    0.11
    0.12    0.47    0.70
    ];
Colors = reshape(Colors, [2 2 3]);

P = [pval, pval; pval, pval];
% Make P symmetric, by copying the upper triangle onto the lower triangle
% PT = P';
% lidx = tril(true(size(P)), -1);
% P(lidx) = PT(lidx);

superbar(Y, 'E', E, 'P', P, 'BarFaceColor', Colors, 'Orientation', 'v', ...
    'ErrorbarStyle', 'I', 'PLineOffset', 0.1, 'PStarShowGT', false);

xlim([0.5 2.5]);
% ylim([0 1]);
set(gca, 'YTick', [0, 0.5, 1])
set(gca, 'XTick', [1, 2]);
ylims = ylim;
expVal = sprintf('%e', pval) ;
keepDigits = expVal(end-2:end) ;
text(1.5, ylims(2)-0.03, ...
    ['$p=$' sprintf(['%0.' num2str(keepDigits) 'f'],pval)], 'interpreter', 'latex', 'horizontalalignment', 'center')


set(gca,'xticklabel',{'UAS-MLCK RNAi', 'control'});
ylabel('probablity of forming three folds')
title('tub>67;tub>15 x UAS-MLCK RNAi')

if ~exist(outdir, 'dir')
    mkdir(outdir)
end
saveas(gcf, fullfile(outdir, 'rnai_60min37C_superbar.pdf')) ;


%% Print contingency table
xx
    