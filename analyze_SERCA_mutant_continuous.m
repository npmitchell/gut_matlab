
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

%202107211653: 67;15 x MLCK RNAi #4 (TRiP 20) 37 for 45min
% 0 --> WT
% 1 --> disrupted but folds
% 2 --> incomplete Fold
% 3 --> missing one Fold
% 4 --> missing two folds
% 5 --> missing all folds
% 6 --> no gut closure
% 7 --> dead before folds
% x --> did not finish/ cannot interpret
close all 
clear
x = 100 ;
scoreLowerLimit = 1 ;
scoreUpperLimit = 6 ;
outdir = '/mnt/data/analysis/sercaR751Q/';
addpath('/mnt/data/code/gut_matlab/addpath_recurse')
addpath_recurse('/mnt/data/code/gut_matlab/plotting/')

%% continuous 37C heatshock for SERCA
% Date, EmbryoID, Genotype, Stage, Result
date = 1 ;
embryoID = 2 ; 
genotype = 3 ;
stage = 4;
score = 5; 
minStage = 13.1 ;
maxStage = 15.2 ;
genoM = 2; 
genoC = 0.3 ;


%% CONTINUOUS 37C heatshock

% Date, EmbryoID, Genotype, Stage, Result
serca_stage_result_continuous = [...
    ... % 202107292013 continuous 37C 2 min /frame mef2G4kSERCA.R751Q and 48YG4kCAAXmCh
202107292013, 1, 2, 14.2, 3 ; ... %e1
202107292013, 2, 2, 14.2, 5; ...
202107292013, 3, 2, 13, 4; ...
202107292013, 4, 0.2, 13, 0; ...
202107292013, 5, 0.2, 13, 1; ...% e5 --> all folds present but middle fold is late (simultaneous) so shape is upset
202107292013, 6, 0.2, 14.1, 0; ...
202107292013, 7, 0.2, 15.12, 0; ...   % e7 
202107292013, 8, 0.2, 14.2, 3; ...  % e8
202107292013, 9, 0.2, 15.1, 0; ... %e9
    ...  % 202107301429 continuous 37C 2 min /frame mef2G4kSERCA.R751Q and 48YG4kCAAXmCh
202107301429, 1, 0.2, 14.1, 4; ... %e1
202107301429, 2, 0.2, 15.1, 0;... % e2
202107301429, 3, 0.2, 15.1, 0; ... %e3 
202107301429, 4, 0.2, 13, 1 ; ... %e4
202107301429, 5, 0.2, 15, 7 ; ...
202107301429, 6, 2, 15.1, 0; ...
202107301429, 7, 0.2, 15.1, 0; ...
202107301429, 8, 2, 15.1, 0 ; ... % e8
202107301429, 9, 2, 13, 3; ... %e9
202107301429, 10, 0.2, 13, 2; ...
202107301429, 11, 0.2, 14.1, 7;...
 ... % 202108030034 Broida 48YG4kxUASCAAXmCh table: date, embryoID, genotype, result
202108030034, 1, 0.2, 13, 5;...
202108030034, 2, 0.2, 15.1, 0 ;...
202108030034, 3, 0.2, 16.1, 0;...
202108030034, 4, 0.2, 15.1, 0;...
202108030034, 5, 0.2, 15.1, 0;...
202108030034, 6, 0.2, 16.1, 0;...
202108030034, 7, 0.2, 16.1, 0;...
202108030034, 8, 0.2, 13, 1;...  % z3
202108030034, 9, 0.2, 13, 5;...
202108030034, 10, 0.2, 13, 4;...
202108030034, 11, 0.2, 14, 1;...
202108030034, 12, 0.2, 16.1, 0;...
202108030034, 13, 0.2, 15.2, 0;...
202108030034, 14, 0.2, 14.2, 4;... 
202108030034, 15, 0.2, 16.1, 0;...
202108030034, 16, 0.2, 15.1, 0 ;...
202108030034, 17, 0.2, 15.1, 0;...   % < -- nice poster 
202108030034, 18, 0.2, 15.1, 0;...
202108030034, 19, 0.2, 15.1, 0;...   % < --- nice poster frame 15a
202108030034, 20, 0.2, 16.1, 0;...
202108030034, 21, 0.2, 13, 3;...
202108030034, 22, 0.2, 15.1, 0;...  % <-- poster
202108030034, 23, 0.2, 16.1, 0;...
202108030034, 24, 0.2, 13, 4 ;...
202108030034, 25, 0.2, 15.1, 2;...  % posterior fold incomplete, shape is odd
202108030034, 26, 0.2, 15.1, 0;...
202108030034, 27, 0.2, 15.1, 0;...
202108030034, 28, 0.2, 13, 0;...
202108030034, 29, 0.2, 15.1, 0;...% nice
];



%% 3.5 hrs 37C
serca_stage_result_3p5hrsOrMore37C = [serca_stage_result_continuous; ...
... % 20210720 -- 3.5 hrs at 37C
20210720, 1, 2, 15.1, 0;...  %e1         
20210720, 2, 2, 16.1, 0;...  %e2                 
20210720, 3, 2, 16.1, 0;...  %e3                 
20210720, 4, 2, 16.1, 0;...  %e4                 
20210720, 5, 2, 15.1, 1;...  %e5                 
20210720, 6, 2, 13, 5 ;...    %e6                 
20210720, 7, 2, 13, 5 ;...    %e7    <---- check what happened             
20210720, 8, 2, 10, x;...    %e8  --> ends at 12         
20210720, 9, 2, 10, x;...    %e9  --> ends at 12          
20210720, 10, 2, 11, x;...    %e10 --> ends at 13  
20210720, 11, 2, 12, x;...    %e11 --> ends at 13    
20210720, 12, 2, 11, 3;...    %e12 --> ends at 14 <-- misses middle fold at 25C hrs later
20210720, 13, 2, 10, x;...    %e13 --> ends at 12            
20210720, 14, 2, 11, x;...    %e14 --> ends at 13
20210720, 15, 2, 13.2, x;...  %e15 -> ends at 15.a
20210720, 16, 2, 12, x;...    %e16 --> ends at 14
...  % 202108021250 -- 3.5 hrs at 37C --> I have the timetrace of the temperature for this one one line per 2.3 seconds
202108021250, 1, 2, 13.2, 5;     % no folds until HS ends, but then folds somewhat
202108021250, 2, 2, 16.1, 0;
202108021250, 3, 2, 16.1, 0;
202108021250, 4, 2, 16.2, 0;
202108021250, 5, 2, 15.1, 3;     % e5
202108021250, 6, 2, 13.2, 5;     % forms all folds after HS
202108021250, 7, 2, 13.1, 5;     % forms posterior fold after HS
202108021250, 8, 2, 13.2, 5; 
202108021250, 9, 2, 15.1, 4; 
202108021250, 10, 2, 14, 4;       % e10
202108021250, 11, 2, 10, x;
202108021250, 12, 2, 13, 5; 
202108021250, 13, 2, 15.1, 3;    % e13
202108021250, 14, 2, 13.2, 5;
202108021250, 15, 2, 14, 2;       % two folds incomplete, really all incomplete
202108021250, 16, 2, 16.1, 0;
202108021250, 17, 2, 13, 3;
202108021250, 18, 2, 14.2, 2;
202108021250, 19, 2, 15.1, 2;
202108021250, 20, 2, 15.1, 4;
202108021250, 21, 2, 15.1, 4; 
202108021250, 22, 2, 15.1, 3;
202108021250, 23, 2, 14, 5;
202108021250, 24, 2, 13, 5;
202108021250, 25, 2, 13, 5;
202108021250, 26, 2, 15.1, 0;
202108021250, 27, 2, 13, 5;
202108021250, 28, 2, 13, 5;
202108021250, 29, 2, 15.1, 0;
202108021250, 30, 2, 13, 4;
202108021250, 31, 2, 13, 5;
202108021250, 32, 2, 14, 5;
202108021250, 33, 2, 15.1, 5;
202108021250, 34, 2, 15.1, 3; % completely missing anterior fold 
202108021250, 35, 2, 15.2, 3;
202108021250, 36, 2, 14.2, 5;
202108021250, 37, 2, 14.1, 5;
202108021250, 38, 2, 14, 3; 
202108021250, 39, 2, 14, 3;
202108021250, 40, 2, 13, 5;
202108021250, 41, 2, 13, 5;
202108021250, 42, 2, 15.1, 0; % "Odd" folds
202108021250, 43, 2, 15.2, 0; 
202108021250, 44, 2, 13, 5;
202108021250, 45, 2, 15.1, 0;
202108021250, 46, 2, 15.1, 5;  % strikingly folds post HS
202108021250, 47, 2, 13, 5;
... % 202108042044_mef2G4kSERCA1_3mpf_7pc488_bioII_37Ccontinuous2044_2p5x10x.lif 
202108042044, 1, 2, 15.1, 1;
202108042044, 2, 2, 15.2, 1 ;
202108042044, 3, 2, 14.1, 3;
202108042044, 4, 2, 14.1, 0;
202108042044, 5, 2, 13.2, 1;
202108042044, 6, 2, 15.1, 1;
202108042044, 7, 2, 15.1, 2;   %<-- revisit
202108042044, 8, 2, 13, 5;
202108042044, 9, 2, 14.2, 0;
202108042044, 10, 2, 13, 1;
202108042044, 11, 2, 16.1, 2 ;  % <--anterior fold did not finish
202108042044, 12, 2, 14.2, 1;
202108042044, 13, 2, 13, 1 ;
202108042044, 14, 2, 13, 2;
202108042044, 15, 2, 15.1, 0 ;
202108042044, 16, 2, 15.1, 2;  %<-- could be 1, anterior fold is incomplete but does curl dorsally
202108042044, 17, 2, 15.2, 0;
202108042044, 18, 2, 14.1, 2;  %<-- anterior fold does not finish but does curl dorsally
202108042044, 19, 2, 13, 1;
202108042044, 20, 2, 13.2, 1;
202108042044, 21, 2, 14.1, 4 ;   %<-- absolutely no anterior fold, posterior fold is much reduced
202108042044, 22, 2, 14.1, 0 ; ...
... % 20210808 Broida mef2G4k x SERCAR751Q table:
202108081302, 1, 2, 14.1, 3; ...
202108081302, 2, 2, 14.1, 3; ... % 13.2 or 14.1 ? Band complete
202108081302, 3, 2, 14.1, 2; ...
202108081302, 4, 2, 13, 5; ...
202108081302, 5, 2, 15.1, 3; ...
202108081302, 6, 2, 13, 4; ...
202108081302, 7, 2, 15.1, 2; ...
202108081302, 8, 2, 14.1, 3; ...
202108081302, 9, 2, 13, 1; ...
202108081302, 10, 2, 14.1, 3; ...  % 13.2 or 14.1? Band complete
202108081302, 11, 2, 15.1, 2; ...
202108081302, 12, 2, 15.1, 3; ...
202108081302, 13, 2, 14.1, 1; ...
202108081302, 14, 2, 15.2, 1; ...
202108081302, 15, 2, 14.1, 0; ...
202108081302, 16, 2, 15.1, 0; ...
202108081302, 17, 2, 13.2, 3 ;...
];

    
    
serca_stage_result_3p5hrsOrMore37C = [serca_stage_result_3p5hrsOrMore37C;...
... % 20210803 bioII mef2 table: date, embryoID, initial HS stage, all normal development (0)
    202108031752, 1, 0.3, 13, 0 ; ... 
    202108031752, 2, 0.3, 14, 0 ; ... 
    202108031752, 3, 0.3, 16.1, 0 ; ... 
    202108031752, 4, 0.3, 12, 0 ; ... 
    202108031752, 5, 0.3, 16.1, 0 ; ... 
    202108031752, 6, 0.3, 16.1, 0 ; ... 
    202108031752, 7, 0.3, 13, 0 ; ... 
    202108031752, 8, 0.3, 16.1, 0 ; ... 
    202108031752, 9, 0.3, 15.1, 0 ; ... 
... % 20210803 bioII tub67;15 table:
    202108031752, 10, 0, 15.1, 0 ; ... 
    202108031752, 11, 0, 15.1, 0 ; ... 
    202108031752, 12, 0, 15.1, 0 ; ... 
    202108031752, 13, 0, 15.1, 0 ; ... 
    202108031752, 14, 0, 11, x ;... <---- not done yet at time of viewing (21:32)
    202108031752, 15, 0, 15.1, 0 ; ... 
    202108031752, 16, 0, 15.2, 0 ; ... 
    202108031752, 17, 0, 13, 0 ; ... 
    202108031752, 18, 0, 15.2, 0 ; ... 
    202108031752, 19, 0, 15.1, 0 ; ... 
    202108031752, 20, 0, 13, 0 ; ... 
    202108031752, 21, 0, 15.1, 0 ; ... 
... % 20210805 Broida mef2;mef2 table:
    202108051250, 1, 0.3, 16.1, 0; ...
    202108051250, 2, 0.3, 15.1, 0; ...
    202108051250, 3, 0.3, 15.1, 2; ...
    202108051250, 4, 0.3, 15.1, 1; ...
    202108051250, 5, 0.3, 13, 3; ...
    202108051250, 6, 0.3, 13, 4; ...
    202108051250, 7, 0.3, 15.2, 1; ...
    202108051250, 8, 0.3, 15.1, 0; ...
    202108051250, 9, 0.3, 13, 4; ...
    202108051250, 10, 0.3, 14, 2; ...
    202108051250, 11, 0.3, 14.1, 2 ;...
    202108051250, 12, 0.3, 15.2, 0; ...
    202108051250, 13, 0.3, 15.1, 0 ;...
    202108051250, 14, 0.3, 15.2, 1 ;...
    202108051250, 15, 0.3, 14.2, 1 ;...
    202108051250, 16, 0.3, 15.1, 0 ;...
    202108051250, 17, 0.3, 15.1, 0 ;...
    202108051250, 18, 0.3, 14, 0 ;...
    202108051250, 19, 0.3, 13, 4 ;...
    202108051250, 20, 0.3, 13.2, 2 ;...
    202108051250, 21, 0.3, 14.1, 0 ;...
    202108051250, 22, 0.3, 15.1, 0 ;...
    202108051250, 23, 0.3, 13.2, 3 ;...
    202108051250, 24, 0.3, 13.2, 1 ;...
...% 202108101608_mef2G4k_6715cntrl_37Ccontinuous_boiII_4mpf_20um_7pc488
202108101608, 23, 0.3, 15.1, 1 ;
202108101608, 24, 0.3, 13, 4 ;
202108101608, 25, 0.3, 14.2, 2;
202108101608, 26, 0.3, 13, 1;
202108101608, 27, 0.3, 15.1, 1;
202108101608, 28, 0.3, 15.1 , 1;
202108101608, 29, 0.3, 13.2, 1;
202108101608, 30, 0.3, 15.1, 0;
202108101608, 31, 0.3, 13, 0;
...% tub67;15 table:
...% 202108101608_mef2G4k_6715cntrl_37Ccontinuous_boiII_4mpf_20um_7pc488
202108101608, 1, 0, 16.1, 1;
202108101608, 2, 0, 15.1, 2 ; % incomplete posterior
202108101608, 3, 0, 15.1, 2;
202108101608, 4, 0, 16.1, 1;
202108101608, 5, 0, 14.2, 3;  % difficult to see, dark, but definitely not normal. At the very least incomplete posterior
202108101608, 6, 0, 15.1, 1 ;
202108101608, 7, 0, 13, 6;
202108101608, 8, 0, 14.1, 1;
202108101608, 9, 0, 13.2,  5;
202108101608, 10, 0, 15.1, 2;
202108101608, 11, 0, 15.2, 1;
202108101608, 12, 0, 15.1, 0; 
202108101608, 13, 0, 14.1, 0; 
202108101608, 14, 0, 15.1, 0; 
202108101608, 15, 0, 15.1, 0;
202108101608, 16, 0, 15.1, 0
202108101608, 17, 0, 13, 0;
202108101608, 18, 0, 14.1, 0;
202108101608, 19, 0, 15.1, 0;
202108101608, 20, 0, 15.1, 0 ;
202108101608, 21, 0, 14.2, 0;
202108101608, 22, 0, 15.1, 0;
];



% 0 --> WT
% 1 --> disrupted but folds
% 2 --> incomplete Fold(s)
% 3 --> missing one Fold
% 4 --> missing two folds
% 5 --> missing all folds
% 6 --> no gut closure
% 7 --> dead before folds
% x --> did not finish/ cannot interpret
results = serca_stage_result_3p5hrsOrMore37C;

%% Quantized stats
for ii = 1:2
    
    close all
    inBin = find(results(:, stage) > minStage & ...
        results(:, stage) < maxStage) ;
    mutantIdx = intersect(inBin, find(results(:, genotype) == genoM)) ;
    
    if ii == 1
        exten = sprintf('_restricted_min%0.1f_max%0.1f', minStage, maxStage) ;
        controlIdx = intersect(inBin, find(results(:, genotype) == genoC)) ;
    else
        exten = sprintf('_allControls_min%0.1f_max%0.1f', minStage, maxStage) ;
        controlIdx = intersect(inBin, find(results(:, genotype) < 1)) ;
    end
    exten = strrep(exten, '.', 'p') ;
    
    % Clean out those that are x or 7
    keepM = find(results(mutantIdx, score) < scoreUpperLimit) ;
    keepC = find(results(controlIdx, score) < scoreUpperLimit) ;
    mutantIdx = mutantIdx(keepM) ;
    controlIdx = controlIdx(keepC) ;

    % N mutant and N control
    nM = length(mutantIdx) ;
    nC = length(controlIdx) ;
    nT = nM + nC ;

    mutant0folds = intersect(mutantIdx, find(results(:, score) == 5)) ;
    mutant1folds = intersect(mutantIdx, find(results(:, score) == 4)) ;
    mutant2folds = intersect(mutantIdx, find(results(:, score) == 3 | results(:, score) == 2)) ;
    mutant3folds = intersect(mutantIdx, find(results(:, score) <= 1)) ;
    control0folds = intersect(controlIdx, find(results(:, score) == 5)) ;
    control1folds = intersect(controlIdx, find(results(:, score) == 4)) ;
    control2folds = intersect(controlIdx, find(results(:, score) == 3 | results(:, score) == 2)) ;
    control3folds = intersect(controlIdx, find(results(:, score) <= 1)) ;

    nM0 = numel(mutant0folds) ;
    nM1 = numel(mutant1folds) ;
    nM2 = numel(mutant2folds) ;
    nM3 = numel(mutant3folds) ;
    nC0 = numel(control0folds) ;
    nC1 = numel(control1folds) ;
    nC2 = numel(control2folds) ;
    nC3 = numel(control3folds) ;
    hisMutant = [nM0, nM1, nM2, nM3] ;
    hisControl = [nC0, nC1, nC2, nC3] ;
    
    nfolds = [0, 1, 2, 3] ;
    colors = define_colors ;
    colors = [
        0.90    0.55    0.55 ;
        0.62    0.76    0.84 ] ;
    bb = barh(nfolds, [ hisMutant', hisControl'], 'histc') ;
    bb(1).FaceColor = colors(1, :);
    bb(2).FaceColor = colors(2, :);
    xlabel('count', 'interpreter', 'latex')
    ylabel('number of folds', 'interpreter', 'latex')
    legend({'UAS-SERCA.R751Q', 'control'})
    saveas(gcf, fullfile(outdir, ['serca_continuous37C_barh' exten '.pdf'])) ;
    saveas(gcf, fullfile(outdir, ['serca_continuous37C_barh' exten '.png'])) ;
    
    boxMutant = [0*ones(nM0, 1); ...
                 1*ones(nM1, 1); ...
                 2*ones(nM2, 1); ...
                 3*ones(nM3, 1)] ;
    boxControl = [0*ones(nC0, 1); ...
                 1*ones(nC1, 1); ...
                 2*ones(nC2, 1); ...
                 3*ones(nC3, 1)] ;
    g1 = repmat({'UAS-SERCA.R751Q'},length(boxMutant),1);
    g2 = repmat({'control'},length(boxControl),1);
    gg = [g1; g2];
    boxplot([boxMutant; boxControl], gg)
    ylabel('number of folds')
    saveas(gcf, fullfile(outdir, ['serca_continuous37C_boxp' exten '.pdf'])) ;
    saveas(gcf, fullfile(outdir, ['serca_continuous37C_boxp' exten '.png'])) ;
   
    % attempt 1--a bit awkward, uses barh and subtightplot
    %     clf
    %     marg_h = 0.1 ;
    %     marg_w = 0.1 ;
    %     ax1 = subtightplot(1,2, 1,[0,0], marg_h, marg_w)
    %     vMutant = hisMutant ./ sum(hisMutant) ;
    %     bm = barh(nfolds,0.5*[vMutant;-vMutant;-vMutant]', 'stacked')
    %     xlim([-0.5, 0.5])
    %     bm(1).FaceColor = colors(1, :);
    %     bm(2).FaceColor = colors(1, :);
    %     bm(3).FaceColor = colors(1, :);
    %     bm(1).EdgeColor = 'none';
    %     bm(2).EdgeColor = 'none';
    %     bm(3).EdgeColor = 'none';
    %     ylabel('number of folds')
    %     xlabel('frequency')
    %     xticks([])
    %     
    %     ax2 = subtightplot(1,2, 2,[0,0], marg_h, marg_w)
    %     vControl = hisControl ./ sum(hisControl) ;
    %     bc = barh(nfolds,0.5*[vControl;-vControl;-vControl]', 'stacked') ;
    %     xlim([-0.5, 0.5])
    %     bc(1).FaceColor = colors(2, :);
    %     bc(2).FaceColor = colors(2, :);
    %     bc(3).FaceColor = colors(2, :);
    %     bc(1).EdgeColor = 'none';
    %     bc(2).EdgeColor = 'none';
    %     bc(3).EdgeColor = 'none';
    %     xticks([])
    %     
    %     linkaxes([ax1, ax2], 'y');
    
    
    %%
    close all
    hf = figure('Position', [100 100 400 400], 'units', 'centimeters');
    
    nMutant = sum(hisMutant) ;
    nControl = sum(hisControl) ;
    vMutant = hisMutant ./ sum(hisMutant) ;
    vControl = hisControl ./ sum(hisControl) ;
    gg = [0,0,0,0, 1, 1, 1, 1];
    % boxyViolinPlot(gg, [nfolds, nfolds], ...
    %     [vMutant, vControl], 1, 'center') ;
    boxyViolinPlot([0,0,0,0], nfolds, ...
        vMutant, 1, 'center', colors(1, :), 'none') ;
    boxyViolinPlot([1,1,1,1], nfolds, ...
        vControl, 1, 'center', colors(2, :), 'none') ;
    
    % Significance
    num = mean(boxMutant) -  mean(boxControl)  ;
    denom = sqrt(nanvar(boxControl) / nControl + ...
        nanvar(boxMutant) / nMutant) ;
    zscoreNumFolds = num / denom ;
    pvalueNumFolds = normcdf(zscoreNumFolds) ;
    
    success = [ mean(boxMutant) ,  mean(boxControl) ] ;
    yerr0 = [std(boxMutant) / sqrt(nMutant), std(boxControl) / sqrt(nControl) ];
    yneg = min(yerr0, success) ;
    ypos = min(yerr0, 3-success) ;
    
    
    % finish plot with labels
    ylim([-0.5, 4])
    set(gca, 'YTick', [0, 1, 2, 3])
    set(gca, 'XTick', [0, 1]);
    expVal = sprintf('%e', pvalueNumFolds) ;
    keepDigits = expVal(end-1:end) ;
    text(0.5, 3.8, ...
        ['$p=$' sprintf(['%0.' num2str(keepDigits) 'f'],pvalueNumFolds)], 'interpreter', 'latex', 'horizontalalignment', 'center')

    set(gca,'xticklabel',{'UAS-SERCA.R751Q', 'control'});
    ylabel('number of folds')
    
    saveas(gcf, fullfile(outdir, ['serca_continuous37C_violin' exten '.pdf'])) ;
    saveas(gcf, fullfile(outdir, ['serca_continuous37C_violin' exten '.png'])) ;
    
    
    %% 
    hold on;
    errorbar([0, 1], [mean(boxMutant), mean(boxControl)], ...
        yneg, ypos, 'LineWidth', 2, 'LineStyle','none', 'color', [0.5, 0.5,0.5])
   
    xlim([-0.5 1.5]);
    ylim([-0.5, 4])
    set(gca, 'YTick', [0, 1, 2, 3])
    set(gca, 'XTick', [0, 1]);
    
    set(gca,'xticklabel',{'UAS-SERCA.R751Q', 'control'});
    ylabel('number of folds')
    % title('mef2-GAL4;klar x UAS-SERCA.R751Q')

    if ~exist(outdir, 'dir')
        mkdir(outdir)
    end
    saveas(gcf, fullfile(outdir, ['serca_continuous37C_violinErr' exten '.pdf'])) ;
    saveas(gcf, fullfile(outdir, ['serca_continuous37C_violinErr' exten '.png'])) ;

    
    
    %% Flip the histogram sideways
    
    close all
    hf = figure('Position', [100 100 400 400], 'units', 'centimeters');
    pos = get(gca, 'Position') ;
    set(gca, 'Position', pos + [0, 0.15, 0, -0.1])
    boxyViolinPlot([0,.2,.4,.6], 0, ...
        0.2, vMutant, 'left', colors(1, :), 'none') ;
    boxyViolinPlot(1 + [0,.2,.4,.6], 0, ...
        0.2, vControl, 'left', colors(2, :), 'none') ;
    
    xticks([0.1, 0.3, 0.5, 0.7, 1.1, 1.3, 1.5, 1.7])
    xticklabels([0,1,2,3, 0,1,2,3])
    yticks([0, 0.5, 1])
    xlabel('number of folds')
    ylabel('frequency')
    xlim([-0.2, 2])
    ylim([0, 1])
    text(0.9, 0.7, ...
        ['$p=$' sprintf(['%0.' num2str(keepDigits) 'f'],pvalueNumFolds)], 'interpreter', 'latex', 'horizontalalignment', 'center')

    text(0.4, 0.6, ['N=', num2str(nMutant)], 'HorizontalAlignment', 'center')
    text(1.4, 0.6, ['N=', num2str(nControl)], 'HorizontalAlignment', 'center')
    text(0.4, -0.2, 'UAS-SERCA.R751Q', 'HorizontalAlignment', 'center')
    text(1.4, -0.2, 'control', 'HorizontalAlignment', 'center')
    saveas(gcf, fullfile(outdir, ['serca_continuous37C_clusterHist' exten '.pdf'])) ;
    saveas(gcf, fullfile(outdir, ['serca_continuous37C_clusterHist' exten '.png'])) ;

end

%% STATS
for ii = 1:2
    
    inBin = find(results(:, stage) > minStage & ...
        results(:, stage) < maxStage) ;
    mutantIdx = intersect(inBin, find(results(:, genotype) == genoM)) ;
    
    if ii == 1
        exten = sprintf('_restricted_min%0.1f_max%0.1f', minStage, maxStage) ;
        controlIdx = intersect(inBin, find(results(:, genotype) == genoC)) ;
    else
        exten = sprintf('_allControls_min%0.1f_max%0.1f', minStage, maxStage) ;
        controlIdx = intersect(inBin, find(results(:, genotype) < 1)) ;
    end
    exten = strrep(exten, '.', 'p') ;
    

    % Clean out those that are x or 7
    keepM = find(results(mutantIdx, score) < scoreUpperLimit) ;
    keepC = find(results(controlIdx, score) < scoreUpperLimit) ;
    mutantIdx = mutantIdx(keepM) ;
    controlIdx = controlIdx(keepC) ;

    % N mutant and N control
    nM = length(mutantIdx) ;
    nC = length(controlIdx) ;
    nT = nM + nC ;

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
        'VariableNames',{'mutant','control'},'RowNames',{'folded','misfolded'}) ;

    [h,pval,stats] = fishertest(xx) ;

    success = [nfoldM / nM, nfoldC/nC] ;
    bar(success)
    hold on;
    yerr0 = sqrt((success .* (1-success)) ./ [nM, nC]) ;
    yneg = min(yerr0, success) ;
    ypos = min(yerr0, 1-success) ;
    errorbar([1,2], success, yneg, ypos , 'LineStyle','none')
    set(gca,'xticklabel',{'UAS-SERCA.R751Q', 'Mef2-GAL4'});
    ylabel('probablity of forming three folds')

    if ~exist(outdir, 'dir')
        mkdir(outdir)
    end
    saveas(gcf, fullfile(outdir, ['serca_continuous37C' exten '.pdf'])) ;
    saveas(gcf, fullfile(outdir, ['serca_continuous37C' exten '.png'])) ;


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

    P = [pval,  pval; pval, pval];
    % Make P symmetric, by copying the upper triangle onto the lower triangle
    % PT = P';
    % lidx = tril(true(size(P)), -1);
    % P(lidx) = PT(lidx);

    superbar(Y, 'E', E, 'P', P, 'BarFaceColor', Colors, 'Orientation', 'v', ...
        'ErrorbarStyle', 'I', 'PLineOffset', 0.1, 'PStarShowGT', false);

    xlim([0.5 2.5]);
    ylims = ylim ;
    if ylims(2) < 1
        ylim([0 1]);
    end
    set(gca, 'YTick', [0, 0.5, 1])
    set(gca, 'XTick', [1, 2]);
    ylims = ylim;
    expVal = sprintf('%e', pval) ;
    keepDigits = expVal(end-1:end) ;
        text(1.5, ylims(2)-0.03, ...
        ['$p=$' sprintf(['%0.' num2str(keepDigits) 'f'],pval)], 'interpreter', 'latex', 'horizontalalignment', 'center')


    set(gca,'xticklabel',{'UAS-SERCA.R751Q', 'control'});
    ylabel('probablity of forming three folds')
    title('mef2-GAL4;klar x UAS-SERCA.R751Q')

    if ~exist(outdir, 'dir')
        mkdir(outdir)
    end
    saveas(gcf, fullfile(outdir, ['serca_continuous37C_superbar' exten '.pdf'])) ;
    saveas(gcf, fullfile(outdir, ['serca_continuous37C_superbar' exten '.png'])) ;


    %% Print contingency table
    xx

    fn = fullfile(outdir, ['serca_continuous37C_contingency' exten '.txt']) ;
    writetable( xx, fn, 'WriteRowNames', true)

    
    
end

error('here')
%% Util -- get temperature trace
temp_fn = '/mnt/data/RNAi/Mef2GAL4klarUASSERCA1/20210802_temperature_trace_until2416ish_startingAround1250.txt' ;
%temp_fn = '/mnt/data/RNAi/Mef2GAL4klarUASSERCA1/20210802_temperature_trace_until1904.txt' ;
temperatures0 = dlmread(temp_fn, ' ', 2, 0);
idx = round(4200:180/2.35:length(temperatures0)) ;
temperatures = temperatures0(idx) ;
timestamp = 0:3:3*(length(temperatures)-1) ;
plot(timestamp/60, temperatures, '.-')
dlmwrite([temp_fn(1:end-4) '_every3min.txt'], temperatures, ' ')

part1 = 'run("Label...", "format=00:00:00 starting=';
part2 = ' interval=0 x=5 y=41 font=36 text=[' ;
part3 = '°C] range=';
fid = fopen([temp_fn(1:end-4) '_every3min_macro.txt'],'w');
for jj = 1 : size( temperatures, 1 )
    starting = num2str(17 + 3 * (jj-1)) ;
    framenum = num2str(jj) ;
    temp = sprintf('%d', round(temperatures(jj))) ;
    fprintf( fid, '%s\n', [part1 starting part2 temp part3 framenum '-' framenum ' use"); ']);
end
fclose(fid);