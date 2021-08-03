
% Here, 
%   13 is bands forming, late 13 is bands formed
%   14a is closing
%   14b is mostly closed but not fully
%   15a is after closure
%   15b is after/during first fold

% 1: mef2G4klar x UAS-SERCA.R751Q
% 0: 48YGAL4 klar x UAS-CAAX mCh

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
x = 100 ;
scoreLowerLimit = 1 ;
scoreUpperLimit = 7 ;


% 202107291353 1 hr 37C 2 min /frame
serca_stage_result_1hr = [...
1, 13,  0; ... %e1
0, 13, 0; ...
1, 14.1, 0; ...
0, 14.1, 0; ...
1, 16.1, 0; ...% e5
0, 16.1, 0; ...
1, 15.1, 3; ...   % e7 disrupted
0, 15.1, 0; ...
1, 17, 0; ... %e9
1, 17.1, 0; ...
1, 17.2, 0; ...
0, 15.1, 0; ...
0, 14.1, 0; ...
1, 13, 0; ...
0, 14, 0; ...
];

%% CONTINUOUS 37C heatshock

% Date, EmbryoID, Genotype, Stage, Result
date = 1 ;
embryoID = 2 ; 
genotype = 3 ;
stage = 4;
score = 5; 
serca_stage_result_continuous = [...
    ... % 202107292013 continuous 37C 2 min /frame
202107292013, 1, 1, 14.2, 3 ; ... %e1
202107292013, 2, 1, 14.2, 5; ...
202107292013, 3, 1, 13, 4; ...
202107292013, 4, 0, 13, 0; ...
202107292013, 5, 0, 13, 1; ...% e5 --> all folds present but middle fold is late (simultaneous) so shape is upset
202107292013, 6, 0, 14.1, 0; ...
202107292013, 7, 0, 15.12, 0; ...   % e7 
202107292013, 8, 0, 14.2, 3; ...  % e8
202107292013, 9, 0, 15.1, 0; ... %e9
    ...  % 202107301429 continuous 37C 2 min /frame
202107301429, 1, 0, 14.1, 4; ... %e1
202107301429, 2, 0, 15.1, 0;... % e2
202107301429, 3, 0, 15.1, 0; ... %e3 
202107301429, 4, 0, 13, 1 ; ... %e4
202107301429, 5, 0, 15, 7 ; ...
202107301429, 6, 1, 15.1, 0; ...
202107301429, 7, 0, 15.1, 0; ...
202107301429, 8, 1, 15.1, 0 ; ... % e8
202107301429, 9, 1, 13, 3; ... %e9
202107301429, 10, 0, 13, 2; ...
202107301429, 11, 0, 14.1, 7;...
];

serca_stage_result_3p5hrs37C = [serca_stage_result_continuous; ...
... % 20210720 -- 3.5 hrs at 37C
20210720, 1, 1, 15.1, 0;...  %e1         
20210720, 2, 1, 16.1, 0;...  %e2                 
20210720, 3, 1, 16.1, 0;...  %e3                 
20210720, 4, 1, 16.1, 0;...  %e4                 
20210720, 5, 1, 15.1, 1;...  %e5                 
20210720, 6, 1, 13, 5 ;...    %e6                 
20210720, 7, 1, 13, 5 ;...    %e7    <---- check what happened             
20210720, 8, 1, 10, x;...    %e8  --> ends at 12         
20210720, 9, 1, 10, x;...    %e9  --> ends at 12          
20210720, 10, 1, 11, x;...    %e10 --> ends at 13  
20210720, 11, 1, 12, x;...    %e11 --> ends at 13    
20210720, 12, 1, 11, 3;...    %e12 --> ends at 14 <-- misses middle fold at 25C hrs later
20210720, 13, 1, 10, x;...    %e13 --> ends at 12            
20210720, 14, 1, 11, x;...    %e14 --> ends at 13
20210720, 15, 1, 13.2, x;...  %e15 -> ends at 15.a
20210720, 16, 1, 12, x;...    %e16 --> ends at 14
...  % 202108021250 -- 3.5 hrs at 37C --> I have the timetrace of the temperature for this one!! one line per 2 seconds
202108021250, 1, 1, 13.2, 5;     % no folds until HS ends, but then folds somewhat
202108021250, 2, 1, 16.1, 0;
202108021250, 3, 1, 16.1, 0;
202108021250, 4, 1, 16.2, 0;
202108021250, 5, 1, 15.1, 3;     % e5
202108021250, 6, 1, 13.2, 5;     % forms all folds after HS
202108021250, 7, 1, 13.1, 5;     % forms posterior fold after HS
202108021250, 8, 1, 13.2, 5; 
202108021250, 9, 1, 15.1, 4; 
202108021250, 10, 1, 14, 4;       % e10
202108021250, 11, 1, 10, x;
202108021250, 12, 1, 13, 5; 
202108021250, 13, 1, 15.1, 3;    % e13
202108021250, 14, 1, 13.2, 5;
202108021250, 15, 1, 14, 2;       % two folds incomplete, really all incomplete
202108021250, 16, 1, 16.1, 0;
202108021250, 17, 1, 13, 3;
202108021250, 18, 1, 14.2, 2;
202108021250, 19, 1, 15.1, 2;
202108021250, 20, 1, 15.1, 4;
202108021250, 21, 1, 15.1, 4; 
202108021250, 22, 1, 15.1, 3;
202108021250, 23, 1, 14, 5;
202108021250, 24, 1, 13, 5;
202108021250, 25, 1, 13, 5;
202108021250, 26, 1, 15.1, 0;
202108021250, 27, 1, 13, 5;
202108021250, 28, 1, 13, 5;
202108021250, 29, 1, 15.1, 0;
202108021250, 30, 1, 13, 4;
202108021250, 31, 1, 13, 5;
202108021250, 32, 1, 14, 5;
202108021250, 33, 1, 15.1, 5;
202108021250, 34, 1, 15.1, 3; % completely missing anterior fold 
202108021250, 35, 1, 15.2, 3;
202108021250, 36, 1, 14.2, 5;
202108021250, 37, 1, 14.1, 5;
202108021250, 38, 1, 14, 3; 
202108021250, 39, 1, 14, 3;
202108021250, 40, 1, 13, 5;
202108021250, 41, 1, 13, 5;
202108021250, 42, 1, 15.1, 0; % "Odd" folds
202108021250, 43, 1, 15.2, 0; 
202108021250, 44, 1, 13, 5;
202108021250, 45, 1, 15.1, 0;
202108021250, 46, 1, 15.1, 5;  % strikingly folds post HS
202108021250, 47, 1, 13, 5;
    ]; ...
    

% 0 --> WT
% 1 --> disrupted but folds
% 2 --> incomplete Fold(s)
% 3 --> missing one Fold
% 4 --> missing two folds
% 5 --> missing all folds
% 6 --> no gut closure
% 7 --> dead before folds
% x --> did not finish/ cannot interpret

inBin = find(serca_stage_result_3p5hrs37C(:, stage) > 13 & ...
    serca_stage_result_3p5hrs37C(:, stage) < 16) ;
mutantIdx = intersect(inBin, find(serca_stage_result_3p5hrs37C(:, genotype) == 1)) ;
controlIdx = intersect(inBin, find(serca_stage_result_3p5hrs37C(:, genotype) == 0)) ;

% Clean out those that are x or 7
keepM = find(serca_stage_result_3p5hrs37C(mutantIdx, score) < scoreUpperLimit) ;
keepC = find(serca_stage_result_3p5hrs37C(controlIdx, score) < scoreUpperLimit) ;
mutantIdx = mutantIdx(keepM) ;
controlIdx = controlIdx(keepC) ;

% N mutant and N control
Nm = length(mutantIdx) ;
Nc = length(controlIdx) ;

badMutant = intersect(mutantIdx, find(serca_stage_result_3p5hrs37C(:, score) > scoreLowerLimit)) ;
goodMutant = intersect(mutantIdx, find(serca_stage_result_3p5hrs37C(:, score) <= scoreLowerLimit)) ;
badControl = intersect(controlIdx, find(serca_stage_result_3p5hrs37C(:, score) > scoreLowerLimit)) ;
goodControl = intersect(controlIdx, find(serca_stage_result_3p5hrs37C(:, score) <= scoreLowerLimit)) ;
fracBadFolds_m = length(badMutant) / Nm 
fracBadFolds_c = length(badControl) / Nc 


% Significance

