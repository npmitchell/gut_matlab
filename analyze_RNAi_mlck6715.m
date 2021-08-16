% analyze_RNAi_results

%202107211653: 67;15 x MLCK RNAi #4 (TRiP 20) 37 for 45min
% 0 --> WT
% 1 --> disrupted but folds
% 2 --> incomplete Fold
% 3 --> missing one Fold
% 4 --> missing two folds
% 5 --> missing all folds
% 6 --> no gut closure
% 7 --> dead
stage_result_mlck = [14, 2;
    14, 3;
    15.2, 1;
    11, 6;
    13, 3;
    12, 5;
    ... % EXCLUDE SINCE DID NOT GET FAR ALONG ENOUGH 12, 5; % <-- 7. haltClosure
    15.3, 0;
    15.2, 0;
    15.1, 4;
    14, 5;  % e11
    15.2, 0;
    14, 4;
    14, 6;
    13, 5;
    13, 7;
    15.1, 4;
    10, 7;
    14, 7;
    11, 7;
    13, 6;
    15.1, 4;
    11, 6;
    13, 4; % <--- 24. how many missing? 2 folds missing (posterior anterior)
    15.1, 2;
    15.1, 7;
    15.2, 1;
    15.1, 1;
    15.2, 0;
    15.1, 2; % <-- 30. check, slowly folds
    13, 6;
    ... % EXCLUDE NOT FAR ENOUGH ALONG 12, 5; % <--- check, rename 30-> 32
    12, 5; % <-- check 33
    15.1, 0;
    12, 7;
    12, 7;
    13, 7;
    15.1, 5;
    15.2, 4; %<-- check seems like one fold  39 (already had one fold!)
    15.2, 3; 
    15.1, 7; % 41
    15.2, 0;
    12, 6;
    15.1, 1;
    14, 4; % 45
    15.1, 2; 
    13, 5;
    15.3, 0;
    ];

stage_result_6715WT = [...
    
] ;
    
    
    
    
    
    
    