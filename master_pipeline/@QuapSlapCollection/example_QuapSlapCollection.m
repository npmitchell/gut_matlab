%% Example usage for QuapSlapCollection()

dirs = {'/mnt/data/48Ygal4UASCAAXmCherry/201902072000_excellent/Time6views_60sec_1p4um_25x_obis1p5_2/data/deconvolved_16bit/', ...
    '/mnt/data/48Ygal4-UAShistRFP/201904031830_great/Time4views_60sec_1p4um_25x_1p0mW_exp0p35_2/data/deconvolved_16bit/', ...
    '/mnt/data/handGAL4klarHandGFPhistGFP/202105072030_1mWGFP/deconvolved_16bit/'} ;

%% Create collection
QS = {} ;
for dd = 1:length(dirs)
    xpdir = dirs{dd} ;
    cd(xpdir) 
    load(fullfile(xpdir, 'xp.mat'), 'xp', 'opts')
    disp('defining QS class instance')
    QS{dd} = QuapSlap(xp, opts) ;    
end
opts = struct() ;
opts.outputDir = '/mnt/data/analysis/QuapSlapCollections/202206/' ;
cd(opts.outputDir)
QSC = QuapSlapCollection(QS, opts) ;

%% align the datasets











