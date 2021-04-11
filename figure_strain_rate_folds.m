%% Average strain rate measurements from multiple QS objects
timeUnits = 'min';
timeInterval = 1;
flipy = true ;

opts0 = struct() ;
opts0.nU = 100 ; 
opts0.nV = 100 ;
opts0.timeUnits = timeUnits ;
opts0.spaceUnits = '$\mu$m' ;
opts0.normalShift = 10 ;
opts0.a_fixed = 2.0 ;
opts0.adjustlow = 1.00 ;                  %  floor for intensity adjustment
opts0.adjusthigh = 99.9 ;                 % ceil for intensity adjustment (clip)
opts0.phiMethod = 'curves3d' ;
opts0.lambda_mesh = 0 ;
opts0.lambda = 0 ;
opts0.lambda_err = 0 ;


meshDirs = {['/mnt/crunch/48Ygal4UASCAAXmCherry/', ...
    '201902072000_excellent/Time6views_60sec_1p4um_25x_obis1p5_2/', ...
    'data/deconvolved_16bit/msls_output'], ...
    };
QSs = {} ;
for qq = 1:length(meshDirs)
    opts = opts0 ;
    opts.meshDir = meshDirs{qq} ;
    opts.flipy = flipy ;
    opts.timeInterval = timeInterval ;
    disp('defining QS')
    QSs{qq} = QuapSlap(xp, opts) ;
    disp('done')
end


opts = struct() ;
opts.outputDir = '/mnt/data/analysis/2021/strainRates/' ;
if ~exist(opts.outputDir, 'dir')
    mkdir(opts.outputDir) ;
end
QSC = QuapSlapCollection(QSs, opts) ;

QSC.collectStrainRates()
QSC.collectStrainRateFolds()
QSC.collectStrainRateLobes()



