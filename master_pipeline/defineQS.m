opts.meshDir = meshDir ;
opts.flipy = flipy ;
opts.timeinterval = timeinterval ;
opts.timeunits = timeunits ;
opts.nV = 100 ;
opts.nU = 100 ;
opts.normalShift = 10 ;
opts.a_fixed = 2.0 ;
% opts.adjustlow = 1.00 ;                  %  floor for intensity adjustment
% opts.adjusthigh = 99.9 ;                 % ceil for intensity adjustment (clip)
opts.adjustlow = 0 ;                  %  floor for intensity adjustment
opts.adjusthigh = 0 ;                 % ceil for intensity adjustment (clip)
opts.phiMethod = 'curves3d' ;
disp('defining QS')
QS = QuapSlap(xp, opts) ;
disp('done')
