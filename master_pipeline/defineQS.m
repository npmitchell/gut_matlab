opts.meshDir = meshDir ;
opts.flipy = flipy ;
opts.timeinterval = timeinterval ;
opts.timeunits = timeunits ;
opts.nV = 100 ;
opts.nU = 100 ;
opts.normalShift = -10 ;
opts.a_fixed = 2.0 ;
disp('defining QS')
QS = QuapSlap(xp, opts) ;
disp('done')
