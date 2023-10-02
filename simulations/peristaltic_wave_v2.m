x = 1 ;
y = 0 ;

ind = 1 ;
rate = 0.01 ;
tps = 0:0.001:1 ;
ntps = length(tps) ;
xx = zeros(ntps, 1) ;
yy = zeros(ntps, 1) ;
for tt = tps

    x = x - rate * y ;
    y = y + rate * x ;

    xx(ind) = x ;
    yy(ind) = y ;
    ind = ind + 1 ;
end

plot(xx, yy, '.')
axis equal