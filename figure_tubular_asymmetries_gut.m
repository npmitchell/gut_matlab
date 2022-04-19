
tmp = load('/mnt/data/analysis/tubular/gut_asymmetries/kinematics_caax.mat')
tmp2 = load('/mnt/data/analysis/tubular/gut_asymmetries/kinematics_histRFP_0p01_0p00_0p01.mat')
d1 = tmp.dorsH2vn ; d2 = tmp2.dorsH2vn ;
v1 = tmp.ventH2vn ; v2 = tmp2.ventH2vn ;
l1 = tmp.leftH2vn ; l2 = tmp2.leftH2vn ;
r1 = tmp.rightH2vn ; r2= tmp2.rightH2vn ;

nU = 100 ;
clf
subplot(1, 2, 1)
plot(linspace(0,1,100), d1-v1) ; hold on;
plot(linspace(0,1,100), d2-v2) ; hold on;
xlabel('position, s/L')
ylabel('DV asymmetry')

subplot(1, 2, 2)

plot(linspace(0,1,100), l1-r1) ; hold on;
plot(linspace(0,1,100), l2-r2) ; hold on;
xlabel('position, s/L')
ylabel('LR asymmetry')