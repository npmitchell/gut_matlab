% Simple (without data listed) analysis


mCNTRL = 4 ;
nCNTRL = 11 ;
dCNTRL = 1 ;

% TRiP 1
m1 = 1 ;
n1 = 17 ;
d1 = 0 ;

% TRiP 2
m2 = 1 ;
n2 = 11 ;
d2 = 1 ;

% TRiP 4
m4 = 21 ;
n4 = 46 ;
d4 = 3 ;

colors = define_colors() ;

clf

% CONTROL
errorbar(1, mCNTRL / (mCNTRL+nCNTRL), sqrt(mCNTRL) / (mCNTRL+nCNTRL), 'Color', colors(1, :)) ;
hold on;
plot(1, mCNTRL/(mCNTRL+nCNTRL), 'o', 'Color', colors(1, :))


% TRiP 1 with 48Y GAL4
errorbar(2, m1 / (m1+n1), sqrt(m1) / (m1+n1), 'Color', colors(2, :)) ;
hold on;
plot(2, m1/(m1+n1), 'o', 'Color', colors(2, :))


% TRiP 2 with Mef2 GAL4
errorbar(3, m2 / (m2+n2), sqrt(m2) / (m2+n2), 'Color', colors(3, :)) ;
hold on;
plot(3, m2/(m2+n2), 'o', 'Color', colors(3, :))


% TRiP 4 with AntpGAL4
errorbar(4, m4 / (m4+n4), sqrt(m4) / (m4+n4), 'Color', colors(4, :)) ;
hold on;
plot(4, m4/(m4+n4), 'o', 'Color', colors(4, :))

xlim([0, 5])
ylim([0, 1])

xticks([])
ylabel('Probabilities for RNAi')
saveas(gcf, '/Users/npmitchell/Desktop/RNAi_week1.png')



