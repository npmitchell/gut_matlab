function buildFateHist(wcell) 
%
% NPMitchell 2020

stages = getStages(wcell) ;
fates = getFates(wcell) ;
wnorm = stages(fates == 0) ;
wdisr = stages(fates == 1) ;
wdead = stages(fates == 2) ;
% indet = stages(fates == 3) ;

kk = 1 ;
stages2do = 1:16 ;
pcnorm = zeros(1, length(stages2do)) ;
pcdisr = zeros(1, length(stages2do)) ;
pcdead = zeros(1, length(stages2do)) ;
nstage = zeros(1, length(stages2do)) ;
for stage = stages2do
    nnorm = sum(wnorm == stage) ;
    ndisr = sum(wdisr == stage) ;
    ndead = sum(wdead == stage) ;
    nthis = double(nnorm + ndisr + ndead) ;
    % fraction in each
    if nthis > 0
        pcnorm(kk) = double(nnorm) / nthis ;
        pcdisr(kk) = double(ndisr) / nthis ;
        pcdead(kk) = double(ndead) / nthis ;
    end
    nstage(kk) = nthis ;
    kk = kk + 1 ;
end
data = [pcnorm; pcdisr; pcdead] ;
labels = {'WT', 'disrupted', 'dies'} ;
bar(data', 'stacked')
% Add numbers
for kk = 1:length(stages2do)
    stage = stages2do(kk) ;
    if nstage(kk) > 0 
        text(stage, 0.5, num2str(nstage(kk)), 'vert','bottom', ...
            'horiz', 'center', 'color', 'w')
    end
end
legend(labels, 'location', 'northeastoutside')
xlabel('stage of heatshock')
ylabel('Probability')
