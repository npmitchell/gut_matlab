%% Supplementary figure with pullback pathline texturepatch overlay
% showing three different timepoints, like 0, 35, and 70 minutes of midgut
% morphogenesis in 3 colors. This shows their overlay is nearly black and
% white (grayscale) in the pullback plane, signaling contrained tissue 
% motion in the parameterization.

%% Pullback pathline texturepatching (PIV pathline) ======================
disp('Create pullback using pullback pathline coords')
tidx0 = tubi.xp.tIdx(tubi.t0) ;
tp2do = tubi.xp.fileMeta.timePoints(tidx0:10:end) ;
tp2do = [tp2do, setdiff(tubi.xp.fileMeta.timePoints, tp2do)] ;
for tt = tp2do
    disp(['PB Pathline texturepatch: NOW PROCESSING TIME POINT ', num2str(tt)]);
    tidx = tubi.xp.tIdx(tt);
    
    % Load the data for the current time point ------------------------
    tubi.setTime(tt) ;
    
    % Establish custom Options for MIP --> choose which pullbacks to use
    pbOptions = struct() ;
    pbOptions.numLayers = [0 0] ; % how many onion layers over which to take MIP
    pbOptions.generate_spsm = false ;
    pbOptions.generate_sp = false ;
    pbOptions.overwrite = false ;
    pbOptions.generate_pivPathline = true ;
    tubi.generateCurrentPullbacks([], [], [], pbOptions) ;
end


%% Then punch the images into Fiji and merge channels cyan/magenta/yellow
% Used either 123 153 183
% or, for revision, I used 123, 158, 193
%