% Compare surface area and volume from different channels

clear 
close all

% Select where figure will go
outdir = '/mnt/data/analysis/2020/' ;
markers = {'caax', 'hrfp', 'la'} ;
labels = {'Membrane', 'Nuclei', 'Actin'}; 
fontsize = 10 ;
markersize = 10 ;
graycolor = 0.8 * [1,1,1] ;
opacity = 0.4 ;
lw_mean = 2 ;
foldText = {'middle', 'fold'} ;
possibleTimes = -100:200 ;  % all possible timestamps in timeUnits
smoothStyle = 'rloess' ;
smoothSpan = 0.2 ;
smoothDegree = 2 ;
cutOffDerivativeEdges = true ;
% For 2 panel figure
axpos = [0.1300    0.5838    0.5439    0.3412] ;
axpos_with_writhe = axpos - [0, 0, 0.07, 0] ;
axpos_derivs = [0.1300    0.1100    axpos(3)   axpos(4)]; 
axpos_derivs_with_writhe = axpos_derivs - [0, 0, 0.07, 0] ; 

% For single panel figures
% axpos = [0.0966    0.1100    0.5756    0.8150 ] ;
% axpos_with_writhe = [0.0966    0.1100    0.48    0.8150 ] ;

if ~exist(outdir, 'dir')
    mkdir(outdir)
end

%% Before running, compute mesh surfacearea and volume for each dataset
% see compute_mesh_surfacearea_volume.m or QuapSlap()

%% Store folding onset for midfold in tf1, anterior fold in tfa, posterior in tfp
% Notes about folding times and LR asymmetry 
% Folds 1,2,3: first TP with self-contact: midfold, antfold, postfold
% First LR symmetry breaking timestep = when compartment 2 moves laterally
% CAAX: 
% - 201902072000_excellent: TP 151, 170, 175, [LR 206] (tps begin 110)
% HRFP:
% - 201901021550_folded_2part: TP 0, 54, 54 [LR 78]
% LifeAct:
% - 201904021800_great: TP 19, 55, 54 [LR 74] (tps begin 1)
tf1_membrane = {151-109, 36 };
tfa_membrane = {178-109, 61};
tfp_membrane = {181-109, 69};
tLRb_membrane = {206-109, 98};
tf1_actin = {19, 0};
tfa_actin = {55, 34};
tfp_actin = {54, 20};
tLRb_actin = {67, };
tf1_nuclei = {7, 52, 77-65 } ; % artificially offset
tfa_nuclei = {40, 82, 106-65} ; % artificially offset
tfp_nuclei = {40, 90, 104-65} ; % artificially offset
tLRb_nuclei = {57, 101, 136-65} ;

% Prepare paths to data
crunch = '/mnt/crunch/' ;
data = '/mnt/data/' ;
% membrane_excellent
caax_root = '48Ygal4UASCAAXmCherry/' ;
caax_paths = {[crunch caax_root '201902072000_excellent/' ...
    'Time6views_60sec_1.4um_25x_obis1.5_2/data/deconvolved_16bit_smallbox/'...
    'msls_output_prnun5_prs1_nu0p00_s0p10_pn2_ps4_l1_l1/'], ...
    [crunch caax_root '201903211930_great/'...
    'Time6views_60sec_1p4um_25x_1p0mW_exp0p150/data/deconvolved_16bit/'...
    'msls_output_prnun5_prs1_nu0p00_s0p10_pn2_ps4_l1_l1/']} ;
caax_nUs = {100, 100, } ;
caax_nVs = {100, 100, } ;
caax_shifts = {10, 10, } ;

% nuclei_folded2part
hist_root = '48Ygal4-UAShistRFP/' ;
hist_paths = {[crunch hist_root '201901021550_folded_2part/'...
    'Time12views_60sec_1.2um_25x_4/data/deconvolved_16bit/'...
    'msls_output_prnu0_prs0_nu0p10_s1p00_pn4_ps4_l1_l1/'], ...
    [crunch hist_root '201904031830_great/Time4views_60sec_1p4um_25x_1p0mW_exp0p35_2/'...
    'data/deconvolved_16bit/msls_output/'], ...
    [crunch hist_root '201903312000_closure_folding_errorduringtwist/'...
    'Time4views_60sec_1p4um_25x_1p0mW_exp0p35_2_folding/data/deconvolved_16bit/'...
    'msls_output_prnun5_prs1_nu0p00_s0p10_pn2_ps4_l1_l1/'],...
    } ;
hist_nUs = {100, 100, 100, } ;
hist_nVs = {100, 100, 100, } ;
hist_shifts = {10, 10, 10, } ;

% actin
la_root = '48YGal4UasLifeActRuby/' ;
la_paths = {[crunch la_root '201904021800_great/Time6views_60sec_1p4um_25x_1p0mW_exp0p150_3/'...
    'data/deconvolved_16bit/msls_output/'], ...
    [data la_root '201907311600_48YGal4UasLifeActRuby_60s_exp0p150_1p0mW_25x_1p4um/'...
    'Time4views_60sec_1p4um_25x_1p0mW_exp0p15/data/deconvolved_16bit/msls_output/']};
la_nUs = {100, 100,} ;
la_nVs = {100, 100,} ;
la_shifts = {10, 10,} ;

paths = {caax_paths, hist_paths, la_paths};
npaths = [length(caax_paths), length(hist_paths), length(la_paths)] ;
nUs = {caax_nUs, hist_nUs, la_nUs} ;
nVs = {caax_nVs, hist_nVs, la_nVs} ;
shifts = {caax_shifts, hist_shifts, la_shifts} ;

% % Build labels
% clear labels
% for ii=1:length(caax_paths)
%     labels{ii} = 'membrane' ;
% end
% for jj=1:length(hrfp_paths)
%     labels{ii + jj} =  'nuclei' ;
% end
% for kk=1:length(la_paths)
%     labels{ii + jj + kk} = 'actin' ;
% end

% Initialize the figure
originalColorOrder = get(groot, 'defaultAxesColorOrder');
originalStyleOrder = get(groot, 'defaultAxesLineStyleOrder');
set(groot,'defaultAxesColorOrder',[0, .4470, .7410; .8500, .3250, .0980],...
      'defaultAxesLineStyleOrder','-|--|:')
linestyle_list = {'-', '--', ':'} ;
% secondax = copyobj(gca, gcf);
hold on;
% color1 = [0, .4470, .7410] ;  % blue
% color2 = [.8500, .3250, .0980] ; % red
% color2 = [0.540,  0.2500, 0.0900] ;  % brown
color1 = [0.4660, 0.6740, 0.1880] ;  % green
color2 = [0.4940, 0.1840, 0.5560] ;  % purple
color3 = [0.5400, 0.2500, 0.0900] ;  % brown
% color4 = [0.8500, 0.3250, 0.0980] ;  % red
color4 = [0.3010, 0.7450, 0.9330] ;  % sky
offys = [0.1, 0.1, 0.5, 0.5] ;
offy2s = [.2, 0.5, 0.7, 0.7] ;
levely = 1.04 ;
textys = 0.5 * offys ;  % [0.01, 1, 0.3, 1] ;
texty2s = 0.7 * offy2s ;

% Build figure in stages: volume, area, length, writhe
for build = 1:4
    close all
    mark_origin = true ;
    ax1 = subplot(2, 1,  1) ;
    if build < 3
        ylim([0.7, 1.5])
    else
        ylim([0.7, 4])
    end
    hold on;
    ax2 = subplot(2, 1, 2);
    if build < 3
        ylim([-0.8, 1.4])
    else
        ylim([-0.8, 4.5])
    end
    hold on;
    
    % dataset index is dmyk. Use this to collate results into a master array
    dmyk = 1;

    % Iterate over each marker
    for mi = 1:length(markers)
        % Obtain the label for this marker
        label = labels{mi} ;
        these_paths = paths{mi} ;

        % Cycle through all datasets of this marker
        for j=1:length(these_paths)
            matdir = these_paths{j} ;
            disp(['seeking data in: ' matdir]) 
            fn = fullfile(matdir, 'surfacearea_volume_stab.mat') ;


            if exist(fn, 'file')
                % Load the surface area and volume from disk
                load(fn, 'aas', 'vvs', 'dt')
                areas = aas ;
                volumes = vvs ;

                % get time offset
                if strcmp(label, 'Membrane')
                    disp('loading membrane tps')
                    t0 = tf1_membrane{j} ;
                    ta = tfa_membrane{j} ;
                    tp = tfp_membrane{j} ;
                    linestyle = linestyle_list{1} ;
                elseif strcmp(label, 'Nuclei')
                    disp('loading nuclei tps')
                    t0 = tf1_nuclei{j} ;
                    ta = tfa_nuclei{j} ;
                    tp = tfp_nuclei{j} ; 
                    linestyle = linestyle_list{2} ;
                elseif strcmp(label, 'Actin')
                    disp('loading actin tps')
                    t0 = tf1_actin{j} ;
                    ta = tfa_actin{j} ;
                    tp = tfp_actin{j} ;            
                    % find which linestyle to use
                    linestyle = linestyle_list{3} ;
                end

                % Plot the data for surface area and volume
                times = 1:dt:dt*length(areas) ;
                times = times - t0 ;

                % Check that no timestamps are out of bounds
                try
                    assert(~any(times<min(possibleTimes)))
                    assert(~any(times>max(possibleTimes)))
                catch
                    error('Please expand the range of the possibleTimes variable')
                end

                % grab time of first/mid fold (t=0)
                [~, t0ind] = min(abs(times)) ;

                % ass is the normed area array (over time)
                ass = areas / areas(t0ind) ;
                vss = volumes / volumes(t0ind) ;

                % Filter the data for derivatives
                windowSize = 7; 
                % sampling = 1:length(times) ;
                % b = (1/windowSize)*ones(1,windowSize);
                % a = 1;
                % asmooth = filter(b, a, ass(sampling)) ;
                % vsmooth = filter(b, a, vss(sampling)) ;
                asmooth = smooth(ass, windowSize) ;
                vsmooth = smooth(vss, windowSize) ;

                % Smooth raw data 
                % asmooth2 = smooth(ass, smoothSpan, smoothStyle, smoothDegree) ;
                % vsmooth2 = smooth(vss, smoothSpan, smoothStyle, smoothDegree);
                asmooth2 = smoothdata(ass, 'rlowess', 5) ;
                vsmooth2 = smoothdata(vss, 'rlowess', 11);

                da = gradient(asmooth) ;
                dv = gradient(vsmooth) ;

                % [envHigh, envLow] = envelope(ass,markersize,'peak');
                % aMean = (envHigh+envLow)/2;
                % [envHigh, envLow] = envelope(vss,markersize,'peak');
                % vMean = (envHigh+envLow)/2;
                % da = gradient(aMean) ;
                % dv = gradient(vMean) ;

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Plot data first
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % figure(fig1)
                axes(ax1)
                % Plot data
                plot(times / 60, vsmooth2, 'Color', [color1 opacity], 'LineStyle', linestyle);
                hold on;
                if build > 1
                    plot(times / 60, asmooth2, 'Color', [color2 opacity], 'LineStyle', linestyle) ;
                end
                % if mark_origin
                %     % % Plot time of first/mid fold (t=0)
                %     % p1 = [times(t0ind), 1 + offys(build) ] ;
                %     % p2 = [times(t0ind), 1] ;
                %     % dp = p2 - p1 ;
                %     % quiver(p1(1), p1(2), dp(1), dp(2), 'k-')
                %     % text(p1(1), p1(2) + textys(build), foldText, 'FontSize' , fontsize, ...
                %     %     'HorizontalAlignment', 'Center', 'Interpreter', 'Latex')
                % end

                % Get indices of anterior and posterior folds
                [~, ia] = min(abs(times + t0 - ta)) ;
                [~, ip] = min(abs(times + t0 - tp)) ;

                % grab time of anterior fold
                if build > 1
                    a1 = [times(ia), ass(ia)] ;
                    plot(a1(1) / 60, a1(2), 'o', 'MarkerSize', markersize, 'Color', graycolor)

                    % grab time of posterior fold
                    p1 = [times(ip), ass(ip)] ;
                    plot(p1(1) / 60, p1(2), 's', 'MarkerSize', markersize, 'Color', graycolor)
                    % p2 = [times(ip), ass(ip)] ;            
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% Collate to total results for averaging
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Find indices of possibleTimes which are represented in times
                [sharedvals,idx] = intersect(possibleTimes, times, 'stable');
                assert(all(sharedvals == times))

                allV(idx, dmyk) = volumes ;
                allVNormed(idx, dmyk) = vss ;

                if build > 1
                    allA(idx, dmyk) = areas ;
                    allANormed(idx, dmyk) = ass ;
                    
                    if ~exist('anteriorFold_time_area', 'var')
                        anteriorFold_time_area = a1 ;
                        posteriorFold_time_area = p1 ;
                    else
                        anteriorFold_time_area(dmyk, :) = a1 ;
                        posteriorFold_time_area(dmyk, :) = p1 ;
                    end
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% Plot the derivatives on other figure
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % figure(fig2);
                axes(ax2) ;
                if cutOffDerivativeEdges
                    tt = times ;
                    tt = tt(windowSize + 1:end-windowSize) ;
                    da = da(windowSize + 1:end-windowSize) ;
                    dv = dv(windowSize + 1:end-windowSize) ;
                else
                    tt = times ;
                    da = da * 100;
                    dv = dv * 100;
                end

                % Plot derivatives
                plot(tt / 60, dv * 60, 'Color', [color1 opacity], 'Linestyle', linestyle);
                hold on;
                if build > 1
                    plot(tt / 60, da * 60, 'Color', [color2 opacity], 'Linestyle', linestyle) ;
                end
                % if mark_origin
                %     % Plot time of first/mid fold (t=0)
                %     p1 = [0, offy2s(build)] ;
                %     p2 = [0, 0] ;
                %     dp = p2 - p1 ;
                %     quiver(p1(1), p1(2), dp(1), dp(2), 'k-')
                %     text(p1(1), p1(2) + texty2s(build), foldText, ...
                %         'HorizontalAlignment', 'Center', 'Interpreter', 'Latex')
                %     mark_origin = false ;
                % end
                
                % Fold timestamps
                if build > 1
                    % Get indices of anterior and posterior folds
                    [~, ia] = min(abs(tt + t0 - ta)) ;
                    [~, ip] = min(abs(tt + t0 - tp)) ;

                    % grab time of anterior fold
                    a1 = [tt(ia), da(ia)] ;
                    af = plot(a1(1) / 60, a1(2) * 60, 'o', 'MarkerSize', markersize, 'Color', graycolor) ;

                    % grab time of posterior fold
                    p1 = [tt(ip), da(ip)] ;
                    pf = plot(p1(1) / 60, p1(2) * 60, 's', 'MarkerSize', markersize, 'Color', graycolor) ;
                end

                %% Now update the dataset index since this one is done
                do_update = true ;
            else
                disp(['Could not find ' fn])
                do_update = false; 
            end

            %% LENGTH & WRITHE
            if build > 2
                % Also load writhe dir for length and writhe
                uvexten = sprintf('nU%04d_nV%04d', nUs{mi}{j}, nVs{mi}{j}) ;
                uvCoordDir = fullfile(matdir, ['gridCoords_' uvexten]) ;
                shiftstr = sprintf('_%03dstep', shifts{mi}{j}) ;
                writheDir = fullfile(uvCoordDir, ['centerline_from_DVhoops' shiftstr], 'writhe') ;
                writhefn = fullfile(writheDir, ['writhe_sphi_' uvexten '_avgpts.mat']) ;
                % length and writhe
                disp(['Seeking writhe & length in ' writhefn])
                if exist(writhefn, 'file')
                    disp('Found writhe file. Loading...')
                    tmp = load(writhefn) ;
                    writhe = tmp.Wr.Levitt ;
                    aplength = tmp.Length_t.lengths ;
                    assert(length(times) == length(aplength))
                    
                    apL = aplength / aplength(t0ind) ;
                    dL = gradient(smoothdata(aplength, 'rlowess', 11), dt) ;
                    dLn = gradient(smoothdata(apL, 'rlowess', 11), dt) ;
                    dWr = gradient(smoothdata(writhe, 'rlowess', 11), dt);
                        
                    %% Plot them
                    % figure(fig1);
                    axes(ax1)
                    plot(times / 60, apL, 'Color', [color3 opacity], ...
                        'Linestyle', linestyle) 
                    % Plot derivatives on figure 2
                    % figure(fig2);
                    axes(ax2) ;
                    plot(times / 60, dLn * 60, 'Color', [color3 opacity], ...
                        'Linestyle', linestyle) 
                    
                    % Plot writhe over time
                    if build > 3
                        % figure(fig1)
                        axes(ax1)
                        yyaxis right
                        plot(times / 60, writhe, 'Color', [color4 opacity], ...
                            'Linestyle', linestyle) 
                        yyaxis left
                        % figure(fig2)
                        axes(ax2)
                        yyaxis right
                        plot(times / 60, dWr * 60, 'Color', [color4 opacity], ...
                            'Linestyle', linestyle) 
                        yyaxis left
                    end
                
                    %% Add to collated list
                    [sharedvals,idx] = intersect(possibleTimes, times, 'stable');
                    assert(all(sharedvals == times))
                    allL(idx, dmyk) = aplength ;
                    allLNormed(idx, dmyk) = apL ;
                    all_dL(idx, dmyk) = dL ;
                    all_dLn(idx, dmyk) = dLn ;
                    allWr(idx, dmyk) = writhe ;
                    all_dWr(idx, dmyk) = dWr ;
                    
                else
                    disp('Writhe file NOT found. Seeking alternate from centerline...')
                end
            end
            
            if do_update 
                dmyk = dmyk + 1;
            end
        end
    end


    %% Data Figure
    % figure(fig1)
    axes(ax1)
    % Label and save figure
    % title('geometric dynamics of the midgut', 'FontSize', fontsize, 'interpreter', 'latex')
    % xlabel('time [hr]', 'FontSize', fontsize, 'interpreter', 'latex')
    if build == 1
        ylabel('$V/V_0$', 'interpreter', 'latex', 'FontSize', fontsize)
    elseif build == 2
        ylabel('$A/A_0$, $V/V_0$', 'interpreter', 'latex', 'FontSize', fontsize)
    elseif build > 2
        ylabel('$L/L_0$, $A/A_0$, and $V/V_0$', ...
            'interpreter', 'latex', 'FontSize', fontsize)
    end
    
    %% Average curves and plot mean surface area and volume
    % mask empty values in master arrays
    emptyID = find(allV == 0);
    goodRow = find(any(allV, 2));
    allVm = allV;
    allVnm = allVNormed;
    allVm(emptyID) = NaN ;
    allVnm(emptyID) = NaN ;
    if build > 1
        allAm = allA;
        allAnm = allANormed;
        allAm(emptyID) = NaN ;
        allAnm(emptyID) = NaN ;
        
        % area of time --> at, normalized area over time --> ant
        at = nanmean(allAm(goodRow, :), 2) ;
        ant = nanmean(allAnm(goodRow, :), 2) ;
    end
    % vomume over time --> vt, normalized volume over time --> vnt
    vt = nanmean(allVm(goodRow, :), 2) ;
    vnt = nanmean(allVnm(goodRow, :), 2) ;
    
    %% Average curves for length and writhe 
    if build > 2
        % mask empty values in master arrays
        emptyID = find(allL == 0);
        goodRow = find(any(allL, 2));
        allLm = allL ;
        allLnm = allLNormed;
        allWrm = allWr ;
        all_dLm = all_dL ;
        all_dLnm = all_dLn ;
        all_dWrm = all_dWr ;
        % masked Length and normed Length
        allLm(emptyID) = NaN ;
        allLnm(emptyID) = NaN ;
        all_dLm(emptyID) = NaN ;
        all_dLnm(emptyID) = NaN ;
        % masked Writhe
        allWrm(emptyID) = NaN ;
        all_dWrm(emptyID) = NaN ;
        % area of time --> at, normalized area over time --> ant
        lent = nanmean(allLm(goodRow, :), 2) ;
        lnt = nanmean(allLnm(goodRow, :), 2) ;
        dlnt = nanmean(all_dLnm(goodRow, :), 2) ;
        wt = nanmean(allWrm(goodRow, :), 2) ;
        dwt = nanmean(all_dWrm(goodRow, :), 2) ;
    end

    %% Add to figure
    % figure(fig1)
    axes(ax1)
    hold on;
    % masked time --> 'timem'
    timem = possibleTimes(goodRow) ;
    
    % Plot volume
    volumemean = plot(timem / 60, vnt, '-', ...
        'color', color1, 'linewidth', lw_mean) ;
    
    % Plot area
    if build > 1
        areamean = plot(timem / 60, ant, '-', ...
            'color', color2, 'linewidth', lw_mean) ;
        % Plot mean fold times
        afoldmean = mean(anteriorFold_time_area(:, 1)) ;
        pfoldmean = mean(posteriorFold_time_area(:, 1)) ;
        [~, aID] = min(abs(timem - afoldmean)) ;
        [~, pID] = min(abs(timem - pfoldmean)) ;

        plot(afoldmean / 60, ant(aID), 'o', 'MarkerSize', markersize, ...
            'Color', color2, 'HandleVisibility', 'off') ;
        plot(pfoldmean / 60, ant(pID), 's', 'MarkerSize', markersize, ...
            'Color', color2, 'HandleVisibility', 'off') ;
    end
    
    % Plot length
    if build > 2
        lengthmean = plot(timem / 60, lnt, '-', ...
            'color', color3, 'linewidth', lw_mean) ;
    end
    
    % Plot writhe
    if build > 3
        yyaxis right
        writhemean = plot(timem / 60, wt, '-', ...
            'color', color4, 'linewidth', lw_mean) ;
        ax = gca;
        ax.YAxis(2).Color = color4 ;        
        ylabel('Writhe, $Wr$', 'interpreter', 'latex')
        yyaxis left
        
    end
    
    % Find t0idx from possible times
    t0idx_possible = find(timem == 0) ;
    
    % axes for the second plot (secondaxes) and the two helping Lines H1 and H2
    hold on 
    if build == 1
        legend([volumemean ], ...
            {'volume'}, ...
            'Location', 'northeastoutside', 'FontSize' , fontsize, 'interpreter', 'latex') ;
    elseif build == 2
        legend([areamean volumemean ], ...
            {'surface area', 'volume'}, ...
            'Location', 'northeastoutside', 'FontSize' , fontsize, 'interpreter', 'latex') ;
    elseif build > 2
        legend([lengthmean areamean volumemean ], ...
            {'length', 'surface area', 'volume'}, ...
            'Location', 'northeastoutside', 'FontSize' , fontsize, 'interpreter', 'latex') ;
    end
    % mf = plot(0, 1, '^', 'MarkerSize', markersize, 'Color', 'k', 'HandleVisibility', 'Off') ;
    ylims = get(ax1, 'ylim') ;
    set(ax1, 'ylim', ylims) ;
    plot([0,0], ylims, 'k--', 'HandleVisibility', 'off') ;
    
    if build < 4
        set(gca, 'position', axpos)
    else
        set(gca, 'position', axpos_with_writhe)
    end
    
    %% Make second legend for axis 1 (main axis)
    % % set(secondax, 'Color', 'none', 'XTick', [], 'YTick', [], 'Box', 'Off')
    % a=axes('position', axpos,'visible','off');
    % delete( get(a, 'Children'))
    % hold on
    % kx = [0, 0] ;
    % ky = [1, 1] ;
    % H1 = plot(kx, ky, '-', 'Color', [0 0 0]);
    % H2 = plot(kx, ky, '--', 'Color', [0 0 0]);
    % H3 = plot(kx, ky, ':', 'Color', [0 0 0]);
    % 
    % hold off
    % if build == 1
    %     legend(a, [H1 H2 H3], ...
    %         {'membrane labeled', 'nuclei labeled', 'actin labeled'}, ...
    %         'Location', 'southeastoutside', 'FontSize' , fontsize, ...
    %         'interpreter', 'latex') ;
    %     set(a, 'position', axpos)
    % else
    %     legend(a, [H1 H2 H3 af pf], ...
    %         {'membrane labeled', 'nuclei labeled', 'actin labeled', ...
    %         'anterior fold', 'posterior fold'}, ...
    %         'Location', 'southeastoutside',  'FontSize' , fontsize, ...
    %         'interpreter', 'latex') ;
    %     set(a, 'position', axpos)
    % end
    
    % Save the figure
    % saveas(gcf, fullfile(outdir, ['area_volume_stab_comparison_build' num2str(build) '.pdf']))
    % saveas(gcf, fullfile(outdir, ['area_volume_stab_comparison_build' num2str(build) '.png']))


    %% Derivatives Figure
    % figure(fig2)
    axes(ax2)
    % Label and save figure
    % title('surface area and volume rate of change', ...
    %     'interpreter', 'latex', 'FontSize', fontsize)
    xlabel('time [hr]', 'FontSize', fontsize, 'interpreter', 'latex')
    if build == 1
        ylabel('rate of change, $\partial_t\tilde{V} $ [hr$^{-1}$]', ...
            'Interpreter', 'Latex', 'FontSize', fontsize)
    elseif build == 2
        ylabel('rate of change, $\partial_t\tilde{A}$, $\partial_t\tilde{V} $ [hr$^{-1}$]', ...
            'Interpreter', 'Latex', 'FontSize', fontsize)
    else 
        ylabel({'rate of change, ', ...
            '$\partial_t\tilde{L}$, $\partial_t\tilde{A}$, $\partial_t\tilde{V} $ [hr$^{-1}$]'}, ...
            'Interpreter', 'Latex', 'FontSize', fontsize)
    end
        
    % Plot mean derivatives
    % Filter the data
    % asmooth2 = smooth(ant, smoothSpan, smoothStyle, smoothDegree) ;
    % vsmooth2 = smooth(vnt, smoothSpan, smoothStyle, smoothDegree) ;
    % This is noisy
    % asmooth2 = smoothdata(ant, 'rlowess', 5) ;
    % vsmooth2 = smoothdata(vnt, 'rlowess', 11);
    
    % This is a better filter for us
    vsmooth2 = smooth(vnt, windowSize) ;
    dv = gradient(vsmooth2) ;
    
    if build > 1
        asmooth2 = smooth(ant, windowSize) ;
        da = gradient(asmooth2) ;
    end
    
    % cut off windowSize points in beginning

    if cutOffDerivativeEdges 
        timed = timem(windowSize+1:end-windowSize) ;
        da = da(windowSize+1:end-windowSize) ;
        dv = dv(windowSize+1:end-windowSize) ;
    else
        timed = timem ;
    end

    hold on;
    % Plot volume
    v2 = plot(timed / 60, dv * 60, '-', 'color', color1, 'linewidth', lw_mean) ;
    
    % Plot area
    if build > 1
        a2 = plot(timed / 60, da * 60, '-', 'color', color2, 'linewidth', lw_mean) ;
        % Plot mean fold times
        [~, aID] = min(abs(timed - afoldmean)) ;
        [~, pID] = min(abs(timed - pfoldmean)) ;
        plot(afoldmean / 60, da(aID) * 60, 'o', 'MarkerSize', markersize, ...
            'Color', color2, 'HandleVisibility', 'off') ;
        plot(pfoldmean / 60, da(pID) * 60, 's', 'MarkerSize', markersize, ...
            'Color', color2, 'HandleVisibility', 'off') ;
    end
    
    % Plot length
    if build > 2
        l2 = plot(timem / 60, dlnt * 60, '-', 'color', color3, 'linewidth', lw_mean) ;
    end
    
    % Plot writhe
    if build > 3
        yyaxis right
        w2 = plot(timem / 60, dwt * 60, '-', 'color', color4, 'linewidth', lw_mean) ;
        ax2.YAxis(2).Color = color4 ;
        ylabel('Writhe change, $\partial_t Wr$ [hr$^{-1}$]', 'interpreter', 'latex')
        yyaxis left
    end
    
    %% LEGEND FOR AXIS 2
    % % axes for the second plot (secondaxes) and the two helping Lines H1 and H2
    % hold on 
    % if build == 1
    %     savlegend = legend(gca, [v2], ...
    %         {'volume rate, $\partial_t V$'}, ...
    %         'Location','northeastoutside', 'FontSize' , fontsize, ...
    %         'interpreter', 'latex') ;
    % elseif build == 2
    %     savlegend = legend(gca, [a2, v2], ...
    %         {'area rate, $\partial_t A$', 'volume rate, $\partial_t V$'}, ...
    %         'Location','northeastoutside', 'FontSize' , fontsize, ...
    %         'interpreter', 'latex') ;
    % elseif build > 2
    %     savlegend = legend(gca, [l2, a2, v2], ...
    %         {'length rate, $\partial_t L$', ...
    %         'area rate, $\partial_t A$', 'volume rate, $\partial_t V$'}, ...
    %         'Location','northeastoutside', 'FontSize' , fontsize, ...
    %         'interpreter', 'latex') ;
    % end

    %% Mark the first fold
    % mf = plot(0, 0, '^', 'MarkerSize', markersize, 'Color', 'k', 'HandleVisibility', 'Off') ;
    ylims = get(ax2, 'ylim') ;
    plot([0 0], ylims, 'k--', 'handlevisibility', 'off')
    set(ax2, 'ylim', ylims) ;
    
    % Set xlimits to match the first panel
    xlims = get(ax1, 'xlim') ;
    set(ax2, 'xlim', xlims) ;
    
    if build < 4
        set(ax2, 'position', axpos_derivs)
    else
        set(ax2, 'position', axpos_derivs_with_writhe)
    end
    
    %% SECOND LEGEND FOR DERIVATE SUBPLOT
    % % set(secondax, 'Color', 'none', 'XTick', [], 'YTick', [], 'Box', 'Off') 
    % a=axes('position', axpos_derivs, 'visible','off');
    % delete( get(a, 'Children'))
    % hold on
    kx = [0, 0] ;
    ky = [1, 1] ;
    H1 = plot(kx, ky, '-', 'Color', [0 0 0]);
    H2 = plot(kx, ky, '--', 'Color', [0 0 0]);
    H3 = plot(kx, ky, ':', 'Color', [0 0 0]);
    % hold off
    if build == 1
        legend(ax2, [H1 H2 H3], ...
            {'membrane labeled', 'nuclei labeled', 'actin labeled'}, ...
            'Location', 'southeastoutside', 'FontSize' , fontsize, ...
            'interpreter', 'latex') ;
        set(ax2, 'position', axpos_derivs) ;
    elseif build < 4
        legend(ax2, [H1 H2 H3 af pf], ...
            {'membrane labeled', 'nuclei labeled', 'actin labeled', ...
            'anterior fold', 'posterior fold'}, ...
            'Location', 'southeastoutside', 'FontSize' , fontsize, ...
            'interpreter', 'latex') ;
        set(ax2, 'position', axpos_derivs) ;
    else
        legend(ax2, [H1 H2 H3 af pf], ...
            {'membrane labeled', 'nuclei labeled', 'actin labeled', ...
            'anterior fold', 'posterior fold'}, ...
            'Location', 'southeastoutside', 'FontSize' , fontsize, ...
            'interpreter', 'latex') ;
        set(ax2, 'position', axpos_derivs_with_writhe) ;
    end

    %% Save the figure
    saveas(gcf, fullfile(outdir, ['area_volume_stab_comparison_derivatives_build' num2str(build) '.pdf']))
    saveas(gcf, fullfile(outdir, ['area_volume_stab_comparison_derivatives_build' num2str(build) '.png']))

    % Reset groot
    set(groot, 'defaultAxesColorOrder', originalColorOrder, ...
        'defaultAxesLineStyleOrder', originalStyleOrder)
end

