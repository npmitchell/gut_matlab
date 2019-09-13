% Compare surface area and volume from different channels

clear 
close all

% Select where figure will go
outdir = '/mnt/data/analysis/' ;
markers = {'caax', 'hrfp', 'la'} ;
labels = {'Membrane', 'Nuclei', 'Actin'}; 

% Before running, run compute_mesh_surfacearea_volume.m
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
tf1_actin = {19, };
tfa_actin = {55, };
tfp_actin = {54, };
tLRb_actin = {67, };
tf1_nuclei = {7, 52, 77-65 } ; % artificially offset
tfa_nuclei = {40, 82, 106-65} ; % artificially offset
tfp_nuclei = {40, 90, 104-65} ; % artificially offset
tLRb_nuclei = {57, 101, 136-65} ;

% Prepare paths to data
rootdir = '/mnt/crunch/' ;
% membrane_excellent
caax_root = [rootdir '48Ygal4UASCAAXmCherry'] ;
caax_paths = {'201902072000_excellent/Time6views_60sec_1.4um_25x_obis1.5_2/data/deconvolved_16bit/msls_output_prnun5_prs1_nu0p00_s0p10_pn2_ps4_l1_l1/', ...
    '201903211930_great/Time6views_60sec_1p4um_25x_1p0mW_exp0p150/data/deconvolved_16bit/msls_output_prnun5_prs1_nu0p00_s0p10_pn2_ps4_l1_l1/', ...
    } ;
% nuclei_folded2part
hrfp_root = [rootdir '48Ygal4-UAShistRFP/'] ;
hrfp_paths = {'201901021550_folded_2part/Time12views_60sec_1.2um_25x_4/data/deconvolved_16bit/msls_output_prnu0_prs0_nu0p10_s1p00_pn4_ps4_l1_l1/', ...
    '201904031830_great/Time4views_60sec_1p4um_25x_1p0mW_exp0p35_2/data/deconvolved_16bit/msls_output_prnun5_prs1_nu0p00_s0p10_pn2_ps4_l1_l1/', ...
    '201903312000_closure_folding_errorduringtwist/Time4views_60sec_1p4um_25x_1p0mW_exp0p35_2_folding/data/deconvolved_16bit/msls_output_prnun5_prs1_nu0p00_s0p10_pn2_ps4_l1_l1/',...
    %'201903312000_closure_folding_errorduringtwist/Time4views_180sec_1p4um_25x_1p0mW_exp0p35_dorsalclosure/data/deconvolved_16bit/msls_output_prnun5_prs1_nu0p00_s0p10_pn2_ps4_l1_l1/', ...
    } ;
% actin
la_root = [rootdir '48YGal4UasLifeActRuby'] ;
la_paths = {'201904021800_great/Time6views_60sec_1p4um_25x_1p0mW_exp0p150_3/data/deconvolved_16bit/msls_output_prnun5_prs1_nu0p00_s0p10_pn2_ps4_l1_l1/', ...
    };
roots = {caax_root, hrfp_root, la_root} ;
paths = {{caax_paths}, {hrfp_paths}, {la_paths}};
npaths = [length(caax_paths), length(hrfp_paths), length(la_paths)] ;

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
figh = figure();
hold on;
fig2 = figure();
% secondax = copyobj(gca, gcf);
hold on;
color1 = [0, .4470, .7410] ;
color2 = [.8500, .3250, .0980] ;
color3 = [0, 0, 0] ;
offy = 0.03 ;
offy2 = .5 ;
levely = 1.04 ;
texty = 0.01 ;
texty2 = 0.1 ;
mark_origin = true ;

% Iterate over each marker
for mi = 1:length(markers)
    % Obtain the label for this marker
    label = labels{mi} ;
    these_paths = paths{mi} ;
    these_paths = these_paths{1} ;
    
    % Cycle through all datasets of this marker
    for j=1:length(these_paths)
        mpath = these_paths{j} ;
        disp(['path: ' mpath])
        matdir = fullfile(roots{mi}, mpath) ;
        disp(['seeking data in: ' matdir]) 
        fn = fullfile(fullfile(roots{mi}, mpath), 'surfacearea_volume_stab.mat') ;
        if exist(fn, 'file')
            % Load the surface area and volume from disk
            load(fn)
            
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
            times = 1:dt:dt*length(aas) ;
            times = times - t0 ;
            
            % grab time of first/mid fold (t=0)
            [~, ind] = min(abs(times)) ;
            
            % aas is the area array (over time)
            % ass is the normed area array (over time)
            ass = aas / aas(ind) ;
            vss = vvs / vvs(ind) ;
            
            % Filter the data
            windowSize = 7; 
            sampling = 1:length(times) ;
            b = (1/windowSize)*ones(1,windowSize);
            a = 1;
            asmooth = filter(b, a, ass(sampling)) ;
            vsmooth = filter(b, a, vss(sampling)) ;
            asmooth2 = smoothdata(ass, 'rlowess', 5) ;
            vsmooth2 = smoothdata(vss, 'rlowess', 11);
            
            da = gradient(asmooth) ;
            dv = gradient(vsmooth) ;
            
            % [envHigh, envLow] = envelope(ass,10,'peak');
            % aMean = (envHigh+envLow)/2;
            % [envHigh, envLow] = envelope(vss,10,'peak');
            % vMean = (envHigh+envLow)/2;
            % da = gradient(aMean) ;
            % dv = gradient(vMean) ;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Plot data first
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            figure(figh)
            if mark_origin
                % Plot data
                ah = plot(times, asmooth2, 'Color', color1, 'LineStyle', linestyle) ;
                vh = plot(times, vsmooth2, 'Color', color2, 'LineStyle', linestyle); 
                % ah = plot(times(sampling), asmooth, 'Color', color1) ;
                % vh = plot(times(sampling), vsmooth, 'Color', color2); 
                % ah = plot(times, aMean, 'Color', color1) ;
                % vh = plot(times, vMean, 'Color', color2); 

                % Plot time of first/mid fold (t=0)
                p1 = [times(ind), ass(ind) + offy] ;
                p2 = [times(ind), ass(ind)] ;
                dp = p2 - p1 ;
                quiver(p1(1), p1(2), dp(1), dp(2), 'k-')
                text(p1(1), p1(2) + texty, 'fold', 'FontSize' , 20, ...
                    'HorizontalAlignment', 'Center', 'Interpreter', 'Latex')
            else
                % Plot data
                plot(times, asmooth2, 'Color', color1, 'LineStyle', linestyle) ;
                plot(times, vsmooth2, 'Color', color2, 'LineStyle', linestyle); 
            end

            % Get indices of anterior and posterior folds
            [~, ia] = min(abs(times + t0 - ta)) ;
            [~, ip] = min(abs(times + t0 - tp)) ;

            % grab time of anterior fold
            p1 = [times(ia), ass(ia)] ;
            plot(p1(1), p1(2), 'o', 'MarkerSize', 10, 'Color', 'k')

            % grab time of posterior fold
            p1 = [times(ip), ass(ip)] ;
            plot(p1(1), p1(2), 's', 'MarkerSize', 10, 'Color', 'k')
            p2 = [times(ip), ass(ip) - 1.2 * offy ] ;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Plot the derivatives on other figure
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            figure(fig2);
            tt = times(sampling) ;
            tt = tt(windowSize + 1:end) ;
            da = da(windowSize + 1:end) * 100;
            dv = dv(windowSize + 1:end) * 100;
                
            if mark_origin
                % Plot data
                a2 = plot(tt, da, 'Color', color1, 'Linestyle', linestyle) ;
                v2 = plot(tt, dv, 'Color', color2, 'Linestyle', linestyle); 

                % Plot time of first/mid fold (t=0)
                p1 = [0, offy2] ;
                p2 = [0, 0] ;
                dp = p2 - p1 ;
                quiver(p1(1), p1(2), dp(1), dp(2), 'k-')
                text(p1(1), p1(2) + texty2, 'fold', ...
                    'HorizontalAlignment', 'Center', 'Interpreter', 'Latex')
                mark_origin = false ;
            else
                % Plot derivatives
                plot(tt, da, 'Color', color1, 'Linestyle', linestyle) ;
                plot(tt, dv, 'Color', color2, 'Linestyle', linestyle); 
            end
            
            % Get indices of anterior and posterior folds
            [~, ia] = min(abs(tt + t0 - ta)) ;
            [~, ip] = min(abs(tt + t0 - tp)) ;

            % grab time of anterior fold
            p1 = [tt(ia), da(ia)] ;
            af = plot(p1(1), p1(2), 'o', 'MarkerSize', 10, 'Color', 'k')

            % grab time of posterior fold
            p1 = [tt(ip), da(ip)] ;
            pf = plot(p1(1), p1(2), 's', 'MarkerSize', 10, 'Color', 'k')
                
        else
            disp(['Could not find ' fn])
        end
    end
end

%% Data Figure
figure(figh)
% Label and save figure
title('Surface Area and Volume of the Midgut', 'FontSize' , 20)
xlabel('Time [min]', 'FontSize' , 20)
ylabel('Normalized Area or Volume', 'FontSize' , 20)

% axes for the second plot (secondaxes) and the two helping Lines H1 and H2
hold on 
savlegend = legend(gca, [ah, vh], {'Surface Area', 'Volume'}, 'Location','northwest', 'FontSize' , 20) ;
% set(secondax, 'Color', 'none', 'XTick', [], 'YTick', [], 'Box', 'Off') 
a=axes('position',get(gca,'position'),'visible','off');
delete( get(a, 'Children'))
hold on
kx = [0, 0] ;
ky = [1, 1] ;
H1 = plot(kx, ky, '-', 'Color', [0 0 0]);
H2 = plot(kx, ky, '--', 'Color', [0 0 0]);
H3 = plot(kx, ky, ':', 'Color', [0 0 0]);
hold off
legend (a, [H1 H2 H3 af pf], {'Membrane', 'Nuclei', 'Actin', 'Anterior Fold', 'Posterior Fold'}, 'Location', 'west', 'FontSize' , 20) ;

% Save the figure
saveas(figh, fullfile(outdir, 'area_volume_stab_comparison.pdf'))
saveas(figh, fullfile(outdir, 'area_volume_stab_comparison.png'))

%% Derivatives Figure
figure(fig2)
% Label and save figure
title('Surface Area and Volume Rate of Change', 'FontSize' , 20)
xlabel('Time [min]', 'FontSize' , 20)
ylabel('Percent change, $\partial A / \partial t$, $\partial V / \partial t$ [min$^{-1}$]', ...
    'Interpreter', 'Latex', 'FontSize' , 20)

% axes for the second plot (secondaxes) and the two helping Lines H1 and H2
hold on 
savlegend = legend(gca, [a2, v2], {'Surface Area', 'Volume'}, 'Location','northwest', 'FontSize' , 20) ;
% set(secondax, 'Color', 'none', 'XTick', [], 'YTick', [], 'Box', 'Off') 
a=axes('position',get(gca,'position'),'visible','off');
delete( get(a, 'Children'))
hold on
kx = [0, 0] ;
ky = [1, 1] ;
H1 = plot(kx, ky, '-', 'Color', [0 0 0]);
H2 = plot(kx, ky, '--', 'Color', [0 0 0]);
H3 = plot(kx, ky, ':', 'Color', [0 0 0]);
hold off
legend (a, [H1 H2 H3 af pf], {'Membrane', 'Nuclei', 'Actin', 'Anterior Fold', 'Posterior Fold'}, 'Location', 'west', 'FontSize' , 20) ;

% Save the figure
saveas(fig2, fullfile(outdir, 'area_volume_stab_comparison_derivatives.pdf'))
saveas(fig2, fullfile(outdir, 'area_volume_stab_comparison_derivatives.png'))

% Reset groot
set(groot, 'defaultAxesColorOrder', originalColorOrder, ...
    'defaultAxesLineStyleOrder', originalStyleOrder)

