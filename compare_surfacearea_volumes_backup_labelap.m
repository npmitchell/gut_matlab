% Compare surface area and volume from different channels

% Select where figure will go
outdir = '/mnt/data/analysis/' ;
markers = {'caax', 'hrfp', 'la'} ;
labels = {'membrane', 'nuclei', 'actin'}; 

% Notes about folding times and LR asymmetry 
% Folds 1,2,3: first TP with self-contact: midfold, antfold, postfold
% First LR symmetry breaking timestep = when compartment 2 moves laterally
% CAAX: 
% - 201902072000_excellent: TP 151, 170, 175, [LR 206] (tps begin 110)
% HRFP:
% - 201901021550_folded_2part: TP 0, 54, 54 [LR 78]
% LifeAct:
% - 201904021800_great: TP 19, 55, 54 [LR 74] (tps begin 1)
tf1_membrane = {151-109, };
tfa_membrane = {178-109, };
tfp_membrane = {181-109, };
tLRb_membrane = {206-109, };
tf1_actin = {19, };
tfa_actin = {55, };
tfp_actin = {54, };
tLRb_actin = {67, };
tf1_nuclei = {20, } ; % artificially offset
tfa_nuclei = {54, } ; % artificially offset
tfp_nuclei = {54, } ; % artificially offset
tLRb_nuclei = {52, } ;

% Prepare paths to data
rootdir = '/mnt/crunch/' ;
% membrane
caax_root = [rootdir '48Ygal4UASCAAXmCherry'] ;
caax_paths = {'201902072000_excellent/Time6views_60sec_1.4um_25x_obis1.5_2/data/deconvolved_16bit/msls_output_prnun5_prs1_nu0p00_s0p10_pn2_ps4_l1_l1/', ...
    '201903211930_great/Time6views_60sec_1p4um_25x_1p0mW_exp0p150/data/', ...
    } ;
% nuclei
hrfp_root = [rootdir '48Ygal4-UAShistRFP/'] ;
hrfp_paths = {'201901021550_folded_2part/Time12views_60sec_1.2um_25x_4/data/deconvolved_16bit/msls_output_prnu0_prs0_nu0p10_s1p00_pn4_ps4_l1_l1/', ...
    '201904031830_great/Time4views_60sec_1p4um_25x_1p0mW_exp0p35_2/data/deconvolved_16bit_4view_thres0p140_pix2_sig15_affinereg_individual_registration/', ...
    '201903312000_closure_folding_errorduringtwist/Time4views_180sec_1p4um_25x_1p0mW_exp0p35_dorsalclosure/data/', ...
    } ;
% actin
la_root = [rootdir '48YGal4UasLifeActRuby'] ;
la_paths = {'201904021800_great/Time6views_60sec_1p4um_25x_1p0mW_exp0p150_3/data/deconvolved_16bit/msls_output_prnun5_prs1_nu0p00_s0p10_pn2_ps4_l1_l1/', ...
    };
roots = {caax_root, hrfp_root, la_root} ;
paths = {caax_paths, hrfp_paths, la_paths};

% Initialize the figure
originalColorOrder = get(groot, 'defaultAxesColorOrder');
originalStyleOrder = get(groot, 'defaultAxesLineStyleOrder');
set(groot,'defaultAxesColorOrder',[0, .4470, .7410; .8500, .3250, .0980],...
      'defaultAxesLineStyleOrder','-|--|:')
figh = figure();
% secondax = copyobj(gca, gcf);
hold on;
color1 = [0, .4470, .7410] ;
color2 = [.8500, .3250, .0980] ;
color3 = [0, 0, 0] ;
offy = 0.03 ;
levely = 1.04 ;
texty = 0.01 ;
mark_origin = true ;

% Iterate over each marker
for mi = 1:length(markers)
    label = labels{mi} ;
    for j=1:length(paths{mi})
        mpaths = paths{mi} ;
        matdir = fullfile(roots{mi}, mpaths{j}) ;
        disp(['seeking data in: ' matdir]) 
        fn = fullfile(fullfile(roots{mi}, mpaths{j}), 'surfacearea_volume.mat') ;
        if exist(fn, 'file')
            % Load the surface area and volume from disk
            load(fn)
            
            % get time offset
            if strcmp(label, 'membrane')
                t0 = tf1_membrane{j} ;
                ta = tfa_membrane{j} ;
                tp = tfp_membrane{j} ;
            elseif strcmp(label, 'nuclei')
                t0 = tf1_nuclei{j} ;
                ta = tfa_nuclei{j} ;
                tp = tfp_nuclei{j} ;
            elseif strcmp(label, 'actin')
                t0 = tf1_actin{j} ;
                ta = tfa_actin{j} ;
                tp = tfp_actin{j} ;
            end
            
            % Plot the data
            times = 1:dt:dt*length(aas) ;
            times = times - t0 ;
            
            % grab time of first/mid fold (t=0)
            [~, ind] = min(abs(times)) ;
            
            ass = aas / aas(ind) ;
            vss = vvs / vvs(ind) ;
            
            if mark_origin
                % Plot data
                ah = plot(times, ass, 'Color', color1) ;
                vh = plot(times, vss, 'Color', color2); 

                % Plot time of first/mid fold (t=0)
                p1 = [times(ind), ass(ind) + offy] ;
                p2 = [times(ind), ass(ind)] ;
                dp = p2 - p1 ;
                quiver(p1(1), p1(2), dp(1), dp(2), 'k-')
                text(p1(1), p1(2) + texty, 'fold', ...
                    'HorizontalAlignment', 'Center', 'Interpreter', 'Latex')
                % plot(times(ind), ass(ind), 'o', 'Color', color3) ;
                % plot(times(ind), vss(ind), 'o', 'Color', color3) ;
                mark_origin = false ;
            else
                % Plot data
                plot(times, ass, 'Color', color1) ;
                plot(times, vss, 'Color', color2); 
            end

            % check if t_antfold = t_postfold
            [~, ia] = min(abs(times + t0 - ta)) ;
            [~, ip] = min(abs(times + t0 - tp)) ;
            if false %ia == ip
                % grab time of ant/post fold (simultaneous)
                p1 = [times(ia), levely] ;
                p2 = [times(ia), ass(ia)] ;
                dp = p2 - p1 ;
                quiver(p1(1), p1(2), dp(1), dp(2), 'k-')
                text(p1(1), p1(2) - texty, 'a,p', ...
                    'HorizontalAlignment', 'Center', 'Interpreter', 'Latex')
            else
                % grab time of anterior fold
                p1 = [times(ia), ass(ia)] ;
                p2 = [times(ia), ass(ia) + 1.2 * offy] ;
                dp = p2 - p1 ;
                quiver(p1(1), p1(2), dp(1), dp(2), 'k-')
                text(p2(1), p2(2) + texty, 'a', ...
                    'HorizontalAlignment', 'Center', 'Interpreter', 'Latex')

                % grab time of posterior fold
                p1 = [times(ip), ass(ip)] ;
                p2 = [times(ip), ass(ip) - 1.2 * offy] ;
                dp = p2 - p1 ;
                quiver(p1(1), p1(2), dp(1), dp(2), 'k-')
                text(p2(1), p2(2) - texty, 'p', ...
                    'HorizontalAlignment', 'Center', 'Interpreter', 'Latex')
            end
        else
            disp(['Could not find ' fn])
        end
    end
end

% Label and save figure
legend(gca, [ah, vh], {'area', 'volume'}, 'Location','northwest')
title('Surface area and volume')
xlabel('Time [min]')
ylabel('Normalized area or volume')

% axes for the second plot (secondaxes) and the two helping Lines H1 and H2
hold on 
%delete( get(secondax, 'Children'))
kx = [0, 0] ;
ky = [1, 1] ;
H1 = plot(kx, ky, '-', 'Color', [0 0 0]);
H2 = plot(kx, ky, '--', 'Color', [0 0 0]);
H3 = plot(kx, ky, ':', 'Color', [0 0 0]);
% set(secondax, 'Color', 'none', 'XTick', [], 'YTick', [], 'Box', 'Off') 
a=axes('position',get(gca,'position'),'visible','off');
legend (a, [H1 H2 H3], {'membrane', 'nuclei', 'actin'}) ;

% Save the figure
saveas(figh, fullfile(outdir, 'area_volume_comparison.pdf'))
saveas(figh, fullfile(outdir, 'area_volume_comparison.png'))

% Reset groot
set(groot, 'defaultAxesColorOrder', originalColorOrder, ...
    'defaultAxesLineStyleOrder', originalStyleOrder)


