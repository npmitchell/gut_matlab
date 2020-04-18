%% Align meshes in time
%
% Isaac Breinyn 2019
%
% Extract information for temporal alignment and scaling of meshes
% Collapses curves using shape parameter s = A / V^(2/3).
%
% To run before: compute_mesh_surfacearea_volume.m
% To run after: plot_aligned_meshes.m
%
% Note: All scaling and alignment is done relative to CAAX Excellent

clear
close all
clc
addpath('/mnt/data/code/gut_matlab/curve_functions/')

%% Prepare paths, Load Data, Calculate Scaling Factors
outdir = '/mnt/data/analysis/SA_volume/' ;
markers = {'caax', 'hrfp', 'la'} ;
labels = {'Membrane', 'Nuclei', 'Actin'};

% Prepare paths to data
rootdir = '/mnt/crunch/' ;

% membrane
caax_root = [rootdir '48Ygal4UASCAAXmCherry'] ;
caax_paths = {'201902072000_excellent/Time6views_60sec_1.4um_25x_obis1.5_2/data/deconvolved_16bit/msls_output_prnun5_prs1_nu0p00_s0p10_pn2_ps4_l1_l1/', ...
    '201903211930_great/Time6views_60sec_1p4um_25x_1p0mW_exp0p150/data/deconvolved_16bit/msls_output_prnun5_prs1_nu0p00_s0p10_pn2_ps4_l1_l1/'} ;

% nuclei
hrfp_root = [rootdir '48Ygal4-UAShistRFP/'] ;
hrfp_paths = {'201901021550_folded_2part/Time12views_60sec_1.2um_25x_4/data/deconvolved_16bit/msls_output_prnu0_prs0_nu0p10_s1p00_pn4_ps4_l1_l1/', ...
    '201904031830_great/Time4views_60sec_1p4um_25x_1p0mW_exp0p35_2/data/deconvolved_16bit/msls_output_prnun5_prs1_nu0p00_s0p10_pn2_ps4_l1_l1/', ...
    '201903312000_closure_folding_errorduringtwist/Time4views_60sec_1p4um_25x_1p0mW_exp0p35_2_folding/data/deconvolved_16bit/msls_output_prnun5_prs1_nu0p00_s0p10_pn2_ps4_l1_l1/',...
    %'201903312000_closure_folding_errorduringtwist/Time4views_180sec_1p4um_25x_1p0mW_exp0p35_dorsalclosure/data/deconvolved_16bit/msls_output_prnun5_prs1_nu0p00_s0p10_pn2_ps4_l1_l1/',...
    } ;

% actin
la_root = [rootdir '48YGal4UasLifeActRuby'] ;
la_paths = {'201904021800_great/Time6views_60sec_1p4um_25x_1p0mW_exp0p150_3/data/deconvolved_16bit/msls_output_prnun5_prs1_nu0p00_s0p10_pn2_ps4_l1_l1/'};

roots = {caax_root, hrfp_root, la_root} ;
paths = {{caax_paths}, {hrfp_paths}, {la_paths}};
npaths = [length(caax_paths), length(hrfp_paths), length(la_paths)] ;

% Initialize the figure
originalColorOrder = get(groot, 'defaultAxesColorOrder');
originalStyleOrder = get(groot, 'defaultAxesLineStyleOrder');
set(groot,'defaultAxesColorOrder',[0, .4470, .7410; .8500, .3250, .0980],...
    'defaultAxesLineStyleOrder','-|--|:')
linestyle_list = {'-', '--', ':'} ;
figh = figure();
hold on;
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


aaa = {} ; % average area cell array
aat = {} ; % average time cell array
sfa = {} ; % scale factor cell array
del = {} ; % delta cell array
Pi = {} ; % dedimensionalized factor cell array
Pavg = {} ; % average dedimensionalized factor cell array

P1 = {} ;
toffmin = {} ;

% Iterate over each marker
for mi = 1:length(markers)
    % Obtain the label for this marker
    label = labels{mi} ;
    these_paths = paths{mi} ;
    these_paths = these_paths{1} ;
    
    % Cycle through all datasets of this marker
    for j=1:length(these_paths)
        mpath = these_paths{j} ;
        %disp(['path: ' mpath])
        matdir = fullfile(roots{mi}, mpath) ;
        %disp(['seeking data in: ' matdir])
        fn = fullfile(fullfile(roots{mi}, mpath), 'surfacearea_volume_stab.mat') ;
        if exist(fn, 'file')
            load(fn) % Load the surface area and volume from disk
            linestyle = linestyle_list{mi} ; % find which linestyle to use
            
            % grab average value of area over time (aa)
            
            aa = .5*(min(aas)+max(aas)); % Ignores asymptotic behavior
            aaa{mi, j} = aa ;
            
            % grab time of aa for this dataset
            
            [~, taa] = min(abs(aa - aas));
            aat{mi, j} = taa ;
            
            % grab time for alignment using dedimensionalized value
            if mi == 1 && j == 1
                Pi = (aas)./(vvs.^(2/3)) ; % dedimensionalized value
                refmdv = .5*(max(Pi)+min(Pi)) ; % avg ddv for ref
                [~, tva] = min(abs(refmdv - Pi)) ;
                tref = tva ;
                toff = 0 ;
                aat{mi, j} = toff ; % save to cell array
            else
                Pi = (aas)./(vvs.^(2/3)) ; % dedimensionalized value
                mdv = .5*(max(Pi)+min(Pi)) ; % avg ddv
                [~, tva] = min(abs(refmdv - Pi)) ;
                toff = tva-tref ; % offset relative to CAAX (no offset)
                aat{mi, j} = toff ; % save to cell array
            end
            P{mi,j} = Pi ;
            
            % grab scaling factor for this dataset
            
            delta = (max(aas)-min(aas)) ; % diff between max and min SA
            del{mi, j} = delta ; % save to cell array
            scalfac = (8.703039e+04)/aa ; % CAAX avg area div by current avg area
            sfa{mi, j} = scalfac ; % save to cell array
            
            % Create SA and vol arrays
            
            ass = aas / aas(tva) ; % normed area array (over time)
            sass = scalfac*aas ; % scaled normed area array (over time)
            vss = vvs / vvs(tva) ; % normed volume array (over time)
            svss = (scalfac*vvs) / vvs(tva) ; % scaled normed volume array (over time)
            
            % Minimization
            if mi == 1 && j == 1
                P{mi,j} = P{mi,j}(30:110);
            else
                P{mi,j} = vertcat((P{mi,j}(2)+zeros(length(P{1,1})-length(P{mi,j}),1)),P{mi,j},(P{mi,j}(end-2)+zeros(1,1)));
            end
            time = zeros(length(P{mi,j}),1);
            for i = 1:length(time)
                time(i) = i;
            end
            P1{mi,j} = horzcat(time, P{mi,j});
            options = optimset('PlotFcns',@optimplotfval);
            toffmin{mi,j} = fminsearch(@(vars)matchCurves(vars, P1{mi,j}, P1{1,1}), [0, 1, 0, 1], options) ;
            
            % Specify what you want to plot!
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            pl = toffmin{mi,j}(4)*(P{mi,j} + toffmin{mi,j}(3)) ; 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Plot the data for surface area and volume
            
            time = toffmin{mi,j}(2)*(time + round(toffmin{mi,j}(1)));
            
            % Filter the data
            
            windowSize = 7;
            sampling = 1:length(time) ;
            b = (1/windowSize)*ones(1,windowSize);
            a = 1;
            asmooth = filter(b, a, pl(sampling)) ;
            asmooth2 = smoothdata(pl, 'rlowess', 5) ;
            
            figure(figh)
            if mark_origin
                ah = plot(time, asmooth2, 'Color', color1, 'LineStyle', linestyle) ; % Plot data
            else
                plot(time, asmooth2, 'Color', color1, 'LineStyle', linestyle) ; % Plot data
            end
        else
            disp(['Could not find ' fn])
        end
    end
end

%% Save area and volume info as .txt
ofn = fullfile(outdir, 'aat_sfa_toffmin.txt') ;
maat = [aat{1,:}, aat{2,:}, aat{3,:}];
msfa = [sfa{1,:}, sfa{2,:}, sfa{3,:}];
mtoffmin = [toffmin{1,:}, toffmin{2,:}, toffmin{3,:}];
dlmwrite(ofn,maat);
dlmwrite(ofn,msfa,'-append', 'delimiter',',','roffset',1);
dlmwrite(ofn,mtoffmin,'-append', 'delimiter',',','roffset',1);

%% Data Figure
figure(figh)
% Label and save figure
title('Alignment of Meshes', 'FontSize' , 20)
xlabel('Time [min]', 'FontSize' , 20)
ylabel('Your label here! $5/mo.', 'FontSize' , 20)

% axes for the second plot (secondaxes) and the two helping Lines H1 and H2
hold on
a=axes('position',get(gca,'position'),'visible','off');
delete( get(a, 'Children'))
hold on
kx = [0, 0] ;
ky = [1, 1] ;
H1 = plot(kx, ky, '-', 'Color', [0 0 0]);
H2 = plot(kx, ky, '--', 'Color', [0 0 0]);
H3 = plot(kx, ky, ':', 'Color', [0 0 0]);
hold off
legend (a, [H1 H2 H3], {'Membrane', 'Nuclei', 'Actin'}, 'Location', 'northeast', 'FontSize' , 10) ;

% Save the figure
%saveas(figh, fullfile(outdir, 'mesh_temporal_alignment.pdf'))
%saveas(figh, fullfile(outdir, 'mesh_temporal_alignment.png'))
