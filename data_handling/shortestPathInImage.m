function [tpath, tpath0] = shortestPathInImage(cij, options)
% Using cij as speed of travel, find shortest path in image
% Note this uses Gabriel Pyre's functions perform_fast_marching
%
% Parameters
% ----------
% 
% Returns
% -------
%
% NPMitchell 2021

%% Options
log_display = false ;
preview = true ;
Woffset = 0. ;
cij_exponent = 1. ; 

if isfield(options, 'logDisplay')
    log_display = options.logDisplay ;
end
if isfield(options, 'preview')
    preview = options.preview ;
end
if isfield(options, 'Woffset')
    Woffset = options.Woffset ;
end
if isfield(options, 'exponent')
    cij_exponent = options.exponent ;
end


%% find the path (rough path via maxima) 
tpath0 = zeros(size(cij, 1), 2) ;
tpath0(:, 1) = 1:size(cij, 1) ;
for tq = 1:size(cij, 1)
    [~, ind] = max(cij(tq, :)) ;        
    tpath0(tq, 2) = ind ;
end

% Use image heatmap as DT to find shortest path
% close all
% imagesc(cij);
% hold on;
% plot(tpath0(:, 2), tpath0(:, 1), 'o') 
ntpguess = min(size(cij)) ;
% xtmp = (1:ntpguess) + tpath0(1, 2) ;
% ytmp = (1:ntpguess) + tpath0(1, 1) ;
% plot(xtmp, ytmp, 'k--')
% xtmp = tpath0(end, 2) - (0:ntpguess-1) ;
% ytmp = tpath0(end, 1) - (0:ntpguess-1) ;
% plot(xtmp, ytmp, 'k--')
% axis equal
% axis tight
% title('Initial maxima')
% pause(2)
% move_on = false ;

% COULD ALLOW FOR MAXIMA TO DEFINE THE PATH
% msg = 'Does path look ok? Enter=yes, Backspace=no, n/Delete=No correspondence' ;
% title(msg)
% disp(msg)
% good_button = false ;
% abort = false ;
% while ~good_button
%     button = waitforbuttonpress() ;
%     if button && strcmp(get(gcf, 'CurrentKey'), 'return')
%         move_on = true ;
%         good_button = true ;
%     elseif button && strcmp(get(gcf, 'CurrentKey'), 'backspace')
%         move_on = false;
%         good_button = true ;
%         % Save this 
%         % title(['Initial correspondence between ' exptIDs{ii} ' and ' exptIDs{jj}])
%         % saveas(gcf, fullfile(corrImOutDir, ['correspondence' ijstr '_initial.png']))
%     elseif button && (strcmp(get(gcf, 'CurrentKey'), 'delete') || strcmp(get(gcf, 'CurrentKey'), 'n') )
%         abort = true ;
%         good_button = true ;
%         move_on = true ;
%     end
% end
tpath = tpath0 ;
move_on = false ;

% Shortest path 
while ~move_on
    close all
    %displays log plot if indicated
    if (log_display == 1)
        imagesc(log(cij))
    else
        imagesc(cij)
    end
    hold on;
    plot(tpath(:, 2), tpath(:, 1), 'ro') 
    plot(tpath(1, 2), tpath(1, 1), 'ks')
    plot(tpath(end, 2), tpath(end, 1), 'k^')
    xtmp = (1:ntpguess) + tpath(1, 2) ;
    ytmp = (1:ntpguess) + tpath(1, 1) ;
    plot(xtmp, ytmp, 'k--')
    xtmp = tpath(end, 2) - (0:ntpguess-1) ;
    ytmp = tpath(end, 1) - (0:ntpguess-1) ;
    plot(xtmp, ytmp, 'k--')
    msg = 'Do endpoints look ok? Enter=yes, Backspace=no/select' ;
    title(msg)
    disp(msg)
    startpt = tpath(1, :) ;
    % Guess start/endpt to be first/last tpath points
    endpt = tpath(end, :) ;
    button = waitforbuttonpress() ;
    if button && strcmp(get(gcf, 'CurrentKey'), 'return')
        move_on = true ;
        disp('set startpt/endpt:')
        disp(['startpt = ', num2str(startpt)])
        disp(['endpt = ', num2str(endpt)])
    elseif button && strcmp(get(gcf, 'CurrentKey'), 'backspace')
        move_on = false;
        msg = 'Start timept is ok? Enter=yes, Backspace=no' ;
        scatter(startpt(2), startpt(1), 100, 'r', 'filled')
        disp(msg)
        title(msg)
        button = waitforbuttonpress() ;
        if button && strcmp(get(gcf, 'CurrentKey'), 'return')
            disp('Great, how about end point?')
        elseif button && strcmp(get(gcf, 'CurrentKey'), 'backspace')
            % New guess for startpoint is maximum along column
            [~, startpt(1, 1)] = max(cij(:, 1)) ;
            scatter(startpt(2), startpt(1), 100, 'r', 'filled')
            msg = 'Ok, how about startpt now? Enter=yes, Backspace=no' ;
            disp(msg)
            title(msg)
            button = waitforbuttonpress() ;
            if button && strcmp(get(gcf, 'CurrentKey'), 'return')
                disp('Great, how about end point?')
            elseif button && strcmp(get(gcf, 'CurrentKey'), 'backspace')
                msg = 'Find the startpoint manually by clicking' ;
                disp(msg)
                title(msg)
                startpt = ginput(1) ;
                startpt = [startpt(2) startpt(1)] ;
            end
        end
        msg = 'End timept is ok? Enter=yes, Backspace=no' ;
        scatter(startpt(2), startpt(1), 100, 'g', 'filled')
        scatter(endpt(2), endpt(1), 100, 'r', 'filled')
        disp(msg)
        title(msg)
        button = waitforbuttonpress() ;
        if button && strcmp(get(gcf, 'CurrentKey'), 'return')
            disp('Great, all done.')
        elseif button && strcmp(get(gcf, 'CurrentKey'), 'backspace')
            % New guess for startpoint is maximum along column
            [~, endpt(1, 1)] = max(cij(:, end)) ;
            endpt(1, 2) = size(cij, 2) ;
            scatter(endpt(2), endpt(1), 100, 'r', 'filled')
            msg = 'Ok, how about endpt now? Enter=yes, Backspace=no' ;
            disp(msg)
            title(msg)
            button = waitforbuttonpress() ;
            if button && strcmp(get(gcf, 'CurrentKey'), 'return')
                disp('Great, all done.')
            elseif button && strcmp(get(gcf, 'CurrentKey'), 'backspace')
                msg = 'Find the endpoint manually by clicking';
                disp(msg)
                title(msg)
                endpt = ginput(1) ;
                endpt = [endpt(2) endpt(1)] ;
            end
        end
    else
        error('bad button press')
    end

    % get shortest path via fastmarching by Gabriel Peyre
    options.propagation_type = 'normal';
    options.Tmax = sum(size(cij))*1.2;
    options.start_points = startpt';
    clf;
    disp('Performing FM');
    options.reduc_factor = 1.0 ;
    options.weight = 1.0;
    %   'W' is the weight matrix (the highest, the slowest the front will move).
    %   'start_points' is a 2 x k array, start_points(:,i) is the ith starting point .
    %   'end_points' is a 2 x 1 array, it is the goal.
    % Define speed of pixels' movement for marching
    % Blur the image a bit to straighten lines, then find curve
    % WW = imgaussfilt((max(cij(:)) - cij) - min(cij(:)), 1) ;
    % WW = imgaussfilt(cij) ;
    WW = imgaussfilt((cij - min(cij(:))).^ cij_exponent + Woffset) ;
    % WW = WW - min(WW(:)) + 1e-2 ;
    % WW = cij - min(cij(:)) + 1e-2;
    % Compute distance transform
    % didn't get slow version to work
    % [DD,S] = perform_front_propagation_2d_slow(WW',...
    %     [startpt(2);startpt(1)], [endpt(2); endpt(1)], 4000, []);
    %
    [DD,S] = perform_fast_marching(WW', [startpt(2)+0.1; startpt(1)+0.1], options) ;

    % Check DD
    % imagesc(DD)
    % button = waitforbuttonpress() ;

    % perform_fmstar_2d(WW', spt, ept, options);    
    disp('Extracting Paths');
    stepsize = 0.1 ;
    thres_dist = 2 ;
    str_options = [stepsize 10000];

    % path extraction
    options.str_options = str_options ;
    options.trim_path = true ;
    options.startpt = [startpt(2), startpt(1)] ;
    options.thres_dist = thres_dist ;
    str_options = options.str_options ;

    % grad = compute_grad(DD);
    % grad = -perform_vf_normalization(grad);
    % Dx = squeeze(grad(:, :, 1)) ;
    % Dy = squeeze(grad(:, :, 2)) ;

    % figure;
    % subplot(1, 2, 1)
    % imagesc(squeeze(grad(:, :, 1))); title('grad(1)')
    % colormap(bwr)
    % caxis([-1,1])
    % subplot(1, 2, 2)
    % imagesc(Dy); title('Dy')
    % colormap(bwr)
    % caxis([-0.01,0.01])
    % 
    % figure;
    % subplot(1, 2, 1)
    % imagesc(squeeze(grad(:, :, 2))); title('grad(2)')
    % colormap(bwr)
    % caxis([-1,1])
    % subplot(1, 2, 2)
    % imagesc(Dx); title('Dx')
    % colormap(bwr)
    % caxis([-1,1])


    % Compute gradient
    [Dy, Dx] = gradient(DD) ;
    Dx = - Dx ;
    Dy = - Dy ;
    normalization = (Dx.^2 + Dy.^2);
    stationary = find(normalization < 1e-6) ;
    normalization(stationary) = 1;
    Dx = Dx .* (1./sqrt(normalization)) ;
    Dy = Dy .* (1./sqrt(normalization)) ;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('Show gradient')
    imagesc(DD'); hold on;
    quiver(Dx', Dy', 1)
    title('DT with gradient')
    ccaxis = caxis ; % get current caxis
    caxis([0, min(3, ccaxis(2))])
    cb = colorbar() ;
    ylabel(cb, 'DT') ;
    axis equal 
    axis tight
    pause(0.001)
    % saveas(gcf, [outfnBase '_DD.png'])

    % Plot with streamlines
    clf
    disp('Saving cijstream')
    [xx, yy ] = meshgrid(1:size(Dx, 1), 1:size(Dx, 2)) ;
    xx2=xx(1:10:size(xx, 1), 1:10:size(xx, 2));
    yy2=yy(1:10:size(yy, 1), 1:10:size(yy, 2));
    cxyz = stream2(xx, yy, Dx', Dy', xx2', yy2') ;
    %displays log plot if indicated
    if (log_display == 1)
        imagesc(log(cij))
    else
        imagesc(cij)
    end
    hold on;
    quiver(xx, yy, Dy', Dx', 0)
    title('cij with gradient')
    plot(startpt(2), startpt(1), 'ko')
    plot(endpt(2), endpt(1), 'ko')
    streamline(cxyz)
    axis equal 
    axis tight
    pause(0.001)
    % saveas(gcf, [outfnBase '_cijstream.png'])

    % path extraction
    % grad = cat(3, Dx, Dy) ;
    % grad = perform_vf_normalization(grad) ;
    % Dx = squeeze(grad(:, :, 1)) ;
    % Dy = squeeze(grad(:, :, 2)) ;
    % works but is uphill
    % path = stream2(Dx, Dy, endpt(1),endpt(2), str_options);
    % messing around here
    path = stream2(xx, yy, Dx', Dy', endpt(2)-1, endpt(1)-1, str_options);
    % works but does not connect to startpt
    % path = stream2(Dy, Dx, endpt(1), endpt(2), str_options) ;
    % path = stream2(-Dx', -Dy', startpt(1), startpt(2), str_options) ;
    path = path{1} ;
    path = [path(:, 2), path(:, 1)] ;

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % disp('Saving DT')
    % figure(1)
    % imagesc(DD') ;
    % title('DT')
    % colorbar()
    % hold on;
    % yellow = [1,1, 0];
    % plot(path(:, 2), path(:, 1), '.-', 'color', yellow)
    % plot(startpt(:, 2), startpt(:, 1), 'ko')
    % plot(endpt(:, 2), endpt(:, 1), 'ko')
    % axis equal 
    % axis tight
    % pause(0.001)
    % % saveas(gcf, [outfnBase '_DT.png'])
    % 
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % disp('Saving Dx')
    % figure(2);
    % imagesc(Dx') ;
    % colormap(bwr)
    % title('Dx')
    % colorbar()
    % %caxis([-1,1])
    % hold on;
    % plot(path(:, 2), path(:, 1), '.-', 'color', yellow)
    % plot(startpt(:, 2), startpt(:, 1), 'ko')
    % plot(endpt(:, 2), endpt(:, 1), 'ko')
    % xlabel(['time, dataset ', exptIDs{jj}])
    % ylabel(['time, dataset ', exptIDs{ii}])
    % title(['\partial_xW for datasets ' exptIDs{ii} ' and ' exptIDs{jj}])
    % axis equal 
    % axis tight
    % pause(0.001)
    % saveas(gcf, [outfnBase '_Dx.png'])
    % 
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % disp('Saving Dy')
    % figure(3);
    % imagesc(Dy') ;
    % title('Dy')
    % colormap(bwr)
    % colorbar()
    % %caxis([-1,1])
    % hold on; 
    % plot(path(:, 2), path(:, 1), '.-', 'color', yellow)
    % plot(startpt(:, 2), startpt(:, 1), 'ko')
    % plot(endpt(:, 2), endpt(:, 1), 'ko')
    % xlabel(['time, dataset ', exptIDs{jj}])
    % ylabel(['time, dataset ', exptIDs{ii}])
    % title(['\partial_yW for datasets ' exptIDs{ii} ' and ' exptIDs{jj}])
    % axis equal 
    % axis tight
    % pause(0.001)
    % saveas(gcf, [outfnBase '_Dy.png'])
    % 
    % disp('Saving cij')
    % figure(4);
    % %displays log plot if indicated
    % if (log_display == 1)
    %     imagesc(log(cij))
    % else
    %     imagesc(cij)
    % end
    % title('cij')
    % colormap(bwr)
    % colorbar()
    % caxis([-1,1])
    % hold on; 
    % plot(path(:, 2), path(:, 1), '.-', 'color', yellow)
    % plot(startpt(:, 2), startpt(:, 1), 'ko')
    % plot(endpt(:, 2), endpt(:, 1), 'ko')
    % axis equal 
    % axis tight
    % pause(0.001)
    % set(gca,'fontsize', 12);
    % % saveas(gcf, [outfnBase '_cij.png'])
    % 
    % % TEST
    % disp('Saving streamlines')
    % xx2=xx(1:10:size(xx, 1), 1:10:size(xx, 2));
    % yy2=yy(1:10:size(yy, 1), 1:10:size(yy, 2));
    % cxyz = stream2(xx, yy, Dx', Dy', xx2', yy2') ;
    % plot(startpt(:, 2), startpt(:, 1), 'ko')
    % plot(endpt(:, 2), endpt(:, 1), 'ko')
    % figure(5);
    % %displays log plot if indicated
    % if (log_display == 1)
    %     imagesc(log(cij))
    % else
    %     imagesc(cij)
    % end 
    % hold on;
    % quiver(Dx', Dy', 0)
    % streamline(cxyz)
    % pause(1)
    % % saveas(gcf, [outfnBase '_streamlines.png'])

    % cpath = compute_geodesic(DD, flipud(ept), opt);
    % Clip the path near the start point
    d2start = vecnorm(path-startpt,2,2) ;
    path = path(d2start > thres_dist, :) ;
    [tmp, ind] = unique(path(:, 1)) ;
    inds = sort(ind) ;
    path = path(inds, :) ;
    cpath = [endpt; path; startpt] ;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Visualize resulting path
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    close all
    %displays log plot if indicated
    if (log_display == 1)
        imagesc(log(cij))
    else
        imagesc(cij)
    end
    hold on;
    plot(startpt(2), startpt(1), 'ro') 
    plot(endpt(2), endpt(1), 'ks')
    plot(cpath(:, 2), cpath(:, 1), '-') ;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Convert cpath to tpath: x timeline is integer, y timeline
    % is float
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    ntimeptsA = round(endpt(1) - startpt(1) + 1) ; 
    xxx = round(startpt(1)):round(endpt(1)) ;
    yyy = interp1(flipud(cpath(:, 1)), ...
        flipud(cpath(:, 2)), xxx, 'pchip') ;
    tpath = [xxx; yyy]' ;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Continue visualization
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot(tpath(:, 2), tpath(:, 1), 'o') ;
    plot(tpath0(:, 2), tpath0(:, 1), 'k.')
    cb = colorbar() ;
    ylabel(cb, 'cij')            
    title('Path ok? Enter=yes, backspace=no/redo')
    button = waitforbuttonpress() ;
    if button && strcmp(get(gcf, 'CurrentKey'), 'return')
        move_on = true ;
    elseif button && strcmp(get(gcf, 'CurrentKey'), 'backspace')
        move_on = false;
        disp('starting over with path detection')
        Woffset = input('Set weight offset = ') ;
        disp(' --> new Woffset')
        if isempty(Woffset)
            Woffset = input('Set weight offset = ') ;
        end

        cij_exponent = input('Set cij potential exponent = ') ;
        disp(' --> new cij exponent')
        if isempty(cij_exponent)
            cij_exponent = input('Set cij potential exponent = ') ;
        end
    end
end
