function skelOut = getSkel(SE, overwrite)
%
%
%
% NPMitchell 2020

if nargin < 2
    overwrite = false ;
end

% Check that all fields are present in SE.thickness
allFieldsPresent = true ;
try
    SE.thickness.skels ;
    SE.thickness.DTs ;
    SE.thickness.skelFull ;
    SE.thickness.skel_ss ;
    SE.thickness.pathlength ;
    if ~isempty(SE.thickness.skels) || ~isempty(SE.thickness.DTs) ...
            || ~isempty(SE.thickness.skelFull) ...
            || ~isempty(SE.thickness.skel_ss) ...
            || SE.thickness.pathlength            
        allFieldsPresent = false ;
    end
catch
    allFieldsPresent = false ;
end


if allFieldsPresent 
    if nargout > 0
        skelOut = SE.skel ;
    end
else
    % Labeled image
    SE.getImL() ;
    % Segment out the membrane/gut cross section
    SE.getBW() ;
    % Landmarks for curve extraction
    SE.getLandmarks() ;
    % Path masks for curve extraction
    SE.getPathMasks(); 

    % Find curve for each pair of landmarks
    skelFn = fullfile(SE.zDir.skel, sprintf('skels_%04d.mat', SE.currentTime)) ;
    if exist(skelFn, 'file') && ~overwrite
        skelOut = load(skelFn, 'skels', 'DT', 'skelFull', ...
            'skel_ss', 'pathlength') ;
        SE.skel = skelOut ;
    else
        options = struct() ;
        skels = cell(0) ;
        skel_ss = cell(0) ;
        skelFull = cell(0) ;
        DTs = cell(0) ;

        % For each string of landmarks
        for qq = 1:length(SE.landmarks)

            % preallocate & unpack landmarks
            skels{qq} = cell(0) ;
            DTs{qq} = cell(0) ;
            lms = SE.landmarks{qq} ;

            % For each pair of adjacent landmarks in a string
            for pp = 1:size(lms.v, 1)

                % The start point is landmark pp
                startpt = lms.v(pp, :) ;
                % The end point is landmark qq+1 except last is lm1
                eptId = mod(pp+1, size(lms.v, 1)) ;
                if eptId == 0
                    eptId = size(lms.v, 1);
                end
                endpt = lms.v(eptId, :) ;

                % Mask bw with this landmark pair's pathMask
                if ~isempty(SE.pathMasks{qq}.mask{pp})
                    bwpp = SE.pathMasks{qq}.mask{pp}' .* SE.bw ;
                    % ensure empty DT and DD for fresh computation
                    options.DT = [] ;
                    options.DD = [] ;
                else
                    bwpp = SE.bw ;
                end

                % Extract the path connecting landmarks
                % imshow(bwpp')
                % title('curve extraction binary')
                % waitfor(gcf)
                [skelpp, DTs{qq}{pp}, DD] = ...
                    extractCurve(bwpp, startpt, endpt, options) ;

                % If no segment/pair-specific masking, clear DT & DD
                if ~isempty(SE.pathMasks{qq}.mask{pp})
                    if ~isfield(options, 'DT')
                        options.DT = DTs{qq}{pp} ;
                    end
                    if ~isfield(options, 'DD')
                        options.DD = DD ;
                    end
                end

                % Collate skeletons into cell
                skels{qq}{pp} = skelpp ;
                
                % Preview the result
                clf
                imshow(SE.bw'); hold on;
                for tmpqq = 1:length(skels)
                    for sii = 1:length(skels{tmpqq})
                        plot(skels{tmpqq}{sii}(:, 1), ...
                            skels{tmpqq}{sii}(:, 2), '-')
                    end
                end
                for tmpid = 1:size(lms.v, 1)
                    plot(lms.v(tmpid, 1), lms.v(tmpid, 2), 'o')    
                end
                pause(1)

            end

            % Resample full curve circuit as one path
            skelFull{qq} = [] ;
            for pp = 1:length(skels{qq})
                skelFull{qq} = cat(1, skelFull{qq}, skels{qq}{pp}) ;
            end
            skel_ss{qq} = curvspace(skelFull{qq}, SE.npts_skel) ;
            pathlength{qq} = ss_from_xyz(skel_ss{qq}) ;
        end 
        
        % Preview the result for this cycle
        clf
        imshow(SE.imL); hold on;
        h = imshow(SE.bw') ;
        set(h, 'AlphaData', 0.1 * double(SE.bw(:, :, 1)));
        for qq = 1:length(skels)
            plot(skel_ss{qq}(:, 1), skel_ss{qq}(:, 2), '.')
        end
        skelFigFn = fullfile(SE.zDir.skel, ...
            sprintf('skels_cycle%04d_%04d.png', qq, SE.currentTime)) ;
        saveas(gcf, skelFigFn)    

        % Save the skeletons
        save(skelFn, 'skels', 'DTs', 'skelFull', 'skel_ss', 'pathlength')

        SE.skel = struct('skels', skels, ...
            'DTs', DTs, 'skelFull', skelFull, ...
            'skel_ss', skel_ss, 'pathlength', pathlength) ;
        
        if nargout > 0
            skelOut = SE.skel ;
        end
    end
end