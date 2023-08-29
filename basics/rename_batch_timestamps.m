%% Rename file timestamps for unpacked data

% datdir = '/mnt/crunch/actb2-mem-cherry_h2afva-gfp/202307111346_4views_1p4um_1ms/' ;
% datdir = '/mnt/crunch/actb2-mem-cherry_h2afva-gfp/202307111346_continue2106_4views_1p4um_1ms/' ;
% datdir = '/mnt/crunch/48YGAL4klarGFPnlsCAAXmCh/202308101028_180s_1p4um_2mW2mW_48YG4knlsGFPCAAXmCh/Time4views_0p25_3p0msexposure/data/'
datdir = '/mnt/crunch/48YGAL4klarLamGFPCAAXmCh/202208122137_180s_1p4um_2mW2mW_48YG4kLamGFPCAAXmCh_4views_0p25_0p5ms/data/' ;
channels = [1,2] ;
angles = 0:45:360 ;
times = 0:199 ;


for ang = angles 
    for ch =channels
        for tt = times
            % Here we change the timestamp of the file
            tout = floor(tt * 0.5) + 2;
           
            % Consider each file that matches the naming convention
            searchstr = ['Time_' sprintf('%06d', tt) '_Angle_' num2str(ang) '_c' num2str(ch) '_ls_1.ome.tif'] ;
            fns = dir(fullfile(datdir, searchstr)) ;

            % Ensure there is only one match
            if length(fns) > 1
                error('More than one match for this filename!')
            elseif length(fns) == 1
                ii = 1 ;
                disp(['Seeking: ' searchstr])
            
    
                % get the filename to move
                fn = fullfile(fns(ii).folder, fns(ii).name)  ;
    
                % Make the output filename
                outname = fullfile(fns(ii).folder, ...
                    ['oTime_' sprintf('%06d', tout) '_Angle_' num2str(ang) '_c' num2str(ch) '_ls_1.ome.tif']) ;
                disp(['renaming ' fns(ii).name ' > ' outname])
                
                % Make sure we don't overwrite anything
                if exist(outname, 'file')
                    error('Output filename already exists! This would overwrite a file!')
                else
        
                    % rename the file here
                    [status ,msg] = movefile(fn, outname) ;
                    if ~status
                        error(['Could not move file:' msg])
                    end
                end
            else
                disp(['Could not find: ' searchstr])
            end

        end
    end
end