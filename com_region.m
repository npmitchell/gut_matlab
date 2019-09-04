function com = com_region(probability_grid, thres, varargin)
    % com = TRAINING_COM(probability_grid, thres, gridxyz)
    % Given a probability cloud, find the center of mass of the largest 
    % connected region of the probability cloud above
    % some threshold, with "mass" proportional to probability
    %
    % INPUTS
    % probability_grid
    % thres : float
    %   Threshold probability used to segment the probability cloud into
    %   connected regions
    % mesh_vertices : N x 3 float array
    %   Positions of the mesh vertices, as N x 3 array
    % varargin: struct
    %   Options struct, with fields xyzgrid and check
    %   options.xyzgrid: 3d float or int array 
    %       xyzgrid positional values matching the probability_grid
    %   options.check: boolean
    %       whether to display
    %
    % OUTPUTS
    % com : 3x1 float 
    %   the position of the cernter of mass of the chunk of probability 
    %   cloud above the supplied threshold thres
    %
    % EXAMPLE USAGE
    % probability_grid = zeros(10,10,10) ;
    % probability_grid(1:5,1:3,1:3) = 1 ; 
    % vertices =  ;
    % match_training_to_vertex(probability_grid, 0.5, vertices)
    
    if isfield(varargin{1}, 'check')
        check = varargin{1}.check ;
    else
        check = false ;
    end
    
    bwcc = bwconncomp(probability_grid > thres) ; 
    npix = cellfun(@numel,bwcc.PixelIdxList);
    [~, indexOfMax] = max(npix); 
    % isolate the largest connected component
    biggest = zeros(size(probability_grid));
    biggest(bwcc.PixelIdxList{indexOfMax}) = 1;
    % Now multiply with probability
    mass = probability_grid .* biggest ;
    % Get center of mass. There are two ways
    if isfield(varargin{1}, 'xyzgrid')
        xyzgrid = varargin{1}.xyzgrid ;
        % Method 1, using mean
        mean_mass = mean(mass(:)) ;
        comX = mean(mass(:) .* xyzgrid(1, :, :, :)) / mean_mass ;
        comY = mean(mass(:) .* xyzgrid(2, :, :, :)) / mean_mass ;
        comZ = mean(mass(:) .* xyzgrid(3, :, :, :)) / mean_mass ;
        com = [ comX comY comZ ] ; 
    else 
        props = regionprops(true(size(mass)), mass, 'WeightedCentroid');
        % Check for flip in xy
        com_tmp = props.WeightedCentroid ;
        com = [ com_tmp(1) com_tmp(2) com_tmp(3) ] ;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Check the com
    if check
        % Show relative to the mesh
        iso = isosurface(probability_grid, 0.5) ;
        patch(iso,'facecolor',[1 0 0],'facealpha',0.1,'edgecolor','none');
        view(3)
        camlight
        hold on
        scatter3(com(1), com(2), com(3))
        axis equal
        hold off
        title('COM and pointcloud. Click any button to continue')
        
        % Show the figure
        disp('Showing com and pointcloud. Click any button to continue')
        drawnow;
        disp(clock)
        try
            waitforbuttonpress
            % Close figure or leave it open
            close(h)
            disp('mouse or key pressed')
        catch
            disp('figure closed')
        end
        disp(clock)
        
        % Show slices
        % tmp1 = round(com(1)) ;
        % tmp2 = round(com(2)) ;
        % imagesc(

        % Plot each section of the intensity data 
        for jj=1:size(mass,1)
            imshow(squeeze(mass(jj,:,:)))
            title(['mass slice ' num2str(jj)])
            pause(0.01)
        end
        for jj=1:size(probability_grid,1)
            imshow(squeeze(probability_grid(jj,:,:)))
            title(['prediction slice' num2str(jj)])
            pause(0.01)
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
return