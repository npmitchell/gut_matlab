function [axs, cbs, meshHandles] = ...
    nFieldsOnSurface(meshes, fields, options)
%[axs, cbs, meshHandles] = nFieldsOnSurface(meshes, fields, options)
%
% Inputs
% ------
% meshes : struct with fields f and v 
%          OR 2x1 cell array of faces and vertices
%          OR triangulation object with fields ConnectivityList and Points
%   The meshes on which to plot the scalar fields. Can be 2D or 3D.
% fields : #vertices x 1 or #faces x 1 float array or length2 cell of
%   magnitude and angle for nematic fields
%   The fields to plot on the surfaces.
% options : optional struct with optional fields
%   clim : numeric or 2x1 numeric
%       colorlimit, by default set to [-clim, clim] if single value
%   clims : numeric, overwrites clim 
%       colorlimit, by default set to [-clim1, clim1] if single value
%   axs : axis instances for each field
%   edgecolor : color specifier
%   cmap : colormap for both axes
%   cmaps : colormap for axis1, overwrites cmap
%   labels : cell of strings
%       titles for each subplot
%   makeCbar : nfields x 1 bool array
%       whether to make colorbar for each axis
%   xyzlims : 3x2 numeric or (nfields)x1 cell array of 3x2 numeric arrays
%       the xyz limits for all panels as [xmin,xmax;ymin,ymax;zmin,zmax]
%   xlim : length 2 numeric array or length nfields cell array of length 2
%       numeric arrays
%       xlim values if xyzlims not provided
%   ylim : length 2 numeric array or length nfields cell array of length 2
%       numeric arrays
%       ylim values if xyzlims not provided
%   zlim : length 2 numeric array or length nfields cell array of length 2
%       numeric arrays
%       zlim values if xyzlims not provided
%   view : length nfields cell array of length 2 arrays
%       viewing angle for each panel
%   
%   
% Returns
% -------
% [axs, cbs, meshHandles] : handles for figure objects
%
% Example usage
% -------------
% labels = {'$\Re \mu$', '$\rho/2$', '$\Re\mu - \rho/2$'} ;
% options.labels = labels ;
% twoScalarFieldsOnSurface(meshes, sf1, sf2, sf3, options)
% sgtitle('deformation at t=0')  
%
% NPMitchell 2020

%% Default options
nfields = length(fields) ;
axisOff = false ;
visible = 'off' ;
if isfield(options, 'visible')
    visible = options.visible ;
end

%% Define axes
set(gcf, 'visible', visible)
if isfield(options, 'axs')
    axs = options.axs ;
else
    if nfields < 4
        for qq = 1:nfields
            axs{qq} = subplot(1, nfields, qq) ;
        end
    elseif nfields == 4
        for qq = 1:nfields
            axs{qq} = subplot(2, 2, qq) ;
        end
    elseif nfields < 7
        for qq = 1:nfields
            axs{qq} = subplot(2, 3, qq) ;
        end
    else
        error('code for default axis arrangement for this number of axes')
    end
end

% Define climits
if isfield(options, 'clim')
    clim = options.clim ;
end
if isfield(options, 'clims')
    clims = options.clims ;
end

% Define xyzlims
if isfield(options, 'xlim')
    if isa(options.xlim, 'cell')
        try
            assert(numel(options.xlim) == nfields)
        catch
            error('number of elements in options.xlim ~= nfields')
        end
        xlimits = options.xlim ;
    else
        xlimit = options.xlim ;
    end
end
if isfield(options, 'ylim')
    if isa(options.ylim, 'cell')
        try
            assert(numel(options.ylim) == nfields)
        catch
            error('number of elements in options.ylim ~= nfields')
        end
        ylimits = options.ylim ;
    else
        ylimit = options.ylim ;
    end
end
if isfield(options, 'zlim')
    if isa(options.zlim, 'cell')
        try
            assert(numel(options.zlim) == nfields)
        catch
            error('number of elements in options.zlim ~= nfields')
        end
        zlimits = options.zlim ;
    else
        zlimit = options.zlim ;
    end
end
if isfield(options, 'xyzlims')
    if isa(options.xyzlims, 'cell')
        try
            assert(numel(options.xyzlim) == nfields)
        catch
            error('number of elements in options.xyzlim ~= nfields')
        end
        xlimits = {options.xyzlims{1}(1, :)} ;
        ylimits = {options.xyzlims{1}(2, :)} ;
        zlimits = {options.xyzlims{1}(3, :)} ;
        for qq = 2:nfields
            xlimits = cat(1, xlimits, options.xyzlims{qq}(1, :)) ;
            ylimits = cat(1, ylimits, options.xyzlims{qq}(2, :)) ;
            zlimits = cat(1, zlimits, options.xyzlims{qq}(3, :)) ;
        end
    else
        xlimit = options.xyzlims(1, :) ;
        ylimit = options.xyzlims(2, :) ;
        zlimit = options.xyzlims(3, :) ;
    end
end

% Unpack view/views
if isfield(options, 'view')
    views = options.view ;
end

% axis off?
if isfield(options, 'axisOff')
    axisOff = options.axisOff ;
end

% Define edgecolor
edgecolor = 'none' ;
if isfield(options, 'edgecolor')
    edgecolor = options.edgecolor ;
end

% Define colormaps
if isfield(options, 'colormap')
    cmap = options.colormap ;
else
    caxis([-1, 1])
    cmap = blueblackred(256) ;
end
if isfield(options, 'colormaps')
    cmaps = options.colormaps ;
end

if isfield(options, 'labels')
    labels = options.labels ;
else
    labels = [] ;
end

if isfield(options, 'makeCbar')
    makeCbar = options.makeCbar ;
else
    makeCbar = ones(nfields, 1) ;
end

if isfield(options, 'masterCbar')
    masterCbar = options.masterCbar ;
    if masterCbar
        if isfield(options, 'masterCbarPosition')
            masterCbarPosition = options.masterCbarPosition ;
        else
            masterCbarPosition = [0.39, 0.05, 0.22, 0.0341];
        end
    end
else
    masterCbar = false ;
end

%% Panels
for qq = 1:nfields
    % Unpack meshes if there are multiple
    if isa(meshes, 'cell')
        if isa(meshes{qq}, 'struct') || ...
                isa(meshes{qq}, 'cell') || ...
                isa(meshes{qq}, 'triangulation')
            mesh = meshes{qq} ;
        else
            mesh = meshes ;
        end
    else
        mesh = meshes ;
    end
    
    set(gcf,'CurrentAxes', axs{qq})
    
    % Check if scalar field
    if ~isa(fields{qq}, 'cell')
        % Scalar field
        if isa(mesh, 'struct')
            if size(mesh.v, 2) == 3
                meshHandles{qq} = trisurf(mesh.f, mesh.v(:, 1), ...
                    mesh.v(:, 2), mesh.v(:, 3), ...
                    'FaceVertexCData', fields{qq}(:), 'edgecolor', edgecolor) ;
            elseif size(mesh.v, 2) == 2
                meshHandles{qq} = trisurf(mesh.f, mesh.v(:, 1), ...
                    mesh.v(:, 2), 0*mesh.v(:, 1), ...
                    'FaceVertexCData', fields{qq}(:), 'edgecolor', edgecolor) ;
            end
        elseif isa(mesh, 'cell')    
            if size(mesh{2}, 2) == 3
                meshHandles{qq} = trisurf(mesh{1}, mesh{2}(:, 1), ...
                    mesh{2}(:, 2), mesh{2}(:, 3), ...
                    'FaceVertexCData', fields{qq}(:), 'edgecolor', edgecolor) ;
            elseif size(mesh{2}, 2) == 2
                meshHandles{qq} = trisurf(mesh{1}, mesh{2}(:, 1), ...
                    mesh{2}(:, 2), 0*mesh{2}(:, 1), ...
                    'FaceVertexCData', fields{qq}(:), 'edgecolor', edgecolor) ;
            end
        elseif isa(mesh, 'triangulation') 
            if size(mesh.Points, 2) == 3
                meshHandles{qq} = trisurf(mesh, ...
                    'FaceVertexCData', fields{qq}(:), 'edgecolor', edgecolor) ;
            elseif size(mesh.Points, 2) == 2
                meshHandles{qq} = trisurf(mesh.ConnectivityList, ...
                    mesh.Points(:, 1), mesh.Points(:, 2), 0*mesh.Points(:, 1), ... 
                    'FaceVertexCData', fields{qq}(:), 'edgecolor', edgecolor) ;
            end    
        end
        axis equal
        
        if makeCbar(qq)
            cbs{qq} = colorbar('location', 'southOutside') ;
        end
        
        if exist('clims', 'var')
            if numel(clims{qq}) == 1
                caxis([-clims{qq}, clims{qq}])
            elseif numel(clims{qq}) == 2
                caxis([clims{qq}(1), clims{qq}(2)])
            end
        elseif exist('clim', 'var')
            if numel(clim) == 1
                caxis([-clim, clim])
            elseif numel(clim) == 2
                caxis([clim(1), clim(2)])
            else
                error('clim provided has > 2 elements')
            end
        else
            caxis([-max(abs(fields{qq})), max(abs(fields{qq}))])
        end
        
        if ~isempty(labels)
            title(labels{qq}, 'Interpreter', 'Latex')   
        end
        if isfield('cmaps', 'var')
            colormap(cmaps{qq})
        else
            colormap(cmap)
        end
    else
        % Nematic field or polar field?
        % todo: code for polar field
        
        % Nematic field
        % Unpack nematic field
        dev = fields{qq}{1} ;
        theta = fields{qq}{2} ;
        
        if exist('clims', 'var')
            if numel(clims{qq}) == 1
                clim_dev = clims{qq} ;
            elseif numel(clims{qq}) == 2
                clim_dev = clims{qq}(2) ;
            end
        elseif exist('clim', 'var')
            if numel(clim) == 1
                clim_dev = clim ;
            elseif numel(clim) == 2
                clim_dev = clim(2) ;
            else
                error('clim provided has > 2 elements')
            end
        else
            clim_dev = max(abs(dev(:))) ;
        end
        caxis([0, clim_dev])
        
        % Intensity from dev and color from the theta
        pm256 = phasemap(256) ;
        indx = max(1, round(mod(2*theta(:), 2*pi)*size(pm256, 1)/(2 * pi))) ;
        colors = pm256(indx, :) ;
        colors = min(dev(:) / clim_dev, 1) .* colors ;
        
        if isa(mesh, 'struct')
            meshHandles{qq} = trisurf(mesh.f, mesh.v(:, 1), ...
                mesh.v(:, 2), mesh.v(:, 3), ...
                'FaceVertexCData', colors, 'edgecolor', edgecolor) ;
        elseif isa(mesh, 'cell')    
            meshHandles{qq} = trisurf(mesh{1}, mesh{2}(:, 1), ...
                mesh{2}(:, 2), mesh{2}(:, 3), ...
                'FaceVertexCData', colors, 'edgecolor', edgecolor) ;
        end
        axis equal
        if ~isempty(labels)
            title(labels{qq}, 'Interpreter', 'Latex')   
        end

        % Colorbar and phasewheel
        if makeCbar(qq)
            colormap(gca, phasemap)
            cbs{qq} = cell(2, 1) ;
            cbs{qq}{1} = phasebar('colormap', phasemap, ...
                'location', [0.82, 0.12, 0.1, 0.135], 'style', 'nematic') ;
            shrink = max(0.6 - 0.1 * (mod(nfields, 3)-2), 0.1) ;
            % axis off
            % view(2)
            cbs{qq}{2} = colorbar('location', 'southOutside') ;
            drawnow
            axpos = get(axs{qq}, 'position') ;
            cbpos = get(cbs{qq}{2}, 'position') ;
            set(cbs{qq}{2}, 'position', [cbpos(1), cbpos(2), cbpos(3)*shrink, cbpos(4)])
            set(axs{qq}, 'position', axpos) ;
            hold on;
        end
        
        % Set colorlimits
        if exist('clims', 'var')
            if numel(clims{qq}) == 1
                caxis([0, clims{qq}])
            elseif numel(clims{qq}) == 2
                caxis([clims{qq}(1), clims{qq}(2)])
            end
        else
            if numel(clim) == 1
                caxis([0, clim])
            elseif numel(clim) == 2
                caxis([clim(1), clim(2)])
            else
                error('clim provided has > 2 elements')
            end
        end
        colormap(gca, gray) 
    end
    
    % Set view from "views"
    if exist('views', 'var')
        if isa(views, 'cell')
            if ~isempty(views{qq})
                view(views{qq}(1), views{qq}(2))
            end
        elseif numel(views) == 2
            view(views(1), views(2))
        else
            error('Number of elements in options.view is not 2 or nfields')
        end
    end
    
    % Are axes off?
    if axisOff
        axis off
    end
    
    % Set xyzlims
    if exist('xlimits', 'var')
        xlim(xlimits{qq}) ;
    elseif exist('xlimit', 'var')
        xlim(xlimit) ;
    end
    if exist('ylimits', 'var')
        ylim(ylimits{qq}) ;
    elseif exist('ylimit', 'var')
        ylim(ylimit) ;
    end
    if exist('zlimits', 'var')
        zlim(zlimits{qq}) ;
    elseif exist('zlimit', 'var')
        zlim(zlimit) ;
    end
end

% Master colorbar
if masterCbar
    axpos = get(gca, 'position') ;
    if exist('cbs', 'var')
        cbs{length(cbs)+1} = colorbar('location', 'southOutside') ;
    else
        cbs{1} = colorbar('location', 'southOutside') ;
    end
    set(gca, 'position', axpos) ;
    set(cbs{length(cbs)}, 'position', masterCbarPosition) ;
end

end