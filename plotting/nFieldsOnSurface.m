function [axs, cbs, meshHandles] = ...
    nFieldsOnSurface(meshes, fields, options)
%nFieldsOnSurface(meshes, fields, options)
%
% Inputs
% ------
% meshes : struct with fields f and v OR 2x1 cell array of faces and vertices
%   the meshes on which to plot the scalar fields
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

% Default options
nfields = length(fields) ;
    
% Define axes
set(gcf, 'visible', 'off')
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

% Define
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

%% Panel 
for qq = 1:nfields
    % Unpack meshes if there are multiple
    if isa(meshes, 'cell')
        if isa(meshes{qq}, 'struct') || isa(meshes{qq}, 'cell')
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
            meshHandles{qq} = trisurf(mesh.f, mesh.v(:, 1), ...
                mesh.v(:, 2), mesh.v(:, 3), ...
                'FaceVertexCData', fields{qq}, 'edgecolor', edgecolor) ;
        elseif isa(mesh, 'cell')    
            meshHandles{qq} = trisurf(mesh{1}, mesh{2}(:, 1), ...
                mesh{2}(:, 2), mesh{2}(:, 3), ...
                'FaceVertexCData', fields{qq}, 'edgecolor', edgecolor) ;
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
            title(labels{3}, 'Interpreter', 'Latex')   
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
end

end