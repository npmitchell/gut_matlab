function [V,F] = remesh_smooth_iterate(V,F, lambda,... 
    tar_length, num_iter, protect_constraints)

% Isotropically remesh the surface
if nargin < 3
    lambda = 0.025 ;
end
if nargin  < 4
    tar_length = 10;
end
if nargin < 5
    num_iter = 5;
end
if nargin < 6
    protect_constraints = false;
end
[F, V, ~, ~] = isotropic_remeshing( F, V, ...
    tar_length, num_iter, protect_constraints );

% Attempt to remove localized mesh spikes by Laplacian relaxation
V = relax_mesh_spikes(F, V, deg2rad(60), pi/2, ...
    'uniform', [], 2, 'implicit', 1000);

% Try to remove self-intersections
[intersects, ~] = mesh_self_intersection_3d(F, V);
intCount = 0;
while intersects

    [V, F] = clean_mesh(V, F, 'MinDist', 0, 'MinArea', 0, ...
        'MinAngle', 0, 'SelfIntersections', 'remove', ...
        'SmallTriangles', 'remove');

    % Peform a another isotropic remeshing
    [F, V, ~, ~] = isotropic_remeshing(F, V, ...
        tar_length, num_iter, protect_constraints);

    % Another round of spike relaxation
    V = relax_mesh_spikes(F, V, deg2rad(60), pi/2, ...
        'uniform', [], 2, 'implicit', 1000);

    [intersects, ~] = mesh_self_intersection_3d(F, V);

    intCount = intCount + 1;
    if intCount > 20
        error('Unable to remove self-intersections');
    end

end

% Another round of isotropic remeshing
[F, V, ~, ~] = isotropic_remeshing(F, V, ...
    tar_length, num_iter, protect_constraints);

% Smooth the entire mesh fixing the boundary
V = laplacian_smooth(V, F, 'cotan', [], lambda, 'implicit', V, 10);

% Peform a final isotropic remeshing
[F, V, ~, ~] = isotropic_remeshing(F, V, ...
    tar_length, num_iter, protect_constraints);

% Mesh quality checks ---------------------------------------------

E = edges(triangulation(F, V));

numBdy = numel(DiscreteRicciFlow.compute_boundaries(F));
if (numBdy ~= 0)
    error( ['Mesh at time point T = %d has %d ' ...
        'boundary components'], t, numBdy );
end

eulerChi = size(F,1) + size(V,1) - size(E,1);
if (eulerChi ~= 2)
    error( ['Mesh at time point T = %d is not ' ...
        'a topological sphere'], t );
end

[intersects, ~] = mesh_self_intersection_3d(F, V);
if intersects
    error( ['Mesh at time point T = %d contains ' ...
        'self-intersections'], t );
end