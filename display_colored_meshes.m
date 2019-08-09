% Show meshes colored by the data interpolated at the grid vertices

% load a mesh 

addpath_recurse('/mnt/data/code/imsaneV1.2.3/generalfunctions/')
addpath_recurse('/mnt/crunch/djcislo/MATLAB/CGAL_Code/')
addpath_recurse('/mnt/data/code/gptoolbox/')
plyfn = 'mesh_apical_000110.ply' ;
mesh = read_ply_mod(plyfn) ;
mesh.v
surf(x,y,z,C)