%% Used Imaris instead!

colors = define_colors(7) ; 
blue   = colors(1, :) ;
red    = colors(2, :) ;
yellow = colors(3, :) ;
purple = colors(4, :) ;
green  = colors(5, :) ;
sky    = colors(6, :) ;
maroon = colors(7, :) ;
lscolor = sky ;

mesh = read_ply_mod('/mnt/data/tubular_test/activeGel/202203221330_region4_ATayarGel/imagesequence/msls_output/mesh_000014.ply') ;
h = trimesh(triangulation(mesh.f, mesh.v), 'EdgeColor', 'none', ...
                    'facecolor', lscolor) ;

axis equal
hold on;
axis off

lighting gouraud    % preferred method for lighting curved surfaces
material dull    % set material to be dull, no specular highlights

% Lighting and view 
view([0, 270])
lgt = camlight('headlight') ;
saveas(gcf, '/mnt/data/tubular_test//activeGel/tubular_fig_202203221330_region4_ATayarGel_000014.png')
            
            
            