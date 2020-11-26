function eulerChar = eulerCharacteristic(mesh)

    % The #Ex2 edge connectivity list of the mesh
    edgeTri = edges( triangulation( mesh.f, mesh.v ) );
    % Check that the input mesh is a topological cylinderChar
    eulerChar = ( size(mesh.v,1) - length(edgeTri) + size(mesh.f,1) ) ;
end