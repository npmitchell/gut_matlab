% Make movie of projected tissue onto mesh surface

%set(0,'defaultFigureColor',[0 0 0])
scrsz = get(0,'ScreenSize');
%set(0,'defaultFigureColor','remove')
%%
patchName = 'cylinder2_index';
eGrids_apical = SOI.embedding.getPatch(patchName).apply();
data = SOI.getField('data_MIP');
colordef black;
%%
for time = 1:2
    %plot3(mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),'.'), axis equal
    %trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3)), axis equal
    %set(0,'defaultFigureColor',[0 0 0],'defaultFigureInvertHardcopy','off')
    %colordef white
    %colordef black;

    cmp = [217,217,215]/255;
    %h = trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),...
    %    'EdgeColor','none','FaceColor',cmp, 'FaceLighting','phong');
    %
    figure('Position',[1 scrsz(4)/1 scrsz(3)/1 scrsz(4)/1],'Color',[0 0 0])
    im = double(data(time).getPatch(patchName).apply{1});
    h = surf(eGrids_apical{3},-eGrids_apical{1},eGrids_apical{2},im);
    axis equal
    shading interp 
    grid off
    colormap gray
    set(gca,'CLim',[7000 20000]);
    view([-180 -85])
    axis off
    axis([-1000 1000 -400 400 -400 400])

    %camva(6.5)
    %camlight top

    set(h,'AmbientStrength',.8)
    %g = light('position',[0  1 0]);
    %set(g,'Color',[1 1 1]/2);
    set(gcf,'Color',[0 0 0]);
    %
    saveas(gcf,sprintf('figs_apical/Time_%06d.png',time))
    close all;
end


%% 

for t = 50
    
    times = max(1,t-2):min(95,t+2);
    v = [];
    for time = times
        outputMesh = fullfile([sprintf('pointCloud_T%06d_mesh.ply',time)]);
        mesh = read_ply_mod(outputMesh);
        v = [v;mesh.v];
    end
    
    clear OBJ
    % swap x and y
    OBJ.vertices = v;
    OBJ.objects(1).type='f';
    OBJ.objects(1).data.vertices=[];
    write_wobj(OBJ, sprintf(fullfile('obj_smooth/test%06d.obj'),time));
    
end


%% 

% rotate the image and record a movie of it.  
%uiopen('/Volumes/Data/fly/embryo/myo-WT-data/BestEstimate_TP32_Residual_Symmetric_12.fig',1)
%set(gca,'Clim',[0 1])
uiopen('/Volumes/Data/fly/embryo/myo-WT-data/Measurement_TP32_Flows_Symmetric_11.fig',1)

axis off
colorbar off
grid off
view([0 0])

dphi = 5;
for phi = 1 : 360/dphi
    
    camorbit(0,dphi)
    pause(.1)
    mov(phi) = getframe;
end

dtheta = 10;
for theta = 1 : 360/dtheta
    camorbit(dtheta,0,'camera')
    pause(.1)
    mov(phi+theta) = getframe;
end

%% adjust size of image; 
sizes = zeros(size(mov,2),3);
for k = 1 : size(mov,2)
    sizes(k,:) = size(mov(k).cdata); 
end
max_sizes = max(sizes);
mov2 = mov; 
for k = 1 : size(mov,2)
    mov2(k).cdata = 205*ones(max_sizes(1),max_sizes(2),max_sizes(3),'uint8');
    
    d_size = floor((max_sizes-sizes(k,:))/2);
    mov2(k).cdata((d_size(1)+1):(d_size(1)+sizes(k,1)),(d_size(2)+1):(d_size(2)+sizes(k,2)),:) = mov(k).cdata;
end