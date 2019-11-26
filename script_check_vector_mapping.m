%% Given vector in 3D wrt triangle, build mapped vector wrt mapped triangle
% NPM 2019 -- for debugging ideas

preview = true ;  % Check intermediate results

%% Compute the tangential velocities in plane
% x123ap represent vertices of a triangle in 3D (for ex that have been 
% obtained by projecting ('p') a triangle onto tangent plane)

% Build example data
x0 = [ 22    22    22    22    22   22    22    22    22    22 ];
y0 = [ 22    38    54    70    86   102   118   134   150   166 ];
x123ap = [  293.7233,  291.5354,  298.6046;
          334.0822,  335.2899,  342.6451;
          380.6488,  380.9907,  374.7066;
          413.5333,  417.3939,  418.6768;
          440.5910,  433.3929,  436.4298;
          456.7345,  450.2074,  459.7222;
          466.7431,  467.8046,  475.2768;
          477.6420,  477.0163,  468.0068;
          459.3940,  461.0635,  465.6446;
          442.1477,  438.4495,  432.2021] ;
y123ap = 1.0e+03 * [1.0484,    1.0440,    1.0475;
            1.0551,    1.0505,    1.0568;
            1.0674,    1.0672,    1.0674;
            1.0574,    1.0554,    1.0593;
            1.0392,    1.0418,    1.0415;
            1.0319,    1.0344,    1.0291;
            1.0384,    1.0388,    1.0362;
            1.0333,    1.0325,    1.0381;
            1.0320,    1.0328,    1.0314;
            1.0263,    1.0217,    1.0282 ];
z123ap = [430.5984,  441.9383, 439.1335;
          431.8409,  440.4366,  432.8782;
          418.4509,  412.2191,  418.9859;
          398.0920,  397.0920,  385.8168;
          358.1608,  350.1777,  359.5054;
          319.7215,  320.9559,  328.4940;
          274.3177,  285.1634,  274.1218;
          233.8080,  224.4000,  233.3765;
          183.9125,  189.2936,  187.1748;
          152.3325,  141.1804,  153.9470 ] ;
normals = [-0.3294,    0.9004,    0.2842;
           -0.2318,    0.8461,    0.4799;
           -0.0152,    0.9988,   -0.0470;
            0.4912,    0.8064,    0.3292;
            0.4517,    0.8832,   -0.1261;
            0.3766,    0.9113,    0.1667;
            0.2425,    0.9682,   -0.0613;
            0.4458,    0.8884,   -0.1094;
            0.1915,    0.9621,   -0.1942;
            0.1011,    0.9052,   -0.4128 ];
pt0 = 1.0e+03 * [0.2965,    1.0473,    0.4372;
                0.3415,    1.0562,    0.4333;
                0.3788,    1.0673,    0.4149;
                0.4177,    1.0583,    0.3897;
                0.4398,    1.0396,    0.3577;
                0.4556,    1.0320,    0.3214;
                0.4714,    1.0374,    0.2774;
                0.4731,    1.0352,    0.2303;
                0.4616,    1.0323,    0.1872;
                0.4391,    1.0253,    0.1494 ];
xtria = [17.4388,   25.9329,   23.1793;
           18.3193,   27.1760,   21.9106;
           24.9450,   20.6358,   23.8249;
           24.7818,   27.7430,   19.8902;
           23.1845,    9.6808,   19.2195;
           24.1444,   13.2477,   27.2213;
           17.0009,   18.7321,   25.8398;
           27.0655,   25.7934,   16.6639;
           22.1661,   20.1502,   25.6980;
           23.4710,   24.7085,   14.2569 ] ;
ytria = [20.3327,   21.0571,   22.8590;
           35.1250,   35.0346,   38.4999;
           53.9381,   55.2250,   51.9133;
           67.4183,   68.1259,   71.0741;
           86.1154,   86.5582,   84.2550;
          102.7900,  100.9031,  100.7149;
          118.4409,  114.9165,  119.6392;
          133.6082,  136.1668,  132.4651;
          151.4253,  149.1347,  149.8995;
          164.6303,  168.5029,  165.8694 ];
v0 = [1.3275,   -3.0032,    7.0538;
       -1.6385,   -1.6555,    4.0994;
        7.2750,   -1.6218,    3.4775;
        5.0790,   -4.0507,    1.9556;
        2.2477,   -3.8687,    0.4626;
        3.8431,   -5.1261,   -0.5771;
       -2.0869,   -1.0925,   -0.7176;
        0.0058,    0.8942,    3.0129;
       -0.7571,    1.0433,   -4.0004;
       -4.8811,   -1.7780,   -3.9824 ]; 
v0t = [0.9532,   -1.9801,    7.3768;
       -1.4191,   -2.4563,    3.6453;
        7.2461,    0.2699,    3.3885;
        5.1418,   -3.9476,    1.9977;
        3.3589,   -1.6962,    0.1525;
        5.0934,   -2.1005,   -0.0235;
       -1.7183,    0.3790,   -0.8108;
       -0.2026,    0.4790,    3.0641;
       -1.0704,   -0.5303,   -3.6828;
       -4.8347,   -1.3626,   -4.1719]; 

%% Create bases of the triangle and decompose vector v0
dx = x123ap - pt0(:, 1) ;
dy = y123ap - pt0(:, 2) ;
dz = z123ap - pt0(:, 3) ;

% dr1 is the vector from pt0 to vertex1 in 3D
r1 = [dx(:, 1), dy(:, 1), dz(:, 1)] ;
r2 = [dx(:, 2), dy(:, 2), dz(:, 2)] ;
r3 = [dx(:, 3), dy(:, 3), dz(:, 3)] ;
% Normalize the vectors to get hats
r1hat = r1 ./ sqrt(sum(r1.^2, 2)) ;
r2hat = r2 ./ sqrt(sum(r2.^2, 2)) ;
r3hat = r3 ./ sqrt(sum(r3.^2, 2)) ;

% check 3d triangle vectors

% dq1 is the vector from (xx, yy) to vertex1 in the plane
q1 = [xtria(:, 1) - x0(:), ytria(:, 1) - y0(:) ];
q2 = [xtria(:, 2) - x0(:), ytria(:, 2) - y0(:) ];
q3 = [xtria(:, 3) - x0(:), ytria(:, 3) - y0(:) ];
% Normalize the vectors
q1hat = q1 ./ sqrt(sum(q1.^2, 2)) ;
q2hat = q2 ./ sqrt(sum(q2.^2, 2)) ;
q3hat = q3 ./ sqrt(sum(q3.^2, 2)) ;
% Build vectors perpendicular to q1hat
q2p = q2hat - dot(q2hat, q1hat, 2) .* q1hat ;
q2phat = q2p ./ sqrt(sum(q2p.^2, 2)) ;
q3p = q3hat - dot(q3hat, q1hat, 2) .* q1hat ;
q3phat = q3p ./ sqrt(sum(q3p.^2, 2)) ;

% check 2d triangle vectors
if preview
    fig = figure('Visible', 'On');
    plot([xtria'; xtria(:, 1)'], [ytria'; ytria(:, 1)'])
    hold on
    plot(x0(:), y0(:), 'k.')
    quiver(x0(:), y0(:), q1(:, 1), q1(:, 2), 0)
    quiver(x0(:), y0(:), q2(:, 1), q2(:, 2), 0)
    quiver(x0(:), y0(:), q3(:, 1), q3(:, 2), 0)
    title('triangles containing flor field sample points')
    xlabel('x [pix]')
    ylabel('y [pix]')
    axis equal
    waitfor(fig)            
end

% Check that the tangential velocity is out of the triangle plane
assert(all(dot(normals, r1, 2) < 1e-1))

% Use dr1,dr2,dr3 to form basis of relative direction for v0t
r2p = r2hat - dot(r2hat, r1hat, 2) .* r1hat ;
r2phat = r2p ./ sqrt(sum(r2p.^2, 2)) ;
r3p = r3hat - dot(r3hat, r1hat, 2) .* r1hat ;
r3phat = r3p ./ sqrt(sum(r3p.^2, 2)) ;
a1 = dot(v0t, r1hat, 2) ;
a2 = dot(v0t, r2phat, 2) ;
a3 = dot(v0t, r3phat, 2) ;

% Velocity in plane is v0t_2d = a1*dq1 + a2*dq2 if conformal. 
% Instead average over different combinations of linearly
% independent vectors:
% v0t_2d = <ai*dqi + aj*dqj>, ij=(1,2), (2,3), (1,3)
v0t_2d_1 = a1 .* q1hat + a2 .* q2phat ;
% Note that v0t_2d_1 is identical to a different choice of basis given
% orthogonalization wrt q1hat
% v0t_2d_1alt = a1 .* q1hat + a3 .* q3phat ;

% Other bases
v0t_2d_2 = a2b .* q2hat + a1b .* q1bhat ;
v0t_2d_3 = a3c .* q3hat + a1c .* q1chat ;

% check that we can reconstruct the velocity 
for j = 1
    % Check the velocity against the vector reconstructed from
    % q2 instead of q1 and q3 instead of q1
    % Build vectors perpendicular to q2hat
    q1b = q1hat - dot(q1hat, q2hat, 2) .* q2hat ;
    q1bhat = q1b ./ sqrt(sum(q1b.^2, 2)) ;
    q3b = q3hat - dot(q3hat, q2hat, 2) .* q2hat ;
    q3bhat = q3b ./ sqrt(sum(q3b.^2, 2)) ;

    % Build vectors perpendicular to q3hat
    q1c = q1hat - dot(q1hat, q3hat, 2) .* q3hat ;
    q1chat = q1c ./ sqrt(sum(q1c.^2, 2)) ;
    q2c = q2hat - dot(q2hat, q3hat, 2) .* q3hat ;
    q2chat = q2c ./ sqrt(sum(q2c.^2, 2)) ;

    % Use dr1_perp,dr2,dr3_perp to form basis of relative direction for v0t
    r1b = r1hat - dot(r1hat, r2hat, 2) .* r2hat ;
    r1bhat = r1b ./ sqrt(sum(r1b.^2, 2)) ;
    r3b = r3hat - dot(r3hat, r2hat, 2) .* r2hat ;
    r3bhat = r3b ./ sqrt(sum(r3b.^2, 2)) ;
    a1b = dot(v0t, r1bhat, 2) ;
    a2b = dot(v0t, r2hat, 2) ;
    a3b = dot(v0t, r3bhat, 2) ;

    % Use dr1_perp,dr2_perp,dr3 to form basis of relative direction for v0t
    r1c = r1hat - dot(r1hat, r3hat, 2) .* r3hat ;
    r1chat = r1c ./ sqrt(sum(r1c.^2, 2)) ;
    r3c = r2hat - dot(r2hat, r3hat, 2) .* r3hat ;
    r3chat = r3c ./ sqrt(sum(r3c.^2, 2)) ;
    a1c = dot(v0t, r1chat, 2) ;
    a2c = dot(v0t, r2hat, 2) ;
    a3c = dot(v0t, r3chat, 2) ;

    close all
    % Plot the velocity in 2d on its containing triangle
    fig = figure;
    triangle = [tria(1, :), tria(1, 1) ] ;
    plot3(m0xy(triangle, 1), m0xy(triangle, 2), 0*m0xy(triangle, 2), '-', 'Color', light_gray)
    hold on
    scatter3(m0xy(tria(1, :), 1), m0xy(tria(1, :), 2), ...
        0*m0xy(tria(1, :), 2), [5,25,50])
    quiver3(x0(j), y0(j), 0, v0t_2d_1(j, 1), v0t_2d_1(j, 2), 0, 0, 'Color', 'k')
    quiver3(x0(j), y0(j), 0, v0t_2d_2(j, 1), v0t_2d_2(j, 2), 0, 0, 'Color', light_gray)
    quiver3(x0(j), y0(j), 0, v0t_2d_3(j, 1), v0t_2d_3(j, 2), 0, 0, 'Color', light_gray)
    v2d_a1 = a1 .* q1hat(j, :) ;
    v2d_a2 = a2 .* q2phat(j, :) ;
    v2d_a3 = a3 .* q3phat(j, :) ;
    quiver3(x0(j), y0(j), 0, v2d_a1(j, 1), v2d_a1(j, 2), 0, 0, 'Color', blue)
    quiver3(x0(j), y0(j), 0, v2d_a2(j, 1), v2d_a2(j, 2), 0, 0, 'Color', orange)
    % quiver3(x0(j), y0(j), 0, v2d_a3(j, 1), v2d_a3(j, 2), 0, 0, 'Color', green)
    % Draw basis vectors
    quiver3(x0(j), y0(j), 0, ...
        q1hat(j, 1), q1hat(j, 2), 0, 0, 'Color', blue, 'Linewidth', 2)
    quiver3(x0(j), y0(j), 0, ...
        q2phat(j, 1), q2phat(j, 2), 0, 0, 'Color', orange, 'Linewidth', 2)
    quiver3(x0(j), y0(j), 0, ...
        q3phat(j, 1), q3phat(j, 2), 0, 0, 'Color', green, 'Linewidth', 2)
    axis equal
    title('2D velocity based on q1')
    print('-dpng','-r500', 'velocity_decomposition2d_a.png')

    % Plot the velocity in 2d on its containing triangle
    figb = figure ;
    triangle = [tria(1, :), tria(1, 1) ] ;
    plot3(m0xy(triangle, 1), m0xy(triangle, 2), 0*m0xy(triangle, 2), '-', 'Color', light_gray)
    hold on
    scatter3(m0xy(tria(1, :), 1), m0xy(tria(1, :), 2), ...
        0*m0xy(tria(1, :), 2), [5,25,50])
    quiver3(x0(j), y0(j), 0, v0t_2d_2(j, 1), v0t_2d_2(j, 2), 0, 0, 'Color', 'k')
    quiver3(x0(j), y0(j), 0, v0t_2d_1(j, 1), v0t_2d_1(j, 2), 0, 0, 'Color', light_gray)
    quiver3(x0(j), y0(j), 0, v0t_2d_3(j, 1), v0t_2d_3(j, 2), 0, 0, 'Color', light_gray)
    v2d_a1 = a1b .* q1bhat(j, :) ;
    v2d_a2 = a2b .* q2hat(j, :) ;
    v2d_a3 = a3b .* q3bhat(j, :) ;
    quiver3(x0(j), y0(j), 0, v2d_a1(j, 1), v2d_a1(j, 2), 0, 0, 'Color', blue)
    quiver3(x0(j), y0(j), 0, v2d_a2(j, 1), v2d_a2(j, 2), 0, 0, 'Color', orange)
    % quiver3(x0(j), y0(j), 0, v2d_a3(j, 1), v2d_a3(j, 2), 0, 0, 'Color', green)
    % Draw basis vectors
    quiver3(x0(j), y0(j), 0, ...
        q1bhat(j, 1), q1bhat(j, 2), 0, 0, 'Color', blue, 'Linewidth', 2)
    quiver3(x0(j), y0(j), 0, ...
        q2hat(j, 1), q2hat(j, 2), 0, 0, 'Color', orange, 'Linewidth', 2)
    quiver3(x0(j), y0(j), 0, ...
        q3bhat(j, 1), q3bhat(j, 2), 0, 0, 'Color', green, 'Linewidth', 2)
    axis equal
    title('2D velocity based on q2')
    print('-dpng','-r500', 'velocity_decomposition2d_b.png')
        
    % Plot the velocity in 2d on its containing triangle
    figc = figure ;
    triangle = [tria(1, :), tria(1, 1) ] ;
    plot3(m0xy(triangle, 1), m0xy(triangle, 2), 0*m0xy(triangle, 2), '-', 'Color', light_gray)
    hold on
    scatter3(m0xy(tria(1, :), 1), m0xy(tria(1, :), 2), ...
        0*m0xy(tria(1, :), 2), [5,25,50])
    quiver3(x0(j), y0(j), 0, v0t_2d_3(j, 1), v0t_2d_3(j, 2), 0, 0, 'Color', 'k')
    quiver3(x0(j), y0(j), 0, v0t_2d_1(j, 1), v0t_2d_1(j, 2), 0, 0, 'Color', light_gray)
    quiver3(x0(j), y0(j), 0, v0t_2d_2(j, 1), v0t_2d_2(j, 2), 0, 0, 'Color', light_gray)
    v2d_a1 = a1c .* q1chat(j, :) ;
    v2d_a2 = a2c .* q2chat(j, :) ;
    v2d_a3 = a3c .* q3hat(j, :) ;
    quiver3(x0(j), y0(j), 0, v2d_a1(j, 1), v2d_a1(j, 2), 0, 0, 'Color', blue)
    % quiver3(x0(j), y0(j), 0, v2d_a2(j, 1), v2d_a2(j, 2), 0, 0, 'Color', orange)
    quiver3(x0(j), y0(j), 0, v2d_a3(j, 1), v2d_a3(j, 2), 0, 0, 'Color', green)
    % Draw basis vectors
    quiver3(x0(j), y0(j), 0, ...
        q1chat(j, 1), q1chat(j, 2), 0, 0, 'Color', blue, 'Linewidth', 2)
    quiver3(x0(j), y0(j), 0, ...
        q2chat(j, 1), q2chat(j, 2), 0, 0, 'Color', orange, 'Linewidth', 2)
    quiver3(x0(j), y0(j), 0, ...
        q3hat(j, 1), q3hat(j, 2), 0, 0, 'Color', green, 'Linewidth', 2)
    axis equal
    title('2D velocity based on q3')
    print('-dpng','-r500', 'velocity_decomposition2d_c.png')


    % Plot the velocity in 3d on its containing triangle
    fig2 = figure;
    meanz = mean(z123ap(1, :)) ;
    hold on
    trngl = [x123ap(1, :)', y123ap(1, :)', z123ap(1, :)'] ;
    trngl4 = [trngl ; trngl(1, :) ] ;
    plot3(trngl4(:, 1), trngl4(:, 2), trngl4(:, 3), '-', 'Color', light_gray)
    scatter3(trngl(:, 1), trngl(:, 2), trngl(:, 3), [5,25,50])
    quiver3(pt0(j, 1), pt0(j, 2), pt0(j, 3), ...
        v0(j, 1), v0(j, 2), v0(j, 3), 0, 'Color', 'k')
    quiver3(pt0(j, 1), pt0(j, 2), pt0(j, 3), ...
        v0t(j, 1), v0t(j, 2), v0t(j, 3), 0, 'Color', 'k')
    quiver3(pt0(j, 1) + v0t(j, 1), ...
        pt0(j, 2) + v0t(j, 2),...
        pt0(j, 3) + v0t(j, 3), ...
        v0n(j) * normals(j, 1), ...
        v0n(j) * normals(j, 2), ...
        v0n(j) * normals(j, 3), 0, 'Color', 'k')
    % draw r1hat, r2hat, r3hat
    quiver3(pt0(j, 1), pt0(j, 2), pt0(j, 3), ...
        r1hat(j, 1), r1hat(j, 2), r1hat(j, 3), 0, 'Color', blue, 'Linewidth', 2)
    quiver3(pt0(j, 1), pt0(j, 2), pt0(j, 3), ...
        r2phat(j, 1), r2phat(j, 2), r2phat(j, 3), 0, 'Color', orange, 'Linewidth', 2)
    quiver3(pt0(j, 1), pt0(j, 2), pt0(j, 3), ...
        r3phat(j, 1), r3phat(j, 2), r3phat(j, 3), 0, 'Color', green, 'LineWidth', 2)
    % Draw construction of v0 from bases
    ar1 = a1(j) * r1hat(j, :) ;
    ar2 = a2(j) * r2phat(j, :) ;
    ar3 = a3(j) * r3phat(j, :) ;
    assert(all(ar2 - ar3 < 1e-7))
    quiver3(pt0(j, 1), pt0(j, 2), pt0(j, 3), ...
        ar1(1), ar1(2), ar1(3), 0, 'Color', sky) ;
    quiver3(pt0(j, 1), pt0(j, 2), pt0(j, 3), ...
        ar2(1), ar2(2), ar2(3), 0, 'Color', yellow) ;
    quiver3(pt0(j, 1), pt0(j, 2), pt0(j, 3), ...
        ar3(1), ar3(2), ar3(3), 0, 'Color', light_green) ;
    axis equal
    xlabel('x')
    ylabel('y')
    zlabel('z')
    title('Velocity decomposition in face')
    print('-dpng','-r500', 'velocity_decomposition3d.png')
end

