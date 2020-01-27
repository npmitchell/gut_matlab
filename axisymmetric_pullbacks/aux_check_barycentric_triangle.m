
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check that baryc gives weights of the three vertices of triangle
% t_contain(i) for pointLocation x0(i), y0(i)
if preview
    close all
    fig = figure('Visible', 'On')
    % triplot(tr0, 'Color', blue, 'LineWidth', 0.0001) 
    hold on
    for j = 1:3
        tritest = tr0.ConnectivityList(fieldfaces(j), :); 
        btest = baryc0(j, :);
        triangle = [tritest, tritest(1) ] ;
        plot(m0xy(triangle, 1), m0xy(triangle, 2), 'g.-')
        plot(x0(j), y0(j), 'o')
        xfind = sum(btest' .* m0xy(tritest, 1)) ;
        yfind = sum(btest' .* m0xy(tritest, 2)) ;
        plot(xfind, yfind, '^', 'MarkerSize', 10)        
        axis equal
    end

    % Draw other connections
    quiver(x0, y0, uu, vv, 0)
    plot(x1, y1, 'ko')
    for j = 1:3
        tritest = tr1.ConnectivityList(t1_contain(j), :); 
        btest = baryc1(j, :);
        triangle = [tritest, tritest(1) ] ;
        plot(m1xy(triangle, 1), m1xy(triangle, 2), 'r.-')
        plot(x1(j), y1(j), 'kx')
        xfind = sum(btest' .* m1xy(tritest, 1)) ;
        yfind = sum(btest' .* m1xy(tritest, 2)) ;
        plot(xfind, yfind, '^', 'MarkerSize', 10)        
        axis equal
    end

    waitfor(fig)
    close all        
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%