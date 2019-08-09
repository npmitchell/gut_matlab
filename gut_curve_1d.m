
% From points extracted from image of gut, extract polynomial fit 
% Choose order of polynomial
order = 4 ;

% Suppose you have a piece of the boundary of a 2D slice of the gut stored 
% as(xx, yy)
yy = [1, 0.9, 0.9, 0.95, 0.6, 0.2, 0.45, 0.7, 0.9, 1.0] ;
xx = linspace(0, 9, length(yy)) ;

% Subtract the mean of the independent axis (axial position of gut)
xx = xx - mean(xx) ;

% Fit the points to an nth order polynomial 
[pp, SS] = polyfit(xx, yy, order) ;

% Evaluate the fit to compare
xfit = linspace(min(xx), max(xx));
yfit = polyval(pp, xfit);

% Plot the result
plot(xx, yy, '.')
hold('on')
plot(xfit, yfit, '-') 
xlabel('position along gut [um]')
ylabel('gut surface')

