
% First run setup.m from imsane directory
addpath('./addpath_recurse/')
addpath('./mesh_handling/')
datadir = '/Users/npmitchell/Dropbox/Soft_Matter/UCSB/gut_morphogenesis/' ;
datadir = [datadir 'data/48Ygal4UasCAAXmCherry/48Ygal4UasCAAXmCherry_20190207200_excellent/'] ;
meshdir = [datadir 'msls_output_prnun5_prs1_nu0p00_s0p10_pn2_ps4_l1_l1/'] ;
meshes = dir([meshdir 'mesh_apical_0*.ply']) ;

for i=1:length(meshes)
    meshfn = meshes(i).name ; 
    mesh = read_ply_mod(meshfn) ;
    
    % Get alignment from AP axis determination
    
    
    
end
    
    
    
% % Fit to a polynomial
% % From points extracted from image of gut, extract polynomial fit 
% % Choose order of polynomial
% order = 4 ;
% 
% % Subtract the mean of the independent axis (axial position of gut)
% xx = xx - mean(xx) ;
% 
% % Fit the points to an nth order polynomial 
% [pp, SS] = polyfit(xx, yy, order) ;
% 
% % Evaluate the fit to compare
% xfit = linspace(min(xx), max(xx));
% yfit = polyval(pp, xfit);
% 
% % Plot the result
% plot(xx, yy, '.')
% hold('on')
% plot(xfit, yfit, '-') 
% xlabel('position along gut [um]')
% ylabel('gut surface')

