%% T1 transition example videos
% Use for figures / SI / videos
Options = struct() ;
QS.visualizeDemoTracks(Options) 


%% Notes: 
param = struct();
param.borderType = 2; % FIXED = 1, FREE = 2
param.fixedType = 1; % ARC_LENGTH = 1, UNIFORM = 2
param.fixedShape = 1; % CIRCLE = 1, SQUARE = 2

% FIXED: BARYCENTRIC = 1, AUTHALIC = 2, CONFORMAL = 3, MEAN = 4
% FREE: LSCM = 1, ARAP = 2
param.paramMethod =  2;

% param.corners = [];
% param.fixedPoints = [];

fprintf('Generating surface parameterization... ');
[ V2D, ~ ] = surface_parameterization( subm.f, subm.v, param );

if param.borderType == 1
    if param.fixedShape == 1
        V2D = 2 .* ( V2D - 0.5 ); % Map disk to unit disk
    end
end

% V2D = conformalParametrization( subm );
fprintf('Done\n');