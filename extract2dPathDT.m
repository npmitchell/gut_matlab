function path = extract2dPathDT(DT, endpt, options)
%
% Find a streamline path in a 2d gradient descent image WW
% Very similar to Gabriel Peyre's compute_geodesic.m
% 
% Npmitchell 2020, based on Gabriel Peyre's code

str_options = options.str_options ;
trim_path = options.trim_path ;
thres_dist = 0 ; 
if isfield(options, 'thres_dist')
    thres_dist = options.thres_dist ;
end

[Dx, Dy] = gradient(DT) ;
Dx = - Dx ;
Dy = - Dy ;

% path extraction
grad = cat(3, Dx, Dy) ;
grad = perform_vf_normalization(grad) ;
Dx = squeeze(grad(:, :, 1)) ;
Dy = squeeze(grad(:, :, 2)) ;
path = stream2(Dx, Dy, endpt(1),endpt(2), str_options);
path = path{1} ;

% Compare to continued options in compute_geodesic.extract_path_2d()
% global grad;
% grad = compute_grad(DT);
% grad = -perform_vf_normalization(grad);
% disp('compute_geodesic.extract_path_2d(): path = ')
% path
% 
% for i=1:length(path)
%     path{i} = path{i}(:,2:-1:1);
% end
% if length(path)==1
%     path = path{1};
% end

if isfield(options, 'start_points')
    start_points = options.startpt ;
else
    start_points = path(end,:);
end
start_points = start_points(:);

if trim_path
    % removing too verbose points
    d = compute_distance_to_points(path', start_points);
    % perform thresholding
    if thres_dist == 0
        T = mmax(d)/300^2 ;
    else
        disp(['Thresholding with ' num2str(thres_dist)])
        plot(d)
        ylim([0, 100])
        T = thres_dist ;
    end
    
    I = find(d<T);
    if not(isempty(I))
        path = path(1:I(1), :);
        path = [path; start_points'];
    else
        path = path';
    end
end

% % complete with a discrete extraction (nasty hack)
% if size(path, 2)~=2 && size(path, 1)==2
%     path = path';
% end
% 
% disp('compute_geodesic.extract_path_2d(): end_points= ')
% end_points
% disp('compute_geodesic.extract_path_2d(): path= ')
% path
% imshow(A')
% path = [path; compute_discrete_geodesic(A, round(path(end,:)))'];