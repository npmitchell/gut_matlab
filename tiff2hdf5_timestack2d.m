% TIFF2HDF5_TIMESTACK2D
% Convert a stack of tiffs to a 3D dataset in hdf5 format

imdir = './' ;
outname = 'conformal_9_timestack.hdf5' ;
imname = 'cmp_1_1_T*.tif' ;

% Obtain the image names
ims = dir(fullfile(imdir, imname)) ;
disp('reading imagestack...')
for i=1:length(ims)
    im = imread(fullfile(ims(i).folder, ims(i).name)) ; 
    if i == 1
       % preallocate size of output dataset
       dset = zeros([length(ims) size(im)]) ;
    end
    dset(i, :, :) = im ;
end

disp('creating hdf5 output...')
if exist(outname, 'file')
    fid = H5F.open(outname);
    fapl = H5F.get_access_plist(fid);
    if H5P.exist(fapl, 'data')
        fprintf('data property exists\n');
        create_fresh = false ;
    else
        fprintf('data property does not exist\n');
        create_fresh = true ;
    end
end

% Add the dataset if the 'data' dataset is not already existent
if create_fresh
    disp('writing to data in hdf5...')
    h5create(outname, '/data', size(dset));
    h5write(outname, '/data', dset);
end

attr_struct = [struct("key", "t", "typeFlags", 8, "resolution", 0, "description", "" ), ...
    struct( "key", "x", "typeFlags", 2, "resolution", 0, "description", "" ), ...
    struct( "key", "y", "typeFlags", 2, "resolution", 0, "description", "" ), ...
    struct( "key", "z", "typeFlags", 2, "resolution", 0, "description", "" ), ...
    struct( "key", "c", "typeFlags", 1, "resolution", 0, "description", "" )] ;
h5writeatt(outname, '/', 'axes', attr_struct) ; 

