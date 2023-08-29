function datXYT = followZPlane(datZXYT, tz, nZplanes, fnout)
% datXYT = followZPlane(dataZXYT, tz, nZplanes, fnout)
%
% Follow the dynamics of a 4D dataset by looking at a particular Z plane at
% at each timepoint. Here we specify which page to look at for each
% timepoint in the argument tz.
% 
% Parameters
% ----------
% dataZXYT : 1 x nZplanes cell array of X x Y x T arrays or string path
%   the data to sample at a given z plane for each timepoint, or the path
%   to the TIFF where this data is stored
% tz : N x 2 array of timepoints and zplanes of interest
%   A set of some timepoints and their associated z planes. If timepoints
%   are not specified, we use a nearest interpolation scheme to identify
%   which z plane would be of interest
% nZplanes : int
%   the number of z planes in the dataset. This is only used if dataZXYT is
%   a string path to the data rather than the data itself. 
% fnout : string path
%   Where to write to disk the output XYT data, where each timepoint is a
%   2D image taken from a particular z plane of the input 4D data.
%
% Returns
% -------
% datXYT : X x Y x T array 
%   XYT image data, where each page datXYT(:, :, n) is a 2D image taken 
%   from a particular z plane of the input 4D data at the nth timepoint.
%
% Example usage:
% tz = [1 20;
% 5 19;
% 8  17;
% 10 16;
% 13 15;
% 18 14;
% 20 13;
% 23 12;
% 26 11;
% 27 10;
% 35 10;
% 38 9;
% 56 9] ;
% nZplanes = 31 ;
% datdir = '/mnt/data/confocal_data/somatic_muscle/202306090954_mef2G4kCAAXmCh_2um8mpf_6t10pc561_GBR/';
% datZXYT = fullfile(datdir, '202306090954_mef2G4kCAAXmCh_2um8mpf_e2.tif') ;
% fnout = fullfile(datdir, '202306090954_mef2G4kCAAXmCh_2um8mpf_e2_followZ.tif') ;
% datXYT = followZPlane(dataZXYT, tz, nZplanes, fnout) ;

% If we passed datZXYT as a string, load the tiff file from that string
if ischar(datZXYT)
    disp('Reading TIFF from filename...')
    datZXYT = readTiff4D(datZXYT, nZplanes) ;
end
ntp = size(datZXYT{1}, 3) ;

% preallocate the data array to be filled. Note we choose uint16 here
datXYT = zeros(size(datZXYT{1}, 1), size(datZXYT{1}, 2), ntp, 'uint16') ;

% Interpolate the lookup table translating time into a Zplane of interest
TZ = interp1(tz(:, 1), tz(:, 2), 1:ntp, "linear") ;

% Fill in the data, following the z plane of interest over time
disp('stuffing datXYT...')
for ii = 1:ntp
    pageID = round(TZ(ii)) ;
    page = datZXYT{pageID}(:, :, ii) ;
    datXYT(:, :, ii) = page ;
end

% Write the single-pane XYT movie to disk as a TIFF
if nargin > 3
    datXYT = reshape(datXYT, [size(datXYT, 1), size(datXYT, 2), 1, 1, ntp]) ;
    % info = imfinfo(dataTXYZ{1});
    % info.BitDepth
    writeTiff5D(datXYT, fnout)
end

return