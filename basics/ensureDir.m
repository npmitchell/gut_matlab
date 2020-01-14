function ensureDir(dirpath)
%ENSUREDIR(dirpath) create a directory if it does not exist already
%
% Parameters
% ----------
% dirpath : str
%   The path to create if it does not exist
%
% 
% NPMitchell 2019

if ~exist(dirpath, 'dir')
    mkdir(dirpath)
end