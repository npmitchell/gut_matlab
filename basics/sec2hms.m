function hms = sec2hms(t, format)
% consider format = '%02d:%02d:%02.f'
% https://stackoverflow.com/questions/12210583/is-there-a-matlab-function-to-convert-elapsed-seconds-to-hhmmss-format/12211011
    if nargin < 2
	format = '%02d:%02d:%02d';
    end
    hours = floor(t / 3600);
    t = t - hours * 3600;
    mins = floor(t / 60);
    secs = t - mins * 60;
    hms = sprintf(format, hours, mins, secs);
end
