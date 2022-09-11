%Noah Mitchell
%change (batch) file names

%make into a function...

directory= './';
fnSearchStr = 'Time*.ome.tif' ; %'*stripe7*.png' ;

existingStrReplace= false;
knownLocReplace = false ;
renameTimeStamps = false ;
divideTimeStampsIntoChannels = true ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if existingStrReplace
    %%%%%%%%%%%%%%%%%%%
    %Replace characters at position of preexisting string in name
    oldchar = '.png' ;
    newchar = '_dist.png' ;
    toRename = dir(fullfile(directory, fnSearchStr));
    nchar = strfind(toRename(1).name, oldchar) ;

    lenstr= length(toRename(1).name);
    fprintf('First new name is:\n  %s\n',[toRename(1).name(1:nchar-1), newchar, toRename(1).name(nchar+length(oldchar):lenstr)])

    for ii=1:length(toRename)
        lenstr= length(toRename(ii).name);
        newname = [toRename(ii).name(1:nchar-1), newchar, toRename(ii).name(nchar+length(oldchar):lenstr)];
        movefile([directory, toRename(ii).name], [directory, newname])
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Replace characters at known position in string
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if knownLocReplace
    %Replace characters at known position in string
    newchar = '.tif' ;
    toRename = dir([directory '/20*']);
    nchar = length(toRename(1).name) ;
    fprintf('First new name is:\n  %s\n',[toRename(ii).name(1:nchar-1), newchar])

    for ii=1:length(toRename)
        lenstr= length(toRename(ii).name);
        newname = [toRename(ii).name(1:nchar-1), newchar] ;
        movefile([directory, toRename(ii).name], [directory, newname])
    end
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
directory = './' ;
fnSearchStr = 'Time*ome.tif';
fns = dir(fullfile(directory, fnSearchStr)) ;

addTime = 1 ;
prepend = 'Time_' ;
postpend = '_Angle_' ;
if renameTimeStamps
    for ii = 1:length(fns)
        fn2rename = fullfile(directory, fns(ii).name) ;
        fn = fns(ii).name ;
        end0 = strfind(fn, postpend) ;
        start0 = strfind(fn, prepend) ;
        tstamp = fn(start0+5:end0-1) ;
        disp(['t=' tstamp])

        % Convert to number, adjust 
        tstampNew = sprintf('1%06d', str2double(tstamp)+addTime) ; 
        newname = [prepend tstampNew fn(end0:end)] ;

        infn = fullfile(directory, fn2rename) ;
        outfn = fullfile(directory, newname) ;
        disp([infn ' -> ' outfn])
        movefile(infn, outfn) ;
    end
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert t=0 to t=0 (ch1), t=1 to t=0 (ch2), t=2 to t=1 (ch1)...
% This is for use after repacking from Vishank with 2 channels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
directory = './' ;
fnSearchStr = 'Time*ome.tif';
fns = dir(fullfile(directory, fnSearchStr)) ;

divTime = 2 ;
prepend = 'Time_' ;
postpend = '_Angle_' ;
if divideTimeStampsIntoChannels
    for ii = 1:length(fns)
        fn2rename = fullfile(directory, fns(ii).name) ;
        fn = fns(ii).name ;
        end0 = strfind(fn, postpend) ;
        start0 = strfind(fn, prepend) ;
        tstamp = fn(start0+5:end0-1) ;
        disp(['t=' tstamp])

        % Convert to number, adjust 
        if mod(str2double(tstamp) / divTime, 1) == 0
            tstampNew = sprintf('%06d', str2double(tstamp)/divTime) ; 
            newname = [prepend tstampNew fn(end0:end)] ;
        else
            tstampNew = sprintf('%06d', floor(str2double(tstamp)/divTime)) ; 
            newname = [prepend tstampNew fn(end0:end)] ;
            newname = strrep(newname, 'c1', 'c2') ;
        end
        
        infn = fullfile(directory, fn2rename) ;
        outfn = fullfile(directory, newname) ;
        disp([infn ' -> ' outfn])
        movefile(infn, outfn) ;
    end
end
