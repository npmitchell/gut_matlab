%Noah Mitchell
%change (batch) file names

%make into a function...

directory=['/Volumes/labshared3-1/noah/2014_06_geoFLAT_SphPatch/2014_06_30_t2_GeoFLATacrylic_05pcplusGeo_4o12pi_L9p4in0p0a_gtb_20pps/imagestack/'];

existingStrReplace=1;
if existingStrReplace==1
    %%%%%%%%%%%%%%%%%%%
    %Replace characters at position of preexisting string in name
    oldchar = 'acrylic_' ;
    newchar = 'acrylic_SphP_' ;
    %oldchar = '2013_12_22_t2_12pc_3DPB_1p0inoff_1sig_coilg_pd_fb_15pps_2013-12-22-154355-'
    %newchar = '2013_12_22_t3_12pc_3DPB_0p8inoff_1sig_coilg_pd_fb_15pps_2013-12-22-154355-'
    toRename = dir([directory '/20*']);
    nchar = strfind(toRename(1).name, oldchar) ;
    
    lenstr= length(toRename(1).name);
    fprintf('First new name is:\n  %s\n',[toRename(1).name(1:nchar-1), newchar, toRename(1).name(nchar+length(oldchar):lenstr)])
    
    for ii=1:length(toRename)
       lenstr= length(toRename(ii).name);
       newname = [toRename(ii).name(1:nchar-1), newchar, toRename(ii).name(nchar+length(oldchar):lenstr)];
       movefile([directory, toRename(ii).name], [directory, newname]) 
    end
end

knownLocReplace=0;
if knownLocReplace==1
    %%%%%%%%%%%%%%%%%%%
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