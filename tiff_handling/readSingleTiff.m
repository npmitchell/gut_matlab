% fileName = sprintf(xp.fileMeta.filenameFormat, xp.currentTime);
% fullFileName = fullfile(xp.fileMeta.dataDir, fileName);

function data = readSingleTiff(fullFileName)
% READSINGLETIFF read a single TIFF file/stack

r = bfGetReader(fullFileName);

r.setSeries(0);

depth = r.getBitsPerPixel;

stackSize = [r.getSizeX(), r.getSizeY(), r.getSizeZ(), r.getSizeT()];

xSize = stackSize(1);
ySize = stackSize(2);
zSize = stackSize(3);

% number of channels
nChannels = r.getSizeC();

nChannelsUsed = 1;

nTimePts = stackSize(4);

% load data in 32 bit

data = zeros([ySize xSize zSize nChannelsUsed], 'single');
            
            for i = 1:r.getImageCount()

                ZCTidx = r.getZCTCoords(i-1) + 1;
                
                % in the fused embryo data coming out of the python script,
                % Z and T are swaped. In general this isn't the case, thus
                % introduce a file metaField swapZT
                
                    zidx = ZCTidx(1);
                    tidx = ZCTidx(3);
             
                cidx = ZCTidx(2);

                % see above: if there is only one timepoint all the planes
                % should be read, if there are multiple timepoints, only
                % the correct time shouldbe read
                if nTimePts == 1 
                    
                    fprintf('.');
                    if rem(i,80) == 0
                        fprintf('\n');
                    end

                    dataCidx = find(1);
                    if ~isempty(dataCidx)
                        
                        data(:,:, zidx, dataCidx) = bfGetPlane(r, i);
                    end
                end
            end
end