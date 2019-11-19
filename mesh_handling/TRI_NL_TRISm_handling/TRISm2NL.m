function [NL] = TRISm2NL(TRISm,varargin)
% wtmi 2010 npmitchell 2013, 2019:
% - aug 2010 (wtmi)
% - june 2011 (wtmi) - varargin is size of NL
% - 2013, 2019 (npmitchell) automatically size NL based on neede rows
% 

if size(varargin,2)==1
    n=varargin{1};
else
    n=15 ;
end

NL = [];
x = max(unique(TRISm));
for kk=1:x
    temp = find(or(TRISm(:,1)==kk,TRISm(:,2)==kk));
    temp = TRISm(temp,:);
    temp = unique(temp(:));
    temp = setdiff(temp,kk)';
    if size(temp,1)>size(temp,2)
        temp = temp';
    end
    if ~isempty(temp);    
        temp =  cat(2,temp,zeros(1,(n-size(temp,2))));
        NL = cat(1,NL,temp);
    else
        temp = zeros(1,n);
        NL = cat(1,NL,temp);
    end
    %temp = temp(1:n);    
end

% Now truncate the unused columns is varargin is not supplied
if size(varargin,2)~=1
    NL = NL(:, 1:find(nonzeros(sum(NL, 1)), 1, 'last' )) ;
end

