% wtmi 2010 modified:
% - aug 2010 (wtmi)
% - june 2011 (wtmi) - varargin is size of NL
% Noah Mitchell commented 2014
%
%Take a 2 column list of particles (identified by their row in some xyz
%list) and give back a neighbor list (NL), in which, for example,
%if particle 1 is bonded to particle 2 3 5 and 9, then the first row of NL
%is
%2 3 5 9 0 0 0 0, 
%where the number of columns (with zeros filling the
%extras) are given by the argument varargin.

function [NL] = TRISm2NL(TRISm,varargin)

if size(varargin,2)==1
    n=varargin{1};
else
    n=6;
end

NL = [];
x = max(unique(TRISm));
for kk=1:x
    temp = or(TRISm(:,1)==kk,TRISm(:,2)==kk);
    temp = TRISm(temp,:);
    temp = unique(temp(:));
    temp = setdiff(temp,kk)';
    if size(temp,1)>size(temp,2)
        temp = temp';
    end
    if ~isempty(temp)
        temp =  cat(2,temp,zeros(1,(n-size(temp,2)))) ; %debugg nm
        NL = cat(1,NL,temp);
    else
        temp = zeros(1,n);
        NL = cat(1,NL,temp);
    end
    %temp = temp(1:n);
    
    
end