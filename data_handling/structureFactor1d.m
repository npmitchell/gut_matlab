function [S,q] = structureFactor1d(L,scatrpos,numBZ)
% STRUCTUREFACTOR1D 
% http://en.wikipedia.org/wiki/Structure_factor
% The structure factor (written in Latex) is defined as
%  S(\vec{q})=\frac{1}{N} \left| \sum_{\ell}^N \exp 
%  \left( i \vec{q}\cdot \vec{R}_{\ell} \right)\right|^2
% 
% inputs are real-space system length (L),
% real-space scatterer positions on integer lattice (scatrpos),
% and number of Brillion zones to plot (numBZ)
% Example use with fully-occupied periodic array:
%     First specify the dimensions,
% L=100; 
%     then specify location of scatterers
% scatrpos=1:L;
%     finally, pass these to this function 
% [S,q] = structureFactor1d(L,scatrpos,1);
%     which returns the structure factor
%     and assoicated momentum-space coordinates
% 
% Example of alternatingly-occupied lattice
% L=100; alt=3; scatrpos=1:alt:L;
% 
% Example of randomly-occupied periodic lattice with filling fraction=0.8, 
% L=10; ff=0.8; scatrpos=rand(L,1); tmp=scatrpos;
% scatrpos=(scatrpos>ff)'.*[1:L]; scatrpos=scatrpos(scatrpos~=0);
% 

x=0;
S=zeros(size(numBZ*L:numBZ*L,1),1);
q=zeros(size(numBZ*L:numBZ*L,1),1);
for n=-numBZ*L:numBZ*L
  x=x+1;
  S(x)=(1/L)*abs(sum(exp(1i*2*pi*n*scatrpos/L)))^2;
  q(x)=2*pi*n/L;
end
plot(q/(2*pi),S/L,'Marker','.','LineWidth',0.5); 
xlabel('q/(2*\pi)'); ylabel('S(q)/L');




