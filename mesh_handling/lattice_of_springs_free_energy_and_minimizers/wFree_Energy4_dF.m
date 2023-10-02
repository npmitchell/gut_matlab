% [F,dF] = wFree_Energy4_dF(xyz,NL,bL0matrix,KL,BLm,KLm)
%
% Parameters
% ----------
% xyz : NxD float array
%   vertex positions
% NL : Nx max#NN int array
%   neighbor list, with zeros for non-neighbors
% bL0matrix : N x max#NN float array
%   rest bond lengths
% KL : #bonds x 1 float array
%   spring constant array
% bL0 : #bonds x 1 float array
%   rest bond lengths for each bond
% kL : #bonds x 1 float array
%   spring constant per bond
% 
% Returns
% -------
% F : free energy
% dF : current forces on each vertex
%
% WTMI, Noah Mitchell 2010-2015
% Calculate the Free energy F and force on each particle dF 
% Note: NL2 is a version of NL with zeros replaced by ones (taken care of by
% KL. In other words, we create artificial coupling of all nodes except 
% node 1 with node 1.
% vk made a slight modification of NL2 (9 Oct 12)
% npm generalized to accept networks with different bond lengths

function [F,dF] = wFree_Energy4_dF(xyz,NL,BL, bL0matrix, KL,bL0,kL)

% Note: to make bLmatrix, do
% bLmatrix = bL2bLmatrix(NL,BL, bL)


temp = (xyz(BL(:,1),:) - xyz(BL(:,2),:));
F = 1/2 * sum(kL .* (sqrt(dot(temp,temp,2)) - bL0).^2); %calc Free energy

% copy NL into NL2
NL2 = NL;
zeros_first_row = NL2(1,:)==0; %vk added to find the column numbers of the first row 
NL2(1, zeros_first_row)=2; %vk added to create artificial coupling of node 1 with node 2 (and not node 1!)
NL2(NL2==0) = 1; %vk create artificial coupling of all nodes except node 1 with node 1


%Each row of dF is the net force on each particle in x,y,z.
%Preallocate dF for speed.
dF = zeros(size(xyz));

%Summing over nearest neighbors, sum the force on each particle
for kk=1:size(NL2,2)
    vec=(xyz-xyz(NL2(:,kk),:)); %Cartesian vectors of NNeighbors
    
    %dis is the scalar distance between each particle and its kk neighbor.
    dis = sqrt(dot(vec,vec,2));
    
    %Debugging
    %fprintf(['size(dis)=',num2str(length(dis)),'\n']) 
    %fprintf(['size(BL(:,kk))=',num2str(length(BL(:,kk))),'\n']) 
    
    %The Lagrangian strain is: (dis-BL(:,kk))./dis  
    % dF is K*displacement*normalvector, where vec/dis is the normal vector
    if size(xyz, 2) == 3
        dF = dF+((KL(:,kk).*(dis-bL0matrix(:,kk))./dis)*[1,1,1]).*vec;  %
    elseif size(xyz, 2) == 2
        dF = dF+((KL(:,kk).*(dis-bL0matrix(:,kk))./dis)*[1,1]).*vec;  %
    end
    
end
