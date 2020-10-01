function [matrix, ia, ib] = uniqueRows(matrix)
% Return a matrix of just the unique rows
% ia : out = in(ia) 
% ib : in = out(ib) 

[matrix,ia,ib]=unique(matrix,'rows','stable') ;
% i=true(size(matrix,1),1);
% i(ia)=false;
% matrix(i,:)=[];