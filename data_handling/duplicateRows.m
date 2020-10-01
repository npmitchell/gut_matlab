function duprows = duplicateRows(matrix)
% Return a matrix of just the duplicate rows
[C,ia,ib]=unique(matrix,'rows','stable') ;
i=true(size(matrix,1),1);
i(ia)=false;
duprows = matrix(i,:);