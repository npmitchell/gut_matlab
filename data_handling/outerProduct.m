function C = outerProduct(A, B)  
% Outer product of two matrices A and B

C = A .* reshape(B, [ones(1, ndims(A)), size(B)]);

% Matlab < R2016b:
% C = bsxfun(@times, A, reshape(B, [ones(1, ndims(A)), size(B)]))

end