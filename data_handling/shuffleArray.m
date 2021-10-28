function arr = shuffleArray(arr)
%arr = shuffleArray(arr)
% Shuffle the indices or rows of an array randomly 
%
%

if length(size(arr)) == 2 && any(size(arr)==1)
    % shuiffle the elements of the array
    arr = arr(randperm(length(v))) ;
else
    % shuffle the rows of the array
    arr(1:size(arr, 1), :) = arr(randperm(size(arr, 1)), :) ;
end