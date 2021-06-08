function arr = deleteRowColumn(arr, i)
% Delete the row and column from @arr that has the @i_th diagonal element.
assert(all(size(arr) >= i))
arr(i, :) = [];
arr(:, i) = [];
end