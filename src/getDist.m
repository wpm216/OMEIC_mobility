function d = getDist(i, j, n, unitCellSize, nImages, nnVecArrays)
% Gets the distance between two indices i and j in the distance array.

% i   = row index of pair whose distance is being calculated
% j   = column index of pair whose distance is being calculated
% n   = number of molecules in the unit cell
% boxSize  = a 3-vector of the unit cell dimensions 
%            in the x, y, and z directions.
% boxCells = a 3-vector of the number of unit cell images
%            in the x, y, and z directions of the simulation box.

% Identify which cell each index belongs to
% unitCellDims = size(u);
% n = unitCellDims(1);
celli = getCell(i, n, nImages);
cellj = getCell(j, n, nImages);

% Identify which index within the unit cell i and j correspond to
% (matlab is origin-1, so we add one to each index)
idxi = mod(i-1, n) + 1;
idxj = mod(j-1, n) + 1;

% Determine which nearest-neighbor unit cell corresponds
% to the (celli, cellj) pair.
celldiff = cellj - celli;
% If moving from celli to cellj requires a translation
% in any direction, the nearest-neighbor cell reflects that.
% nearestNeighbor = (celldiff ~= 0) .* sign(celldiff);
nearestNeighbor = sign(celldiff);
translations = celldiff - nearestNeighbor;

assert(all((celli + translations + nearestNeighbor) == cellj), "messed up cells")

% Calculate distance, which is equal to the nearest-neighbor distance
% plus any additional translations needed.
nnIdx = 1 + sum((nearestNeighbor + 1) .* [1 3 9]);
nnVec = reshape(nnVecArrays(idxi, idxj, :, nnIdx), [1, 3]);
d = norm(nnVec + translations .* unitCellSize);

end