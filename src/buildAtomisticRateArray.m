function [K, fullDistArr] = buildAtomisticRateArray(unitSize, fullSize, unitPositions, nnVecArrays, ...
                            betaValue, islandRadius, leftIslandLoc, rightIslandLoc, pdbname)

% Generates the array of inter-particle tunneling rates
% for a pss-rich matrix populated by rod-like particles.

% Parameters:
% -----------
% unitSize = box size of unit cell (Angstroms)
% fullSize = size of the full cell to be built (in Angstroms)
% unitPositions = positions of all clusters in the unit cell
% nnVecArrays = n x n x 3 x 27 array of pairwise cluster nn distance vectors
% betaValue = hopping decay coefficient. Alessandro set it to 0.6.
% islandRadius = radius of islands (Angstroms)
% leftIslandLoc = x,y,z coordinates of left island 
% rightIslandLoc = x,y,z coordinates of right island 
% pdbname = name of pdb file to write for this system. 
%           if bool(pdbname) == 0, then no pdb is written.

%                                              %
% Calculate all of the inter-cluster distances % 
%                                              % 

% number of images needed to fill the box, rounding up
nImages = ceil(fullSize./unitSize);

% initialize the distance array
nClustersPerImage = size(nnVecArrays, 1);
nClustersTotal = prod(nImages) * nClustersPerImage;
clusterDistArr = zeros(nClustersTotal, nClustersTotal); % we'll add the islands at the end

% we want to keep track of what box each cluster is in 
% and its index within the unit cell
clusterBoxIndices = zeros(nClustersTotal, 3); 
clusterUnitCellIndex = zeros(nClustersTotal, 1);

% populate the upper diagonal half of the distance array
% element-by-element
for i = 1:nClustersTotal
    clusterBoxIndices(i, :) = getCell(i, nClustersPerImage, nImages);
    clusterUnitCellIndex(i) = mod(i-1, nClustersPerImage) + 1;
    for j = i+1:nClustersTotal
        clusterDistArr(i, j) = getDist(i, j, nClustersPerImage, unitSize, nImages, nnVecArrays);
        clusterDistArr(j, i) = clusterDistArr(i, j); % it's a symmetric matrix
    end    
end


%                                                             %
% Get positions and cluster-island distances for each cluster % 
%                                                             % 

clusterPositions = cell(nClustersTotal, 1);
someIslandOverlapIdxs = [];
fullIslandOverlapIdxs = [];
leftIslandMinDists = zeros(nClustersTotal, 1);
rightIslandMinDists = zeros(nClustersTotal, 1);
for i = 1:nClustersTotal
    
    % Get the positions of this cluster
    idx = clusterUnitCellIndex(i);
    pos = unitPositions(:, :, idx);
    zeroRows = all(pos == 0, 2); % remove padding of dummy [0,0,0] coordinates
    pos = pos(~zeroRows, :);
    shift = (clusterBoxIndices(i, :)-1) .* unitSize; % translate into correct box
    clusterPositions{i} = pos + repmat(shift, size(pos, 1), 1);
    
    % Get minimum island-cluster distance
    leftIslandCDists = sqrt(sum((leftIslandLoc - clusterPositions{i}).^2, 2));
    rightIslandCDists = sqrt(sum((rightIslandLoc - clusterPositions{i}).^2, 2));
    leftIslandMinDists(i) = min(leftIslandCDists) - islandRadius;
    rightIslandMinDists(i) = min(rightIslandCDists) - islandRadius;
    
    % Figure out if it overlaps with either island
    % First, check for partial overlap. These won't be deleted, but their
    % island distance will be set to zero.
    leftIslandOverlap = any(leftIslandCDists < islandRadius);
    rightIslandOverlap = any(rightIslandCDists < islandRadius);
    someIslandOverlap = leftIslandOverlap || rightIslandOverlap;
    if someIslandOverlap
       someIslandOverlapIdxs = [someIslandOverlapIdxs i]; %#ok<AGROW>
       
       if leftIslandOverlap
          leftIslandMinDists(i) = 0;
       elseif rightIslandOverlap
           rightIslandMinDists(i) = 0;
       else
           error("Logic error in island overlap.")
       end
       
    end
    
    % Second, check for complete overlap with the island. These will be
    % deleted.
    fullIslandOverlap = all(leftIslandCDists < islandRadius) || all(rightIslandCDists < islandRadius);
    if fullIslandOverlap
       fullIslandOverlapIdxs = [fullIslandOverlapIdxs i]; %#ok<AGROW>
    end
end


%                              %
% Add cluster-island distances %
%                              %

% The first row and column are devoted to the left island,
% whereas the last row and column are for the right island.
fullDistArr = zeros(nClustersTotal + 2, nClustersTotal + 2);
fullDistArr(2:nClustersTotal+1, 2:nClustersTotal+1) = clusterDistArr;
for i = 1:nClustersTotal
    fullDistArr(1, i+1) = leftIslandMinDists(i);
    fullDistArr(i+1, 1) = leftIslandMinDists(i);
    fullDistArr(nClustersTotal+2, i+1) = rightIslandMinDists(i);
    fullDistArr(i+1, nClustersTotal+2) = rightIslandMinDists(i);
end

%Add inter-island distance
interIslandDist = norm(leftIslandLoc - rightIslandLoc) - 2 * islandRadius;
fullDistArr(1, nClustersTotal + 2) = max(interIslandDist, 0);
fullDistArr(nClustersTotal + 2, 1) = max(interIslandDist, 0);


%                     %
% Write pdb of system %
%                     %
if pdbname
    writeAtomisticPDB(pdbname, nClustersTotal, clusterPositions, fullSize, ...
                leftIslandLoc, rightIslandLoc, someIslandOverlapIdxs)
end

%                                          %
% Remove clusters overlapping with islands % 
%                                          %

% Shift idxs because we added a row for the left island
fullIslandOverlapIdxs = fullIslandOverlapIdxs + 1;
fullIslandOverlapIdxs = sort(fullIslandOverlapIdxs, 'descend');
for i = 1:size(fullIslandOverlapIdxs, 2)
    idx = fullIslandOverlapIdxs(i);
    fullDistArr = deleteRowColumn(fullDistArr, idx);
end


%                             %
% Calculate the hopping rates %
%                             %

K = exp(-betaValue * fullDistArr);
%subtract off-diagonal (set it to zero)
K = K - diag((diag(K))); 
% assert(all(all(diag(K) == 0)), "Diagonal elements of K should be zero.")

% Steady-state condition - no accumulation of charge. 
% Diagonal elements of islands will be further edited in the sscompute method. 
for i = 1:size(K, 2)
    t=sum(K(:,i));
    K(i,i)=-t;
end

end