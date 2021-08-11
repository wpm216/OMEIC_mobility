function [K, fullDistArr] = buildAtomisticRateArray(unitSize, fullSize, unitPositions, nnVecArrays, ...
                            betaValue, grainRadius, leftgrainLoc, rightgrainLoc, pdbname)

% Generates the array of inter-particle tunneling rates
% for a pss-rich matrix populated by rod-like particles.

% Parameters:
% -----------
% unitSize = box size of unit cell (Angstroms)
% fullSize = size of the full cell to be built (in Angstroms)
% unitPositions = positions of all clusters in the unit cell
% nnVecArrays = n x n x 3 x 27 array of pairwise cluster nn distance vectors
% betaValue = hopping decay coefficient. Alessandro set it to 0.6.
% grainRadius = radius of grains (Angstroms)
% leftgrainLoc = x,y,z coordinates of left grain 
% rightgrainLoc = x,y,z coordinates of right grain 
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
clusterDistArr = zeros(nClustersTotal, nClustersTotal); % we'll add the grains at the end

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
% Get positions and cluster-grain distances for each cluster % 
%                                                             % 

clusterPositions = cell(nClustersTotal, 1);
somegrainOverlapIdxs = [];
fullgrainOverlapIdxs = [];
leftgrainMinDists = zeros(nClustersTotal, 1);
rightgrainMinDists = zeros(nClustersTotal, 1);
for i = 1:nClustersTotal
    
    % Get the positions of this cluster
    idx = clusterUnitCellIndex(i);
    pos = unitPositions(:, :, idx);
    zeroRows = all(pos == 0, 2); % remove padding of dummy [0,0,0] coordinates
    pos = pos(~zeroRows, :);
    shift = (clusterBoxIndices(i, :)-1) .* unitSize; % translate into correct box
    clusterPositions{i} = pos + repmat(shift, size(pos, 1), 1);
    
    % Get minimum grain-cluster distance
    leftgrainCDists = sqrt(sum((leftgrainLoc - clusterPositions{i}).^2, 2));
    rightgrainCDists = sqrt(sum((rightgrainLoc - clusterPositions{i}).^2, 2));
    leftgrainMinDists(i) = min(leftgrainCDists) - grainRadius;
    rightgrainMinDists(i) = min(rightgrainCDists) - grainRadius;
    
    % Figure out if it overlaps with either grain
    % First, check for partial overlap. These won't be deleted, but their
    % grain distance will be set to zero.
    leftgrainOverlap = any(leftgrainCDists < grainRadius);
    rightgrainOverlap = any(rightgrainCDists < grainRadius);
    somegrainOverlap = leftgrainOverlap || rightgrainOverlap;
    if somegrainOverlap
       somegrainOverlapIdxs = [somegrainOverlapIdxs i]; %#ok<AGROW>
       
       if leftgrainOverlap
          leftgrainMinDists(i) = 0;
       elseif rightgrainOverlap
           rightgrainMinDists(i) = 0;
       else
           error("Logic error in grain overlap.")
       end
       
    end
    
    % Second, check for complete overlap with the grain. These will be
    % deleted.
    fullgrainOverlap = all(leftgrainCDists < grainRadius) || all(rightgrainCDists < grainRadius);
    if fullgrainOverlap
       fullgrainOverlapIdxs = [fullgrainOverlapIdxs i]; %#ok<AGROW>
    end
end


%                              %
% Add cluster-grain distances %
%                              %

% The first row and column are devoted to the left grain,
% whereas the last row and column are for the right grain.
fullDistArr = zeros(nClustersTotal + 2, nClustersTotal + 2);
fullDistArr(2:nClustersTotal+1, 2:nClustersTotal+1) = clusterDistArr;
for i = 1:nClustersTotal
    fullDistArr(1, i+1) = leftgrainMinDists(i);
    fullDistArr(i+1, 1) = leftgrainMinDists(i);
    fullDistArr(nClustersTotal+2, i+1) = rightgrainMinDists(i);
    fullDistArr(i+1, nClustersTotal+2) = rightgrainMinDists(i);
end

%Add inter-grain distance
intergrainDist = norm(leftgrainLoc - rightgrainLoc) - 2 * grainRadius;
fullDistArr(1, nClustersTotal + 2) = max(intergrainDist, 0);
fullDistArr(nClustersTotal + 2, 1) = max(intergrainDist, 0);


%                     %
% Write pdb of system %
%                     %
if pdbname
    writeAtomisticPDB(pdbname, nClustersTotal, clusterPositions, fullSize, ...
                leftgrainLoc, rightgrainLoc, somegrainOverlapIdxs, [])
end

%                                          %
% Remove clusters overlapping with grains % 
%                                          %

% Shift idxs because we added a row for the left grain
fullgrainOverlapIdxs = fullgrainOverlapIdxs + 1;
fullgrainOverlapIdxs = sort(fullgrainOverlapIdxs, 'descend');
for i = 1:size(fullgrainOverlapIdxs, 2)
    idx = fullgrainOverlapIdxs(i);
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
% Diagonal elements of grains will be further edited in the sscompute method. 
for i = 1:size(K, 2)
    t=sum(K(:,i));
    K(i,i)=-t;
end

end