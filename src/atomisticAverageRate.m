function [S, VarS, xBoxSize, v] = atomisticAverageRate(datfiles, boxLength, ...
                                   betaValue, dirname, islandRadius)

% datfiles: a 1-D array of matlab objects, each containing:
%   d.l: box dimensions of (cubic) unit cell.
%   d.r: positions of clusters in unit cell. 
%   d.v: shortest pairwise inter-cluser vectors for all 27 nearest-neighbor
%           images (has dimensions [size(d.r, 2), size(d.r, 2), 3, 27])
%   d.d: shortest pairwise inter-cluster distances for all nn images 
%           (not used)
% boxLength: number of images in the direction of charge carrier flow.
%             We set number of images in the other two directions to 5.
% beta: tunneling decay coefficient. alessandro set to 0.6

nfiles = size(datfiles, 2);
S=zeros(nfiles * 3,1);
t=zeros(nfiles * 3,1);
for j=1:nfiles 
    
    % Load data
    data = datfiles(j);
    unitSize = data.l;
    unitPositions = data.r;
    nnVecArrays = data.v;
    
    % Make boxes in x, y, and z directions
    xBoxSize = [boxLength, 5, 5] .* unitSize;
    xLeftIsland = [0, 5/2, 5/2].* unitSize;
    xRightIsland = [boxLength, 5/2, 5/2] .* unitSize;
    
    yBoxSize = [5, boxLength, 5] .* unitSize;
    yLeftIsland = [5/2, 0, 5/2].* unitSize;
    yRightIsland = [5/2, boxLength, 5/2] .* unitSize;
    
    zBoxSize = [5, 5, boxLength] .* unitSize;
    zLeftIsland = [5/2, 5/2, 0] .* unitSize;
    zRightIsland = [5/2, 5/2, boxLength] .* unitSize;
    
    boxes = [xBoxSize; yBoxSize; zBoxSize];
    
    islands = [xLeftIsland, xRightIsland;
               yLeftIsland, yRightIsland;
               zLeftIsland, zRightIsland];
    
    for k=1:3
        
        fprintf('\t Replicate %d --', (j-1)*3 + k)
        tstart = tic;
        
        if dirname
            pdbname = sprintf('%s/c%d_L%d_r%d.pdb', dirname, j, boxLength, k);
        else
            pdbname = 0;
        end

        % Run the model (x direction)
        [K, ~] = buildAtomisticRateArray(unitSize, boxes(k, :), unitPositions, nnVecArrays, ...
                            betaValue, islandRadius, islands(k, 1:3), islands(k, 4:6), pdbname);
        v = sscompute(K);  % compute the steady state rate 
        S(3*(j-1)+k) = log(v);
        t(3*(j-1)+k) = log(v);
        
        telapsed = toc(tstart);
        fprintf(' finished in %8.3f seconds\n', telapsed);
    
    end
end
S=S/log(10); 
VarS=std(t)/log(10);% average steady state rate in Log