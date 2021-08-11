function [atm_lengths, atm_distances] = shortestPathAtomistic(beta, boxLength)

% This function finds the fastest hopping path through each atomistic 
% configuration of the PSS-rich matrix.
% We calculate all of the inter-particle hopping times (1/rate)
% and use Dijkstra's algorithm to find the fastest path through tie film.

% Load all simulation data.
d = '../data/distance_matrices';
[r1to2, r1to5, r1to10, r1to20] = loadAtomisticData(d);
simulations = [r1to2, r1to5, r1to10, r1to20].';

% ------ Start doing the loops ------ %
atm_costs = zeros(size(simulations, 1), 15);
atm_routes = zeros(size(simulations, 1), 15, 125);
atm_lengths = zeros(size(simulations, 1), 15);
atm_distances = zeros(size(simulations, 1), 15, 125);

% Loop over PEDOT:PSS weight ratios
for i=1:size(simulations, 1)
    % Loop over five replicates of each weight ratio
    for j=1:5
        
        % Load data for a single replicate
        data = simulations(i, j);
        unitSize = data.l;
        unitPositions = data.r;
        nnVecArrays = data.v;

        % Make boxes in x, y, and z directions
        xBoxSize = [boxLength, 5, 5] .* unitSize;
        xLeftgrain = [0, 5/2, 5/2].* unitSize;
        xRightgrain = [boxLength, 5/2, 5/2] .* unitSize;

        yBoxSize = [5, boxLength, 5] .* unitSize;
        yLeftgrain = [5/2, 0, 5/2].* unitSize;
        yRightgrain = [5/2, boxLength, 5/2] .* unitSize;

        zBoxSize = [5, 5, boxLength] .* unitSize;
        zLeftgrain = [5/2, 5/2, 0] .* unitSize;
        zRightgrain = [5/2, 5/2, boxLength] .* unitSize;

        boxes = [xBoxSize; yBoxSize; zBoxSize];

        grains = [xLeftgrain, xRightgrain;
                   yLeftgrain, yRightgrain;
                   zLeftgrain, zRightgrain];

        % Do calculations in x,y,z directions
        for k=1:3
            idx = (j-1)*3+k;
            
            fprintf("Atomistic cells: concentration #%d, replicate #%d/15\n", i, (j-1)*3+k)
            
            % Find the inter-particle distance array, @dmat
            [~, dmat] = buildAtomisticRateArray(unitSize, boxes(k, :), unitPositions, nnVecArrays, ...
                                beta, 0, grains(k, 1:3), grains(k, 4:6), 0);

            % Find the quickest path between the two infinitesimal grains.
            % The matrix of pairwise particle travel times is calculated
            % as 1/(pairwise particle hopping rates).
            % The diagonal elements are set to zero to remove self-hopping.
            tmat = exp(beta * dmat) - diag(ones(1, size(dmat, 1))); 
            % @cost is the total travel time of the shortest path.
            % @route is a list of molecle indices the path takes.
            [cost, route] = dijkstra(tmat, 1, size(dmat, 1));
            
            atm_costs(i, idx) = cost;
            atm_routes(i, idx, 1:length(route)) = route;
            atm_lengths(i, idx) = length(route);

            % find the length of each hop in the quickest path
            for x=1:length(route)-1
                % disallow hops below 5 Angstroms
                d = dmat(route(x), route(x+1));
                if d > 5
                    atm_distances(i, idx, x) = d;
                else
                    atm_distances(i, idx, x) = 0;
                    atm_lengths(i, idx) = atm_lengths(i, idx) - 1;
                end
            end
            
%             % write a pdb for troubleshooting
%             pdbname = sprintf("pointlike_pdbs/ATM_n%d_rep%d.pdb", i, idx);
%             nClustersTotal = 125 * size(nnVecArrays, 1);
%             writeAtomisticPDB(pdbname, nClustersTotal, pos, boxes(k, :), ...
%                  grains(k, 1:3), grains(k, 4:6), [], route)

        end
    end
end

