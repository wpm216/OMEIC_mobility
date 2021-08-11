function [rod_lengths, rod_distances] = shortestPathRodlike(pedotsPerCell, ...
                        cellsOneDim, nreps, unitLength, beta)

% This function finds the fastest hopping path through rod-like 
% representations of the PSS-rich matrix.
% We calculate all of the inter-particle hopping times (1/rate)
% and use Dijkstra's algorithm to find the fastest path through tie film.


% ------ Start doing the loops ------ %
rod_costs = zeros(length(pedotsPerCell), length(nreps));
rod_routes = zeros(length(pedotsPerCell), length(nreps), 125);
rod_lengths = zeros(length(pedotsPerCell), length(nreps));
rod_distances = zeros(length(pedotsPerCell), length(nreps), 125);

for i=1:length(pedotsPerCell)
    for j=1:nreps
        fprintf('Rodlike: n = %d, rep %d\n', pedotsPerCell(i), j);        
       
        % make simulation box
        l = unitLength*cellsOneDim;
        box = [l l l];
        cells = cellsOneDim^3;

        nPedots = zeros(1, 24);
        nPedots(18) = round(pedotsPerCell(i) * cells);
        nMols = sum(nPedots);
        
        % Find the inter-particle distance array, @dmat
        [~, dmat] = buildRodlikeRateArray(box, nPedots, beta, 0, "");
        
        % Find the quickest path between the two infinitesimal grains.
        % The matrix of pairwise particle travel times is calculated
        % as 1/(pairwise particle hopping rates).
        % The diagonal elements are set to zero to remove self-hopping.
        tmat = exp(beta * dmat) - diag(ones(1, nMols + 2));
        % @cost is the total travel time of the shortest path.
        % @route is a list of molecle indices the path takes.
        [cost, route] = dijkstra(tmat, 1, nMols+2);
        
        rod_costs(i, j) = cost;
        rod_routes(i, j, 1:length(route)) = route;
        rod_lengths(i, j) = length(route);
        
        % find the length of each hop in the quickest path
        for k=1:length(route)-1
            % disallow hops below 5 Angstroms
            d = dmat(route(k), route(k+1));
            if d > 5
                rod_distances(i, j, k) = d;
            else
                rod_distances(i, j, k) = 0;
                rod_lengths(i, j) = rod_lengths(i, j) - 1;
            end             
        end
        
%         % make a pdb file for troubleshooting
%         pdbname = sprintf('pointlike_pdbs/ROD_n%d_rep%d.pdb', pedotsPerCell(i), j);
%         writeRodlikePDB(pdbname, coords, box, ...
%                                [0 l/2 l/2], [l l/2 l/2], route);

    end
end


    