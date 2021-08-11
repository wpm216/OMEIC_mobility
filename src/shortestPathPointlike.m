function [pnt_lengths, pnt_distances] = shortestPathPointlike(pedotsPerCell, ...
                        cellsOneDim, nreps, unitLength, beta)

% This function finds the fastest hopping path through point-like 
% representations of the PSS-rich matrix.
% We calculate all of the inter-particle hopping times (1/rate)
% and use Dijkstra's algorithm to find the fastest path through tie film.

% Initialize containers. 
pnt_costs = zeros(length(pedotsPerCell), length(nreps));
pnt_routes = zeros(length(pedotsPerCell), length(nreps), 125); % fastest paths will be shorter than 125 steps
pnt_lengths = zeros(length(pedotsPerCell), length(nreps));
pnt_distances = zeros(length(pedotsPerCell), length(nreps), 125);

for i=1:length(pedotsPerCell)
    for j=1:nreps
        fprintf('Pointlike: n = %d, rep %d\n', pedotsPerCell(i), j);
       
        % make simulation box
        l = unitLength*cellsOneDim;
        box = [l l l];
        cells = cellsOneDim^3;

        % this replaces each N=18 pedot with a pointlike particle.
        nMols = round(pedotsPerCell(i) * cells);
        pedotCoords = rand(nMols, 3) .* box;
        
        % add grains
        coords = [0 l/2 l/2 ; pedotCoords; l l/2 l/2];

        % calculate inter-particle distances
        dmat = zeros(nMols + 2, nMols + 2);
        for x1=1:nMols+2
            for x2=x1:nMols+2
                v = coords(x1, :) - coords(x2, :);
                dmat(x1, x2) = norm(v); 
                dmat(x2, x1) = dmat(x1, x2); 
            end
        end
        
        % Find the quickest path between the two infinitesimal grains.
        % The matrix of pairwise particle travel times is calculated
        % as 1/(pairwise particle hopping rates).
        % The diagonal elements are set to zero to remove self-hopping.
        tmat = exp(beta * dmat) - diag(ones(1, nMols + 2));
        % @cost is the total travel time of the shortest path.
        % @route is a list of molecle indices the path takes.
        [cost, route] = dijkstra(tmat, 1, nMols+2);
        
        pnt_costs(i, j) = cost;
        pnt_routes(i, j, 1:length(route)) = route;
        pnt_lengths(i, j) = length(route);
        
        % find the length of each hop in the quickest path
        for k=1:length(route)-1
            d = dmat(route(k), route(k+1));
            % disallow hops below 5 Angstroms
            if d > 5
                pnt_distances(i, j, k) = d;
            else
                pnt_distances(i, j, k) = 0;
                pnt_lengths(i, j) = pnt_lengths(i, j) - 1;
            end
        end
        
%         % make a pdb file for troubleshooting
%         pdbname = sprintf('pointlike_pdbs/n%d_rep%d.pdb', pedotsPerCell(i), j);
%         writePointlikePDB(pdbname, coords(2:nMols+1, :), box, ...
%                                coords(1, :), coords(nMols+2, :), route);

    end
end


    