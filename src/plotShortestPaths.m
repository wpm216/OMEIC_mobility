% This function generates SI Extended Data Figure 8b, which plots the
% average hopping distance vs. PEDOT to PSS weight ratio. See SI Section S5
% for a methodological overview.

addpath('dijkstra') % add dijkstra's algorithm implementation to path.

%%%%%%%%%%%

% Simulation parameters for atomistic, rodlike, and pointlike
% representations.

pedotsPerCell = [20.6 9.8 5.2 2.7 1.65]; % number of rodlike and pointlike particles per cell.
                                         % added sig figs for curve
                                         % smoothness
unitLength = 85; % Angstroms
cellsOneDim = 5; % number of cells along one direction of the box.
% together, the box we place particles in is (5 * 85)^3 Angstroms^3.
% this is approximately the size of the atomistic simulation boxes.
% there are pedotsPerCell(i) * cellsOneDim^3 particles in the box.

nreps = 50; % number of replicates for rodlike and pointlike particles at each concentration
            % (each element of pedotsPerCell). there are 15 replicates for 
            % the atomistic simulations.
            
beta = 0.3; % tunneling attentuation coefficient (1/Angstroms)

%%%%%%%%%%%

% Get the fastest paths (and corresponding hopping distances) for 
% atomistic, rodlike, and pointlike representations of the PSS-rich matrix.
[atm_lengths, atm_distances] = shortestPathAtomistic(beta, cellsOneDim);
[rod_lengths, rod_distances] = shortestPathRodlike(pedotsPerCell, ...
                        cellsOneDim, nreps, unitLength, beta);
[pnt_lengths, pnt_distances] = shortestPathPointlike(pedotsPerCell, ...
                        cellsOneDim, nreps, unitLength, beta);

%%%%%%%%%%%

% Plot results

figure;
hold on;

% Atomistic data
x = [4.6 8.1 13.53 23.85];
atm_mean_dists = mean(sum(atm_distances, 3) ./ atm_lengths, 2);
atm_var_dists = var(sum(atm_distances, 3) ./ atm_lengths, 0, 2);
atm_err_dists = sqrt(atm_var_dists) * 2.131;
errorbar(x, atm_mean_dists, atm_err_dists);

% Rodlike data
x = [4.6 8.1 13.53 23.85 37.5];
rod_mean_dists = mean(sum(rod_distances, 3) ./ rod_lengths, 2);
rod_var_dists = var(sum(rod_distances, 3) ./ rod_lengths, 0, 2);
rod_err_dists = sqrt(rod_var_dists) * 2.01;
errorbar(x, rod_mean_dists, rod_err_dists);

% Pointlike data
x = [4.6 8.1 13.53 23.85 37.5];
pnt_mean_dists = mean(sum(pnt_distances, 3) ./ pnt_lengths, 2);
pnt_var_dists = var(sum(pnt_distances, 3) ./ pnt_lengths, 0, 2);
pnt_err_dists = sqrt(pnt_var_dists) * 2.01;
errorbar(x, pnt_mean_dists, pnt_err_dists);

% Power law fit (rate_predicted = p(2) * concentration ^ p(1))
p = polyfit(log(pedotsPerCell), log(pnt_mean_dists), 1);
y_predicted = exp(p(1) * log(pedotsPerCell) + p(2));
plot(x, y_predicted, 'k--')

% Add plot labels and legend
xlabel("PEDOT to PSS weight ratio (1:X)")
ylabel("Average tunneling distance (Angstroms)")
legend("Atomistic", "Rodlike", "Pointlike", "Power Law Fit")
