function rateTableAtomistic = productionAtomistic(betaValue, grainRadius, d)

% This script calculates the charge transport rate between two PEDOT-rich 
% grains embedded in a PSS matrix with atomistic PEDOT particles.
% The CT rate is calculated as a function of PEDOT to PSS weight ratio,
%  represented in the distance matrices in r1to2 (PEDOT to PSS weight ratio
%  of 1 to 2 in the matrix), r1to5, r1to10, and r1to20. 
%
% Multiple box lengths are used to calculate the CT rate vs. distance.
%  The unit cell from atomistic MD simulation is tesselated to fill 
%  the box used in this CT calculation.

% Load atomistic simulation data.
[r1to2, r1to5, r1to10, r1to20] = loadAtomisticData(d); % r ~ weight ratio
simulations = [r1to2, r1to5, r1to10, r1to20].';

% Initialize containers.
rateTableAtomistic = zeros(4, 5, 15); % zeros(concentrations, boxLengths, replicates)
% rateTableAtomistic_var = zeros(4, 5);

% Loop through PEDOT to PSS weight ratios.
for i=1:4 % 1:4
    for boxLength=3:7 % 3:7 is  240 to 560 angstroms, approximately
        fprintf('Starting simulation %d, box length %d\n', i, boxLength)
        dirname = sprintf('pdbs_conc_%d', i);
        mkdir(dirname);
        
        % Calculate average CT rate between grains for the given
        %  weight ratio and box length.
        [S, ~, ~, ~] = atomisticAverageRate(simulations(i,:), ...
                                boxLength, betaValue, 0, grainRadius);
        rateTableAtomistic(i, boxLength - 2, :) = S;
%         rateTableAtomistic_var(i, boxLength - 2) = VarS;
    end
end
