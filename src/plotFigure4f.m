
% Generate computational traces for Figure 4f:
%   - atomistic mobility, beta = 0.3
%   - rodlike mobility, beta = 0.3
%   - pointlike mobility, beta = 0.3
% Traces are lotted against approximate experimental results 
%  (see publication for exact data).

%%%%%%%%%%

% General parameters for atomistic, rodlike, and pointlike simulations.
beta = 0.3; % 1/Angstroms

% Grain and box parameters
grainRadius = 100; % Angstroms
unitCellLength = 85; % Angstroms
boxLengths = [3 4 5 6 7]; % multiples of lengths of unit cell
nreps = [25 50 500 500]; % number of replicates to do for each (pedotsPerCell, boxLengths) pair
                       % (rodlike and pointlike)

% PEDOT-focused parameters
pedotsPerCell = [21 10 5 3]; % number of N=6 and N=12 (each) pedot lengths in cell (atomistic)
pedotLengthDistribution = zeros(1, 24); % number of PEDOTs of length=k in each cell 
                                        % (scaling factor applied to elements of pedotsPerCell(i))
pedotLengthDistribution(18) = 1;

% Sampling parameters
grain_configurations = "../data/grain_MC_sims_varying_dilution/grain_0.675";
atomistic_d_mats = '../data/distance_matrices';
n_bootstrap_samples = 50;

% Miscellaneous parameters
write_pdb = 0; % flag to write PDBs (1 = true, 0 = false)

%%%%%%%%%%

% Run simulations.

% for atomistic samples. estimate error with mean of samples.
rateTableAtomistic_beta_03 = productionAtomistic(beta, grainRadius, atomistic_d_mats);
atomistic_mobilities = zeros(15, 4);
for i=1:15
    [~, y] = filmMobility(rateTableAtomistic_beta_03(:, :, i), grain_configurations, 0, []);
    atomistic_mobilities(i, :) = y;
end

% for rodlike samples. do bootstrapping in filmMobility to estimate error.
rateTableRodlike_beta_03 = productionRodlike(pedotsPerCell, boxLengths, ...
    unitCellLength, pedotLengthDistribution, nreps, beta, grainRadius, write_pdb);
rodlike_mobilities = zeros(n_bootstrap_samples, 4);
for i=1:n_bootstrap_samples
    [~, y] = filmMobility(rateTableRodlike_beta_03, grain_configurations, 1, nreps);
    rodlike_mobilities(i, :) = y;
end

% for pointlike samples. do bootstrapping in filmMobility to estimate error.
rateTablePointlike_beta_03 = productionPointlike(pedotsPerCell, beta, ... 
                            nreps, grainRadius, boxLengths, unitCellLength);
pointlike_mobilities = zeros(n_bootstrap_samples, 4);
for i=1:n_bootstrap_samples
    [~, y] = filmMobility(rateTablePointlike_beta_03, grain_configurations, 1, nreps);
    pointlike_mobilities(i, :) = y;
end

% Plot data.

figure;
hold on;

% Experiment.
x = [1.7 5.1 1.7*5 1.7*10 1.7* 15]; % omitting last point for clarity
y = [0.4 0.11 -0.3 -0.95 -2.25];    % omitting last point for clarity
plot(x, y, '.b', 'Markersize', 30)

% Simulation.
x = [4.6 8.1 13.53 23.85];

% means
mean_atm = mean(atomistic_mobilities, 1);  % atomistic, beta = 0.3
mean_rod = mean(rodlike_mobilities, 1);    % rodlike, N=18, beta = 0.3
mean_pnt = mean(pointlike_mobilities, 1);  % pointlike, beta = 0.3

% 95% CI on mean
err_atm = sqrt(var(atomistic_mobilities)/15) * 1.753;
err_rod = sqrt(var(rodlike_mobilities)/50) *  2.0086;
err_pnt = sqrt(var(pointlike_mobilities)/50) *  2.0086;

% do plot
t1 = errorbar(x, mean_atm - mean_atm(1) - 0.1, err_atm, '-o');
set(t1, 'MarkerFaceColor', get(t1,'Color'));
t2 = errorbar(x, mean_rod - mean_rod(1) - 0.1, err_rod, '-o');
set(t2, 'MarkerFaceColor', get(t2,'Color'));
t3 = errorbar(x, mean_pnt - mean_pnt(1) - 0.1, err_pnt, '-o');
set(t3, 'MarkerFaceColor', get(t3,'Color'));

legend('Experiment', 'Atomistic', 'Rod-like ($N = 18$)', 'Point-like', ...
                'Interpreter', 'latex')

xlabel('1:X (PEDOT to PSS weight ratio)')

ylabel('log_{10}(\mu) (a.u)')
pbaspect([1 1 1])
hold off
