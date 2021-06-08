% Generate computational traces for Figure 4f:
%   - atomistic mobility, beta = 0.3
%   - rodlike mobility, beta = 0.3
%   - pointlike mobility, beta = 0.3
% Traces are lotted against approximate experimental results 
%  (see publication for exact data).

d = "../data/island_MC_sims_varying_dilution/island_0.675";

% for atomistic samples. estimate error with mean of samples.
beta03_rates = zeros(15, 4);
for i=1:15
    [~, y] = filmMobility(rateTableAtomistic_beta03(:, :, i), d, 0, []);
    beta03_rates(i, :) = y;
end

% for rodlike samples. do bootstrapping in wesMobility to estimate error.
n = 50;
n18_rates = zeros(n, 4);
rateArrays = zeros(4, 5, n);
for i=1:n
    [~, y] = filmMobility(n18total_b03, d, 1, [50 50 250 250]);
    n18_rates(i, :) = y;
end

% for pointlike samples. do bootstrapping in wesMobility to estimate error.
n = 50;
pl_rates = zeros(n, 4);
for i=1:n
    [~, y] = filmMobility(rateTablePointlike_alessandro_beta03_n18_rd2, ...
                            d, 1, [50 50 50 50]);
    pl_rates(i, :) = y;
end

figure;
hold on;

% experiment
x = [1.7 5.1 1.7*5 1.7*10 1.7* 15]; % omitting last point for clarity
y = [0.4 0.11 -0.3 -0.95 -2.25];    % omitting last point for clarity
plot(x, y, '.b', 'Markersize', 30)

x = [4.6 8.1 13.53 23.85];

% means
b03m = mean(beta03_rates, 1); % atomistic, beta = 0.3
n18m = mean(n18_rates, 1);    % rodlike, N=18, beta = 0.3
plm = mean(pl_rates, 1);      % pointlike, beta = 0.3


% 95% CI on mean
b03e = sqrt(var(beta03_rates)/15) * 1.753;
n18e = sqrt(var(n18_rates)/50) *  2.0086;
ple = sqrt(var(pl_rates)/50) *  2.0086;

% do plot
t1 = errorbar(x, b03m - b03m(1) - 0.1, b03e, '-o');
set(t1, 'MarkerFaceColor', get(t1,'Color'));
t2 = errorbar(x, n18m - n18m(1) - 0.1, n18e, '-o');
set(t2, 'MarkerFaceColor', get(t2,'Color'));
t3 = errorbar(x, mean_pl_y - mean_pl_y(1) - 0.1, ple, '-o');
set(t3, 'MarkerFaceColor', get(t3,'Color'));

legend('Experiment', 'Atomistic', 'Rod-like ($N = 18$)', 'Point-like', ...
                'Interpreter', 'latex')

xlabel('1:X (PEDOT to PSS weight ratio)')

ylabel('log_{10}(\mu) (a.u)')
pbaspect([1 1 1])
hold off
