% Generate computational traces for Figure 4e:
%   - atomistic mobility at beta = 0.3
%   - atomistic mobility at beta = 0.6
%   - atomistic mobility at beta = 1.0
% Traces are plotted against approximate experimental results.

% Sampling parameters
n_samples = 15; % we take the average of each atomistic replicate.
d = "../data/grain_MC_sims_varying_dilution/grain_0.675";

% Generate mobility curves from atomistic data.
% First, get inter-grain CT rates vs beta.
% beta = tunnelling attenuation coefficient (= {0.3, 0.6, 1.0}).
grainRadius = 100; % Angstroms
d_mats = '../data/distance_matrices';
rateTableAtomistic_beta03 = productionAtomistic(0.3, grainRadius, d_mats);
rateTableAtomistic_beta06 = productionAtomistic(0.6, grainRadius, d_mats);
rateTableAtomistic_beta1  = productionAtomistic(1.0, grainRadius, d_mats);

% Next, calculate mobility through entire film mediated by inter-grain CT.
beta03_rates = zeros(n_samples, 4);
beta06_rates = zeros(n_samples, 4);
beta1_rates = zeros(n_samples, 4);
for i=1:n_samples
    [x, y] = filmMobility(rateTableAtomistic_beta03(:, :, i), d, 0, []);
    beta03_rates(i, :) = y;
    [~, y] = filmMobility(rateTableAtomistic_beta06(:, :, i), d, 0, []);
    beta06_rates(i, :) = y;
    [~, y] = filmMobility(rateTableAtomistic_beta1( :, :, i), d, 0, []);
    beta1_rates(i, :) = y;
end


% Plot data.
figure;
hold on;

% Experimental mobility.
x_exp = [1.7 5.1 1.7*5 1.7*10 1.7*15 37.5];
y_exp = [0.4 0.11 -0.3 -0.95 -2.25 -3.71];
plot(x_exp, y_exp, '.b', 'Markersize', 30)

% Simulated mobility (atomistic representation of PSS-rich matrix).
shift = 0.1; % for visual clarity
% Mean mobility vs. beta.
b03m = mean(beta03_rates, 1);
b06m = mean(beta06_rates, 1);
b1m  = mean(beta1_rates, 1);
% 95% CI on mean (two-sided T distribution)
b03e = sqrt(var(beta03_rates)/15) * 1.753;
b06e = sqrt(var(beta06_rates)/15) * 1.753;
b1e  = sqrt(var(beta1_rates)/15)  * 1.753;
% do plot
t1 = errorbar(x, b03m - b03m(1) + shift, b03e, '-o');
set(t1, 'MarkerFaceColor', get(t1,'Color'));
t2 = errorbar(x, b06m - b06m(1) + shift, b06e, '-o');
set(t2, 'MarkerFaceColor', get(t2,'Color'));
t3 = errorbar(x, b1m  - b1m(1)  + shift, b1e, '-o');
set(t3, 'MarkerFaceColor', get(t3,'Color'));

legend('Experiment', 'Atomistic ($\beta = 0.3$ \AA$^{-1}$)', 'Atomistic ($\beta = 0.6$ \AA$^{-1}$)', ...
                        'Atomistic ($\beta = 1.0$ \AA$^{-1}$)', 'Interpreter','latex')
ylabel('log_{10}(\mu) (a.u)')
xlabel('1:X (PEDOT to PSS weight ratio)')
pbaspect([1 1 1])
hold off
