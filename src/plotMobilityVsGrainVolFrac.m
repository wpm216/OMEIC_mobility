% SI panel: atomistic results (beta = 0.3) with different grain volume
% fraction assumptions (0.55, 0.675, 0.8 initial volume fractions)

 
% File paths and parameters
grain_configurations = "../data/grain_MC_sims_varying_dilution/grain_0.675";
atomistic_d_mats = '../data/distance_matrices';
beta = 0.3;
d0 = "../data/grain_MC_sims_varying_dilution/%s";

% Gather kss data and apply to different grain concentrations
rateTableAtomistic_beta03 = productionAtomistic(beta, grainRadius, atomistic_d_mats);
atm_rates_v_0550 = zeros(15, 4);
atm_rates_v_0675 = zeros(15, 4);
atm_rates_v_0800 = zeros(15, 4);
for i=1:15
    d = sprintf(d0, "grain_0.55");
    [x, y] = filmMobility(rateTableAtomistic_beta03(:, :, i), d, 0, []);
    atm_rates_v_0550(i, :) = y;
    d = sprintf(d0, "grain_0.675");
    [~, y] = filmMobility(rateTableAtomistic_beta03(:, :, i), d, 0, []);
    atm_rates_v_0675(i, :) = y;
    d = sprintf(d0, "grain_0.8");
    [~, y] = filmMobility(rateTableAtomistic_beta03(:, :, i), d, 0, []);
    atm_rates_v_0800(i, :) = y;
end

figure;
hold on;

% experiment
x_exp = [1.7 5.1 1.7*5 1.7*10 1.7*15 37.5];
y_exp = [0.4 0.11 -0.3 -0.95 -2.25 -3.71];
plot(x_exp, y_exp, 'k--', 'Markersize', 30)

% x-axis data for simulation
x0550 = [5.4, 10.1, 17.6, 32.5];
x0675 = [4.6, 8.1, 13.5, 23.9];
x0800 = [3.8, 6.2, 9.6, 15.8];

% means
m0550 = mean(atm_rates_v_0550, 1);
m0675 = mean(atm_rates_v_0675, 1);
m0800 = mean(atm_rates_v_0800, 1);

% 95% CI on mean
e0550 = sqrt(var(atm_rates_v_0550)/15) * 1.753;
e0675 = sqrt(var(atm_rates_v_0675)/15) * 1.753;
e0800 = sqrt(var(atm_rates_v_0800)/15) * 1.753;

% do plot
shift = 0.1;
t1 = errorbar(x0550, m0550 - m0550(1) + 0.1, e0550, '-o');
set(t1, 'MarkerFaceColor', get(t1,'Color'));
t2 = errorbar(x0675, m0675 - m0675(1) + 0.15, e0675, '-o');
set(t2, 'MarkerFaceColor', get(t2,'Color'));
t3 = errorbar(x0800, m0800 - m0800(1) + 0.2, e0800, '-o');
set(t3, 'MarkerFaceColor', get(t3,'Color'));


legend('Experiment', '$v_0 = 0.550$',  '$v_0 = 0.675$', ...  
                        '$v_0 = 0.800$', 'Interpreter', 'latex')

xlabel('1:X (PEDOT to PSS weight ratio)')

ylabel('log_{10}(\mu) (a.u)')
pbaspect([1 1 1])
hold off
