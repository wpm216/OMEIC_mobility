% This script calculates the charge transport rate between two PEDOT-rich 
% grains embedded in a PSS  matrix with point-like PEDOT particles.
% The CT rate is calculated as a function of concentration (modulated by 
%  @molars) and distance (modulated by @lengths).

r=100.0;  % radius of dots in angstrom

molars = [21 10 5 3]/3.683/100; % effective concentration of particles,
                                % equal to atomistic and rodlike values.
lengths = [3 4 5 6 7] * 85;
nreps = 50;
beta = 0.3;

rate_table = zeros(length(molars), length(lengths), nreps);

for j=1:length(molars)
    molar=molars(j);
    for i=1:length(lengths)   
        l=lengths(i);      
        
        % Set simulation parameters.
        V = l * 5*85 * 5*85 ;    % Volume of total box
        V1 = V - 4/3*pi*r^3;     % Volume occupied by the grains
        nn=round(V1/1667*molar); % Number of molecules
        nPerCell = nn/V1 * (85^3);
        
        fprintf('nPerCell = %8.3f, l = %d, ', nPerCell, l);
        start = tic;
        rate_table(j,i,:) = pointlikeRate(l,nn,nreps,beta).';
        finish = toc(start);
        fprintf('t_per_rep = %8.3f\n', finish/nreps);
    end
end
rate_table = rate_table / log(10);