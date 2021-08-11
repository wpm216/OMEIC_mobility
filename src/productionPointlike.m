function rateTablePointlike = productionPointlike(pedotsPerCell, beta, ... 
                            nreps, grainRadius, boxLengths, unitSize)

% This script calculates the charge transport rate between two PEDOT-rich 
% grains embedded in a PSS  matrix with point-like PEDOT particles.
% The CT rate is calculated as a function of concentration (modulated by 
%  @molars) and distance (modulated by @lengths).


molars = pedotsPerCell/3.683/100; % effective concentration of particles,
                                  % equal to atomistic and rodlike values.
lengths = boxLengths * unitSize;


rateTablePointlike = zeros(length(molars), length(lengths), max(nreps));

for j=1:length(molars)
    molar=molars(j);
    for i=1:length(lengths)   
        l=lengths(i);      
        
        % Set simulation parameters.
        V = l * 5*85 * 5*85 ;    % Volume of total box
        V1 = V - 4/3*pi*grainRadius^3;     % Volume occupied by the grains
        nn=round(V1/1667*molar); % Number of molecules
        nPerCell = nn/V1 * (85^3);
        
        fprintf('nPerCell = %8.3f, l = %d, ', nPerCell, l);
        start = tic;
        rateTablePointlike(j,i,1:nreps(j)) = pointlikeRate(l,nn,nreps(j),beta).';
        finish = toc(start);
        fprintf('t_per_rep = %8.3f\n', finish/nreps(j));
    end
end
rateTablePointlike = rateTablePointlike / log(10);
