function rateTableRodlike = productionRodlike(pedotsPerCell, boxLengths, ...
unitLength, pedotLengthDistribution, nreps, beta, grainRadius, write_pdb)

% This script calculates the charge transport rate between two PEDOT-rich 
% grains embedded in a PSS  matrix with rod-like PEDOT particles.
% The CT rate is calculated as a function of concentration (modulated by 
%  @pedotsPerCell) and distance (modulated by @boxLengths).
% The length of the rods can be adjusted on lines 34-41.



% ------ Start doing the loops ------ %
rateTableRodlike = zeros(length(pedotsPerCell), size(boxLengths, 2), max(nreps));

for i=1:length(pedotsPerCell)
    
    for j=1:length(boxLengths)
        
        box = [unitLength*boxLengths(j) unitLength*5 unitLength*5];
        cells = 5 * 5 * boxLengths(j);
        % we subtract "7" to account for the grain volume in the box,
        % to keep the concentration correct.
        nPedots = pedotLengthDistribution .* pedotsPerCell(i) * cells - ...
                    7 * (pedotLengthDistribution > 0);

 
        fprintf("--> Rodlike: %d per cell, %d x %d x %d box:\n", pedotsPerCell(i), ...
                    boxLengths(j), 5, 5)
        
        tstart = tic;
        for k=1:nreps(i)
            
            if write_pdb == 1
                pdbname = sprintf("rodlike_pdbs/n%d_c%d_l%d_r%d.pdb", n, i, j, k);
            else
                pdbname = "";
            end
                
            [K, ~] = buildRodlikeRateArray(box, nPedots, beta, grainRadius, pdbname);
            v = sscompute(K);  % compute the steady state rate 
            rateTableRodlike(i, j, k) = log(v);

        end
        telapsed = toc(tstart);
        fprintf('     finished %d reps at %8.3f seconds per rep (%8.3f total) \n', ...
          nreps(i), telapsed/nreps(i), telapsed);
    end
end

% scale from base e to base 10
rateTableRodlike=rateTableRodlike/log(10); 
