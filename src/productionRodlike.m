% This script calculates the charge transport rate between two PEDOT-rich 
% grains embedded in a PSS  matrix with rod-like PEDOT particles.
% The CT rate is calculated as a function of concentration (modulated by 
%  @pedotsPerCell) and distance (modulated by @boxLengths).
% The length of the rods can be adjusted on lines 34-41.

% We loop over the vectors in the order they're defined here.
pedotsPerCell = [21 10 5 3]; % number of N=6 and N=12 (each) pedot lengths in cell
boxLengths = [3 4 5 6 7]; % multiples of lengths of unit cell
nreps = 50; % number of replicates to do

% Parameters
unitLength = 85; % Angstroms
beta = 0.3;
islandRadius = 100;
write_pdb = 0;

% ------ Start doing the loops ------ %
S = zeros(length(pedotsPerCell), size(boxLengths, 2), nreps);
t = zeros(length(pedotsPerCell), size(boxLengths, 2), nreps);

for i=1:length(pedotsPerCell)
    
    for j=1:length(boxLengths)
        
        box = [unitLength*boxLengths(j) unitLength*5 unitLength*5];
        cells = 5 * 5 * boxLengths(j);
        
        % vector of number of molecules each length, corresponding to the
        %   index position (e.g. [0, 0, 1, 0, 0, 2] has one trimer and two
        %   hexamers).
        nPedots = zeros(1, 24);
        
        % For rods of length N=6:
%         nPedots(6) = round((pedotsPerCell(i) * cells-7)*3);

        % For rods of length N=12:
%         nPedots(12) = round((pedotsPerCell(i) * cells-7)*1.5);

        % For rods of length N=18:
        nPedots(18) = round(pedotsPerCell(i) * cells  - 7); 
 
        fprintf("--> Rodlike: %d per cell, %d x %d x %d box:\n", pedotsPerCell(i), ...
                    boxLengths(j), 5, 5)
        
        tstart = tic;
        for k=1:nreps
            
            if write_pdb == 1
                pdbname = sprintf("rodlike_pdbs/n%d_c%d_l%d_r%d.pdb", n, i, j, k);
            else
                pdbname = "";
            end
                
            K = buildRodlikeRateArray(box, nPedots, beta, islandRadius, pdbname);
            v = sscompute(K);  % compute the steady state rate 
            S(i, j, k) = log(v);
            t(i, j, k) = log(v);
        end
        telapsed = toc(tstart);
        fprintf('     finished %d reps at %8.3f seconds per rep (%8.3f total) \n', ...
          nreps, telapsed/nreps, telapsed);
    end
end
     
S=S/log(10); 
VarS=std(t)/log(10);% average steady state rate in Log
