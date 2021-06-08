function [r1to2, r1to5, r1to10, r1to20] = loadAtomisticData(d)

% Load array of inter-PEDOT distances calculated from atomistic simulation.

r1to2 = [load(sprintf('%s/1to2_rep1.mat', d))
              load(sprintf('%s/1to2_rep2', d))
              load(sprintf('%s/1to2_rep3', d))
              load(sprintf('%s/1to2_rep4', d))
              load(sprintf('%s/1to2_rep5', d))];

r1to5 = [load(sprintf('%s/1to5_rep1.mat', d))
              load(sprintf('%s/1to5_rep2.mat', d))
              load(sprintf('%s/1to5_rep3.mat', d))
              load(sprintf('%s/1to5_rep4.mat', d))
              load(sprintf('%s/1to5_rep5.mat', d))];
          
r1to10 = [load(sprintf('%s/1to10_rep1.mat', d))
              load(sprintf('%s/1to10_rep2.mat', d))
              load(sprintf('%s/1to10_rep3.mat', d))
              load(sprintf('%s/1to10_rep4.mat', d))
              load(sprintf('%s/1to10_rep5.mat', d))];
          
r1to20 = [load(sprintf('%s/1to20_rep1.mat', d))
              load(sprintf('%s/1to20_rep2.mat', d))
              load(sprintf('%s/1to20_rep3.mat', d))
              load(sprintf('%s/1to20_rep4.mat', d))
              load(sprintf('%s/1to20_rep5.mat', d))];


end