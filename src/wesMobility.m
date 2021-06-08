function [x,y] = wesMobility(rateTable, islandDir, bootstrap, rowZMaxs)

% this script calculates the mobility using inter-grain CT rates.
% rateTable comes from productionAtomistic, productionRodlike, or 
%  productionPointlike.

if bootstrap == 1
    rateArray = bootstrapRateTable(rateTable, rowZMaxs).';
else
    rateArray = mean(rateTable, 3, 'omitnan').';
end
    
    
islandRadius = 100; % Angstroms
boxLength = 2263 + 4 * islandRadius; % same as Alessandro

% unitCellSizes = [85.5660 84.2844 83.6185 83.4340]; % for atomistic
unitCellSizes = [85 85 85 85]; % for rod-like and point-like

% steps to get the mobility for a given concentration:
% 1. for all simulation replicates:
%       get the distance matrix
%       linearly interpolate the rate matrix element-by-element
%       use the sscompute function to get the rate

mobilityAtomistic_test = zeros(4, 24);
for i = 1:4
    
    % rateTable is a matrix with Log(rates) for different distances and concentration
    % the following command generates a fine grid of rates as a function of the
    % distance between particle at a given concentration

    d0=200; dd=5; dmax=4000; %range of distances and gridsize
    ndis=round(dmax-d0)/dd;
    rate_dist=zeros(ndis,1);
    s = unitCellSizes(i);
    for k=(1:ndis) 
        d=d0+(k-1)*dd;       
        % The table isn't meant to be interpolated in the y-direction:
        %   we simply look up the correct row.
        rate_dist(k)=10^linterpol(rateArray,3*s,s,1,1,d,i); 
    end
    
    for j = 1:24
        
        % Load the coordinates of the particles
        fname = sprintf("%s/fort.1%d%02d", islandDir, i, j);
        coord = load(fname);
        
        % Calculate the distance array for the particles 
        % Start by adding first and last 'particles' representing box ends
        X=[-2*islandRadius, 0, 0; coord; boxLength + 2*islandRadius, 0, 0]; 

        % Calculate the overall rate.
        K=interGrainCTRate(X,rate_dist,d0,dd);
        v=sscompute(K); 
        mobilityAtomistic_test(i,j)=v;
        
    end
end

x = [4.60 8.10 13.53 23.85]; % from analysis of film dilution 
y = mean(log10(mobilityAtomistic_test), 2);

