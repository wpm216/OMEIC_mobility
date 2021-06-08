function writeAtomisticPDB(pdbname, nClustersTotal, clusterPositions, fullSize, ...
                    leftIslandLoc, rightIslandLoc, islandOverlapIdxs)

% Write a pdb file from CT simulations with atomistic PEDOT particles.

fid = fopen(pdbname, 'wt');
n = 0; % number of atoms total in system
fprintf(fid, "TITLE     SIMULATION CELL\n");
fprintf(fid, "CRYST1 %8.3f %8.3f %8.3f %6.2f %6.2f %6.2f P 1           1\n", ...
                fullSize, 90, 90, 90);
for i = 1:nClustersTotal
    
    positions = clusterPositions{i};
    
    overlap = any(islandOverlapIdxs == i);
    if overlap
        resname = "DEL";
    else
        resname = "MOL";
    end
    
    for j = 1:size(positions, 1)
        n = n + 1;
        nprinted = mod(n, 10000);
        fprintf(fid, 'ATOM  %5d %4s %3s %1s %3d    %8.3f%8.3f%8.3f  1.00  1.00\n', ...
                    nprinted, "C", resname, "X", mod(i, 1000), positions(j, :));
    end
end
% Add islands
fprintf(fid, 'ATOM  %5d %4s %3s %1s %3d    %8.3f%8.3f%8.3f  1.00  1.00\n', ...
            mod(n+1, 10000), "I", "ISL", "X", mod(i+1, 1000), leftIslandLoc);
fprintf(fid, 'ATOM  %5d %4s %3s %1s %3d    %8.3f%8.3f%8.3f  1.00  1.00\n', ...
            mod(n+2, 10000), "I", "ISL", "X", mod(i+2, 1000), rightIslandLoc);
fclose(fid);

end