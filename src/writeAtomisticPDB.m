function writeAtomisticPDB(pdbname, nClustersTotal, clusterPositions, fullSize, ...
                    leftgrainLoc, rightgrainLoc, grainOverlapIdxs, [])

% Write a pdb file from CT simulations with atomistic PEDOT particles.

fid = fopen(pdbname, 'wt');
n = 0; % number of atoms total in system
fprintf(fid, "TITLE     SIMULATION CELL\n");
fprintf(fid, "CRYST1 %8.3f %8.3f %8.3f %6.2f %6.2f %6.2f P 1           1\n", ...
                fullSize, 90, 90, 90);
            
% Add left grain
fprintf(fid, 'ATOM  %5d %4s %3s %1s %3d    %8.3f%8.3f%8.3f  1.00  1.00\n', ...
            1, "I", "ISL", "X", 1, leftgrainLoc);
        
n = n + 1;

for i = 1:nClustersTotal
    
    positions = clusterPositions{i};
    
    overlap = any(grainOverlapIdxs == i);
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

% Add right grain
fprintf(fid, 'ATOM  %5d %4s %3s %1s %3d    %8.3f%8.3f%8.3f  1.00  1.00\n', ...
            mod(n+1, 10000), "I", "ISL", "X", mod(i+1, 1000), rightgrainLoc);


% add dummy particles corresponding to a given tunneling route, if desired
nAtoms = n+2;
if length(route) > 1
    for j = 2:length(route)-1
        n = n + 1;
        cluster = clusterPositions(route(j)-1);
        mean_pos = mean(cluster{1}, 1);
        fprintf(fid, 'ATOM  %5d %4s %3s %1s %3d    %8.3f%8.3f%8.3f  1.00  1.00\n', ...
            mod(n+1, 10000), "I", "TNL", "X", mod(i+2, 1000), mean_pos);
    end  
end
fclose(fid);


% VMD doesn't read CONECT records, so we write here an executable
% for the VMD tkconsole to add bonds between all atoms and
% a hopping path between the rods
fid = fopen(sprintf("%s.tk", pdbname), 'wt');
for i = nAtoms-1:n
    fprintf(fid, 'topo addbond %d %d\n', i-1, i);
end
fprintf(fid, 'topo addbond %d %d\n', i, 0);
fclose(fid);

end