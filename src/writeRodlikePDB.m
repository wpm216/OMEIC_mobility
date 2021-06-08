function writeRodlikePDB(pdbname, rodCoords, box, leftIslandLoc, rightIslandLoc)

% Write a pdb file from CT simulations with rod-like PEDOT particles.

fid = fopen(pdbname, 'wt');
resname = "ROD";
n = 1; % number of atoms total in system
fprintf(fid, "TITLE     SIMULATION CELL\n");
fprintf(fid, "CRYST1 %8.3f %8.3f %8.3f %6.2f %6.2f %6.2f P 1           1\n", ...
                box, 90, 90, 90);
            
for i = 1:size(rodCoords, 1)
    
    fprintf(fid, 'ATOM  %5d %4s %3s %1s %3d    %8.3f%8.3f%8.3f  1.00  1.00\n', ...
                    n, "C", resname, "X", mod(i, 1000), rodCoords(i, 1:3));
    fprintf(fid, 'ATOM  %5d %4s %3s %1s %3d    %8.3f%8.3f%8.3f  1.00  1.00\n', ...
                    n+1, "C", resname, "X", mod(i, 1000), rodCoords(i, 4:6));
    n = n+2;
end

% Add islands
fprintf(fid, 'ATOM  %5d %4s %3s %1s %3d    %8.3f%8.3f%8.3f  1.00  1.00\n', ...
            mod(n, 10000), "I", "ISL", "X", mod(i+1, 1000), leftIslandLoc);
fprintf(fid, 'ATOM  %5d %4s %3s %1s %3d    %8.3f%8.3f%8.3f  1.00  1.00\n', ...
            mod(n+1, 10000), "I", "ISL", "X", mod(i+2, 1000), rightIslandLoc);

fprintf(fid, "END\n");
fclose(fid);

% VMD doesn't read CONECT records, so we write here an executable
% for the VMD tkconsole to add bonds between all atoms

fid = fopen(sprintf("%s.tk", pdbname), 'wt');

for i = 1:size(rodCoords, 1)
    fprintf(fid, 'topo addbond %d %d\n', 2*i-2, 2*i-1);
end

fclose(fid);

end