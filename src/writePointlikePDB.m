function writePointlikePDB(pdbname, pntCoords, box, leftgrainLoc, rightgrainLoc, path)

fid = fopen(pdbname, 'wt');
resname = "PNT";
n = 1; % number of atoms total in system
fprintf(fid, "TITLE     SIMULATION CELL\n");
fprintf(fid, "CRYST1 %8.3f %8.3f %8.3f %6.2f %6.2f %6.2f P 1           1\n", ...
                box, 90, 90, 90);
            
% Add left grain
fprintf(fid, 'ATOM  %5d %4s %3s %1s %3d    %8.3f%8.3f%8.3f  1.00  1.00\n', ...
            mod(n, 10000), "I", "ISL", "X", n, leftgrainLoc);

n = n+1;
for i = 1:size(pntCoords, 1)
    
    fprintf(fid, 'ATOM  %5d %4s %3s %1s %3d    %8.3f%8.3f%8.3f  1.00  1.00\n', ...
                    mod(n, 10000), "C", resname, "X", mod(n, 1000), pntCoords(i, :));
    n = n+1;
end

% Add right grain
fprintf(fid, 'ATOM  %5d %4s %3s %1s %3d    %8.3f%8.3f%8.3f  1.00  1.00\n', ...
            mod(n, 10000), "I", "ISL", "X", mod(n, 1000), rightgrainLoc);


fprintf(fid, "END\n");
fclose(fid);

% VMD doesn't read CONECT records, so we write here an executable
% for the VMD tkconsole to print out a hopping path between the atoms
fid = fopen(sprintf("%s.tk", pdbname), 'wt');
for i = 1:length(path)-1
    fprintf(fid, 'topo addbond %d %d\n', path(i)-1, path(i+1)-1);
end
fclose(fid);

end