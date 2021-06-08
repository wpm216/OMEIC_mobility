function a=getCell(idx, n, d)
% Get the cell number, a, of a given index in the distance array.
% Assumes the distance array is ordered first by increasing x,
%  then by increasing y, then by increasing z. 
%  I.e., if we draw an analogy with an odometer, the x-index is in
%  the 'ones' place, y is in the 'tens', and z is in the 'hundreds.'
%
% idx : index in the distance array (1-origin). An N-dimensional array 
%        will have N indices associated with each element.
% n   : number of molecules in one copy of the system.
% d   : number of periodic images in each direction in the box, 
%        passed as a list ([nx,ny,nz]).
% a   : cell number, returned as a list ([ax, ay, az]).
%        for each ai, 1 <= ai <= ni.

% assert((1 <= idx) && (idx <= n * prod(d)), 'Index out of bounds.');

idx0 = idx - 1; % covert to 0-origin for modulo operations
cellIdx = floor(idx0/n); % each cell has n molecules; this gets the cell index.
az = floor(cellIdx/(d(1) * d(2))); 
rem = cellIdx - az * d(1) * d(2);
ay = floor(rem/d(1));
rem = rem - ay * d(1);
ax = rem;

a = [ax+1, ay+1, az+1]; % convert from 0-origin to 1-origin
% assert(all(a <= d), 'Error in index to cell mapping.');
end

