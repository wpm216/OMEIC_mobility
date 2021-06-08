function K = buildRodlikeRateArray(box, rodLens, betaValue, ...
                                            islandRadius, pdbname)

% Generates the array of inter-particle tunneling rates
% for a pss-rich matrix populated by rod-like particles.

% Parameters:
% -----------
% box = 3-vector of upper corner coordinates (rectangular prism cell)
% nMols = vector of number of molecules each length, corresponding to the
%           index position (e.g. [0, 0, 1, 0, 0, 2] has one trimer and two
%           hexamers).
% betaValue = tunneling attenuation coefficient in inverse angstroms
% islandRadius = radius of islands (angstroms)
% pdbname = (optional) location to write pdb file of system, if desired.

function rod = makeRod(rodLength, box)
    % Generate a rigid rod with random first endpoint in box and random
    % orientation.
    
    % Get first endpoint of rod
    pos1 = rand(1,3) .* box;
    % Get vector of rod backbone (random, normalized vector on unit
    % sphere)
    angle = randn(1,3);
    angle = angle / norm(angle);
    % Get second endpoint of rod
    pos2 = pos1 + angle * rodLength;
    rod = [pos1 pos2];
end

function isValid = isRodValid(rod, lftIsl, rgtIsl, islRad)
    % Check if the rod configuration is valid: 
    %    neither endpoint is resides in either island.
    % (Figuring out if there's no overlap in the middle of the rod feels like an
    % edge case which would be expensive to tease out without little
    % upside).
    p1 = rod(1:3);
    p2 = rod(4:6);
%     p1inside = all(p1 < box) && all(p1 > [0 0 0]); 
%     p2inside = all(p2 < box) && all(p2 > [0 0 0]); 
%     
%     lftIsl = [0, box(2)/2, box(3)/2];
%     rgtIsl = [box(1), box(2)/2, box(3)/2];
    p1overlap = (norm(p1 - lftIsl) < islRad) || (norm(p1 - rgtIsl) < islRad);
    p2overlap = (norm(p2 - lftIsl) < islRad) || (norm(p2 - rgtIsl) < islRad);
  
%     isValid = p1inside && p2inside && ~p1overlap && ~p2overlap;
    isValid = ~p1overlap && ~p2overlap;
end

function dist = rodIslandDistance(rod, islandLoc, islandRadius)
    % The closest point of the rod to the sphere is also the closest to the
    % center of the sphere. So, we find the distance to the center and then
    % subtract the islandRadius to get the distance to the surface.
    
    p1 = rod(1:3);
    p2 = rod(4:6);
    il = islandLoc;
    
    if dot(il-p2, p2-p1) > 0
        % p2 is the closest point to the island
        dist = norm(il - p2);
    elseif dot(il-p1, p1-p2) > 0
        % p1 is the closest point to the island
        dist = norm(il - p1);
    else
        % the closest point is somewhere in the middle of the segment
        dist = norm(cross(p2-p1, islandLoc-p1)/norm(p2-p1)) - islandRadius;
        % If distance is negative, just take its absolute value (it's
        % slightly inside the island; we approximate it as being outside).
        dist = abs(dist);
    end
end

%
% Scripting starts here.
%

leftIslandLoc = [0, box(2)/2, box(3)/2];
rightIslandLoc = [box(1), box(2)/2, box(3)/2];

nRods = sum(rodLens);
rodCoords = zeros(nRods, 6); % coordinates are stored [x1 y1 z1 x2 y2 z2] for each rod
x = 1;
for i=1:length(rodLens)
    % Get rod length
    rodLength = 3.88 * i; % Angstroms
    for j=1:rodLens(i)
        valid = 0;
        while ~valid
            rod = makeRod(rodLength, box);
            valid = isRodValid(rod, leftIslandLoc, rightIslandLoc, islandRadius);
        end
        % If we're here, we generated a valid rod
        rodCoords(x, :) = rod;
        x = x+1;
    end
end

if ~(pdbname=="")
    writeRodlikePDB(pdbname, rodCoords, box, leftIslandLoc, rightIslandLoc)
end

distanceArray = zeros(nRods + 2, nRods + 2);
% Do island-island distances
distanceArray(1, nRods+2) = box(1) - 2*islandRadius;
distanceArray(nRods+2, 1) = box(1) - 2*islandRadius;
% Do rod-island distances
for i=1:nRods
    rod = rodCoords(i, :);
    leftIslandDist = rodIslandDistance(rod, leftIslandLoc, islandRadius);
    rightIslandDist = rodIslandDistance(rod, rightIslandLoc, islandRadius);
    distanceArray(1, i+1) = leftIslandDist;
    distanceArray(i+1, 1) = leftIslandDist;
    distanceArray(nRods+2, i+1) = rightIslandDist;
    distanceArray(i+1, nRods+2) = rightIslandDist;
end
% Do rod-rod distances
for i=1:nRods
    for j=i+1:nRods
        rod1 = rodCoords(i, :);
        rod2 = rodCoords(j, :);
        dist = rodlikeDistanceNonrobust(rod1, rod2);
        distanceArray(i+1, j+1) = dist;
        distanceArray(j+1, i+1) = dist;
    end
end


%                             %
% Calculate the hopping rates %
%                             %

K = exp(-betaValue * distanceArray);
%subtract off-diagonal (set it to zero)
K = K - diag((diag(K))); 
assert(all(all(diag(K) == 0)), "Diagonal elements of K should be zero.")

% Steady-state condition - no accumulation of charge.
% Diagonal elements of islands will be further edited in the sscompute method. 
for i = 1:size(K, 2)
    t=sum(K(:,i));
    K(i,i)=-t;
end

end
