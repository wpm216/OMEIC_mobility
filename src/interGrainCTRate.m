function K=interGrainCTRate(X,rate_dist,d0,dd)
% generate the matrix of the rates for a set of particles whose 
% particles are in X. "Particle" 1 and N are planes at x=0 and x=xL
% the rates between particles depends on the distance as in the vector
% rate_dist(i), with point i corresponding to distance d0+i*dd.

N=size(X,1);
K=zeros(N,N); % initialize rate matrix
dmat = zeros(N, N);

% rate between source and drain and the other particles
for i=2:N-1
    d=norm(X(i,1)-X(1,1));      % x distance from particle 1 
    d=max(d0+dd,d);                   % distance forced to be larger than d0
   
      K(1,i)= rate_dist( round((d-d0)/dd) ); % and corresponding rate
      K(1,i)=K(1,i)*20^(rand);
      K(i,1)=K(1,i);              
    
    d=norm(X(i,1)-X(N,1));      % distance from particle N
    d=max(d0+dd,d);                   % distance forced to be larger than d0
    
    
      K(N,i)= rate_dist( round((d-d0)/dd) ); % and corresponding rate
      K(i,N)=K(N,i);              
    % no direct source drain rate    
end

% rates between particles 2-N
for i=2:N-1
    for j=i+1:N-1
        d=norm(X(i,:)-X(j,:));  % distance between particles 
        d=max(d0+dd,d);            % distance forced to be larger than d0
        dmat(i, j) = d;
        dmat(j, i) = d;
           K(i,j)= rate_dist( round((d-d0)/dd) );
           K(j,i)= K(i,j);
        
    end
end

% diagonal elements
for i=1:N
    t=sum(K(:,i));
    K(i,i)=-t;
end
return


        

    
    
