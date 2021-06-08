function K=pointlikeRateArray(X, beta)
% generate the matrix of the rates from coordinates in X
% all points except first and last are connected with
% rate proportional to exp(-d*beta) where d is their distance
% first and last points are connected with the others with rates
% proportional to exp(-(d-r)*beta). r is the radius of source
% and drain particles, beta is the tunneling attenuation factor
% as they don't change much they are set within the function
r=100.;
N=size(X,1);
K=zeros(N,N); % initialize rate

% rate between source and drain and the other particles
for i=[2:N-1] 
    d=norm(X(i,:)-X(1,:));      % distance from particle 1
    K(1,i)=exp(-(d-r)*beta);
    K(i,1)=exp(-(d-r)*beta);    % and corresponding rate
    d=norm(X(i,:)-X(N,:));      % distance from particle N
    K(N,i)=exp(-(d-r)*beta);
    K(i,N)=exp(-(d-r)*beta);    % and corresponding rate
    % no direct source drain rate (may be added)    
end
% rates between particles 2-N
for i=[2:N-1]
    for j=[i+1:N-1]
        d=norm(X(i,:)-X(j,:));  % distance between particles
        K(i,j)=exp(-d*beta);
        K(j,i)=K(i,j);
    end
end

% rate between source and drain 
i=1; j=N; 
d=norm(X(i,:)-X(j,:));  % distance between particles 
K(i,j)=exp(-d*beta);
K(j,i)=K(i,j);
   

% diagonal elements
for i =[1:N]
    t=sum(K(:,i));
    K(i,i)=-t;
end
return


        

    
    
