function [S]=pointlikeRate(L,n,nrep,beta)

% Generates the array of inter-particle tunneling rates
% for a pss-rich matrix populated by point-like particles
% and calculates the effective CT rate between them.
% Does this @nrep times.

S=zeros(nrep,1);
for j=1:nrep
    X=buildPointlikeDistanceArray(L,n); % generate the distance array
    K=pointlikeRateArray(X, beta);    % compute the rate between particles
    v=sscompute(K);  % compute the steady state rate 
    S(j)=log(v);
end
