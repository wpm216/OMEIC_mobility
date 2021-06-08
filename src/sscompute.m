function kss=sscompute(K)
%compute the steady state rate from a matrix of first order kinetics
%rate K, assuming the 1 and N and a source and a drain site;
N=size(K,1);
A=zeros(N,1); 
A(1)=-1;
kmax=max(max(K));          % the largest rate there is
K(N,N)=K(N,N)-kmax;     % drain is drained fast
x=linsolve(K,A);
kss=1/x(1);             % this is the steady state rate
return