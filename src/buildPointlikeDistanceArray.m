function X=buildPointlikeDistanceArray(L,n)
% this function generates the geometry of a model of two dots of radius r
% separated by a medium containing redox molecules
% L = length of box
% n = number of point particles (grains are added on top of this)

r=100.0;             % radius of dots in angstrom
X=zeros(n+2,3);
X(1,1)=0; X(1,2)=0; X(1,3)=0;        % first and last special point
X(n+2,1)=L; X(n+2,2)=0; X(n+1,3)=0;
for i=2:n+1
    w=0;
    while w==0
      x=rand()*L;        % generate coordinates in the total volume  
      y=(rand()*5-2.5)*85;
      z=(rand()*5-2.5)*85;
      d1=x^2+y^2+z^2;  % distance from first dot
      d2=(L-x)^2+y^2+z^2; % distance from second dot
      if d1 > r^2 && d2 > r^2 
          X(i,1)=x;
          X(i,2)=y;
          X(i,3)=z;
          w=1;
      end
    end
end
