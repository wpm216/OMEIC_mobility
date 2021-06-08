function a=linterpol(A,x0,dx,y0,dy,xp,yp)
% Linearly interpolate A, a matrix representing a 2D function.

% A(i,j) corresponds to x=x0+(i-1)*dx and y=y0+(j-1)*dy
% it returns the value of the interpolated function at the coordinates xp, yp
% it returns values (extrapolated) if out of bound values are requested
i1=floor((xp-x0)/dx)+1;  
j1=floor((yp-y0)/dy)+1;   % need an action for out of bound

i1=min(i1,size(A,1)-1); % these 4 lines perform extrapolations if interpolation not possible
j1=min(j1,size(A,2)-1);
i1=max(1,i1);
j1=max(1,j1);


X=(xp-(x0+(i1-1)*dx)) / dx ;
Y=(yp-(y0+(j1-1)*dy)) / dy ;

a=A(i1,j1)*(1-X)*(1-Y) + A(i1+1,j1)*X*(1-Y);
a=a + A(i1,j1+1)*(1-X)*Y + A(i1+1,j1+1)*X*Y; 

return 