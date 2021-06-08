      program daslkj
      integer NUMBER(4)
c this version make a box of decreasing concentration
      real L, ll
      real C(2048,3)
c spheres are quite rigid with this DD parameters
      DD=35.

c with this selection of of parameters a box of 864
c particles is formed and places in fort 10
      is=7
      L=2263.
      ll=L/is

      n=1574
      do i=1,n
       read (99,*)(C(i,j),j=1,3)
c      C(i,1)=rand()*L
c      C(i,2)=rand()*L
c      C(i,3)=rand()*L
      enddo

c step size is 0.1*ll
      step=0.100*ll

c perform 10^6 MC moves for equilibration 
      write (*,*)' performing equilibration '
      icr=0.
c     MCmoves=200000
c     MCmoves=0
      do imc=1,MCmoves 
       k=int(rand()*n)+1 
       X=C(k,1)+(rand()-0.5)*step
       Y=C(k,2)+(rand()-0.5)*step
       Z=C(k,3)+(rand()-0.5)*step
c      do this only if you don't go out of the box
c------
       if (X.le.L.and.Y.le.L.and.Z.le.L.and.X.gt.0.and.
     x     Y.gt.0..and.Z.gt.0.) then

c       compute the change in energy due to this move
        DE=0.
        do i=1,n
         if (i.ne.k) then
            R2=DI2(C(i,1),C(k,1),C(i,2),C(k,2),C(i,3),C(k,3),L)  
            E2=Energy(R2,DD*imc*1.0E-6)
            R2=DI2(C(i,1),X,C(i,2),Y,C(i,3),Z,L)  
            E1=Energy(R2,DD*imc*1.0E-6)
            DE=DE+E1-E2
         endif 
        enddo

c       accept the move with MC rule
c       write (11,*)(C(k,j),j=1,3),'->',X,Y,Z,DE
        if (min(1.0,exp(-DE)).gt.rand()) then
         C(k,1)=X
         C(k,2)=Y
         C(k,3)=Z
         icr=icr+1
C        write (11,*)'accept'
c       else
c        write (11,*)'REJECT'
        endif
       
       endif
c------
      enddo

     
      do i=1,n 
       write (12,10)'  ',(C(i,j),j=1,3)
      enddo


      NUMBER(1)=1574
      NUMBER(2)=1058
      NUMBER(3)=717
      NUMBER(4)=451
cCCCC study 9 systems reducing the number of particles 
c     by 94 particles every time (except the first one) 
cHHHHHHHHHHHHHHHH
      do in=1,4 
        n=NUMBER(in) 
        write (*,*)n
cHHHHHHHHHHHHHHHH
 

c       PERFORM several MC simulations on this system
cXXXXXXX
        do ik=0,24
cXXXXXXX

        icr=0.
        MCmoves=100000
c       MCmoves=15000
        do imc=1,MCmoves 
         k=int(rand()*n)+1 
         X=C(k,1)+(rand()-0.5)*step
         Y=C(k,2)+(rand()-0.5)*step
         Z=C(k,3)+(rand()-0.5)*step
c        do this only if you don't go out of the box
c--------
         if (X.le.L.and.Y.le.L.and.Z.le.L.and.X.gt.0.and.
     x       Y.gt.0..and.Z.gt.0.) then
     
c         compute the change in energy due to this move
          DE=0.
          do i=1,n
           if (i.ne.k) then
             R2=DI2(C(i,1),C(k,1),C(i,2),C(k,2),C(i,3),C(k,3),L)  
             E2=Energy(R2,DD)
             R2=DI2(C(i,1),X,C(i,2),Y,C(i,3),Z,L)  
             E1=Energy(R2,DD)
             DE=DE+E1-E2
           endif 
          enddo
     
c         accept the move with MC rule
c         write (11,*)(C(k,j),j=1,3),'->',X,Y,Z,DE
          if (min(1.0,exp(-DE)).gt.rand()) then
           C(k,1)=X
           C(k,2)=Y
           C(k,3)=Z
           icr=icr+1
          endif
         
         endif
c--------
        enddo 
c       compute total energy
        ETOT=0.
        do i=1,n
         do j=i+1,n
          R2=DI2(C(i,1),C(j,1),C(i,2),C(j,2),C(i,3),C(j,3),L)  
          ETOT=ETOT+Energy(R2,DD)
         enddo
        enddo
        write (*,*)n,ik,icr,MCmoves,ETOT
     
c       output the configuration at that point
        do i=1,n 
         write (1000+in*100+ik,10)'  ',(C(i,j),j=1,3)
        enddo

cXXXXXXX
        enddo     
cXXXXXXX


 

cHHHHHHHHHHHHHHHH
      enddo       
cHHHHHHHHHHHHHHHH



   10 format (A2,3F12.4)
      stop 
      end

      real function Energy(R2,DD)
      r=sqrt(R2)
      Energy= 10 * exp (-(r-220.)/DD)
      return
      end  


      real function DI2(X1,X2,Y1,Y2,Z1,Z2,S)
      Dx=X1-X2
      if (Dx.gt.S/2)  Dx=Dx-S
      if (Dx.lt.-S/2) Dx=Dx+S
      Dy=Y1-Y2
      if (Dy.gt.S/2)  Dy=Dy-S
      if (Dy.lt.-S/2) Dy=Dy+S
      Dz=Z1-Z2
      if (Dz.gt.S/2)  Dz=Dz-S
      if (Dz.lt.-S/2) Dz=Dz+S
      DI2=Dx**2+Dy**2+Dz**2
      return
      end




      
