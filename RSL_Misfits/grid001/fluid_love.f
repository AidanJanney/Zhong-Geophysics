      program fluid_lovenumber
c this subroutine computes the fluid love number for a given density profile

      parameter (nlvl=200)
      implicit double precision (a-h,o-z)
      dimension r(nlvl),rho(nlvl),vp(nlvl),vs(nlvl),vfsq(nlvl)
     1 ,g(nlvl),visc(nlvl),a(2,2),b(2,2),c(2,2),aa(2),x(4)

      pi=3.1415926536d0
      con=pi*6.673d-11
c read in the density profile, we only need r(i) and rho(i) here
c      open (1,file='model.paul',status='unknown')
c      open (1,file='model.wei',status='unknown')
      open (1,file='model.vm5a',status='unknown')
c      open (1, file='model.Lambeck',status='unknown')
c      open (1, file='model.vm7',status='unknown')

      read (1,100) temp,n,(r(i),rho(i),vp(i),vs(i),g(i)
     1 ,visc(i), i=1,n)
 100  format(A20,I10/(5F10.2,E12.4))
 1021 format(5A4,I10/(4F10.0,F10.5,E12.4))

c      do 5 i=1,n
c      if (r(i).lt.3485500. ) then
c      rho(i)=9895.
c      rho(i)=10957.
c      endif
c 5    continue
      do 92 i=1,n
c      print *, r(i), rho(i)
      if (visc(i).ge.1.e30) visc(i)=1.e26
      if (visc(i).lt.1.e-30) visc(i)=0.
   92 continue

c change start - make a uniform core
      do 1191, i=1,n
c      vp(i)=vp(i)*1.e3
c turn it on for structures from citcomsve runs
c      if (r(i) .lt. 3485500.) then
      if (i.lt.67) then
c      rho(i)= 10005.4 !9825.0 !10005.4
c      vp(i) = vp(i)*1.e3
      endif
c make a mantle with constant density and shear modulus
      if (i.ge.67) then
c      rho(i)= 4604.4 !4400.0 !4604.4
c      vp(i)=vp(i)*1.e3
c      vp(i)=2.*sqrt(1.4305e11/rho(i))
c      vs(i)=sqrt(1.4305e11/rho(i))
c      visc(i)=1.e21
      endif

c set viscosity for lith, upper mantle and lower mantle       
      if (i.ge.67.and.i.lt.134) then
      visc(i)=1.0e22 !3.16e21 !0.2e22 !5.e22
c      rho(i)=4400.
c      vs(i)=sqrt(1.4305e11/rho(i))
      endif

      if (i.ge.134.and.i.lt.158) then
      visc(i)=1.0e19 !1.5e20 !1.25e20
c      rho(i)=3300.
c      vs(i)=sqrt(1.4305e11/rho(i))
      endif

      if (i.ge.158) then
      visc(i)=1.e26
c      rho(i)=3300
      endif

1191  continue
c change end





c a test case for Archie's model: I skipped the 1st entry only to make it consistent with prem
      go to 4
      n=6
      r(2)=3503500.
      r(3)=5700.e3
      r(4)=5960.e3
      r(5)=6250.e3
      r(6)=6370.e3
      rho(2)=9900.
      rho(3)=4617.
      rho(4)=4227.
      rho(5)=4047.
      rho(6)=4047.
c      rho(6)=4047.
c      rho(7)=4047.
 4    continue
c compute g(i)
      if (.TRUE.) then
      G(1)=0.
      F=0.
c      DO 7001 I=2,N
c7001  VFSQ=1.e-7
      DO 10 I=2,N
      J=I-1
      RJI=R(I)-R(J)
      if (RJI .eq. 0) go to 15
      DEL=RJI/3.
      X(1)=R(J)
      X(2)=X(1)+DEL
      X(4)=R(I)
      X(3)=X(4)-DEL
      DO 14 K=1,4
      DEL=(X(K)-R(J))/RJI
      RO=RHO(J)+DEL*(RHO(I)-RHO(J))
   14 X(K)=RO*X(K)*X(K)
      F=F+.125*RJI*(X(1)+X(4)+3.*(X(2)+X(3)))
   15 G(I)=4.*CON*F/(R(I)*R(I))
c   15 G(I)=9.8

      print *, i, r(i), rho(i), g(i)
   10 continue
      endif

      if(.false.) then
      do 2 i=2,n
      if (i.eq.2) then
      g(i)=con*4./3.*rho(i)*r(i)
c      g(i)=9.8
      print *, i, r(i), rho(i), g(i)
      else
      g(i)=g(i-1)*(r(i-1)/r(i))**2
     2 +4./3.*con*rho(i)*(r(i)-r(i-1)**3/r(i)**2)
c      g(i)=9.8
      print *, i, r(i), rho(i), g(i)
      endif
 2    continue
      endif
 
c compute the deg-2 body tide fluid love number: fkf
c at the first layer: aa = [i=2] = [C0*(r/a0)**l 0*(a0/r)**(l+1)] = [C0 0], C0 is the unknown
c the equations we use is: b*[i] = a*[i-1] -> [i] = binv*a*[i-1]
      
      l=2
      aa(1)=1.
      aa(2)=0.

      do 3 i=3,n
      b(1,1)=(r(i-1)/r(i))**l
      b(1,2)=(r(i)/r(i-1))**(l+1)
      b(2,1)=l/r(i)*(r(i-1)/r(i))**(l-1)
      b(2,2)=-(l+1)/r(i)*(r(i)/r(i-1))**(l+2)
      call inv2(b,c)

      a(1,1)=1.
      a(1,2)=1.
      a(2,1)=l/r(i-1)+4.*con/g(i-1)*(rho(i)-rho(i-1))
      a(2,2)=-(l+1)/r(i-1)+4.*con/g(i-1)*(rho(i)-rho(i-1))
     
      temp1=a(1,1)*aa(1)+a(1,2)*aa(2)
      aa(2)=a(2,1)*aa(1)+a(2,2)*aa(2)
      aa(1)=temp1

      temp1=c(1,1)*aa(1)+c(1,2)*aa(2)
      aa(2)=c(2,1)*aa(1)+c(2,2)*aa(2)
      aa(1)=temp1
 3    continue
      c0=(2*float(l)+1)/r(n)/((2*l+1)/r(n)*aa(1)
     2 -4*con*rho(n)/g(n)*(aa(1)+aa(2)))
      fkf=c0*(aa(1)+aa(2))-1.

      rhobar=9.82/con*3./4./r(n)
      print *, 'g(n) = ', g(n), 'rhobar = ', rhobar
     2 ,'fkf = ',fkf, ' fkf*(1+delta) = ', fkf*(1+.008)
      end

      subroutine inv2(a,c)
c invert a 2x2 matrix
      implicit real*8 (a-h,o-z)
      dimension a(2,2),c(2,2)
      det=a(1,1)*a(2,2)-a(2,1)*a(1,2)
      do 565 i=1,2
      do 565 j=1,2
 565  c(i,j)=0.
      c(1,1)=a(2,2)/det
      c(1,2)=-a(1,2)/det
      c(2,1)=-a(2,1)/det
      c(2,2)=a(1,1)/det
      return
      end





