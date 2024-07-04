      program stokes
c
c   FOR ICE5 only.  Because this program assume there are 96 time
c indices
c
c  finds time rates of change of Stokes' coefficients due to
c   PGR.  Does this for any viscosity model you want (read
c   in from file model.(imodel), with imodel read in from
c in.stokes.ice5).
c   And, for any version of the ice-5 model (reads in the Ylm
coefficients of the
c   ice model at any epoch, from file "name" - read in from
c in.stokes.ice5).
c  It can do this for lots of input viscosity/ice models - the
c   number of models considered is nauto, read in from in.stokes.ice5.
c
      parameter(lmax=101,nmodels=36)
      character*2 numb(nmodels)
      character*40 name,name1,namev
      double precision fkf

      open (88,file='in.stokes.ice5.PW',status='unknown')

c  nauto is the number of viscosity/ice models considered in this run.
c  new notes: here nauto is always set to be 1, because the love number is solved in a separate routine
c  and we don't use any viscosity models in this program.
      read (88,*) nauto

      do 10 n=1,nauto

c  do the sea level equation iteration if use_sle = 1, and num_iter is the number of iterations
      read (88,*) use_sle, num_iter, polar_wander
      print *,use_sle, 'number of iteration = ', num_iter
     2 ,' polar wander = ',polar_wander

c  name is the file with the ylm coeffs of the ice model (out.ylm...)
c  lmaxylm is the maximum no of l's in this input ice file
      read (88,*) name,lmaxylm
      print *,'file39',name
      open (39,file=name,
     2  status='unknown')

c  lmaxb is the maximum no of l's we will include
      read (88,*) lmaxb
      print *,lmaxb
      if (lmaxb.gt.lmax) then
      print *,'error #16',lmax,lmaxb
      stop
      endif

c  name1 is the output file containing the Stokes' coeff rates
      read (88,*) name1
      print *,name1
      open (13,file=name1,status='unknown')

c read in the fluid love number
      read (88,*) fkf

      nepochs=50 !55 for W12 !114 for Lambeck global  !50 for ice6g
      do 1 nobs=1,nepochs
C      do 1 nobs=nepochs,nepochs
      call love(lmaxb,nobs)
      if (use_sle.eq.1) then
      call solve_sle(lmaxb,lmaxylm,num_iter,fkf,polar_wander)
      else
c  turn off the sle iteration, read in the original ice load
      close (39)
C      open (49,file=name,status='unknown')
      open(49,file='final_load.ice6g',status='unknown')
      endif
      call cuplift (lmaxb,lmaxylm,use_sle,fkf,polar_wander,nobs)
      close(13)
1     continue

      close(50)
10    continue

      stop
      end

      subroutine love(lmaxb,nobs)
c
c  this subroutine uses the Love nos found by the compressible code, to
c   find effective love numbers for the ice-5 epochs
c  Present time is 0
c
c new note - this subroutine determines the precision of the final result, here we use e17.10 to save pot and u

      parameter(lmax=101,nepochs=50)
      parameter (ncoll=200)
      implicit doubleprecision (a-h,o-z)
      common /cetm/ rh(0:lmax,ncoll),rk(0:lmax,ncoll),rl(0:lmax,ncoll),
     2  s(0:lmax,ncoll)
     2  ,eh(0:lmax),ek(0:lmax),el(0:lmax),nn(0:lmax)
      common /renorm/ stokyr
      dimension 
     2   u(nepochs+1,0:lmax),
     2   pot(nepochs+1,0:lmax)
     2  ,u_nr(nepochs+1,0:lmax), pot_nr(nepochs+1,0:lmax)
     2  ,v(nepochs+1,0:lmax),v_nr(nepochs+1,0:lmax)
     3  ,tepoch(0:nepochs)
c save non-rate geoid and uplift in u_nr and pot_nr for testing
      open (22,file='out.love',status='unknown')

      pi=3.14159265
      dtr=pi/180.
c
      stokyr=1000.*365.25*86400.
      stoyr=365.25*86400.
c
c  the following finds the epochs of the model
      open (66,file='years.reverse',status='unknown')
      do 1542 n=1,nepochs
      read (66,*) tepoch(n)
      tepoch(n)=-tepoch(n)*stokyr
c      tepoch(n)=-tepoch(n)*1.e3
1542  continue
      close(66)

      tmp = tepoch(nobs)
      do 1543 n=1,nepochs
      tepoch(n)=tepoch(n)-tmp
1543  continue
c
c  the following is a relic from ice4.  it assumes there is a
c   start time for the de-glaciation, 90 kyr before the start time
c   of ice5.  The ice thickness change will be the difference between
c   present day and the start of ice5.  So, the ice change will be
c   small, and will be spread out between about 210 and 122 kyrs ago.
c   So I assume that including it will notmake much difference.  And
c   if I include it, I don't have to figure out just how the program
c   would have to be changed if it were not included.
c      tepoch(0)=tepoch(1)-90.*stokyr
c      tepoch(0)=tepoch(1)-90000.

      call read_love (lmaxb)


c change start - here I include lmaxb+1, to store the rotation love numbers for degree (2,1) term
      do 934 mm=1,nobs+1
      do 10 n=0,lmaxb+1
      pot(mm,n)=0.
      u(mm,n)=0.
      v(mm,n)=0.
      pot_nr(mm,n)=0.
      u_nr(mm,n)=0.
      v_nr(mm,n)=0.
      if (nn(n).gt.0) then
c
      ntop=nn(n)
      do 11 k=1,ntop
      if (mm.eq.1) then
      factor=-exp(tepoch(mm)*s(n,k))
      else if (mm.eq.nobs+1) then
      factor=-1.
      else
      ddum=tepoch(mm)*s(n,k)
      ddum1=tepoch(mm-1)*s(n,k)
      factor=
     2  (exp(ddum)-exp(ddum1))/s(n,k)
      endif
      if (nobs.eq.1) factor=0.
      pot(mm,n)=pot(mm,n)+rk(n,k)*factor
      u(mm,n)=u(mm,n)+rh(n,k)*factor
      v(mm,n)=v(mm,n)+rl(n,k)*factor
      pot_nr(mm,n)=pot_nr(mm,n)+rk(n,k)*factor/(-s(n,k))
      u_nr(mm,n)=u_nr(mm,n)+rh(n,k)*factor/(-s(n,k))
      v_nr(mm,n)=v_nr(mm,n)+rl(n,k)*factor/(-s(n,k))
11    continue
      endif
c  the following converts to yr**-1
      if (mm.ne.nobs+1 .and. mm.ne.1) then
      pot_nr(mm,n)=pot_nr(mm,n)/(tepoch(mm-1)-tepoch(mm))
      u_nr(mm,n)=u_nr(mm,n)/(tepoch(mm-1)-tepoch(mm))
      v_nr(mm,n)=v_nr(mm,n)/(tepoch(mm-1)-tepoch(mm))
      pot(mm,n)=pot(mm,n)/(tepoch(mm-1)-tepoch(mm))*stoyr
      u(mm,n)=u(mm,n)/(tepoch(mm-1)-tepoch(mm))*stoyr
      v(mm,n)=v(mm,n)/(tepoch(mm-1)-tepoch(mm))*stoyr
      endif
10    continue
934   continue
      do 20 n=0,lmaxb+1
      do 20 mm=1,nobs+1
      write (22,320) n,mm,pot(mm,n),u(mm,n),v(mm,n)
     2 ,pot_nr(mm,n),u_nr(mm,n),v_nr(mm,n)
320   format(i3,1x,i3,1x,e17.10,1x,e17.10,1x,e17.10,1x,e17.10
     2 ,1x,e17.10,1x,e17.10)
20    continue
      close (22)

      return
      end

      subroutine read_love (lmaxb)
      parameter (lmax=101)
      parameter (ncoll=200)
      implicit real*8(a-h,o-z)
      common /cetm/ rh(0:lmax,ncoll),rk(0:lmax,ncoll),
     2 rl(0:lmax,ncoll),s(0:lmax,ncoll)
     2 ,eh(0:lmax),ek(0:lmax),el(0:lmax),nn(0:lmax)

      open (999,file='out.dahlen.love',status='unknown')
      open (994,file='out.dahlen.coll',status='unknown')

      nn(0)=0

      do 998 l=1,lmaxb
      read (994,*)
      read (994,*) ll,ns
      if (ll.ne.l) then
      print *, 'l error! #1',l,ll
      stop
      endif
c      if (ns.ne.20) then
c      print *, 'ns error! #1',ns
c      stop
c      endif
      nn(l)=ns
      read (999,997) ll,eh(l),ek(l),el(l)
997   format(12X,I3,17X,F12.8,17X,F12.8,17X,F12.8,27X)
      if (ll.ne.l) then
      print *, 'l error #2',l,ll
      endif
      do 996 n=1,ns
      read (999,*)
      read (994,995) m,s(l,n),rh(l,n),rk(l,n),rl(l,n)
995   format(I5,18X,4(1X,e17.10))
c make s go to 1/yrs
c      s(l,n)=s(l,n)*86400.*365.25
      if (m.ne.n) then
      print *, 'n error #2',n,m
      endif
      rh(l,n)=-rh(l,n)
      rl(l,n)=-rl(l,n)
996   continue
c change eh
c      eh(l)=-eh(l)
998   continue
      close(999)
      close(994)

c degree (2,1) term - read in rotation love numbers
c      go to 2000
      open (999,file='out.btide',status='unknown')
      open (994,file='out.rot',status='unknown')
      read (999,997) ll,eh(lmaxb+1),ek(lmaxb+1),el(lmaxb+1)
      if (ll.ne.2) then
      print *, 'not degree 2 term!',ll
      stop
      endif
      read (994,*)
      read (994,*) ll,ns
      nn(lmaxb+1)=ns
      if (ll.ne.2) then
      print *, 'not degree 2 term!',ll
      stop
      endif
      do 1999 n=1,ns
      read (994,995) m,s(lmaxb+1,n),rh(lmaxb+1,n),rk(lmaxb+1,n)
     2 ,rl(lmaxb+1,n)

      if (m.ne.n) then
      print *, 'n error #3',n,m
      stop
      endif
1999   continue
      close(999)
      close(994)
2000  continue

      open(993,file='out.elastic.tmp',status='unknown')
      do 992 ll=0,lmaxb+1
      write (993,997) ll,eh(ll),ek(ll),el(ll)
992   continue
      close(993)
      return
      end

      

      subroutine cuplift(lmaxb,lmaxylm,use_sle,fkf,polar_wander,nobs)
c  This program finds the geoid rates, using the effective 
c   ICE-3g Love numbers found by subroutine love
c     nth and nphi are the number of theta and phi values in
c     the grid.  lmax is the maximum no. of l values computed
      parameter (nth=181,nphi=361,
     2           lmaxa=101,
     2           lmax=101,ltot=(lmax+1)**2,nepochs=50)
      character*40 filen
      dimension th(nth),phi(nphi),alm(2,0:lmaxa,0:lmaxa)
     4  ,u(nepochs+1,0:lmax)
     4  ,pot(nepochs+1,0:lmax)
     4  ,u_nr(nepochs+1,0:lmax),pot_nr(nepochs+1,0:lmax)
     4  ,v(nepochs+1,0:lmax),v_nr(nepochs+1,0:lmax)
     4  ,geoid(2,0:lmax,0:lmax),uplift(2,0:lmax,0:lmax)
     4  ,geoid_nr(2,0:lmax,0:lmax),uplift_nr(2,0:lmax,0:lmax)
     4  ,hori(2,0:lmax,0:lmax),hori_nr(2,0:lmax,0:lmax)
     4  ,eh(0:lmax),ek(0:lmax),el(0:lmax)
c      common /model/ rmu(10),rho(10),a(10)
c     2 ,smu(10),visc(10),rhs(3),gt(10),nlvl
      double precision fkf

      open(3993,file='out.elastic.tmp',status='unknown')
      do 3992 n=0,lmaxb+1
      read (3993,3997) nl,eh(n),ek(n),el(n)
3992   continue
3997   format(12X,I3,17X,F12.8,17X,F12.8,17X,F12.8,27X)
      close(3993)


c change start - check geoid_nr and uplift_nr
c      open (41,file='non_rate',status='unknown')
      open (42,file='check_sum',status='unknown')
      open (40,file='out.love',status='unknown')

      if (use_sle.eq.1) then
      open (49,file='final_load.ice6g',status='unknown')
      print *, 'loading the final dynamic ocean load'
      endif
c change start - here I include lmaxb+1 to store the rotation love numbers
      do 340 n=0,lmaxb+1
      do 340 mm=1,nobs+1
      read (40,320) nz,mmz,pot(mm,n),u(mm,n),v(mm,n)
     2 ,pot_nr(mm,n),u_nr(mm,n),v_nr(mm,n)
320   format(i3,1x,i3,1x,e17.10,1x,e17.10,1x,e17.10,1x,e17.10
     2 ,1x,e17.10,1x,e17.10)
      if ((nz.ne.n).or.(mmz.ne.mm)) then
      print *,'error #51:',nz,n,mmz,mm
      stop
      endif
340   continue

      pi=3.14159265
      fpi=4.*pi
      dtr=pi/180.
c  dth and dphi are the theta and phi intervals (between grid pts)
c     th and phi are arrays with the values of theta and phi at
c     the gridpts
      dth=180./float(nth-1)*dtr
      dphi=360./float(nphi-1)*dtr
      do 10 i=1,nth
      th(i)=float(i-1)*dth
10    continue
      do 11 j=1,nphi
      phi(j)=float(j-1)*dphi-180.*dtr
11    continue

      coefu=4.*pi*6.673e-8*6.371e8/984.0
c      coefu=4.*pi*6.673e-8*6.370e8/980.0
c citcom.canada: 951.154150978950
c951.195851596772 citcom's vm2
c 982.390126455720 prem
c982.872037937798
c 950.17,982.872037937798
c change g value for different model!
c      do 869 l=0,lmaxb
c      do 869 m=0,l
c      do 869 j=1,2
c      geoid(j,l,m)=0.
c      uplift(j,l,m)=0.
c      geoid_nr(j,l,m)=0.
c      uplift_nr(j,l,m)=0.
c869   continue
c

      do 869 l=0,lmaxb
      do 869 m=0,l
      do 869 j=1,2
      geoid(j,l,m)=0.
      uplift(j,l,m)=0.
      hori(j,l,m)=0.
      geoid_nr(j,l,m)=0.
      uplift_nr(j,l,m)=0.
      hori_nr(j,l,m)=0.
869   continue

      write(filen,3610) 'PW.ice6g.',nobs+100
3610  format(a,i3)
      print *, filen
      open(41,file=filen,status='unknown')
      close(49)
c close the ice load file, and then read from the 1st line again.
      if (use_sle.eq.1) then
      open (49,file='final_load.ice6g',status='unknown')
      else
c      open(49,file='ice6g+ocn.lmax100.txt',status='unknown')
      open (49,file='final_load.ice6g',status='unknown')
      endif

      do 569 nint=1,nobs
      do 950 i=1,3
950   read (49,150)
150   format(1x)
c
c  Now, read in the ice mass Ylm coefs
c  
      do 312 l=0,lmaxylm
      do 312 m=0,l
      read (49,2903) l1,m1,alm(1,l,m),alm(2,l,m)
2903   format(3x,i3,3x,i3,1x,6(e14.5,2x))
      if ((l1.ne.l).or.(m1.ne.m)) then
      print *, 'error #2',l1,l,m1,m
      stop
      endif
      if (l.le.lmaxb) then
      alm(1,l,m)=alm(1,l,m)*100.*.9174
      alm(2,l,m)=alm(2,l,m)*100.*.9174
c     geruo's ice density .9174
      endif
312   continue

c  Now, find the geoid coefs
c change start - here I include lmaxb+1
      do 30 l=0,lmaxb+1
      if (polar_wander.lt.0.5.and.l.eq.lmaxb+1) then
      go to 30
      endif
      if (nint.ne.nobs) then
c  in the following, I use two love numbers, because an individual
c    love number, u(nint), should really be multiplied by: (the load at
c    the start of interval nint) minus (the load at the end of the
c    interval).  u(nint) is bounded by alm(nint-1) and alm(nint)
      ee=-pot(nint,l)+pot(nint+1,l)
      eu=-u(nint,l)+u(nint+1,l)
      ev=-v(nint,l)+v(nint+1,l)
      ee_nr=-pot_nr(nint,l)+pot_nr(nint+1,l)
      eu_nr=-u_nr(nint,l)+u_nr(nint+1,l)
      ev_nr=-v_nr(nint,l)+v_nr(nint+1,l)
      else
c  the load at the last time (nepochs): I change it so it is consistent
c  with the python code
c      ee=-pot(nint,l)
c      eu=-u(nint,l)
c      ev=-v(nint,l)
c  add present day elastic response (the last epoch)
      if (l.ne.lmaxb+1) then
      ee_nr=-pot_nr(nint,l)+pot_nr(nint+1,l)+1.+ek(l)
      eu_nr=-u_nr(nint,l)+u_nr(nint+1,l)+eh(l)
      ev_nr=-v_nr(nint,l)+v_nr(nint+1,l)+el(l)
      ee=-pot(nint,l)+(1.+ek(l))/500.
      eu=-u(nint,l)+eh(l)/500.
      ev=-v(nint,l)+el(l)/500.
      else
      ee_nr=-pot_nr(nint,l)+pot_nr(nint+1,l)+ek(l)
      eu_nr=-u_nr(nint,l)+u_nr(nint+1,l)+eh(l)
      ev_nr=-v_nr(nint,l)+v_nr(nint+1,l)+el(l)
      ee=-pot(nint,l)+(ek(l))/500.
      eu=-u(nint,l)+eh(l)/500.
      ev=-v(nint,l)+el(l)/500.

c degree (2,1) is a little bit different only because the sign of the input data is different, and also we don't need to add 1 for the geoid correction term
      endif
      endif
      if (nint.eq.nobs-1) then
      if (l.ne.lmaxb+1) then
      ee=ee-(1.+ek(l))/500.
      eu=eu-eh(l)/500.
      ev=ev-el(l)/500.
      else
      ee=ee-(ek(l))/500.
      eu=eu-eh(l)/500.
      ev=ev-el(l)/500.
      endif
      endif
c 500. at above is 500 years

c change start - update the (2,1) term using rotation love numbers
      if (l.eq.lmaxb+1) then
      do 931 j=1,2
      geoid(j,2,1)=geoid(j,2,1)+
     2 alm(j,2,1)*ee/float(2*2+1)*coefu
      uplift(j,2,1)=uplift(j,2,1)+
     2 alm(j,2,1)*eu/float(2*2+1)*coefu
      hori(j,2,1)  =  hori(j,2,1)+
     2 alm(j,2,1)*ev/float(2*2+1)*coefu

      geoid_nr(j,2,1)=geoid_nr(j,2,1)+
     2 alm(j,2,1)*ee_nr/float(2*2+1)*coefu
      uplift_nr(j,2,1)=uplift_nr(j,2,1)+
     2 alm(j,2,1)*eu_nr/float(2*2+1)*coefu
      hori_nr(j,2,1)  =  hori_nr(j,2,1)+
     2 alm(j,2,1)*ev_nr/float(2*2+1)*coefu

931   continue
      write(42,6699) nint, alm(1,2,1)*ee_nr/float(2*2+1)*coefu,
     2 alm(2,2,1)*ee_nr/float(2*2+1)*coefu
6699  format(3x,i3,2(3x,e17.10))
      go to 30
      endif
c change end
      do 333 m=0,l
      do 31 j=1,2
      geoid(j,l,m)=geoid(j,l,m)+
     2  alm(j,l,m)*ee/float(l*2+1)*coefu
      uplift(j,l,m)=uplift(j,l,m)+
     2  alm(j,l,m)*eu/float(l*2+1)*coefu
      hori(j,l,m)  =  hori(j,l,m)+
     2  alm(j,l,m)*ev/float(l*2+1)*coefu

      geoid_nr(j,l,m)=geoid_nr(j,l,m)+
     2   alm(j,l,m)*ee_nr/float(l*2+1)*coefu
      uplift_nr(j,l,m)=uplift_nr(j,l,m)+
     2   alm(j,l,m)*eu_nr/float(l*2+1)*coefu
      hori_nr(j,l,m)  =  hori_nr(j,l,m)+
     2   alm(j,l,m)*ev_nr/float(l*2+1)*coefu

31    continue
333    continue
30    continue
c  
569   continue
c      close(49) ! modified by kaixuan 7/10/2019

c add in the centrif potential
      if (polar_wander.gt.0.5) then
      geoid_nr(1,2,1)=geoid_nr(1,2,1)*(1.+1./fkf)
      geoid_nr(2,2,1)=geoid_nr(2,2,1)*(1.+1./fkf)
      endif


      do 3612 l=0,lmaxb
      do 3612 m=0,l
      write (41,2903) l,m,geoid_nr(1,l,m),geoid_nr(2,l,m)
     2 ,uplift_nr(1,l,m),uplift_nr(2,l,m)
     2 ,hori_nr(1,l,m),hori_nr(2,l,m)
3612  continue
      close(41)


c new note - the convention is consistent with the one from geo.py
c      coef0=1./6.371e8/sqrt(4.*3.14159)
c      coefm=1./6.371e8/sqrt(2.*3.14159)
c      do 7899 l=0,lmaxb
c      do 7899 m=0,l
c      if (m.eq.0) then
c      clm_nr=geoid_nr(1,l,m)*coef0
c      slm_nr=geoid_nr(2,l,m)*coef0
c      clm1_nr=uplift_nr(1,l,m)*coef0
c      slm1_nr=uplift_nr(2,l,m)*coef0
c      clm2_nr=hori_nr(1,l,m)*coef0
c      slm2_nr=hori_nr(2,l,m)*coef0

c      write (41,2903) l,m,clm_nr,slm_nr,clm1_nr,slm1_nr,clm2_nr,slm2_nr
c      else
c      clm_nr=geoid_nr(1,l,m)*coefm*(-1.)**float(m)
c      slm_nr=-geoid_nr(2,l,m)*coefm*(-1.)**float(m)
c      clm1_nr=uplift_nr(1,l,m)*coefm*(-1.)**float(m)
c      slm1_nr=-uplift_nr(2,l,m)*coefm*(-1.)**float(m)
c      clm2_nr=hori_nr(1,l,m)*coefm*(-1.)**float(m)
c      slm2_nr=-hori_nr(2,l,m)*coefm*(-1.)**float(m)
c      write (41,2903) l,m,clm_nr,slm_nr,clm1_nr,slm1_nr,clm2_nr,slm2_nr
c      endif
cc change start - check the geoid and uplift (not rate)
cc      write (41,2903) l,m,geoid_nr(1,l,m),geoid_nr(2,l,m)
cc     2 ,uplift_nr(1,l,m),uplift_nr(2,l,m)
c7899   continue
c      close(41)



c      print *, alm(1,2,0)*(1.+ek(2))/float(2*2+1)*coefu,
c     2  alm(1,2,0)*eh(2)/float(2*2+1)*coefu
c new note - the convention is consistent with the one from geo.py
      coef0=1./6.371e8/sqrt(4.*3.14159)
      coefm=1./6.371e8/sqrt(2.*3.14159)
      do 789 l=0,lmaxb
      do 789 m=0,l
      if (m.eq.0) then
      clm=geoid(1,l,m)*coef0
      slm=geoid(2,l,m)*coef0
      clm1=uplift(1,l,m)*coef0
      slm1=uplift(2,l,m)*coef0
      clm2=hori(1,l,m)*coef0
      slm2=hori(2,l,m)*coef0

      write (13,2903) l,m,clm,slm,clm1,slm1,clm2,slm2
      else
      clm=geoid(1,l,m)*coefm*(-1.)**float(m)
      slm=-geoid(2,l,m)*coefm*(-1.)**float(m)
      clm1=uplift(1,l,m)*coefm*(-1.)**float(m)
      slm1=-uplift(2,l,m)*coefm*(-1.)**float(m)
      clm2=hori(1,l,m)*coefm*(-1.)**float(m)
      slm2=-hori(2,l,m)*coefm*(-1.)**float(m)
      write (13,2903) l,m,clm,slm,clm1,slm1,clm2,slm2
      endif
c change start - check the geoid and uplift (not rate)
c      write (41,2903) l,m,geoid_nr(1,l,m),geoid_nr(2,l,m)
c     2 ,uplift_nr(1,l,m),uplift_nr(2,l,m)
789   continue
c  
      close (40)
      return
      end



      subroutine solve_sle(lmaxb,lmaxylm,num_iter,fkf,polar_wander)
c
c  this subroutine uses the Love nos found by the compressible code, to
c   find the ocean load
c  the output is a load file containing (ice_load+ocean_load)
c  Present time is 0
c

      parameter(nlat=181,nlon=361,lmaxa=101,lmax=101,nepochs=50)
      parameter (ncoll=200)
      implicit doubleprecision (a-h,o-z)
      common /cetm/ rh(0:lmax,ncoll),rk(0:lmax,ncoll),rl(0:lmax,ncoll),
     2  s(0:lmax,ncoll)
     2  ,eh(0:lmax),ek(0:lmax),el(0:lmax),nn(0:lmax)
      common /renorm/ stokyr
      common /ijlm/ rNU_ij(1:nlat,1:nlon),rNU_lm(2,0:lmax,0:lmax)
     2 ,rLO_lm(2,0:lmax,0:lmax),rI
      dimension
     2  u_nr(nepochs+1,0:lmax), pot_nr(nepochs+1,0:lmax),Ao(1:nepochs)
     2  ,v_nr(nepochs+1,0:lmax)
     3  ,tepoch(0:nepochs),almt(2,0:lmaxa,0:lmaxa,1:nepochs)
     4  ,geoid_nr(2,0:lmax,0:lmax),uplift_nr(2,0:lmax,0:lmax)
     4  ,hori_nr(2,0:lmax,0:lmax)
     5  ,th(nlat),phi(nlon)
     6  ,ocn_xyt(1:nlat,1:nlon,1:nepochs)
     7  ,almt0(2,0:lmaxa,0:lmaxa,1:nepochs)


      pi=3.14159265
      fpi=4.*pi
      dtr=pi/180.
c      fkf=0.9422d0
c      fkf=0.940954d0


      stokyr=1000.*365.25*86400.
      stoyr=365.25*86400.

      coefu=4.*pi*6.673e-8*6.371e8/984.0
c citcom.canada: 951.154150978950
c950.7 incomp
c951.195851596772
c982.390126455720
c982.872037937798
c 950.17, 982.872037937798
c change g value for different model!

c check geoid and uplift file
      open (41,file='geo+uplift',status='unknown')
      open (42,file='chk_sum',status='unknown')

c  the following finds the epochs of the model
      open (66,file='years.reverse',status='unknown')
      do 1542 n=1,nepochs
      read (66,*) tepoch(n)
      tepoch(n)=-tepoch(n)*stokyr
1542  continue
      close(66)
c      tepoch(0)=tepoch(1)-90*stokyr     !changed by kaixuan 2/3/2020

      call read_love (lmaxb)
      ek(0)=-1.
c read in the ice+(static ocean load)
      do 5690 nint=1,nepochs
      do 950 i=1,3
950   read (39,150)
150   format(1x)
      do 312 l=0,lmaxylm
      do 312 m=0,l
      read (39,2903) l1,m1,almt(1,l,m,nint),almt(2,l,m,nint)
2903   format(3x,i3,3x,i3,1x,6(e14.5,2x))
      if ((l1.ne.l).or.(m1.ne.m)) then
      print *, 'error #2',l1,l,m1,m
      stop
      endif
      if (l.le.lmaxb) then
      almt(1,l,m,nint)=almt(1,l,m,nint)*100.*.9174
      almt(2,l,m,nint)=almt(2,l,m,nint)*100.*.9174
      almt0(1,l,m,nint)=almt(1,l,m,nint)
      almt0(2,l,m,nint)=almt(2,l,m,nint)
      endif
c     geruo's ice density .9174
312   continue
5690  continue
      close(39)

c read in the land/ocean mask and ocean area for each epoch
      open(44,file='ocn6g_xyt.txt',status='unknown')
      do 4 k=1,nepochs
      read(44,*)
      do 4 i=1,nlat
      do 4 j=1,nlon
      read(44,*) tlat,tlon,ocn_xyt(i,j,k)
 4    continue
      close(44)
      open(44, file='ocn6g_lm.txt',status='unknown')
      do 7 k=1,nepochs
      read(44,*)
      read(44,*) l1,m1,Ao(k)
      do 7 i=1,5150
      read(44,*)
 7    continue
      close(44)

c here starts the iterations
      do 9 iter=1,num_iter
      print *, 'iter = ', iter, '...'
c loop through the observation time (same as load times here)
      do 6 nobs=1,nepochs
c      do 6 nobs=10,10
c sum viscous residues
      do 934 mm=1,nobs+1
      do 910 n=0,lmaxb+1
      pot_nr(mm,n)=0.
      u_nr(mm,n)=0.
      if (nn(n).gt.0) then
      ntop=nn(n)
      do 911 k=1,ntop
      if (mm.eq.1) then
      factor=0.
      else if (mm.eq.nobs+1) then
      factor=-1.
      else
      ddum=(tepoch(mm)-tepoch(nobs))*s(n,k)
      ddum1=(tepoch(mm-1)-tepoch(nobs))*s(n,k)
      factor=
     2  (exp(ddum)-exp(ddum1))/s(n,k)
      endif
      pot_nr(mm,n)=pot_nr(mm,n)+rk(n,k)*factor/(-s(n,k))
      u_nr(mm,n)=u_nr(mm,n)+rh(n,k)*factor/(-s(n,k))
911    continue
      endif
c      if (mm.ne.nobs+1) then
      if(mm.ne.nobs+1 .and. mm.ne.1) then
      pot_nr(mm,n)=pot_nr(mm,n)/(tepoch(mm-1)-tepoch(mm))
      u_nr(mm,n)=u_nr(mm,n)/(tepoch(mm-1)-tepoch(mm))
      endif
910   continue
934   continue

      do 869 l=0,lmaxb
      do 869 m=0,l
      do 869 j=1,2
      geoid_nr(j,l,m)=0.
      uplift_nr(j,l,m)=0.
      rNU_lm(j,l,m)=0.
869   continue

c combine the residues with the loading history
      do 569 nint=1,nobs
      do 30 l=0,lmaxb+1
      if (polar_wander.lt.0.5.and.l.eq.lmaxb+1) then
      go to 30
      endif
      if (nint.ne.nobs) then
      ee_nr=-pot_nr(nint,l)+pot_nr(nint+1,l)
      eu_nr=-u_nr(nint,l)+u_nr(nint+1,l)
      else
c  the load at the last time (nepochs): I change it so it is consistent
c  with the python code
c  add elastic response for the last epoch
      if (l.ne.lmaxb+1) then
      ee_nr=-pot_nr(nint,l)+pot_nr(nint+1,l)+1.+ek(l)
      eu_nr=-u_nr(nint,l)+u_nr(nint+1,l)+eh(l)
      else
      ee_nr=-pot_nr(nint,l)+pot_nr(nint+1,l)+ek(l)
      eu_nr=-u_nr(nint,l)+u_nr(nint+1,l)+eh(l)
      endif
      endif
c update the (2,1) term using rotation love numbers
      if (l.eq.lmaxb+1) then
      do 931 j=1,2
      geoid_nr(j,2,1)=geoid_nr(j,2,1)+
     2 almt(j,2,1,nint)*ee_nr/float(2*2+1)*coefu
      uplift_nr(j,2,1)=uplift_nr(j,2,1)+
     2 almt(j,2,1,nint)*eu_nr/float(2*2+1)*coefu
931   continue
      go to 30
      endif
      do 333 m=0,l
      do 31 j=1,2
      geoid_nr(j,l,m)=geoid_nr(j,l,m)+
     2   almt(j,l,m,nint)*ee_nr/float(l*2+1)*coefu
      uplift_nr(j,l,m)=uplift_nr(j,l,m)+
     2   almt(j,l,m,nint)*eu_nr/float(l*2+1)*coefu
31    continue
333    continue
30    continue
569   continue

c add in the centrif potential
      if (polar_wander.gt.0.5) then
      geoid_nr(1,2,1)=geoid_nr(1,2,1)*(1.+1./fkf)
      geoid_nr(2,2,1)=geoid_nr(2,2,1)*(1.+1./fkf)
      endif

c find NU_lm
      do 5 l=0,lmaxb
      do 5 m=0,l
      rNU_lm(1,l,m)=(geoid_nr(1,l,m)-uplift_nr(1,l,m))/100.
      rNU_lm(2,l,m)=(geoid_nr(2,l,m)-uplift_nr(2,l,m))/100.
  5   continue

c transform NU_lm to spacial domain:NU_ij
      call ylmspacial(lmaxb,1,1)
      do 2 i=1,nlat
      do 2 j=1,nlon
      rNU_ij(i,j)=rNU_ij(i,j)*ocn_xyt(i,j,nobs)
  2   continue

c find rI= (NU*O)_00
      call ylmspacial(lmaxb,-1,0)
      print *, 't = ', nobs, 'I(t) = ', rI/Ao(nobs)
c      call ylm_gauleg(rNU_ij,rLO_lm,lmaxb,0)
c      rI=rLO_lm(1,0,0)

c find (NU-rI/Ao)*O
      do 8 i=1,nlat
      do 8 j=1,nlon
      rNU_ij(i,j)=rNU_ij(i,j)-rI/Ao(nobs)*ocn_xyt(i,j,nobs)
  8   continue

c transform back to lm domain to get dynamic ocean load
c      call ylm_gauleg(rNU_ij,rLO_lm,lmaxb,1)
      call ylmspacial(lmaxb,-1,1)
      rLO_lm(1,0,0)=0.

c update the total load with almt0 + rLO_lm
      do 3 l=0,lmaxb
      do 3 m=0,l
      almt(1,l,m,nobs)=almt0(1,l,m,nobs)+rLO_lm(1,l,m)*100.*.9174
      almt(2,l,m,nobs)=almt0(2,l,m,nobs)+rLO_lm(2,l,m)*100.*.9174
c     geruo's ice density .9174
  3   continue

  6   continue
  9   continue

c output the updated load
      open (43,file='final_load.ice6g',status='unknown')
      do 5691 nint=1,nepochs
      write(43,*)
      write(43,*)
      write(43,*) 't=',nint
      do 313 l=0,lmaxylm
      do 313 m=0,l
      write (43,2903) l,m
     2 ,almt(1,l,m,nint)/(100*.9174),almt(2,l,m,nint)/(100*.9174)
c geruo's ice density .9174
313   continue
5691  continue
      close(43)
      close(41)
      close(42)

      return
      end

      subroutine ylmspacial(lmaxb,flag,deg0)
      parameter (nlat=181,nlon=361)
      parameter (lmax=101)
      implicit doubleprecision (a-h,o-z)
      external factl,plgndr,plegendre
      integer flag,deg0
      complex*16 c_alm(0:lmax,0:lmax),rexp(nlon,0:lmax),h
      dimension ylm(nlat,0:lmax,0:lmax),rlat(nlat)
     2 ,rlon(nlon),alat(nlat),alon(nlon)
      common /ijlm/ rNU_ij(1:nlat,1:nlon),rNU_lm(2,0:lmax,0:lmax)
     2 ,rLO_lm(2,0:lmax,0:lmax),rI

      pi=3.14159265
      dtr=pi/180.

      dlat=180./float(nlat-1)
      do 1720 l=0,lmaxb
      do 1720 m=0,l
c      clm=sqrt((2.*l+1.)/4./pi*factl(l-m)/factl(l+m))
      do 1720 i=1,nlat
      rlat(i)=dlat*float(i-1)
      cs=cos(dtr*rlat(i))
c      ylm(i,l,m)=plgndr(l,m,cs)*clm
      ylm(i,l,m)=plegendre(l,m,cs)
1720  continue
      dlon=360./float(nlon-1)
      do 1721 m=0,lmaxb
      do 1721 j=1,nlon
      rlon(j)=dlon*float(j-1)
      cphi=cos(float(m)*rlon(j)*dtr)
      sphi=sin(float(m)*rlon(j)*dtr)
      if (flag.eq.1) then
      rexp(j,m)=cmplx(cphi,sphi)
      if (m.ne.0) rexp(j,m)=2.*rexp(j,m)
      else
      rexp(j,m)=cmplx(cphi,-sphi)
      endif
1721  continue

      if (flag.eq.1) then
      do 1722 l=0,lmaxb
      do 1722 m=0,l
      c_alm(l,m)=cmplx(rNU_lm(1,l,m),rNU_lm(2,l,m))
1722  continue
      open (100,file='globalfield',status='unknown')
      do 720 i=1,nlat
      do 720 j=1,nlon
      rNU_ij(i,j)=0.
      h=0.
      do 916 l=0,lmaxb
      do 916 m=0,l
      h=h+ylm(i,l,m)*c_alm(l,m)*rexp(j,m)
916   continue
      hh=real(h)
      rNU_ij(i,j)=hh
c      write (100,*) rlat(i),rlon(j),rNU_ij(i,j)
720   continue

      else
c extended Simpson's rule:
      do 72 i=4,nlat-3
      alat(i)=1.
  72  continue
      do 73 i=4,nlon-3
      alon(i)=1.
  73  continue
      alat(1)=3./8.
      alat(2)=7./6.
      alat(3)=23./24.
      alat(nlat-2)=23./24.
      alat(nlat-1)=7./6.
      alat(nlat)=3./8.
      do 74 i=1,3
      alon(i)=alat(i)
      alon(nlon-i+1)=alat(nlat-i+1)      
  74  continue

      do 917 l=0,lmaxb*deg0
      do 917 m=0,l
      h=0.
      do 721 i=1,nlat
      do 721 j=1,nlon
      h=h+ylm(i,l,m)*rNU_ij(i,j)*rexp(j,m)*sin(dtr*rlat(i))
     2 *alat(i)*alon(j)
721   continue
      h=h*(2.*pi**2)/float(nlat-1)/float(nlon-1)
      if (deg0.eq.0) then
      rI=real(h)
      go to 917
      endif
      rLO_lm(1,l,m)=real(h)
      rLO_lm(2,l,m)=imag(h)
917   continue
      endif
      end

      subroutine ylm_gauleg(fij,flm,lmaxb,deg0)
c this subroutine finds the spherical harmonic expansion coefficients up to degree lmaxb
c for a given spacial field based on guass-legendre quadrature formula
c the field is evenly spaced (1deg x 1deg) on a sphere
      parameter (nlat=181,nlon=361)
      parameter (lmax=101)
      implicit double precision (a-h,o-z)
      external factl,plgndr
      integer deg0
      dimension fij(1:nlat,1:nlon),flm(2,0:lmax,0:lmax)
     2 ,gp_lat(2*(nlat-1)),gp_lon(2*(nlon-1)),rlat(nlat),rlon(nlon)
     3 ,gp_f(2*(nlat-1),2*(nlon-1))
     4 ,ylm(2*(nlat-1),0:lmax,0:lmax)

      pi=3.14159265
      dtr=pi/180.
      dlat=180./float(nlat-1)
      dlon=360./float(nlon-1)
      c1=(1.+1./sqrt(3.))/2.
      c2=(1.-1./sqrt(3.))/2.

c find the spacial grids
      do 1 i=1,nlat
      rlat(i)=dlat*float(i-1)
 1    continue
      do 2 i=1,nlon
      rlon(i)=dlon*float(i-1)
 2    continue

c find the guassian points
      do 3 i=1,nlat-1
      gp_lat(2*i-1)=(c1*rlat(i)+c2*rlat(i+1))*dtr
      gp_lat(2*i)=(c2*rlat(i)+c1*rlat(i+1))*dtr
  3   continue
      do 4 j=1,nlon-1
      gp_lon(2*j-1)=(c1*rlon(j)+c2*rlon(j+1))*dtr
      gp_lon(2*j)=(c2*rlon(j)+c1*rlon(j+1))*dtr
  4   continue
c loop through each element and interpolate the field at gaussian points
      do 5 i=1,(nlat-1)
      do 5 j=1,(nlon-1)
      gp_f(2*i-1,2*j-1)=c1*(c1*fij(i,j)+c2*fij(i+1,j))
     2 +c2*(c1*fij(i,j+1)+c2*fij(i+1,j+1))
      gp_f(2*i,2*j-1)  =c1*(c2*fij(i,j)+c1*fij(i+1,j))
     2 +c2*(c2*fij(i,j+1)+c1*fij(i+1,j+1))
      gp_f(2*i,2*j)    =c2*(c2*fij(i,j)+c1*fij(i+1,j))
     2 +c1*(c2*fij(i,j+1)+c1*fij(i+1,j+1))
      gp_f(2*i-1,2*j)  =c2*(c1*fij(i,j)+c2*fij(i+1,j))
     2 +c1*(c1*fij(i,j+1)+c2*fij(i+1,j+1))
  5   continue

c sum the quadrature for each (l,m)
      area=(pi**2)/float(nlat-1)/float(nlon-1)/2.
      do 6 l=0,lmaxb*deg0
      do 6 m=0,l
      flm(1,l,m)=0.
      flm(2,l,m)=0.
      do 7 i=1,2*(nlat-1)
      cs=cos(gp_lat(i))
      ylm(i,l,m)=plgndr(l,m,cs)
     2 *sqrt((2.*l+1.)/4./pi*factl(l-m)/factl(l+m))
      sint=sin(gp_lat(i))
      do 7 j=1,2*(nlon-1)
      c_cos=cos(float(m)*gp_lon(j))
      s_sin=-sin(float(m)*gp_lon(j))
      flm(1,l,m)=flm(1,l,m)+gp_f(i,j)*ylm(i,l,m)*c_cos*sint
      flm(2,l,m)=flm(2,l,m)+gp_f(i,j)*ylm(i,l,m)*s_sin*sint
  7   continue
      flm(1,l,m)=flm(1,l,m)*area
      flm(2,l,m)=flm(2,l,m)*area
  6   continue
      end


      function factl(m)
      implicit doubleprecision (a-h,o-z)
      factl=1.
      if (m.eq.0) return
      do 10 i=1,m
      factl=float(i)*factl
10    continue
      return
      end
c
      FUNCTION PLGNDR(L,M,X)
      implicit doubleprecision (a-h,o-z)
      IF(M.LT.0.OR.M.GT.L.OR.ABS(X).GT.1.) print *, 'bad arguments'
      PMM=1.
      IF(M.GT.0) THEN
        SOMX2=SQRT((1.-X)*(1.+X))
        FACT=1.
        DO 11 I=1,M
          PMM=-PMM*FACT*SOMX2
          FACT=FACT+2.
11      CONTINUE
      ENDIF
      IF(L.EQ.M) THEN
        PLGNDR=PMM
      ELSE
        PMMP1=X*(2*M+1)*PMM
        IF(L.EQ.M+1) THEN
          PLGNDR=PMMP1
        ELSE
          DO 12 LL=M+2,L
            PLL=(X*(2*LL-1)*PMMP1-(LL+M-1)*PMM)/(LL-M)
            PMM=PMMP1
            PMMP1=PLL
12        CONTINUE
          PLGNDR=PLL
        ENDIF
      ENDIF
      RETURN
      END

      function plegendre(l,m,x)
      implicit doubleprecision (a-h,o-z)
      PI=3.141592653589793
      IF(M.LT.0.OR.M.GT.L.OR.ABS(X).GT.1.) print *, 'bad arguments'
      pmm=1.
      if(m.gt.0) then
        omx2=(1.-x)*(1.+x)
        fact=1.
        DO 11 I=1,M
          pmm=pmm*omx2*fact/(fact+1.)
          fact=fact+2.
11      CONTINUE
      ENDIF
      pmm=sqrt((2.*m+1.)*pmm/(4.*PI))
      if (mod(m,2).eq.1) pmm=-pmm
      IF(L.EQ.M) THEN
        plegendre=pmm
      ELSE
        PMMP1=X*sqrt(2.*M+3.)*PMM
        IF(L.EQ.M+1) THEN
          plegendre=PMMP1
        ELSE
          oldfact=sqrt(2.*m+3.)
          DO 12 LL=M+2,L
            fact=sqrt((4.*ll**2-1.)/(ll**2-m**2))
            PLL=(X*PMMP1-PMM/oldfact)*fact
            oldfact=fact
            PMM=PMMP1
            PMMP1=PLL
12        CONTINUE
          plegendre=PLL
        ENDIF
      ENDIF
      RETURN
      END





