      PROGRAM MAIN
c      use mpmodule
C   PROGRAM TO COMPUTE LOAD LOVE NUMBERS for a viscoelastic earth,
c    using the collocation method
c 9999999999999999999
c  11/30/99:
c  As far as I can tell, the solution vector is ordered as follows:
c  (1)=U
c  (2)=V
c  (3)=P
c  (4)=Q
c  (5)=PHI
c  (6)=g1E
c  And these variables are related to those in my Admont notes (the
c    vaariables denoted as "*_me" below) as follows:
c  U = U_me
c  -sqrt[2*l*(l+1)]*V = V_me
c  P = P_me
c  -sqrt[2*l*(l+1)]*Q = Q_me
c  PHI = PHI_me
c  g1E - (l+1)/r*PHI = g1E_me
c 9999999999999999999
      parameter(nnodes=200,numbs=1000)
      implicit double precision(a-h,o-z)
c      implicit type (mp_real) (a-h,o-z)
      COMMON CON,R,RHO,VP,VS,VFSQ,FLAM,FMU,G,RHOBAR,W,WSQ,
     1DEL1,DEL2,QS,QF,Q,EPS,SFL3,FORD,VERTNO,PI,FL,
     2FL1,FL2,FL3,N1,N2,N,JCOM,L
      COMMON/E/NDIS,KDIS
      COMMON/F/C,NIR,NER
      COMMON/P/CAPW,EA,EAH,ELL
      COMMON/C/STEP
      COMMON/J/NPT
      COMMON/R2/R2
      COMMON viscUpper, viscLower !read viscosity from files
      COMMON  coll1, coll2, coll3  !read collocation tuning para from files
      DIMENSION C(nnodes),FLAM(nnodes),FMU(nnodes),G(nnodes),
     2  R(nnodes),R2(nnodes),RHO(nnodes)
     1,VFSQ(nnodes),VP(nnodes),VS(nnodes),BULK(nnodes),
     2  ELL(nnodes),STEP(8),KDIS(28),MOD
     2EL(5),LRD(36),wp(nnodes),ws(nnodes)
     3 ,gmu(nnodes),glam(nnodes),visc(nnodes),wisc(nnodes)
     4 ,svalue(numbs)
     5 ,soln(0:numbs,3),temp(numbs),zmatrix(numbs,numbs)
     6 ,indx(numbs), r_test(6), rho_test(6), smu_test(6)
     7 ,visc_test(6), check_sum(3), check_s0(3)
     8 ,soln_rot(0:numbs,3),temp_rot(numbs),z_trunc(numbs,numbs)
c     9 ,rho1(nnodes),vs1(nnodes)
c  1066A
c     open (10, file='model.prem',status='unknown')
c      open (10, file='model.low',status='unknown')
      open (10, file='model.vm5a',status='unknown')
c       open (10, file='model.paul',status='unknown')
c       open (10, file='model.wei',status='unknown')
c      open (10, file='model.Lambeck',status='unknown')
c      open (10, file='model.vm7',status='unknown')
c      open (10, file='model.vm2',status='unknown')
c model.citcom.all the citcom model
c      open (10, file='model.vm2',status='unknown')
c      open (15, file='viscosity.model',status='unknown')
      open (11, file='in.dahlen.viscoelastic',status='unknown')
      open (16, file='out.dahlen.love',status='unknown')
      open (20, file='out.dahlen.coll',status='unknown')
      open (21, file='check_s0',status='unknown')
      open (22, file='chk_det', status='unknown')
      open (23, file='out.btide',status='unknown')
      open (24, file='out.rot',status='unknown')
c     read viscosity and collocation tuning para from files
      open (100, file='in.para.collocation', status='unknown')
      read (100,*) viscUpper
      read (100,*) viscLower
      print *, viscUpper, viscLower
c   Read collocation para from file
      read (100,*) coll1
      read (100,*) coll2
      read (100,*) coll3

      print *, coll1 
      print *, coll2 
      print *, coll3
      close(100)

c      fkf=0.9422d0
c      fkf=0.940954d0
      read (11,*) fkf
      PI=3.1415926536d0
      CON=PI*6.673D-11
      CAPW=7.292116D-05
      EA=(1.)/(2.982582)
      EPS=1.D-08
      NPT=1
      T3=2./3.
c    1 READ(10,102) (MODEL(J),J=1,5),N,(R(I),RHO(I),VP(I),VS(I),VFSQ(I),
c     1 I=1,N)
      read(10,1021) (MODEL(J),J=1,5),n,(r(i),rho(i),vp(i),vs(i),tmp
     1 ,visc(i),i=1,n)
 1021 format(5A4,I10/(4F10.0,F10.5,E12.4))

      do 92 i=1,n
      if (visc(i).ge.1.e30) visc(i)=1.e26
      if (visc(i).lt.1.e-30) visc(i)=0.
   92 continue

c change start - make a uniform core
      do 1191, i=1,n
c      vp(i)=vp(i)*1.e3
c turn it on for structures from citcomsve runs
c      if (r(i) .lt. 3485500.) then
      if (i.lt.67) then
c      rho(i)= 10005.4 !10008.0 !9825.0 !10005.4
c      vp(i) = vp(i)*1.e3
      endif
c make a mantle with constant density and shear modulus
      if (i.ge.67) then
c      rho(i)= 4604.4 !4444.8 !4400.0 !4604.4
c      vp(i)=vp(i)*1.e3
c      vp(i)=2.*sqrt(1.4305e11/rho(i))
c      vs(i)=sqrt(1.4305e11/rho(i))
c      visc(i)=2.e21
      endif


c      vs1(i)=sqrt(rho(i)/rho1(i))*vs(i)
c      vs(i)=vs1(i)
c      rho(i)=rho1(i)
c      print *, vs(i)**2*rho(i)
c      print *, rho(i)
c      print *, r(i), visc(i)
c set viscosity for lith, upper mantle and lower mantle       
      if (i.ge.67.and.i.lt.134) then
      visc(i)= viscLower !3.16e22 !1.0e22!3.16e21 !0.2e22 !5.0e22
c      rho(i)=4400.
c      vs(i)=sqrt(1.4305e11/rho(i))
      endif

      if (i.ge.134.and.i.lt.158) then !(1)149 for wei (2) 162 for model.Lambeck
      visc(i)= viscUpper !1.5e20 !1.25e20
c      rho(i)=3300.
c      vs(i)=sqrt(1.4305e11/rho(i))
      endif

c      if (i.ge.149.and.i.lt.159) then  ! wei's viscosity with a weak layer
c      visc(i)= 3.0e19 !1.5e20 !1.25e20
c      rho(i)=3300.
c      vs(i)=sqrt(1.4305e11/rho(i))
c      endif

      if (i.ge.158) then
      visc(i)=1.e26
c      print *,visc(i)
c      rho(i)=3300
      endif

1191  continue
c change end


cc  read in the viscosity model
c      do 91 i=1,n
c      read (15,1020) zr,visc(i)
c1020  format(f10.0,2x,f7.4)
c      if (zr.ne.r(i)) then
c      print *,'error #1',i zr,r(i)
c      stop
c      endif
cc  Change viscosity to SI units (Pa-sec)
c      visc(i)=visc(i)*1.e21
cc     print *,visc(i)
c91    continue
  102 FORMAT(5A4,I10/(4F10.2,E10.3))

c ------------- change  ----------------- start  -----------------

      go to 9995
      r_test(1)=3503500.
      r_test(2)=5201000.
      r_test(3)=5701000.
      r_test(4)=5971000.
      r_test(5)=6251000.
      r_test(6)=6371000.
      
      rho_test(1)=9900.
      rho_test(2)=4617.
      rho_test(3)=4617.
      rho_test(4)=4227.
      rho_test(5)=4047.
      rho_test(6)=4047.

      smu_test(1)=0.
      smu_test(2)=2.254e11
      smu_test(3)=2.254e11
      smu_test(4)=1.227e11
      smu_test(5)=0.7587e11
      smu_test(6)=0.4558e11

      visc_test(1)=0.
      visc_test(2)=3.6e21
      visc_test(3)=0.9e21
      visc_test(4)=0.9e21
      visc_test(5)=0.9e21
      visc_test(6)=0.
      
      k=1
      do 9999 i=1,n
      vp(i)=vp(i)*1.e3
      if (r(i).le.r_test(k)) then
      rho(i)=rho_test(k)
      visc(i)=visc_test(k)
      if (i.gt.67) then
      vs(i)=sqrt(smu_test(k)/rho_test(k))
      endif
      else
      k=k+1
      if (r(i-1).eq.r(i-2)) then
      rho(i-1)=rho_test(k)
      visc(i-1)=visc_test(k)
      vs(i-1)=sqrt(smu_test(k)/rho_test(k))
      endif
      rho(i)=rho_test(k)
      visc(i)=visc_test(k)
      vs(i)=sqrt(smu_test(k)/rho_test(k))
      endif
 9999 continue
 9995 continue

      READ(10,103) RHOBAR,SIGRHO,Z,SIGZ
  103 FORMAT(F10.0,E10.3,F10.5,E10.3)

      READ(10,104) NDIS,(KDIS(KK),KK=1,NDIS)
  104 FORMAT(16I5)

      CALL DG
      print *, 'first rhobar = ', rhobar
      NIR=1
      NER=N
      RN=R(N)
c      READ(10,103) RHOBAR,SIGRHO,Z,SIGZ
c  103 FORMAT(F10.0,E10.3,F10.5,E10.3)
c      rhobar=5355.13
      RABOHR=RHOBAR
      ZED=.5*Z
      DO 620 I=1,N
      R(I)=R(I)/RN
      R2(I)=R(I)*R(I)
  620 RHO(I)=RHO(I)/RHOBAR
      CALL PREINT
  970 DO 971 I=1,N
  971 C(I)=RHO(I)
      CALL GLINT(Q1)
      
      DO 972 I=1,N
  972 C(I)=RHO(I)*R(I)*R(I)
      CALL GLINT(Q2)
      A1=26.25*(ZED-Q2)-6.25+18.75*Q1
      A2=43.75*(Q2-ZED)+8.75-26.25*Q1
c      DO 973 I=1,N
c  973 RHO(I)=RHO(I)-A1-A2*R(I)*R(I)
      DO 974 I=1,N
      R(I)=R(I)*RN
  974 RHO(I)=RHO(I)*RHOBAR
      CALL DG
      print *, 'does DG change rhobar?', rhobar
      print *, 'old rhobar = ', rabohr
c      RHOBAR=RABOHR

c change end

      GN=CON*RHOBAR*RN
      print *,'surface g =',gn*4./3.
      V=SQRT (GN*RN)
      WN=V/RN
      CAPW=CAPW/WN
C   CON IS VISCOUS PARAMETER.  ABOVE, IT WAS USED IN SUBROUTINE DG
      CON=1.
      DO 7 I=1,N
      R(I)=R(I)/RN
      VP(I)=VP(I)/V
      VS(I)=VS(I)/V
      RHO(I)=RHO(I)/RHOBAR
      gMU(I)=RHO(I)*VS(I)*VS(I)
      gLAM(I)=RHO(I)*VP(I)*VP(I)-2.*gMU(I)
      BULK(I)=gLAM(I)+T3*gMU(I)
c  I think I should now normalize visc, since it will be used below
c      only after dividing into a normalized mu
      wisc(i)=visc(i)/(v*v*rhobar)
c     write (99,*) i,gmu(i),glam(i),visc(i)
c      print *, 'surface g = ', g(N)
    7 G(I)=G(I)/GN
      CALL ELLIP
      DO 701 I=1,N
      R(I)=R(I)*RN
      RHO(I)=RHO(I)*RHOBAR
      VP(I)=VP(I)*V
      VS(I)=VS(I)*V
      print *, i, g(i)*gn
  701 G(I)=G(I)*GN
   30 I=1
    2 LINES=1
      NP=NP+1
c      WRITE(16,902) NP,(MODEL(J),J=1,5),N
  902 FORMAT(1H1,116X,I3/10X,10HMODEL NAME,10X,6HLEVELS//5X,5A4,I10///4X
     1,5HLEVEL,2X,6HRADIUS,11X,3HRHO,13X,2HVP,13X,2HVS,11X,7HGRAVITY,9X,
     211HELLIPTICITY/)
c    3 WRITE(16,9021) I,R(I),RHO(I),VP(I),VS(I),G(I),ELL(I)
c change start-end remove the output
    3 continue
 9021 FORMAT(4X,I3,F11.0,3F15.0,2E18.4)
c      IF(I-N) 4,6,6
      IF(I-N .LT. 0) then 
      GOTO 4
      ELSE IF (I-N .GE. 0) then 
      GOTO 6
      ENDIF
    4 I=I+1
c      IF(LINES-50)5,2,2
      IF(LINES-50 .LT. 0) THEN
      GOTO 5
      ELSE IF(LINES-50 .GE. 0) THEN
      GOTO 2
      ENDIF
    5 LINES=LINES+1
      GO TO 3
    6 CONTINUE
c      WRITE(16,903) RHOBAR,Z
  903 FORMAT(//10X,8HRHOBAR= ,F10.4,10X,3HZ= ,F10.8)
      EAR=1./EA
      EAHR=1./EAH
c      WRITE(16,9031) EAR,EAHR
 9031 FORMAT(  10X,23HELLIPTICITY = ONE OVER ,F7.3,/10X,35HHYDROSTATIC E
     1LLIPTICITY = ONE OVER ,F7.3)
      DO 702 I=1,N
      R(I)=R(I)/RN
      RHO(I)=RHO(I)/RHOBAR
c  change start
      wP(I)=wP(I)/V
      wS(I)=wS(I)/V
c      vp(i)=vp(i)/v
c      vs(i)=vs(i)/v
c  change end
  702 G(I)=G(I)/GN
      CALL STEPS(EPS,VERTNO)
      W=0.
      WSQ=W*W
c      WRITE(16,706)
  706 FORMAT(1H1,//////)

c  read in the min (should be >0) and max l values 
c      to consider, and the spacing
c      between the l values (probably should be 1)
      read (11,*) l1,l2,nl
      if (nl.le.0) then
      print *,'error #11',nl
      stop
      endif
c      write (20,310) 
310   format(1x,'elastic model = ?'/'LM viscosity = ?'
     2  /'UM viscosity = ?'/'Lithospheric thickness (km) = ?')
c      write (20,300) l1,l2,nl
c      write (60,300) l1,l2,nl
300   format(1x,'l1, l2, are:',2i5,2x,'l spacing=',i5)

c  read in the min and max of the log of the period range 
c    that we find solns in;
c  Also, read in the no of steps we divide the log of the period 
c    range into.
      read (11,*) plog1, plog2, ns
      if (ns.gt.numbs) then
      print *,'error #10: ns too large',ns
      stop
      endif
c      write (20,302) plog1,plog2,ns
c      write (60,302) plog1,plog2,ns
302   format(1x,'log of the period range (in years): from ',f11.2,2x,
     2 ' to ' f11.2,4x,'no of steps:',i4)

      do 911 i=1,n
      wp(i)=vp(i)
      ws(i)=vs(i)
911   continue


c  find the solutions for each l
      DO 705 L=l1,l2,nl
      print *,'about to start l=',l
      write (20,304) l,ns
c      write (60,305) l
304   format(1x,'******************************************'/
     2 3x,i5,2x,i5,4x'(value of l, and number of s values)')
c305   format(1x,'******************************************'/
c     2 3x,i5,4x'(value of l)')
      FL=L
      FL1=FL+1.
      FL2=FL+FL1
      FL3=FL*FL1
      SFL3=SQRT(FL3)
c
c   loop through the s-values, finding solutions
      dlog=(plog2-plog1)/float(ns)

c change start: use a different set of s values of degree 1 and 2
      if (l.eq.1) then
      dlog=(coll1-plog1)/float(ns)
      else if (l.gt.1.and.l.le.2) then
      dlog=(coll2-plog1)/float(ns)
      else if (l.gt.2.and.l.lt.11) then
      dlog=(coll3-plog1)/float(ns)

c for 4610(shijie's density) for mantle density :7.4; 9.9; 7.0
c for g=9.8 case:
c 7.1;9.95;9.45  2.0, 6.6 10  or 8.235 for l>3
c for 100 km lithosphere with 10e25 PaS kf = 1.1016615908859579 ( 2.0  6.2  20):
c 8.175 or 9.75 for l=1; 9.8 for l=2; 7.0 others 
c for one layer viscosity model with kf=1.008: 8.175 8.95 7.0
c for paul with vm5a viscosity: 8.0  9.585  7.
c for paul: 8.0   6.488   7.0
c for 2-layer density
c  7.5,12.8,7.0
c for incom 4400 4745 density
c     8.5 for l=1, 10.0 for l=2, 7.0 for l=3-10
c for fkf = 0.940954d0, VM5a:
c 8.5 for l=1, 8.588 for l=2, 7.0 for l=3-10
c for fkf = 0.935759d0 , VM5a
c 8.6 for l=2, others are the same
c ====
c when i include l love number, i set the three parameters as 
c 8.5 for l=1~2, 7 for l=3~11
c which is different from the runs i did in March, where the parameters were
c 10 for l=1, 5 for l=2~7
      endif
c change end
      do 50 nmodes=0,ns
c     print *,'about to start nmodes=',nmodes
      plog=plog1+dlog*float(nmodes-1)
      per=10**(plog)*86400.*365.25
      s=1./per
      if (nmodes.ne.0) svalue(nmodes)=s
      do 954 n1=1,n
      if ((visc(n1).gt.0.).and.(nmodes.ne.0)) then
      vs(n1)=ws(n1)*sqrt(s/(s+gmu(n1)/visc(n1)))
      vp(n1)=wp(n1)*sqrt(s/(s+gmu(n1)/visc(n1)))*
     2  sqrt(s+gmu(n1)/visc(n1)*(2./3.*gmu(n1)+glam(n1))/
     3  (glam(n1)+2.*gmu(n1)))
      fmu(n1)=gmu(n1)*s/(s+gmu(n1)/wisc(n1))
      flam(n1)=
     2  glam(n1)*(s+gmu(n1)/wisc(n1)*(2.*gmu(n1)/3./glam(n1) + 1.))
     3  /(s+gmu(n1)/wisc(n1))
c     write (99,*) n1,r(n1)*rn,fmu(n1)/gmu(n1)-1.,
c    2   flam(n1)/glam(n1)-1.,wisc(n1)
      else
      fmu(n1)=gmu(n1)
      flam(n1)=glam(n1)
      vp(n1)=wp(n1)
      vs(n1)=ws(n1)
      endif
954   continue

c check det - change start
      write (22,7999) s
7999  format(e20.10)
c change end

      CALL LOVENO(0,HPR,FKPR,FLPR)
c      write (20,987) 10**plog,hpr,fkpr,flpr
c987   format(e12.4,3(2x,e16.8))
c  covert h love nos to displacement (in cm per cgs forcing potential)
      soln(nmodes,1)=-hpr
c      soln(nmodes,1)=-hpr/(gn*4./3.)/100.
c  soln(nmodes,2) is the potential (in cgs per cgs forcing potential)
      soln(nmodes,2)=fkpr
c  covert l love nos to displacement (in cm per cgs forcing potential)
c   Except when I compare this with the output of ../post.coll.f,
c     I find a factor of 2 difference.  I'm not sure why; but
c     rather than figure it out at the moment, I'm just going
c     to ignore horizontal displacements.   So, comment out the
c     next line, and seet kmax=2, rather than =3.
c      soln(nmodes,3)=flpr/(gn*4./3.)/100.*sqrt(l*(l+1.)/2.)
      soln(nmodes,3)=-flpr
      kmax=3
c      kmax=2
      WRITE(16,703) L,HPR,FKPR,FLPR,s
  703 FORMAT(10X,2HL=,I3,10X,7HHPRIME=,F12.8,10X,7HKPRIME=,F12.8,10X,7HL
     1PRIME=,F12.8,10X,e17.10)
c i deleted ////// at the end of 703

c change start - rotation love numbers (degree 2,1 terms), find the body tide love numbers
      if (l.eq.2) then
      call loveno(1,ht,fkt,flt)
      soln_rot(nmodes,1)=-ht
      soln_rot(nmodes,2)=fkt
      soln_rot(nmodes,3)=-flt
c      write(23,2703) l,ht,fkt,flt,s
c 2703 format(10x,2hl=,I3,10x,7hhbtide=,f12.8,10x,7hkbtide=,f12.8,10x,7hl
c     1btide=,f12.8,10x,e17.10)
      endif 
c change end

50    CONTINUE

c truncate if kt stops monotonically rising
      if (l.eq.2) then
      ns_trunc = ns
      do 2250 nmodes=0,ns
      if (nmodes.gt.1) then
      if (soln_rot(nmodes,2).lt.soln_rot(nmodes-1,2)) then
      ns_trunc = nmodes-1
      print *, 'truncating rotation love number s-space at index'
     1,ns_trunc
      go to 2251
      endif
      endif
2250  continue
2251  continue

c find rotation love numbers
      do 2252 nmodes=0,ns_trunc
      ht = -soln_rot(nmodes,1)
      fkt = soln_rot(nmodes,2)
      flt= -soln_rot(nmodes,3)
      fkpr = soln(nmodes,2)
      hr = ht*(1.+fkpr)/(fkf-fkt)
      fkr = fkt*(1.+fkpr)/(fkf-fkt)
      flr = flt*(1.+fkpr)/(fkf-fkt)
      write(23,2704) l,hr,fkr,flr,svalue(nmodes)
 2704 format(10x,2hl=,I3,10x,7hhrotat=,f12.8,10x,7hkrotat=,f12.8,10x
     1,7hlrotat=,f12.8
     2,10x,e17.10)
      soln_rot(nmodes,1) = hr
      soln_rot(nmodes,2) = fkr
      soln_rot(nmodes,3) = flr
2252  continue
     
c find residues using collocation technique
      do 2253 ni=1,ns_trunc
      taui=1./svalue(ni)
      do 2253 nj=1,ns_trunc
      tauj=1./svalue(nj)
      z_trunc(ni,nj)=taui*tauj/(taui+tauj)
2253  continue
      call zludcmp(z_trunc,ns_trunc,numbs,indx,d)
      do 2254 k=1,kmax
      check_sum(k)=0.
      check_s0(k)=soln_rot(ns_trunc,k)-soln_rot(0,k)
      do 2255 ni=1,ns_trunc
      temp_rot(ni)=soln_rot(ni,k)-soln_rot(0,k)
2255  continue
      call zlubksb(z_trunc,ns_trunc,numbs,indx,temp_rot)
      do 2256 ni=1,ns_trunc
      soln_rot(ni,k)=temp_rot(ni)
c Mark's test - change start
      check_sum(k)=check_sum(k)+soln_rot(ni,k)/svalue(ni)
c change end
2256  continue
c Mark's test
      print *, '2,1 terms'
      print *, '@@', check_sum(k)/check_s0(k), soln_rot(0,k)
2254  continue
      write (24,2259) l,ns_trunc
2259  format(1x,'******************************************'/
     2 3x,i5,2x,i5,4x'(value of l, and number of s values)')
      do 2257 ni=1,ns_trunc
      period=1/svalue(ni)
      period=period/86400./365.25
      write (24,2258) ni,period,svalue(ni),(soln_rot(ni,k),k=1,kmax)
2258  format(i5,8(1x,e17.10))
2257  continue
      endif


c
c  Now, find the residue values for this collocation technique
c
c  First, construct the matrix of relaxation times.
c
      do 500 ni=1,ns
      taui=1./svalue(ni)
      do 500 nj=1,ns
      tauj=1./svalue(nj)
      zmatrix(ni,nj)=taui*tauj/(taui+tauj)
500   continue

c  Invert that matrix
      call zludcmp(zmatrix,ns,numbs,indx,d)
c  Subtract the elastic solution from the s-dependent soln,
c     and solve for the residues  
      do 5001 k=1,kmax
      check_sum(k)=0.
      check_s0(k)=soln(ns,k)-soln(0,k)
      do 5002 ni=1,ns
      temp(ni)=soln(ni,k)-soln(0,k)
5002  continue
      call zlubksb(zmatrix,ns,numbs,indx,temp)
      do 5004 ni=1,ns
      soln(ni,k)=temp(ni)
c Mark's test - change start
      check_sum(k)=check_sum(k)+soln(ni,k)/svalue(ni)
c change end
5004  continue
c Mark's test
      print *, '@@', check_sum(k)/check_s0(k), soln(0,k)
      write (21,9001) l, check_sum(k)/check_s0(k)*1e2,
     1check_s0(k)+soln(0,k),check_s0(k),check_sum(k)
9001  format(i3,4(4x,f12.8))
5001  continue
      do 5003 ni=1,ns
      period=1/svalue(ni)
      period=period/86400./365.25
      write (20,2004) ni,period,svalue(ni),(soln(ni,k),k=1,kmax)
2004  format(i5,8(1x,e17.10))
5003  continue
662   continue
c
  705 CONTINUE
      STOP
      END

      SUBROUTINE LOVENO(i_bodytide,HPR,FKPR,FLPR)
C     FOR USE WITH AN EARTH MODEL WITHOUT A SURFICIAL OCEANIC LAYER.
      parameter(nnodes=200)
      implicit double precision(a-h,o-z)
      integer i_bodytide
c      type (mp_real) 
      COMMON CON,R,RHO,VP,VS,VFSQ,FLAM,FMU,G,RHOBAR,W,WSQ,
     1DEL1,DEL2,QS,QF,Q,EPS,SFL3,FORD,VERTNO,PI,FL,
     2FL1,FL2,FL3,N1,N2,N,JCOM,L
      COMMON/E/NDIS,KDIS
      DIMENSION R(nnodes),RHO(nnodes),VP(nnodes),VS(nnodes),
     2  VFSQ(nnodes),FLAM(nnodes),FMU(
     1nnodes),G(nnodes),PHIC(nnodes),PHIPRC(nnodes),AS(6,3),
     2  B(3,3),BINV(3,3),AR(3),
     2ASA(6),P(nnodes),S(nnodes),KDIS(28)
      NIC=KDIS(1)
      NICP=NIC+1
      NOC=KDIS(2)
      NOCP=NOC+1
      IN=2
      IF(L.GT.20) IN=15
c change start - test for starting solution
c      CALL SPSJFG(N,AS)
c      DO 9990 K=1,3
c      Do 9990 J=1,6
c 9990 print *, AS(J,K)
c change end
      CALL SPSJFG(IN,AS)
      DO 1 K=1,3
      AS(2,K)=AS(2,K)/SFL3
      AS(5,K)=-AS(5,K)
    1 AS(6,K)=-AS(6,K)
      NS=IN+1
      Y=R(IN)
c  test for incompressible (homogeneous) earth - change start
c      go to 9994
c  test for compressible homogeneous fluid earth - change start
      go to 9989
      F1=AS(3,1)-RHO(NS)*AS(5,1)-RHO(NS)*G(IN)*AS(1,1)
      F2=AS(3,2)-RHO(NS)*AS(5,2)-RHO(NS)*G(IN)*AS(1,2)
      F3=AS(3,3)-RHO(NS)*AS(5,3)-RHO(NS)*G(IN)*AS(1,3)
      Q1=AS(4,1)
      Q2=AS(4,2)
      Q3=AS(4,3)
      BIGA=-(F3/F2-Q3/Q2)/(F1/F2-Q1/Q2)
      BIGB=-(F3/F1-Q3/Q1)/(F2/F1-Q2/Q1)
      F=AS(5,3)+BIGB*AS(5,2)+BIGA*AS(5,1)
      FPR=AS(6,3)-(FL1/R(IN))*AS(5,3)-4.*RHO(NS)*AS(1,3)+BIGB*(AS(6
     1,2)-(FL1/R(IN))*AS(5,2)-4.*RHO(NS)*AS(1,2))+BIGA*(AS(6,1)-(FL
     21/R(IN))*AS(5,1)-4.*RHO(NS)*AS(1,1))
      PHIC(IN)=F
      PHIPRC(IN)=FPR
      DO 9930 I=NS,N
      X=Y
      Y=R(I)
      CALL CMTR(X,Y,F,FPR)
      PHIC(I)=F
 9930 PHIPRC(I)=FPR
      ROC=RHO(N)
      GC=G(N)
      FACC=PHIPRC(N)/PHIC(N)+FL1/R(N)
      DO 995 J=1,6
      DO 995 K=1,3
  995 AS(J,K)=0.
      AS(1,1)=1.
      AS(3,1)=ROC*GC
      AS(6,1)=4.*ROC
      AS(2,2)=1.
      AS(5,3)=1.
      AS(3,3)=ROC
      AS(6,3)=FACC
      go to 9988
 9989 continue
c  change end
      DO 4 I=NS,NIC
      X=Y
      Y=R(I)
    4 CALL SMTR(X,Y,AS)
      F1=AS(3,1)-RHO(NICP)*AS(5,1)-RHO(NICP)*G(NIC)*AS(1,1)
      F2=AS(3,2)-RHO(NICP)*AS(5,2)-RHO(NICP)*G(NIC)*AS(1,2)
      F3=AS(3,3)-RHO(NICP)*AS(5,3)-RHO(NICP)*G(NIC)*AS(1,3)
      Q1=AS(4,1)
      Q2=AS(4,2)
      Q3=AS(4,3)
      BIGA=-(F3/F2-Q3/Q2)/(F1/F2-Q1/Q2)
      BIGB=-(F3/F1-Q3/Q1)/(F2/F1-Q2/Q1)
      F=AS(5,3)+BIGB*AS(5,2)+BIGA*AS(5,1)
      FPR=AS(6,3)-(FL1/R(NIC))*AS(5,3)-4.*RHO(NICP)*AS(1,3)+BIGB*(AS(6
     1,2)-(FL1/R(NIC))*AS(5,2)-4.*RHO(NICP)*AS(1,2))+BIGA*(AS(6,1)-(FL
     21/R(NIC))*AS(5,1)-4.*RHO(NICP)*AS(1,1))
      PHIC(NICP)=F
      PHIPRC(NICP)=FPR
      NS=NICP+1
      Y=R(NICP)
      DO 30 I=NS,NOC
      X=Y
      Y=R(I)
      CALL CMTR(X,Y,F,FPR)
      PHIC(I)=F
   30 PHIPRC(I)=FPR
      ROC=RHO(NOC)
      GC=G(NOC)
      FACC=PHIPRC(NOC)/PHIC(NOC)+FL1/R(NOC)
      DO 5 J=1,6
      DO 5 K=1,3
    5 AS(J,K)=0.
      AS(1,1)=1.
      AS(3,1)=ROC*GC
      AS(6,1)=4.*ROC
      AS(2,2)=1.
      AS(5,3)=1.
      AS(3,3)=ROC
      AS(6,3)=FACC
c change start - degree 1 term
      if (l.eq.1) then
      as(1,2)=1.
      as(2,2)=1.
      as(5,2)=-gc
c      as(1,1)=3*as(1,1)-as(1,2)
c      as(2,1)=3*as(2,1)-as(2,2)
c      as(5,1)=3*as(5,1)-as(5,2)
c      print *, as(6,3),(2*l+1)/r(noc)
      endif
c change end
      NS=NOCP+1
      Y=R(NOCP)
 9994 continue
c      print *, 'R(NOC) = ', y
      DO 6 I=NS,N
      if (l.eq.1) then
      as(1,2)=1.
      as(2,2)=1.
      as(5,2)=-g(i-1)
      endif
      X=Y
      Y=R(I)
    6 CALL SMTR(X,Y,AS)
 9988 continue
c change start - just in case..
      if (l.eq.1) then
      as(1,2)=1.
      as(2,2)=1.
      as(5,2)=-g(n)
      endif
c change end
      DO 7 K=1,3
      B(1,K)=AS(3,K)
      B(2,K)=AS(4,K)
    7 B(3,K)=AS(6,K)
      if (l.eq.1) then
      call inv22(b,binv)
      xt=as(5,1)/g(n)
      yt=as(5,3)/g(n)
      binv(2,1)=binv(1,1)*xt+binv(3,1)*yt
      binv(2,2)=binv(1,2)*xt+binv(3,2)*yt
      else
      CALL MATADJ(B,BINV)
      DET=B(1,1)*BINV(1,1)+B(1,2)*BINV(2,1)+B(1,3)*BINV(3,1)
c change start - check det
      write (22,7998) DET
7998  format(e20.10)
      Do 7996 K=1,3
7996  write (22,7997) B(K,1), B(K,2), B(K,3)
7997  format(3(e20.10,2X))
c change end
      DO 48 J=1,3
      DO 48 K=1,3
   48 BINV(J,K)=BINV(J,K)/DET
      endif
c change start - add body tide option
      if (i_bodytide.eq.0) then
      DO 47 J=1,3
c  The following is used for load love nos
   47 AR(J)=-G(N)*BINV(J,1)-4.*BINV(J,3)
      else
      do 472 j=1,3
c  The following is used for body tide love nos
  472 AR(J)=-4.*BINV(J,3)
      endif
c change end
      DO 41 J=1,6
   41 ASA(J)=0.
      DO 42 J=1,6
      DO 42 K=1,3
   42 ASA(J)=ASA(J)+AS(J,K)*AR(K)
      FKPR=-1.-(FL2/(4.*R(N)))*ASA(5)
      HPR=G(N)*(FL2/(4.*R(N)))*ASA(1)
      FLPR=G(N)*(FL2/(4.*R(N)))*ASA(2)
c      flpr=g(n)*(fl2/(4.*r(n)))*asa(2)/sqrt(fl3/2.)
c      print *, 'g(n) = ', g(n)
      RETURN
      END

      subroutine inv22(a,c)
c invert a 2x2 matrix
      implicit double precision (a-h,o-z)
      dimension a(3,3),c(3,3)
      det=a(1,1)*a(2,3)-a(2,1)*a(1,3)
      do 565 i=1,3
      do 565 j=1,3
 565  c(i,j)=0.
      c(1,1)=a(2,3)/det
      c(3,2)=a(1,1)/det
      c(3,1)=-a(2,1)/det
      c(1,2)=-a(1,3)/det
      return
      end

      SUBROUTINE DG
C     UNALTERED DPFG VERSION
      parameter(nnodes=200)
      implicit double precision(a-h,o-z)
      COMMON CON,R,RHO,VP,VS,VFSQ,FLAM,FMU,G,RHOBAR,W,WSQ,
     1DEL1,DEL2,QS,QF,Q,EPS,SFL3,FORD,VERTNO,PI,FL,
     2FL1,FL2,FL3,N1,N2,N,JCOM,L
      COMMON/E/NDIS,KDIS
      DIMENSION R(nnodes),RHO(nnodes),VP(nnodes),VS(nnodes),
     2  VFSQ(nnodes),FLAM(nnodes),FMU(
     1nnodes),G(nnodes),X(4),KDIS(28)

      G(1)=0.
      F=0.
      ii=1
c      DO 7001 I=2,N
c7001  VFSQ=1.e-7
      DO 10 I=2,N
c change start - compute g
c      if (i.eq.2) then
c      g(i)=con*4./3.*rho(i)*r(i)
c      else
c      g(i)=g(i-1)*(r(i-1)/r(i))**2
c     2 +4./3.*con*rho(i)*(r(i)-r(i-1)**3/r(i)**2)
c      print *, g(i)
c      endif
c      f=g(i)*(r(i)**2)/4./con
c      go to 10
c change end
      J=I-1
      RJI=R(I)-R(J)
c     test for density change - change start
c      go to 6
c     change end
C    1 IF(RHO(I)) 6,6,2
    1 IF(RHO(I) .LE. 0) THEN
      GOTO 6
      ELSE IF (RHO(I) .GT. 0) THEN
      GOTO 2
      ENDIF
C    2 IF(RJI) 5,5,3
    2 IF(RJI .LE. 0) THEN
      GOTO 5
      ELSE IF(RJI .GT. 0) THEN
      GOTO 3
      ENDIF
    3 DEL=RJI/3.
      X(1)=R(J)
      X(2)=X(1)+DEL
      X(4)=R(I)
      X(3)=X(4)-DEL
      DO 4 K=1,4
      DEL=(X(K)-R(J))/RJI
      RO=RHO(J)+DEL*(RHO(I)-RHO(J))
    4 X(K)=RO*X(K)*X(K)
      F=F+.125*RJI*(X(1)+X(4)+3.*(X(2)+X(3)))
    5 G(I)=4.*CON*F/(R(I)*R(I))
c   set g as constant 9.8 for one layer mantle model
C    5 G(I)=9.8
c test for density change - change start-end
      GO TO 10
      if (i.eq.(kdis(ii)+1)) then
      if (ii.lt.ndis) ii = ii+1
      go to 10
      else
      go to 6      
      endif
c      if (i.lt.34) go to 6
c      if ((i.gt.34) .and. (i.lt.67)) go to 6
c      if ((i.gt.67) .and. (i.lt.136)) go to 6
c      if (i.gt.136) go to 6
c      if (((i.gt.66) .and. (i.lt.136))) go to 6
c      if (i.gt.136) then
c      go to 6
c      else
      go to 10
c      endif
c change end
    6 IF(I.GT.2) GO TO 8
    7 ROP=0.
      GP=4.*CON*RHO(1)/3.
      GO TO 9
    8 C=VP(J)*VP(J)-4.*VS(J)*VS(J)/3.
      ROP=-RHO(J)*(VFSQ(J)/G(J)+G(J)/C)
      GP=4.*CON*RHO(J)-2.*G(J)/R(J)
    9 RHO(I)=RHO(J)
      G(I)=G(J)
      SR=.5*ROP
      SG=.5*GP
      ER=ROP
      EG=GP
c      print *, 1,RJI*SR
      RHO(I)=RHO(I)+RJI*SR
      G(I)=G(I)+RJI*SG
      RIJ=.5*(R(I)+R(J))
      V=.5*(VP(I)+VP(J))
      S=.5*(VS(I)+VS(J))
      C=V*V-4.*S*S/3.
      VF=.5*(VFSQ(I)+VFSQ(J))
      ROP=-RHO(I)*(VF/G(I)+G(I)/C)
      GP=4.*CON*RHO(I)-2.*G(I)/RIJ
      A=.292893218814
      SR=A*(ROP-ER)
      SG=A*(GP-EG)
      ER=ER+3.*SR-A*ROP
      EG=EG+3.*SG-A*GP
c      print *, 2,SR*RJI
      RHO(I)=RHO(I)+RJI*SR
      G(I)=G(I)+RJI*SG
      ROP=-RHO(I)*(VF/G(I)+G(I)/C)
      GP=4.*CON*RHO(I)-2.*G(I)/RIJ
      A=1.707106781186
      SR=A*(ROP-ER)
      SG=A*(GP-EG)
      ER=ER+3.*SR-A*ROP
      EG=EG+3.*SG-A*GP
c      print *, 3,SR*RJI, ER*RJI

      RHO(I)=RHO(I)+RJI*SR
      G(I)=G(I)+RJI*SG
      C=VP(I)*VP(I)-4.*VS(I)*VS(I)/3.
      ROP=-RHO(I)*(VFSQ(I)/G(I)+G(I)/C)
      GP=4.*CON*RHO(I)-2.*G(I)/R(I)
      A=1./6.
      SR=A*(ROP-2.*ER)
      SG=A*(GP-2.*EG)
      ER=ER+3.*SR-.5*ROP
      EG=EG+3.*SG-.5*GP
      RHO(I)=RHO(I)+RJI*(SR-ER/3.)
      G(I)=G(I)+RJI*(SG-EG/3.)
      F=G(I)*R(I)*R(I)/(4.*CON)
   10 CONTINUE
      RHOBAR=3.*F/(R(N)**3)

c      do 111  i=1,N
c      print *,G(i)
c  111 CONTINUE

      RETURN
      END
      SUBROUTINE ELLIP
C     FOR USE IN LOVE NUMBER PROGRAMS.
C     DOES NOT DEFINE DELTA M SUB ELL AND DOES NOT CALL DERIV.
      parameter(nnodes=200)
C      implicit real*16(a-h,o-z)
      implicit double precision(a-h,o-z)
      COMMON CON,R,RHO,VP,VS,VFSQ,FLAM,FMU,G,RHOBAR,W,WSQ,
     1DEL1,DEL2,QS,QF,Q,EPS,SFL3,FORD,VERTNO,PI,FL,
     2FL1,FL2,FL3,N1,N2,N,JCOM,L
      COMMON/F/C,NIR,NER
      COMMON/P/CAPW,EA,EAH,ELL
      COMMON/MP/FMU1,BULK1,RHO1,G1
      COMMON/R2/R2
      DIMENSION R(nnodes),RHO(nnodes),VP(nnodes),VS(nnodes),
     2  VFSQ(nnodes),FLAM(nnodes),FMU(
     1nnodes),G(nnodes),C(nnodes),R2(nnodes),P(nnodes),
     2  T(nnodes),ETA(nnodes),ELL(nnodes)
      NIR=1
      DO 1 I=1,N
    1 C(I)=RHO(I)
      NER=0
      DO 2 I=1,N
      NER=NER+1
      CALL GLINT(Q)
    2 P(I)=Q
      DO 3 I=1,N
    3 C(I)=RHO(I)*R2(I)
      NER=0
      DO 4 I=1,N
      NER=NER+1
      CALL GLINT(Q)
    4 T(I)=Q
      DO 5 I=2,N
      ZED=(2.*T(I))/(3.*R2(I)*P(I))
    5 ETA(I)=(25./4.)*(1.-1.5*ZED)*(1.-1.5*ZED)-1.
      ETA(1)=0.
      FM=.75*CAPW*CAPW
      EAH=(5.*FM)/(2.*(ETA(N)+2.))
      DO 6 I=2,N
    6 C(I)=ETA(I)/R(I)
      C(1)=0.
      NER=N
      CALL GLINT(Q)
      NER=0
      DO 7 I=1,N
      NER=NER+1
      CALL GLINT(Q1)
      ARG=Q1-Q
    7 ELL(I)=EAH*EXP(ARG)
      RETURN
      END
      SUBROUTINE GLINT(Q)
C     ALTERED DPFG VERSION  WORKS FOR ARBITRARY UPPER AND LOWER LIMITS
      parameter(nnodes=200)
C      implicit real*16(a-h,o-z)
      implicit double precision(a-h,o-z)
      COMMON CON,R,RHO,VP,VS,VFSQ,FLAM,FMU,G,RHOBAR,W,WSQ,
     1DEL1,DEL2,QS,QF,O,EPS,SFL3,FORD,VERTNO,PI,FL,
     2FL1,FL2,FL3,N1,N2,N,JCOM,L
      COMMON/F/F,NIR,NER
      COMMON/INT/B(nnodes)
      COMMON/R2/R2(nnodes)
      DIMENSION R(nnodes),RHO(nnodes),VP(nnodes),VS(nnodes),
     2  VFSQ(nnodes),FLAM(nnodes),FMU(
     1nnodes),G(nnodes),F(nnodes),A(nnodes),C(nnodes)
      Q=0.
      IF(NER.EQ.NIR) RETURN
      NIRP=NIR+1
      NERM=NER-1
      IF(NERM.LT.NIRP) GO TO 3
      DO 1 I=NIRP,NERM
    1 Q=F(I)*B(I)+Q
    3 Q=Q+F(NIR)*A(NIR)+F(NER)*C(NER)
      RETURN
      ENTRY PREINT
      QTR=.25
      C6=1./6.
      C12=1./12.
      Y=0.
      YSQ=0.
      B(1)=0.
      DO 2 I=2,N
      IM1=I-1
      X=Y
      Y=R(I)
      H=Y-X
      XSQ=YSQ
      YSQ=R2(I)
      XY=C6*X*Y
      AI=(QTR*XSQ+XY+C12*YSQ)*H
      A(I)=AI
      B(IM1)=AI+B(IM1)
      C(I)=(QTR*YSQ+XY+C12*XSQ)*H
    2 B(I)=C(I)
      A(1)=B(1)
      RETURN
      END
      SUBROUTINE SHANKS
c - change start
c I exchange label 6 and 8, to make it order 8 when IN(i.e. 'I' in this
c subroutine)=8.
c - change end

C     UNALTERED DPFG VERSION
C      implicit real*16(a-h,o-z)
      implicit double precision(a-h,o-z)
      COMMON/I/B,C,DX,I
      DIMENSION B(78),C(12)
      C(1)=0.
      GO TO (1,2,3,4,5,6,7,8),I
    1 B(1)=DX
      RETURN
    2 C(2)=DX
      B(1)=DX
      B(2)=.5*DX
      B(3)=B(2)
      RETURN
    3 C(2)=.5*DX
      C(3)=DX
      B(1)=C(2)
      B(2)=-DX
      B(3)=DX+DX
      B(4)=.1666666667*DX
      B(5)=.6666666667*DX
      B(6)=B(4)
      RETURN
    4 C(2)=.01*DX
      C(3)=.6*DX
      C(4)=DX
      B(1)=C(2)
      B(2)=-17.46122448979*DX
      B(3)=18.06122448979*DX
      B(4)=59.69127516778*DX
      B(5)=-60.53065635308*DX
      B(6)=1.839381185303*DX
      B(7)=-2.5555555556*DX
      B(8)=2.853392683901*DX
      B(9)=.5767419962335*DX
      B(10)=.1254208754209*DX
      RETURN
    5 C(2)=.25*DX
      C(3)=C(2)
      C(4)=.5*DX
      C(5)=.75*DX
      C(6)=DX
      B(1)=C(2)
      B(2)=.125*DX
      B(3)=.125*DX
      B(4)=0.
      B(5)=-C(4)
      B(6)=C(6)
      B(7)=.1875*DX
      B(8)=0.
      B(9)=0.
      B(10)=.5625*DX
      B(11)=-.42857142857*DX
      B(12)=.285714285714*DX
      B(13)=1.7142857143*DX
      B(14)=-B(13)
      B(15)=1.14285714286*DX
      B(16)=.077777777778*DX
      B(17)=0.
      B(18)=.355555555556*DX
      B(19)=.13333333333*DX
      B(20)=B(18)
      B(21)=B(16)
      I=6
      RETURN
    8 C(2)=.1111111111*DX
      C(3)=.1666666667*DX
      C(4)=.3333333333*DX
      C(5)=.5*DX
      C(6)=.6666666667*DX
      C(7)=.8333333333*DX
      C(8)=DX
      B(1)=C(2)
      B(2)=.04166666667*DX
      B(3)=.125*DX
      B(4)=.1666666667*DX
      B(5)=-.5*DX
      B(6)=.6666666667*DX
      B(7)=-.625*DX
      B(8)=3.375*DX
      B(9)=-3.*DX
      B(10)=.75*DX
      B(11)=24.55555556*DX
      B(12)=-109.*DX
      B(13)=96.33333333*DX
      B(14)=-11.33333333*DX
      B(15)=.1111111111*DX
      B(16)=-3.8125*DX
      B(17)=14.125*DX
      B(18)=-9.833333333*DX
      B(19)=-1.375*DX
      B(20)=1.666666667*DX
      B(21)=.0625*DX
      B(22)=8.731707317*DX
      B(23)=-25.35365854*DX
      B(24)=12.21951219*DX
      B(25)=10.17073171*DX
      B(26)=-5.536585366*DX
      B(27)=-.1097560976*DX
      B(28)=.8780487805*DX
      B(29)=.04880952381*DX
      B(30)=0.
      B(31)=.2571428571*DX
      B(32)=.03214285714*DX
      B(33)=.3238095238*DX
      B(34)=B(32)
      B(35)=B(31)
      B(36)=B(29)
      I=8
      RETURN
    7 C(2)=.2222222222*DX
      C(3)=.3333333333*DX
      C(4)=.5*DX
      C(5)=.1666666667*DX
      C(6)=.8888888889*DX
      C(7)=.1111111111*DX
      C(8)=.8333333333*DX
      C(9)=DX
      B(1)=C(2)
      B(2)=.08333333333*DX
      B(3)=.25*DX
      B(4)=.125*DX
      B(5)=0.
      B(6)=.375*DX
      B(7)=.1064814814*DX
      B(8)=0.
      B(9)=.09722222222*DX
      B(10)=-.03703703704*DX
      B(11)=-5.673525377*DX
      B(12)=0.
      B(13)=-18.63374486*DX
      B(14)=7.22085048*DX
      B(15)=17.97530864*DX
      B(16)=.693329904*DX
      B(17)=0.
      B(18)=1.991769547*DX
      B(19)=-.7105624143*DX
      B(20)=-1.874643875*DX
      B(21)=.01121794872*DX
      B(22)=-.5634259259*DX
      B(23)=0.
      B(24)=-2.013888889*DX
      B(25)=1.261073318*DX
      B(26)=1.851282051*DX
      B(27)=.05951726845*DX
      B(28)=.2387755102*DX
      B(29)=.09356936416*DX
      B(30)=0.
      B(31)=-.4855491329*DX
      B(32)=-.08092485549*DX
      B(33)=2.761227212*DX
      B(34)=-.3964976497*DX
      B(35)=-1.852251794*DX
      B(36)=.9604268564*DX
      B(37)=.05148809524*DX
      B(38)=0.
      B(39)=0.
      B(40)=.3587949466*DX
      B(41)=.2967032967*DX
      B(42)=-.02758886522*DX
      B(43)=B(42)
      B(44)=B(41)
      B(45)=B(37)
      I=9
      RETURN
    6 C( 2)= .11111111111*DX
      C( 3)= .16666666667*DX
      C( 4)= .25*DX
      C( 5)= .1*DX
      C( 6)= C(3)
      C( 7)= .5*DX
      C( 8)= .666666666667*DX
      C( 9)= .33333333333*DX
      C(10)= .83333333333*DX
      C(11)= C(10)
      C(12)= DX
      B( 1)= C(2)
      B( 2)= .041666666667*DX
      B( 3)= .125*DX
      B( 4)= .0625*DX
      B( 5)= 0.
      B( 6)= .1875*DX
      B( 7)= .058*DX
      B( 8)= 0.
      B( 9)= .066*DX
      B(10)= -.024*DX
      B(11)= .033950617284*DX
      B(12)= 0.
      B(13)= 0.
      B(14)= .0041152263374*DX
      B(15)= .12860082305*DX
      B(16)= -.58333333333*DX
      B(17)= 0.
      B(18)= 0.
      B(19)= 2.1111111111*DX
      B(20)= 3.4722222222*DX
      B(21)= -4.5*DX
      B(22)= -.12345678901*DX
      B(23)= 0.
      B(24)= 0.
      B(25)= -.1316872428*DX
      B(26)= .51440329218*DX
      B(27)= 0.
      B(28)= .40740740741*DX
      B(29)= 3.6265432099*DX
      B(30)= 0.
      B(31)= 0.
      B(32)= -10.666666667*DX
      B(33)= -19.290123457*DX
      B(34)= 26.*DX
      B(35)= .74691358025*DX
      B(36)= -.083333333333*DX
      B(37)= .90432098765*DX
      B(38)= 0.
      B(39)= 0.
      B(40)= -2.6296296296*DX
      B(41)= -4.2438271605*DX
      B(42)= 5.6666666667*DX
      B(43)= -.36419753086*DX
      B(44)= .5*DX
      B(45)= DX
      B(46)= .80432098765*DX
      B(47)= 0.
      B(48)= 0.
      B(49)= -2.6296296296*DX
      B(50)= -4.2438271605*DX
      B(51)= 6.1666666667*DX
      B(52)= .63580246914*DX
      B(53)= 0.
      B(54)= 0.
      B(55)= .1*DX
      B(56)= -1.9410569106*DX
      B(57)= 0.
      B(58)= 0.
      B(59)= 6.9376693767*DX
      B(60)= 11.009485095*DX
      B(61)= -14.926829268*DX
      B(62)= .085365853659*DX
      B(63)= -.16463414634*DX
      B(64)= -.43902439024*DX
      B(65)= -.29268292683*DX
      B(66)= .73170731707*DX
      B(67)= .04880952381*DX
      B(68)= 0.
      B(69)= 0.
      B(70)= 0.
      B(71)= 0.
      B(72)= .25714285714*DX
      B(73)= .32380952381*DX
      B(74)= .032142857143*DX
      B(75)= B(74)
      B(76)= .042857142857*DX
      B(77)= .21428571429*DX
      B(78)=B(67)
      I=12
      RETURN
      END
      SUBROUTINE STEPS(EPS,VERTNO)
C     UNALTERED DPFG VERSION
C      implicit real*16(a-h,o-z)
      implicit double precision(a-h,o-z)
      COMMON/C/STEP
      DIMENSION STEP(8)
c      PS=qLOG(EPS)
      PS=DLOG(EPS)
      VERTNO=-PS
      FAC=1.
      DO 2 N=1,8
      FN=N+1
      FAC=FAC*FN
c      X=(QLOG(FAC)+PS)/FN
      X=(DLOG(FAC)+PS)/FN
      X=EXP(X)
      S=X
      DO 1 I=1,N
    1 S=X*EXP(-S/FN)
    2 STEP(N)=S
      RETURN
      END
      SUBROUTINE MATADJ(C,CA)
c      implicit real*16(a-h,o-z)
      implicit double precision(a-h,o-z)
      DIMENSION C(3,3),CA(3,3)
      CA(1,1)=C(2,2)*C(3,3)-C(3,2)*C(2,3)
      CA(1,2)=C(3,2)*C(1,3)-C(1,2)*C(3,3)
      CA(1,3)=C(1,2)*C(2,3)-C(2,2)*C(1,3)
      CA(2,1)=C(3,1)*C(2,3)-C(2,1)*C(3,3)
      CA(2,2)=C(1,1)*C(3,3)-C(3,1)*C(1,3)
      CA(2,3)=C(2,1)*C(1,3)-C(1,1)*C(2,3)
      CA(3,1)=C(2,1)*C(3,2)-C(3,1)*C(2,2)
      CA(3,2)=C(3,1)*C(1,2)-C(1,1)*C(3,2)
      CA(3,3)=C(1,1)*C(2,2)-C(2,1)*C(1,2)
      RETURN
      END
      SUBROUTINE SPSJFG(I,A)
C     UNALTERED DPFG VERSION.
      parameter(nnodes=200)
C      implicit real*16(a-h,o-z)
      implicit double precision(a-h,o-z)
      COMMON CON,R,RHO,VP,VS,VFSQ,FLAM,FMU,G,RHOBAR,W,WSQ,
     1DEL1,DEL2,QS,QF,Q,EPS,SFL3,FORD,VERTNO,PI,FL,
     2FL1,FL2,FL3,N1,N2,N,JCOM,L
      COMMON/J/KG
      DIMENSION R(nnodes),RHO(nnodes),VP(nnodes),VS(nnodes),
     2  VFSQ(nnodes),FLAM(nnodes),FMU(
     1nnodes),G(nnodes),A(6,3)
      X=R(I)
      RO=RHO(I)
      GR=G(I)
      P=VP(I)
      VPSQ=P*P
      P=VS(I)
      VSSQ=P*P
      ZETA=4.*CON*RO
      XI=GR/X
      XSQ=X*X
      ALFSQ=(WSQ+ZETA+XI)/VPSQ
      BETASQ=WSQ/VSSQ
      GAMSQ=4.*FL3*XI*XI/(VPSQ*VSSQ)
      DELSQ=SQRT((BETASQ-ALFSQ)*(BETASQ-ALFSQ)+GAMSQ)
      FKSQ=.5*(ALFSQ+BETASQ+DELSQ)
      QSQ=FKSQ-DELSQ
      XL=X**L
      XLP1=XL*X
      XLP2=XLP1*X
      XLM1=XL/X
      FKXSQ=FKSQ*XSQ
      QXSQ=QSQ*XSQ
      QK=QXSQ*FKXSQ
      QPK=QXSQ+FKXSQ
      FLU=RO*VPSQ
      FU=RO*VSSQ
      D0=1.
      H0=-VPSQ/(FL1*VSSQ)
      K=1
    5 C=2.
      B=FL2+2.
      F=2.
      C2=1./(C*B)
      D1=C2*(FL3*XI*H0/VPSQ-ALFSQ*D0)*XSQ
      H1=C2*(XI*D0/VSSQ-BETASQ*H0)*XSQ
      U=C2*(FL3*H0+(FL+F)*D0)
      V=C2*(D0+(FL1+F)*H0)
      P=C2*D0
      S=(FL2+F)*P-U
      H=H0+H1
      D=D0+D1
    6 C1=C2
      C=C+2.
      B=B+2.
      F=F+2.
      C2=1./(C*B)
      UN=C2*(FL3*H1+(FL+F)*D1)
      VN=C2*(D1+(FL1+F)*H1)
      PN=C2*D1
      SN=(FL2+F)*PN-UN
      D2=-C2*(D1*ALFSQ-H1*FL3*XI/VPSQ)*XSQ
      H2=-C2*(H1*BETASQ-D1*XI/VSSQ)*XSQ
      D=D+D2
      H=H+H2
      D0=D1
      D1=D2
      H0=H1
      H1=H2
      U=U+UN
      V=V+VN
      P=P+PN
      S=S+SN
      TE=ABS (D2/D)
c      IF(TE-EPS) 7,6,6
      IF(TE-EPS .LT. 0) THEN
      GOTO 7
      ELSE IF(TE-EPS .GE. 0) THEN
      GOTO 6
      ENDIF
    7 TE=ABS(H2/H)
C      IF(TE-EPS) 8,6,6
      IF(TE-EPS .LT. 0) THEN
      GOTO 8
      ELSE IF(TE-EPS .GE. 0) THEN
      GOTO 6
      ENDIF
    8 C=C+2.
      B=B+2.
      F=F+2.
      C2=1./(C*B)
      UN=C2*(FL3*H1+(FL+F)*D1)
      VN=C2*(D1+(FL1+F)*H1)
      PN=C2*D1
      SN=(FL2+F)*PN-UN
      A(1,K)=(U+UN)*XLP1
      A(2,K)=(V+VN)*XLP1
      A(3,K)=FLU*D*XL+2.*FU*(FL3*A(2,K)-2.*A(1,K))/X
      A(4,K)=FU*(H*XL+2.*(A(1,K)-A(2,K))/X)
      A(5,K)=ZETA*(P+PN)*XLP2
      A(6,K)=ZETA*(S+SN)*XLP1
      GO TO (9,10),K
    9 K=2
      D0=0.
      H0=-1.
      GO TO 5
   10 A(1,3)=XLM1*FL
      A(2,3)=XLM1
      A(3,3)=2.*FU*(FL3*A(2,3)-2.*A(1,3))/X
      A(4,3)=2.*FU*(A(1,3)-A(2,3))/X
      B=XI*FL-WSQ
      A(5,3)=B*XL
      A(6,3)=(FL2*B-ZETA*FL)*XLM1
      A(5,2)=A(5,2)+FL1*VSSQ*XL
      A(6,2)=A(6,2)+FL2*FL1*VSSQ*XLM1
      IS=I
      DO 20 I=1,3,2
      A(2,I)=SFL3*A(2,I)
   20 A(4,I)=SFL3*A(4,I)
      DO 21 I=1,5,2
   21 A(I,2)=A(I,2)/SFL3
      A(6,2)=A(6,2)/SFL3
      I=IS
      IF(KG.EQ.1) RETURN
      DO 30 J=1,6
   30 A(J,3)=0.
      A(5,1)=0.
      A(5,2)=0.
      A(6,1)=0.
      A(6,2)=0.
      RETURN
      END

      SUBROUTINE CCOEF(X,I,A)
      parameter(nnodes=200)
C      implicit real*16(a-h,o-z)
      implicit double precision(a-h,o-z)
      COMMON CON,R,RHO,VP,VS,VFSQ,FLAM,FMU,G,RHOBAR,W,WSQ,
     1DEL1,DEL2,QS,QF,Q,EPS,SFL3,FORD,VERTNO,PI,FL,
     2FL1,FL2,FL3,N1,N2,N,JCOM,L
      DIMENSION R(nnodes),RHO(nnodes),VP(nnodes),VS(nnodes),
     2  VFSQ(nnodes),FLAM(nnodes),FMU(
     1nnodes),G(nnodes),A(2,2)
      J=I
    1 J=J-1
C      IF(R(J)-R(I)) 2,1,1
      IF(R(J)-R(I) .LT. 0) THEN
      GOTO 2
      ELSE IF(R(J)-R(I) .GE. 0) THEN
      GOTO 1
      ENDIF
    2 DEL=(X-R(J))/(R(I)-R(J))
      RO=RHO(J)+DEL*(RHO(I)-RHO(J))
      GR=G(J)+DEL*(G(I)-G(J))
      ROPR=(RHO(I)-RHO(J))/(R(I)-R(J))
      Z=1./X
      ZSQ=Z*Z
      A(1,1)=0.
      A(1,2)=1.
      A(2,1)=FL3*ZSQ+4.*CON*ROPR/GR
      A(2,2)=-2.*Z
      RETURN
      END

      SUBROUTINE CMTR(X,Y,U,V)
      parameter(nnodes=200)
C      implicit real*16(a-h,o-z)
      implicit double precision(a-h,o-z)
      COMMON CON,R,RHO,VP,VS,VFSQ,FLAM,FMU,G,RHOBAR,W,WSQ,
     1DEL1,DEL2,QS,QF,Q,EPS,SFL3,FORD,VERTNO,PI,FL,
     2FL1,FL2,FL3,N1,N2,N,JCOM,L
      COMMON/C/STEP
      COMMON/I/B,C,DX,IN
      COMMON/E/NDIS,KDIS
      DIMENSION R(nnodes),RHO(nnodes),VP(nnodes),VS(nnodes),
     2  VFSQ(nnodes),FLAM(nnodes),FMU(
     1nnodes),G(nnodes),F(2),A(2,2),S(2),C(12),B(78),H(12,2),
     2  STEP(8),KDIS(28)
C      IF(Y-X) 1,1,2
      IF(Y-X .LE. 0) THEN
      GOTO 1
      ELSE IF(Y-X .GT. 0) THEN
      GOTO 2
      ENDIF
    1 RETURN
    2 YS=Y
      XS=X
      NOC=KDIS(2)
      F(1)=U
      F(2)=V
      I=2
C    3 IF(Y-R(I)) 5,5,4
    3 IF(Y-R(I) .LE. 0) THEN
      GOTO 5
      ELSE IF(Y-R(I) .GT. 0) THEN
      GOTO 4
      ENDIF
    4 I=I+1
C      IF(I-NOC) 3,5,5
      IF(I-NOC .LT. 0) THEN
      GOTO 3
      ELSE IF(I-NOC .GE. 0) THEN
      GOTO 5
      ENDIF
    5 RO=RHO(I)
      FLAMB=FLAM(I)
    6 QSQ=ABS(4.*CON*RO*RO/FLAMB-FL3/(X*X))
      Q=SQRT(QSQ)+1./X
      DX=STEP(8)/Q
      Y=X+DX
C      IF(Y-YS) 11,11,10
      IF(Y-YS .LE. 0) THEN
      GOTO 11
      ELSE IF(Y-YS .GT. 0) THEN
      GOTO 10
      ENDIF
   10 Y=YS
   11 DX=Y-X
      DS=Q*DX
      DO 13 J=1,7
      IF(DS.LE.STEP(J)) GO TO 12
      GO TO 13
   12 IN=J
      GO TO 14
   13 CONTINUE
      IN=8
   14 CALL SHANKS
      S(1)=F(1)
      S(2)=F(2)
      DO 25 NI=1,IN
      Z=X+C(NI)
      CALL CCOEF(Z,I,A)
      H(NI,1)=A(1,1)*F(1)+A(1,2)*F(2)
      H(NI,2)=A(2,1)*F(1)+A(2,2)*F(2)
      F(1)=S(1)
      F(2)=S(2)
      DO 25 M=1,NI
      K1=M+NI*(NI-1)/2
      F(1)=F(1)+B(K1)*H(M,1)
   25 F(2)=F(2)+B(K1)*H(M,2)
      X=Y
c      IF(Y-YS) 6,26,26
      IF(Y-YS .LT. 0) THEN
      GOTO 6
      ELSE IF(Y-YS .GE. 0) THEN
      GOTO 26
      ENDIF
   26 X=XS
      Y=YS
      U=F(1)
      V=F(2)
      RETURN
      END
      SUBROUTINE SCOEF(X,I,C)
C     CONVENTION USED IS DELSQ PHI = 4 PI G RHO.
      parameter(nnodes=200)
C      implicit real*16(a-h,o-z)
      implicit double precision(a-h,o-z)
      COMMON CON,R,RHO,VP,VS,VFSQ,FLAM,FMU,G,RHOBAR,W,WSQ,
     1DEL1,DEL2,QS,QF,Q,EPS,SFL3,FORD,VERTNO,PI,FL,
     2FL1,FL2,FL3,N1,N2,N,JCOM,L
      DIMENSION R(nnodes),RHO(nnodes),VP(nnodes),VS(nnodes),
     2  VFSQ(nnodes),FLAM(nnodes),FMU(
     1nnodes),G(nnodes),A(6,6),B(6,6),C(6,6)
      J=I
    1 J=J-1
C      IF(R(J)-R(I)) 2,1,1
      IF(R(J)-R(I) .LT. 0) THEN
      GOTO 2
      ELSE IF(R(J)-R(I) .GE. 0) THEN
      GOTO 1
      ENDIF
    2 DEL=(X-R(J))/(R(I)-R(J))
      RO=RHO(J)+DEL*(RHO(I)-RHO(J))
      FU=FMU(J)+DEL*(FMU(I)-FMU(J))
      FLU=FLAM(J)+DEL*(FLAM(I)-FLAM(J))
      GR=G(J)+DEL*(G(I)-G(J))
      DO 3 K=1,6
      DO 3 J=1,6
      A(J,K)=0.
      B(J,K)=0.
    3 C(J,K)=0.
      D=1./(FLU+FU+FU)
      E=(3.*FLU+FU+FU)*D*FU
      ZETA=4.*CON*RO
      A(3,1)=4.*E
      A(4,1)=-2.*E
      A(3,2)=A(4,1)*FL3
      A(4,2)=2.*FU*(2.*(FLU+FU)*FL3*D-1.)
      B(1,1)=-2.*FLU*D
      B(2,1)=-1.
      B(3,1)=-4.*RO*GR
      B(4,1)=RO*GR
      B(6,1)=ZETA*FL1
      B(1,2)=FLU*D*FL3
      B(2,2)=1.
      B(3,2)=B(4,1)*FL3
      B(6,2)=-ZETA*FL3
      B(3,3)=-4.*FU*D
      B(4,3)=-FLU*D
      B(3,4)=FL3
      B(4,4)=-3.
      B(3,5)=RO*FL1
      B(4,5)=-RO
      B(5,5)=-FL1
      B(6,6)=FL-1.
      C(3,1)=-RO*WSQ
      C(5,1)=ZETA
      C(4,2)=C(3,1)
      C(1,3)=D
      C(2,4)=1./FU
      C(3,6)=-RO
      C(5,6)=1.
      Z=1./X
      ZSQ=Z*Z
      C(1,1)=Z*B(1,1)
      C(2,1)=Z*B(2,1)
      C(3,1)=ZSQ*A(3,1)+Z*B(3,1)+C(3,1)
      C(4,1)=ZSQ*A(4,1)+Z*B(4,1)
      C(6,1)=Z*B(6,1)
      C(1,2)=Z*B(1,2)
      C(2,2)=Z*B(2,2)
      C(3,2)=ZSQ*A(3,2)+Z*B(3,2)
      C(4,2)=ZSQ*A(4,2)+C(4,2)
      C(6,2)=Z*B(6,2)
      C(3,3)=Z*B(3,3)
      C(4,3)=Z*B(4,3)
      C(3,4)=Z*B(3,4)
      C(4,4)=Z*B(4,4)
      C(3,5)=Z*B(3,5)
      C(4,5)=Z*B(4,5)
      C(5,5)=Z*B(5,5)
      C(6,6)=Z*B(6,6)
      C(3,5)=-C(3,5)
      C(3,6)=-C(3,6)
      C(4,5)=-C(4,5)
      C(5,1)=-C(5,1)
      C(6,1)=-C(6,1)
      C(6,2)=-C(6,2)
      RETURN
      END

      SUBROUTINE SMTR(X,Y,F)
      parameter(nnodes=200)
C      implicit real*16(a-h,o-z)
      implicit double precision(a-h,o-z)
      COMMON CON,R,RHO,VP,VS,VFSQ,FLAM,FMU,G,RHOBAR,W,WSQ,
     1DEL1,DEL2,QS,QF,Q,EPS,SFL3,FORD,VERTNO,PI,FL,
     2FL1,FL2,FL3,N1,N2,N,JCOM,L
      COMMON/C/STEP
      COMMON/I/B,C,DX,IN
      COMMON/J/KG
      DIMENSION R(nnodes),RHO(nnodes),VP(nnodes),VS(nnodes),
     2  VFSQ(nnodes),FLAM(nnodes),FMU(
     1nnodes),G(nnodes),F(6,3),A(6,6),S(6,3),C(12),
     2  B(78),H(12,6,3),STEP(8)
C      IF(Y-X) 1,1,2
      IF(Y-X .LE. 0) THEN
      GOTO 1
      ELSE IF(Y-X .GT. 0) THEN
      GOTO 2
      ENDIF
    1 RETURN
    2 YS=Y
      XS=X
      KK=4-KG
      JJ=12+KG*(2*KG-8)
      I=2
c    3 IF(Y-R(I)) 5,5,4
    3 IF(Y-R(I) .LE. 0) THEN
      GOTO 5
      ELSE IF(Y-R(I) .GT. 0) THEN
      GOTO 4
      ENDIF
    4 I=I+1
C      IF(I-N) 3,5,5
      IF(I-N .LT. 0) THEN
      GOTO 3
      ELSE IF(I-N .GE. 0) THEN
      GOTO 5
      ENDIF
    5 V=VP(I)
      VPSQ=V*V
      VSSQ=VS(I)*VS(I)
      ZETA=4.*CON*RHO(I)
      RI=R(I)
      XI=G(I)/RI
      ALFSQ=(WSQ+ZETA+XI)/VPSQ
      BETASQ=WSQ/VSSQ
      GAMSQ=4.*FL3*XI*XI/(VSSQ*VPSQ)
      DELSQ=SQRT((BETASQ-ALFSQ)*(BETASQ-ALFSQ)+GAMSQ)
      FKSQ=.5*(ALFSQ+BETASQ+DELSQ)
      QSQ=FKSQ-DELSQ
      SFL3=SQRT(FL3)
    6 Q=SFL3/X
      QS=SQRT(ABS(FKSQ-FL3/(X*X)))+1./X
      QF=SQRT(ABS(QSQ-FL3/(X*X)))+1./X
      GO TO (61,7,60),KG
   60 Q=QS+QF
      GO TO 9
C   61 IF(Q-QF) 7,8,8
   61 IF(Q-QF .LT. 0) THEN
      GOTO 7
      ELSE IF(Q-QF .GE. 0) THEN
      GOTO 8
      ENDIF 
    7 Q=QF
C    8 IF(Q-QS) 81,9,9
    8 IF(Q-QS .LT. 0) THEN
      GOTO 81
      ELSE IF(Q-QS .GE. 0) THEN
      GOTO 9
      ENDIF
   81 Q=QS
    9 DX=STEP(8)/Q
      Y=X+DX
C      IF(Y-YS) 11,11,10
      IF(Y-YS .LE. 0) THEN
      GOTO 11
      ELSE IF(Y-YS .GT. 0) THEN
      GOTO 10
      ENDIF
   10 Y=YS
   11 DX=Y-X
      DS=Q*DX
      DO 13 J=1,7
      IF(DS.LE.STEP(J)) GO TO 12
      GO TO 13
   12 IN=J
      GO TO 14
   13 CONTINUE
      IN=8
   14 CALL SHANKS
      DO 25 J=1,JJ
      DO 25 K=1,KK
   25 S(J,K)=F(J,K)
      DO 27 NI=1,IN
      Z=X+C(NI)
      CALL SCOEF(Z,I,A)
      DO 26 J=1,JJ
      DO 26 K=1,KK
      H(NI,J,K)=0.
      DO 26 M=1,JJ
   26 H(NI,J,K)=H(NI,J,K)+A(J,M)*F(M,K)
      DO 27 J=1,JJ
      DO 27 K=1,KK
      F(J,K)=S(J,K)
      DO 27 M=1,NI
      K1=M+NI*(NI-1)/2
   27 F(J,K)=F(J,K)+B(K1)*H(M,J,K)
      X=Y
C      IF(Y-YS) 6,29,29
      IF(Y-YS .LT. 0) THEN
      GOTO 6
      ELSE IF(Y-YS .GE. 0) THEN
      GOTO 29
      ENDIF
   29 X=XS
      Y=YS
      RETURN
      END


      SUBROUTINE zlubksb(a,n,np,indx,b)
c      implicit real*16 (a-h,o-z)
      implicit double precision(a-h,o-z)
      INTEGER n,np,indx(n)
      dimension a(np,np),b(n)
      INTEGER i,ii,j,ll
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14    continue
      return
      END

      SUBROUTINE zludcmp(a,n,np,indx,d)
c      implicit real*16 (a-h,o-z)
      implicit double precision(a-h,o-z)
      INTEGER n,np,indx(n),NMAX
      dimension a(np,np)
      PARAMETER (NMAX=500,TINY=1.0e-20)
      INTEGER i,imax,j,k
      dimension vv(NMAX)
      d=1.
      do 12 i=1,n
        aamax=0.
        do 11 j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11      continue
c        if (aamax.eq.0.) pause 'singular matrix in ludcmp'
        if (aamax.eq.0) print *, 'singular matrix in ludcmp'
        vv(i)=1./aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0.
        do 16 i=j,n
          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.)a(j,j)=TINY
        if(j.ne.n)then
          dum=1./a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      return
      END
