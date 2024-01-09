      SUBROUTINE SECEX(X)
      IMPLICIT REAL*8(A-H,O-Z)
c      parameter(nl=340)
      COMMON/SECX/CSEC(20),DCSDX(20),CISEC(20),DCISDX
      DIMENSION C(20),CI(20)
      XE=X
C     IF(X.GT.1.) XE=0.99
CF 1016
      IF(X.lT.0.) XE=1.e-4
      CALL EXSEC(XE,C,CI)
      DO I=1,20
         CSEC(I)=C(I)
      enddo
      DO I=1,20
         CISEC(I)=CI(I)
      enddo
c!!! no sec ion!!!!
      DO I=1,20
         CISEC(I)=1.d-20*cisec(i)
      enddo
      RETURN
      END

      SUBROUTINE SEC(ESEC,XEL,HEFF,SEL)
      IMPLICIT REAL*8(A-H,O-Z)
C      ************************************************************
C      ******
C      THIS ROUTINE CALCULATES THE HEATING EFFICIENCY AND NUMBER OF
C      ADDITIONAL IONIZATIONS DUE TO SECONDARY ELECTRONS, WITH
C      ENERGY ESEC. THE FITS ARE INTERPOLATIONS TO MONTE-CARLO CALC-
C      CULATIONS.
C      ******
C      ************************************************************
      Y=LOG10(XEL)
      E=ESEC
      HEFF=1.
      SEL=1.
      IF(E.LT.10.2) GOTO 10
      IF(XEL.GE.1.) GOTO 10
      IF(Y.LT.-0.5) GOTO 11
      HEFF=81.*Y/E**2.-10.3*Y/E+1.+.2214*Y
      GOTO 10
  11  IF(Y.LT.-1.) GOTO 12
      HEFF=(28.6+138.2*Y)/E**2.+(-4.0-18.3*Y)/E+1.1104+.4422*Y
      GOTO 10
  12  IF(Y.LT.-1.5) GOTO 13
      HEFF=(-169.4-59.8*Y)/E**2.+(17.5+3.2*Y)/E+.9584+.2902*Y
      GOTO 10
  13  IF(Y.LT.-2.0) GOTO 14
      HEFF=(-149.7-46.6*Y)/E**2.+(13.9+0.8*Y)/E+1.1015+.3856*Y
      GOTO 10
  14  IF(Y.LT.-2.5) GOTO 15
      HEFF=(-386.3-164.9*Y)/E**2.+(41.6+14.6*Y)/E+.7027+.1862*Y
      GOTO 10
  15  IF(Y.LT.-3.0) GOTO 16
      HEFF=(138.3+44.9*Y)/E**2.+(-3.4-3.4*Y)/E+.4612+.0896*Y
      GOTO 10
  16  IF(Y.LT.-4.0) GOTO 17
      HEFF=(-3.9-2.5*Y)/E**2.+(13.3+2.2*Y)/E+.3313+.0463*Y
      GOTO 10
  17  HEFF=6.1/E**2.+4.5/E+.1461
  10  IF(HEFF.GT.1.) HEFF=1.
      RETURN
      END


      SUBROUTINE EXSEC(X,CEX,CION)
C
C     1 = HEI, 2 = O I,  3 = O II,   4 = CA I, 5 = CA II, 6 = MG I
C     7 = C I, 8 = NA I, 9 = MG II, 10 = S I, 11 = SI I  12 = NE I
C    13 = NE II 14 = AR I 15 = AR II 16 = FE I 17 = FE II18 = H I
C     ALL RATES SHOULD BE DIVIDED BY DENS**2
C
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 ION
      include 'param'
c      parameter(nl=340)
      COMMON/NION/IONq
      COMMON/TPAR/RIN,DRQ,R,TDAY
      COMMON/GAMPAR/GAMLUM,GAELH,GAHE
      COMMON/HYDROS/HSCALE,DEN1,XEL,AMEAN
      common/ind/ik
      common/ionx/xion(md,14,27)
      common/abl/abn(15)
      DIMENSION CEX(20),CION(20),ION(20)
      dimension axsi(20)
      DATA PI/3.1415926E0/,ELCH/1.60219E-12/,AMU/1.660531E-24/
      DO K=1,20
         CION(K)=1.e-30
         cex(K)=1.e-30
         ion(k)=-100.
      ENDDO
      gaelh=0.
c      gahe=0.
      
      he=1.
      hyex=0.
c  hydrogen or helium dominated
      if(abn(1).gt.0.1.or.abn(2).gt.0.1) then
         heat=1.-toto_ex_ion
         abhi=abn(1)*xion(ik,1,1)
         abhei=abn(2)*xion(ik,2,1)
         call f_nonthermal(abhi,abhei,x,f_heat,
     &     f_ion_he,f_ex_he,f_ion_h,f_ex_h)         
         GAELH=f_heat*GAHE
         he=f_heat
         cion(1)=f_ion_he*gahe/(elch*24.5)
         cion(18)=f_ion_h*gahe/(elch*13.6)
         if(ionq.eq.5) then
c  use approximate excitation rates from KF92
            cex(2)=0.82*gahe*f_ex_h/(elch*10.2)
            cex(3)=0.12*gahe*f_ex_h/(elch*12.09)
            cex(4)=0.05*gahe*f_ex_h/(elch*12.75)
            cex(5)=0.03*gahe*f_ex_h/(elch*13.056)
         endif
         if(ionq.eq.3) then
            cex(2)=gahe*f_ex_he/(elch*19.2)
         endif
         goto 6543
      endif
c                         !!      IF(X.GT.1.) GOTO 6543
c      IF(X.GT.0.98) GOTO 6543
c      IF(X.GT.0.999) GOTO 6543
      if(x.lt.0.98) then
         XL=LOG10(X)
         iext=0
      else
         xi=0.98
         xl=log10(xi)
         iext=1
      endif
c     HYDROGEN DOMINATED
      ESEC=1.E3
c     CALL SEC(ESEC,X,HEFF,SEL)
c     IF(ABn(2).GT.0.3) GOTO 60
c     OXYGEN DOMINATED
      ION(16)=-1.523
      ION(17)=-1.523
      IF(XL.GT.-1.00000) GOTO 50
      HE= 0.85012+0.27247*XL+0.02703*XL**2
      ION( 1)=-2.44824-0.09560*XL-0.01888*XL**2
      ION( 2)=-2.05749-0.24897*XL-0.04278*XL**2
      ION( 3)=-2.46878-0.15229*XL-0.03719*XL**2
      ION( 4)=-1.66678-0.13942*XL-0.01920*XL**2
      ION( 5)=-2.25347-0.07530*XL-0.01100*XL**2
      ION( 6)=-1.83600-0.13622*XL-0.02602*XL**2
      ION( 7)=-1.80086-0.10583*XL-0.01849*XL**2
      ION( 8)=-1.80164-0.09341*XL-0.00618*XL**2
      ION( 9)=-2.36815-0.03717*XL-0.00207*XL**2
      ION(10)=-1.62796-0.06683*XL-0.00966*XL**2
      ION(11)=-1.62177-0.10160*XL-0.01563*XL**2
      ION(12)=-1.93059-0.03954*XL-0.00591*XL**2
      ION(14)=-1.56809-0.01799*XL+0.00091*XL**2
      GOTO  60
 50   IF(XL.GT.-0.30103) GOTO 55
      HE= 0.82860+0.21277*XL-0.01115*XL**2
      ION( 1)=-2.57233-0.21134*XL-0.01054*XL**2
      ION( 2)=-2.18307-0.43588*XL-0.10411*XL**2
      ION( 3)=-2.53715-0.26787*XL-0.08441*XL**2
      ION( 4)=-1.74414-0.25944*XL-0.06185*XL**2
      ION( 5)=-2.33686-0.03342*XL+0.11426*XL**2
      ION( 6)=-1.95195-0.34338*XL-0.11724*XL**2
      ION( 7)=-1.99154-0.51203*XL-0.23401*XL**2
      ION( 8)=-1.96099-0.41032*XL-0.16374*XL**2
      ION( 9)=-2.52601-0.28823*XL-0.09526*XL**2
      ION(10)=-1.82242-0.41406*XL-0.16242*XL**2
      ION(11)=-1.75820-0.30483*XL-0.08243*XL**2
      ION(12)=-2.08786-0.17835*XL+0.01254*XL**2
      ION(14)=-1.73749-0.13983*XL+0.04846*XL**2
      GOTO  60
 55   IF(XL.GT. 0.00004) GOTO 56
      HE= 0.83097+0.19072*XL-0.11060*XL**2
      ION( 1)=-2.62139-0.55947*XL-0.62560*XL**2
      ION( 2)=-2.23986-0.77837*XL-0.61517*XL**2
      ION( 3)=-2.58118-0.55265*XL-0.54455*XL**2
      ION( 4)=-1.82371-0.62313*XL-0.39196*XL**2
      ION( 5)=-2.45797-0.50263*XL-0.10790*XL**2
      ION( 6)=-1.98587-0.47955*XL-0.19531*XL**2
      ION( 7)=-1.99122-0.71637*XL-0.91636*XL**2
      ION( 8)=-1.98230-0.54215*XL-0.36654*XL**2
      ION( 9)=-2.59339-0.28247*XL+0.66738*XL**2
      ION(10)=-1.83154-0.49657*XL-0.33597*XL**2
      ION(11)=-1.80315-0.28671*XL+0.47372*XL**2
      ION(12)=-2.15078-0.52567*XL-0.44686*XL**2
      ION(14)=-1.82305-0.46815*XL-0.09795*XL**2
      GOTO  60
 56   HE= 0.83097+0.35645*XL+0.19409*XL**2
      ION( 1)=-2.62141-0.08512*XL-0.39237*XL**2
      ION( 2)=-2.23989-0.01840*XL-0.47092*XL**2
      ION( 3)=-2.58120-0.09856*XL-1.10769*XL**2
      ION( 4)=-1.82374+0.00649*XL-0.86178*XL**2
      ION( 5)=-2.45799-0.11213*XL-0.27513*XL**2
      ION( 6)=-1.98588-0.16311*XL-0.29836*XL**2
      ION( 7)=-1.99125-0.01754*XL-0.60476*XL**2
      ION( 8)=-1.98231-0.15042*XL-0.18665*XL**2
      ION( 9)=-2.59340-0.10417*XL-0.68904*XL**2
      ION(10)=-1.83156-0.00896*XL-0.82896*XL**2
      ION(11)=-1.80315-0.12404*XL-0.50833*XL**2
      ION(12)=-2.15079-0.19971*XL-0.15258*XL**2
      ION(14)=-1.82306-0.18302*XL-0.30904*XL**2
 60   CONTINUE
      IF(ABn(2).LT.0.3) GOTO 120
C     HE DOMINATED
      IF(XL.GT.-2.00000) GOTO 31
      HE=+0.76606+0.22527*XL+0.02330*XL**2
      ION( 1)= -1.98021 -0.16043*XL -0.01839*XL**2
      ION( 2)= -1.65366 -0.24070*XL -0.03785*XL**2
      ION( 3)= -1.66609 +0.03123*XL +0.00749*XL**2
      ION( 4)= -1.46773 -0.27453*XL -0.02831*XL**2
      
      ION( 6)= -1.59638 -0.24246*XL -0.03757*XL**2
      ION( 7)= -1.25043 +0.01179*XL +0.00925*XL**2
      
c     ION( 9)= -2.46719 -0.41284*XL -0.05766*XL**2
      ION(10)= -0.75351 +0.21784*XL +0.03691*XL**2
      ION(11)= -1.30055 -0.10029*XL -0.00713*XL**2
      ION(12)= -1.50910 -0.03468*XL -0.00268*XL**2
      ION(14)= -1.61882 -0.36051*XL -0.05624*XL**2
      goto 120
 31   IF(XL.GT.-1.00000) GOTO 32
      HE=+0.99317+0.44060*XL+0.07419*XL**2
      ION( 1)= -2.20200 -0.37718*XL -0.07132*XL**2
      ION( 2)= -1.31674 +0.09051*XL +0.04353*XL**2
      ION( 3)= -1.81405 +0.01724*XL +0.03748*XL**2
      ION( 4)= -1.29513 -0.15595*XL -0.01217*XL**2
      
      ION( 6)= -0.56756 +1.01425*XL +0.33358*XL**2
      ION( 7)= -1.52973 -0.24367*XL -0.04865*XL**2
      
c     ION( 9)= -2.69182 -0.37885*XL +0.01549*XL**2
      ION(10)= -0.39997 +0.96367*XL +0.32144*XL**2
      ION(11)= -2.52918 -1.96409*XL -0.63188*XL**2
      ION(12)= -1.78415 -0.34629*XL -0.08972*XL**2
      ION(14)= -2.13126 -1.14584*XL -0.32079*XL**2
      goto 120
 32   IF(XL.GT.-0.30103) GOTO 33
      HE=+1.00812+0.42881*XL+0.04744*XL**2
      ION( 1)= -2.34930 -0.65321*XL -0.20005*XL**2
      ION( 2)= -1.78821 -0.76416*XL -0.33967*XL**2
      ION( 3)= -2.14954 -0.73874*XL -0.38300*XL**2
      ION( 4)= -2.30073 -2.35916*XL -1.20978*XL**2
      ION( 6)= -1.79739 -0.76700*XL -0.21785*XL**2
      ION( 7)= -1.73076 -0.73382*XL -0.33778*XL**2
c     ION( 9)= -0.60010 +4.95778*XL +3.26041*XL**2
      ION(10)= -2.15632 -2.52493*XL -1.41081*XL**2
      ION(11)= -1.23438 +0.14190*XL +0.17932*XL**2
      ION(12)= -1.87810 -0.58727*XL -0.23676*XL**2
      ION(14)= -1.65906 -0.97742*XL -0.62458*XL**2
      goto 120
 33   IF(XL.GT.-0.04576) GOTO 120
      HE=+0.99332+0.34986*XL-0.05151*XL**2
      ION( 1)= -2.34542 -0.65243*XL -0.24031*XL**2
      ION( 2)= -1.73438 -0.46256*XL +0.06821*XL**2
      ION( 3)= -2.13966 -0.33522*XL +0.84837*XL**2
      ION( 4)= -1.55621 -0.89530*XL -4.56283*XL**2
      ION( 6)= -1.59713 -1.47803*XL -4.78974*XL**2
      ION( 7)= -1.71968 -0.91049*XL -1.04686*XL**2
c     ION( 9)= -2.60086 -4.70150*XL -6.74825*XL**2
      ION(10)= -1.44071 -1.10230*XL -4.58188*XL**2
      ION(11)= -1.46024 -0.56482*XL +0.32409*XL**2
      ION(12)= -1.91839 -0.97600*XL -1.08341*XL**2
      ION(14)= -1.56207 -1.24450*XL -2.58211*XL**2
c                         !!!!!
 120  IF(ABn(1).LT.1.e-30) GOTO 122
c     HYDROGEN DOMINATED PLASMA (SHULL AND VAN STEENBERG 85)
      IF(X.LT.1.D0) HE=0.9971*(1.-(1.-X**0.2663)**1.3163)
      IF(X.LT.1.D0) HIION=0.3908*(1.-X**0.4092)**1.7592
c                         !!      IF(X.GE.1.) HE=1.
      IF(X.GE.1.) HIION=1.E-10
      ION(18)=-LOG10(13.6/HIION)
c     CORRECT FOR HE ABUNDANCE 0.1 OF H IN THE S&S FIT
      HEION=1.E-10
      IF(X.LT.1.D0) HEION=0.554*(1.-X**.4614)**1.6660
      ION(1)=-LOG10(24.58/HEION)
 122  CONTINUE
      ION(13)=ION(3)
      ION(15)=ION(3)
c     HEATING EFFICIENCY NOT CORRECTED FOR IRON 
      call ionpot(x,axsi)
      sumh=0.
      do i=1,18
         sumh=axsi(i)/10**(-ion(i))+sumh
      enddo
      if(x.le.1.) then
         hyex=abn(1)*xion(ik,1,1)*.4766*(1.-x**.2735)**1.5221
c     hyex=axsi(18)*hyex/13.6
      else
         hyex=0.
      endif
 8278 format(' s ',1pe12.4,10e12.4)
      sumt=sumh+he+hyex
      do i=1,18
         ion(i)=ion(i)-log10(sumt)
      enddo
      he=he/sumt
      hyex=0.
      if(x.le.1.) then
         hyex=.4766*(1.-x**.2735)**1.5221/sumt
      endif
      sumh=0.
      do i=1,18
         sumh=axsi(i)/10**(-ion(i))+sumh
      enddo
      DO  I=1,20
         CION(I)=GAHE*10.**ION(I)/ELCH
      enddo
c                         !n-t
      DO  I=1,20
         if(ionq.eq.5) then
            hyex=.4766*(1.-x**.2735)**1.5221/sumt
            if(i.eq.2) then
c     n=2
               Cex(I)=0.82*GAHE*hyex/(ELCH*10.2)
            elseif(i.eq.3) then
c     n=3
               Cex(I)=0.12*GAHE*hyex/(ELCH*12.08)
            elseif(i.eq.4) then
c     n=4
               Cex(I)=0.05*GAHE*hyex/(ELCH*12.75)
            elseif(i.eq.5) then
c     n=5
               Cex(I)=0.03*GAHE*hyex/(ELCH*13.056)
            else
               cex(i)=0.
            endif
         else
            cex(i)=0.
         endif
      enddo
 6543 CONTINUE
      GAELH=HE*GAHE
      RETURN
      END
      
      
      subroutine f_nonthermal(abhi,abhei,xe,f_heat_2,
     &     f_ion_he,f_ex_he,f_ion_h,f_ex_h)         
      implicit real*8(a-h,o-z)
      if(xe.gt.0.99) then
         xel=0.99
      elseif(xe.lt.1.e-4) then
         xel=1.e-4
      else
         xel=xe
      endif
      xs=log10(xel)-log10(1.-xel)+0.00001
      a=8.4000E-01  
      b=1.4200E+00
      f_heat= 6.53639E-01 + 5.15011E-01*xs + 1.66968E-01*xs**2 +
     &     2.37903E-02*xs**3 + 1.20906E-03*xs**4
      f_heat=f_heat*(-atan(a*xs)+b)

      if(xel.lt.1.e-10) then
         f_ion_he=0.
      else
         a=8.4000E-01  
         b=1.1600E+00 
         f_ion_he= 1.02592E-02 -2.69490E-03*xs + 5.30982E-03*xs**2 +
     &        2.67128E-03*xs**3 + 3.32947E-04*xs**4
      endif
      f_ion_he=f_ion_he*(1.-xel)*(-atan(a*xs)+b)
 
      if(xel.lt.1.e-4) then
          a=3.4000E-01  
         b=2.0000E-02
         f_ex_he= 2.91160E-03 +  6.11349E-04*xs + 5.56751E-05*xs**2 
      else
         a=9.4000E-01  
         b=1.2000E+00  
         f_ex_he= 2.38400E-03 -4.69763E-04*xs + 1.34543E-03*xs**2+
     &        6.47551E-04*xs**3 + 7.93871E-05*xs**4
      endif
      f_ex_he=f_ex_he*(1.-xel)*(-atan(a*xs)+b)

      if(xel.lt.1.e-4) then
         a=3.6000E-01  
         b=2.0000E-02
         f_ion_h= 5.91124E-01+  5.79456E-02*xs + 1.84557E-03*xs**2
      else
         a=7.8000E-01  
         b=1.1200E+00  
         f_ion_h= 5.98357E-02 -1.86917E-02*xs +  3.41991E-02*xs**2 +
     &        1.50335E-02*xs**3 + 1.67892E-03*xs**4
      endif
      f_ion_h=f_ion_h*(1.-xel)*(-atan(a*xs)+b)

      if(xel.lt.1.e-4) then
         a=4.0000E-01  
         b=2.0000E-02
         f_ex_h= 5.75337E-01+  2.47695E-02*xs -1.26653E-03*xs**2
      else
         a=8.4000E-01  
         b=1.1600E+00  
         f_ex_h= 6.63175E-02 -1.80824E-02*xs + 3.95371E-02*xs**2 +
     &        1.64258E-02*xs**3 +  1.80579E-03*xs**4
      endif
      f_ex_h=f_ex_h*(1.-xel)*(-atan(a*xs)+b)

      f_heat_2=1.-f_ion_he-f_ex_he-f_ion_h-f_ex_h 
      eff_ion_h=0.9*(1.-xel)*13.6/f_ion_h
      eff_ex_h=0.9*(1.-xel)*10.2/f_ex_h
      eff_ion_he=0.1*(1.-xel)*24.5/f_ion_he
      eff_ex_he=0.1*(1.-xel)*19.2/f_ex_he
      f_ion_h=abhi*13.6/eff_ion_h
      f_ex_h=abhi*10.2/eff_ex_h
      f_ion_he=abhei*24.5/eff_ion_he
      f_ex_he=abhei*19.2/eff_ex_he
      f_heat_2=1.-f_ion_he-f_ex_he-f_ion_h-f_ex_h 
      return
      end





      SUBROUTINE EXSEC_old(X,CEX,CION)
C
C     1 = HEI, 2 = O I,  3 = O II,   4 = CA I, 5 = CA II, 6 = MG I
C     7 = I, 8 = NA I, 9 = MG II, 10 = S I, 11 = SI I  12 = NE I
C    13 = NE II 14 = AR I 15 = AR II 16 = FE I 17 = FE II18 = H I
C     ALL RATES SHOULD BE DIVIDED BY DENS**2
C
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 ION
c      parameter(nl=340)
      COMMON/TPAR/RIN,DRQ,R,TDAY
      COMMON/GAMPAR/GAMLUM,GAELH,GAHE
      COMMON/HYDROS/HSCALE,DEN1,XEL,AMEAN
      common/abl/abn(15)
      DIMENSION CEX(20),CION(20),ION(20)
      dimension axsi(20)
      DATA PI/3.1415926E0/,ELCH/1.60219E-12/,AMU/1.660531E-24/
      DO K=1,20
         CION(K)=1.e-30
         ion(k)=-100.
      ENDDO
      gaelh=0.
c      gahe=0.
      he=1.
      IF(X.GT.1.) GOTO 6543
      XL=LOG10(X)
c     HYDROGEN DOMINATED
      ESEC=1.E3
c      CALL SEC(ESEC,X,HEFF,SEL)
c      IF(ABn(2).GT.0.3) GOTO 60
c     OXYGEN DOMINATED
      ION(16)=-1.523
      ION(17)=-1.523
      IF(XL.GT.-1.00000) GOTO 50
      HE= 0.85012+0.27247*XL+0.02703*XL**2
      ION( 1)=-2.44824-0.09560*XL-0.01888*XL**2
      ION( 2)=-2.05749-0.24897*XL-0.04278*XL**2
      ION( 3)=-2.46878-0.15229*XL-0.03719*XL**2
      ION( 4)=-1.66678-0.13942*XL-0.01920*XL**2
      ION( 5)=-2.25347-0.07530*XL-0.01100*XL**2
      ION( 6)=-1.83600-0.13622*XL-0.02602*XL**2
      ION( 7)=-1.80086-0.10583*XL-0.01849*XL**2
      ION( 8)=-1.80164-0.09341*XL-0.00618*XL**2
      ION( 9)=-2.36815-0.03717*XL-0.00207*XL**2
      ION(10)=-1.62796-0.06683*XL-0.00966*XL**2
      ION(11)=-1.62177-0.10160*XL-0.01563*XL**2
      ION(12)=-1.93059-0.03954*XL-0.00591*XL**2
      ION(14)=-1.56809-0.01799*XL+0.00091*XL**2
      GOTO  60
50      IF(XL.GT.-0.30103) GOTO 55
      HE= 0.82860+0.21277*XL-0.01115*XL**2
      ION( 1)=-2.57233-0.21134*XL-0.01054*XL**2
      ION( 2)=-2.18307-0.43588*XL-0.10411*XL**2
      ION( 3)=-2.53715-0.26787*XL-0.08441*XL**2
      ION( 4)=-1.74414-0.25944*XL-0.06185*XL**2
      ION( 5)=-2.33686-0.03342*XL+0.11426*XL**2
      ION( 6)=-1.95195-0.34338*XL-0.11724*XL**2
      ION( 7)=-1.99154-0.51203*XL-0.23401*XL**2
      ION( 8)=-1.96099-0.41032*XL-0.16374*XL**2
      ION( 9)=-2.52601-0.28823*XL-0.09526*XL**2
      ION(10)=-1.82242-0.41406*XL-0.16242*XL**2
      ION(11)=-1.75820-0.30483*XL-0.08243*XL**2
      ION(12)=-2.08786-0.17835*XL+0.01254*XL**2
      ION(14)=-1.73749-0.13983*XL+0.04846*XL**2
      GOTO  60
55      IF(XL.GT. 0.00004) GOTO 56
      HE= 0.83097+0.19072*XL-0.11060*XL**2
      ION( 1)=-2.62139-0.55947*XL-0.62560*XL**2
      ION( 2)=-2.23986-0.77837*XL-0.61517*XL**2
      ION( 3)=-2.58118-0.55265*XL-0.54455*XL**2
      ION( 4)=-1.82371-0.62313*XL-0.39196*XL**2
      ION( 5)=-2.45797-0.50263*XL-0.10790*XL**2
      ION( 6)=-1.98587-0.47955*XL-0.19531*XL**2
      ION( 7)=-1.99122-0.71637*XL-0.91636*XL**2
      ION( 8)=-1.98230-0.54215*XL-0.36654*XL**2
      ION( 9)=-2.59339-0.28247*XL+0.66738*XL**2
      ION(10)=-1.83154-0.49657*XL-0.33597*XL**2
      ION(11)=-1.80315-0.28671*XL+0.47372*XL**2
      ION(12)=-2.15078-0.52567*XL-0.44686*XL**2
      ION(14)=-1.82305-0.46815*XL-0.09795*XL**2
      GOTO  60
56    HE= 0.83097+0.35645*XL+0.19409*XL**2
      ION( 1)=-2.62141-0.08512*XL-0.39237*XL**2
      ION( 2)=-2.23989-0.01840*XL-0.47092*XL**2
      ION( 3)=-2.58120-0.09856*XL-1.10769*XL**2
      ION( 4)=-1.82374+0.00649*XL-0.86178*XL**2
      ION( 5)=-2.45799-0.11213*XL-0.27513*XL**2
      ION( 6)=-1.98588-0.16311*XL-0.29836*XL**2
      ION( 7)=-1.99125-0.01754*XL-0.60476*XL**2
      ION( 8)=-1.98231-0.15042*XL-0.18665*XL**2
      ION( 9)=-2.59340-0.10417*XL-0.68904*XL**2
      ION(10)=-1.83156-0.00896*XL-0.82896*XL**2
      ION(11)=-1.80315-0.12404*XL-0.50833*XL**2
      ION(12)=-2.15079-0.19971*XL-0.15258*XL**2
      ION(14)=-1.82306-0.18302*XL-0.30904*XL**2
60    CONTINUE
      IF(ABn(2).LT.0.3) GOTO 120
C     HE DOMINATED
      IF(XL.GT.-2.00000) GOTO 31
      HE=+0.76606+0.22527*XL+0.02330*XL**2
      ION( 1)= -1.98021 -0.16043*XL -0.01839*XL**2
      ION( 2)= -1.65366 -0.24070*XL -0.03785*XL**2
      ION( 3)= -1.66609 +0.03123*XL +0.00749*XL**2
      ION( 4)= -1.46773 -0.27453*XL -0.02831*XL**2

      ION( 6)= -1.59638 -0.24246*XL -0.03757*XL**2
      ION( 7)= -1.25043 +0.01179*XL +0.00925*XL**2

c      ION( 9)= -2.46719 -0.41284*XL -0.05766*XL**2
      ION(10)= -0.75351 +0.21784*XL +0.03691*XL**2
      ION(11)= -1.30055 -0.10029*XL -0.00713*XL**2
      ION(12)= -1.50910 -0.03468*XL -0.00268*XL**2
      ION(14)= -1.61882 -0.36051*XL -0.05624*XL**2
      goto 120
31    IF(XL.GT.-1.00000) GOTO 32
      HE=+0.99317+0.44060*XL+0.07419*XL**2
      ION( 1)= -2.20200 -0.37718*XL -0.07132*XL**2
      ION( 2)= -1.31674 +0.09051*XL +0.04353*XL**2
      ION( 3)= -1.81405 +0.01724*XL +0.03748*XL**2
      ION( 4)= -1.29513 -0.15595*XL -0.01217*XL**2

      ION( 6)= -0.56756 +1.01425*XL +0.33358*XL**2
      ION( 7)= -1.52973 -0.24367*XL -0.04865*XL**2

c      ION( 9)= -2.69182 -0.37885*XL +0.01549*XL**2
      ION(10)= -0.39997 +0.96367*XL +0.32144*XL**2
      ION(11)= -2.52918 -1.96409*XL -0.63188*XL**2
      ION(12)= -1.78415 -0.34629*XL -0.08972*XL**2
      ION(14)= -2.13126 -1.14584*XL -0.32079*XL**2
      goto 120
32    IF(XL.GT.-0.30103) GOTO 33
      HE=+1.00812+0.42881*XL+0.04744*XL**2
      ION( 1)= -2.34930 -0.65321*XL -0.20005*XL**2
      ION( 2)= -1.78821 -0.76416*XL -0.33967*XL**2
      ION( 3)= -2.14954 -0.73874*XL -0.38300*XL**2
      ION( 4)= -2.30073 -2.35916*XL -1.20978*XL**2
      ION( 6)= -1.79739 -0.76700*XL -0.21785*XL**2
      ION( 7)= -1.73076 -0.73382*XL -0.33778*XL**2
c      ION( 9)= -0.60010 +4.95778*XL +3.26041*XL**2
      ION(10)= -2.15632 -2.52493*XL -1.41081*XL**2
      ION(11)= -1.23438 +0.14190*XL +0.17932*XL**2
      ION(12)= -1.87810 -0.58727*XL -0.23676*XL**2
      ION(14)= -1.65906 -0.97742*XL -0.62458*XL**2
      goto 120
33    IF(XL.GT.-0.04576) GOTO 120
      HE=+0.99332+0.34986*XL-0.05151*XL**2
      ION( 1)= -2.34542 -0.65243*XL -0.24031*XL**2
      ION( 2)= -1.73438 -0.46256*XL +0.06821*XL**2
      ION( 3)= -2.13966 -0.33522*XL +0.84837*XL**2
      ION( 4)= -1.55621 -0.89530*XL -4.56283*XL**2
      ION( 6)= -1.59713 -1.47803*XL -4.78974*XL**2
      ION( 7)= -1.71968 -0.91049*XL -1.04686*XL**2
c      ION( 9)= -2.60086 -4.70150*XL -6.74825*XL**2
      ION(10)= -1.44071 -1.10230*XL -4.58188*XL**2
      ION(11)= -1.46024 -0.56482*XL +0.32409*XL**2
      ION(12)= -1.91839 -0.97600*XL -1.08341*XL**2
      ION(14)= -1.56207 -1.24450*XL -2.58211*XL**2
c!!!!!
120   IF(ABn(1).LT.1.e-30) GOTO 122
C     HYDROGEN DOMINATED PLASMA (SHULL AND VAN STEENBERG 85)
c!!      IF(X.LT.1.D0) HE=0.9971*(1.-(1.-X**0.2663)**1.3163)
      IF(X.LT.1.D0) HIION=0.3908*(1.-X**0.4092)**1.7592
c!!      IF(X.GE.1.) HE=1.
      IF(X.GE.1.) HIION=1.E-10
      ION(18)=-LOG10(13.6/HIION)
C     CORRECT FOR HE ABUNDANCE 0.1 OF H IN THE S&S FIT
      HEION=1.E-10
      IF(X.LT.1.D0) HEION=0.0554*(1.-X**.4614)**1.6660+1.
c!!      ION(1)=-LOG10(24.58/HEION)
122   CONTINUE
      ION(13)=ION(3)
      ION(15)=ION(3)
C     HEATING EFFICIENCY NOT CORRECTED FOR IRON 
      call ionpot(x,axsi)
      sumh=0.
      do i=1,18
         sumh=axsi(i)/10**(-ion(i))+sumh
      enddo
      hyex=.4766*(1.-x**.2735)**1.5221
c      hyex=axsi(18)*hyex/13.6
      sumt=sumh+he
      do i=1,18
         ion(i)=ion(i)-log10(sumt)
      enddo
      he=he/sumt
      sumh=0.
      do i=1,18
         sumh=axsi(i)/10**(-ion(i))+sumh
      enddo
      DO I=1,20
         CION(I)=GAHE*10.**ION(I)/ELCH
      enddo
6543  CONTINUE
      GAELH=HE*GAHE
      RETURN
      END

      SUBROUTINE ionpot(x,axsi)
      IMPLICIT REAL*8(A-H,O-Z)
      include 'param'
      common/ind/ik
      common/ionx/xion(md,14,27)
      common/abl/abn(15)
      dimension axsi(20)
      axsi(1)=abn(2)*xion(ik,2,1)*24.5876
      axsi(2)=abn(5)*xion(ik,5,1)*13.614
      axsi(3)=abn(5)*xion(ik,5,2)*35.108
      axsi(4)=abn(13)*xion(ik,13,1)*6.09
      axsi(5)=abn(13)*xion(ik,13,2)*11.8
      axsi(6)=abn(8)*xion(ik,8,1)*7.6
      axsi(7)=abn(3)*xion(ik,3,1)*11.264
      axsi(8)=abn(7)*xion(ik,7,1)*5.12
      axsi(9)=abn(8)*xion(ik,8,2)*15.0
      axsi(10)=abn(11)*xion(ik,11,1)*10.3
      axsi(11)=abn(10)*xion(ik,10,1)*8.1
      axsi(12)=abn(6)*xion(ik,6,1)*21.559
      axsi(13)=abn(6)*xion(ik,6,2)*41.07
      axsi(14)=abn(12)*xion(ik,12,1)*15.7
      axsi(15)=abn(12)*xion(ik,12,2)*27.8
      axsi(16)=abn(14)*xion(ik,14,1)*7.8
      axsi(17)=abn(14)*xion(ik,14,2)*16.2
      axsi(18)=abn(1)*xion(ik,1,1)*13.595
      axsi(19)=0.
      axsi(20)=0.
      return
      end


      SUBROUTINE ESCAPE(T0,T0T,WL,A21,VTERM,BE,DBEDTA)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/QSOM/qso,IAGN,ISTAT,ISOBOL
C     EVALUATE ESCAPE PROBABILITY IN BOTH DIRECTIONS
      CALL ESCAP(T0,WL,A21,VTERM,BE1,DBEDTA1)
      BE2=0.
      DBEDTA2=0.
c!!!!      IF(ISTAT.EQ.1) CALL ESCAP(T0T,WL,A21,VTERM,BE2,DBEDTA2)
      BE=BE1+BE2
      DBEDTA=DBEDTA1+DBEDTA2
      RETURN
      END

      SUBROUTINE ESCAP(T00,WL,A21,VTERM,BE,DBEDTA)
      IMPLICIT REAL*8(A-H,O-Z)
C 
C     ESCAPE PROBABILTY
C
      COMMON/QSOM/qso,IAGN,ISTAT,ISOBOL
      COMMON/TPAR/RIN,DRQ,R,TDAY
      common/deltar/drscgs
      COMMON/NION/IONq
      IF(ISTAT.EQ.1) THEN
         t0=t00
         goto 234
C     ESCAPE PROB. FROM FERLAND
c!!
         IF(VTERM.LE.0.) VTERM=1.285E6        
         A=7.96E-10*WL*A21/VTERM
         AT=A*T0
         IF(AT.GT.1.) THEN
            B=1.6 + 3./((2.*A)**.12*(1.+AT))
            DBDT = -3.*A/((2.*A)**.12 * (1. + AT)**2)
         ELSEIF(AT.GT.0.) THEN
            B=1.6 + 3./((2.*A)**.12*(1.+1./SQRT(AT)))
            DBDT = 3./(2.*SQRT(AT)*T0*(2.*A)**.12*(1.+1./SQRT(AT))**2)
         ELSE
            B=1.6
            DBDT=0.
         ENDIF
         IF(B.GT.5) THEN
            B=5.
            DBDT=0.
         ENDIF
         BE = 1./(1. + B * T0)
         DBEDTA = - BE**2*(DBDT*T0 + B)
c     The escape probability used for partial redistr.
c     used in cloudy
 234     tau=t00
         a = 0.6451+TAU 
         b1 = 1.+7.3E-6*TAU
         b = 0.47+1.08/b1
c     be = 1. / ( (0.6451+TAU) * (0.47+1.08/(1.+7.3E-6*TAU) ) )
         be = 1. / ( a * b )
C!!     ONLY ESCAPE FROM ONE SIDE OF THE SLAB
         BE=BE/2.
c     dbedta=tau*be*(-1.+0.6451*1.08*7.3e-6*be/(1.+7.3e-6*tau))
         dbedta=be*(-1./a +7.884e-6/(b*b1**2))
         DBEDTA=DBEDTA/2.
         GOTO 123
C     ESCAPE PROB. FROM FERLAND AND NETZER
         IF(T0.LT.1) THEN
            Q=1.11*T0**0.826
            BE=1./(1.+Q)
            DQDT=0.
            IF(T0.GT.0.) DQDT=0.826*Q/T0
            DBEDTA=-DQDT/(1.+Q)**2
         ELSE
            Q=1.11*T0**1.071/(1.+(LOG10(T0)/5.)**5)
            BE=1./(1.+Q)
            DQDT=1.071*Q/T0-6.26E-4*Q**2*(LOG10(T0))**4/T0**2.071
            DBEDTA=-DQDT/(1.+Q)**2
         ENDIF
C     ESCAPE PROB. FROM KWAN AND KROLIK
         GOTO 123
         IF(T0.LT.1.E-5) THEN
            BE=1.-T0
            DBEDTA=-1.
         ELSEIF(T0.LE.1.) THEN
            BE=(1.-EXP(-2.*T0))/(2.*T0)
            DBEDTA=-BE/T0+EXP(-2.*T0)/T0
         ELSE
C     DAMPING CONSTANT
C     A VALID ONLY FOR H AND T=1.E4K. SCALES AS 1./SQRT(T)
c!!   wl in A? see also calling subroutines!
            A=7.96E-10*WL*A21/VTERM
            SQPI=1.77245
            TC=(SQPI*(-LOG(A)+3.102))/A
            BE=1./(SQPI*T0*(1.2+0.5*SQRT(LOG(T0))/(1.+T0/TC)))
            DBEDTA=-BE*(1./T0+0.5*SQPI*SQRT(LOG(T0))*BE*
     &           (1./(2.*LOG(T0))-1./(1.+T0/TC))/(1.+T0/TC))
         ENDIF
C!!   ONLY ESCAPE FROM ONE SIDE OF THE SLAB
         BE=BE/2.
         DBEDTA=DBEDTA/2.
      ELSE
C     SOBOLEV FOR V = R
         drsobolev=vterm*tday*8.64e4
         if(drsobolev.gt.abs(drscgs)) then
            t0=abs(drscgs)*t00/drsobolev
         else
            t0=t00
         endif

c!!!! skip non-sobolev here. do it when calc. tau

         t0 = t00

         IF(T0.LT.1.E-5) THEN
            BE=1.-T0/2.+T0**2/6.
            DBEDTA=-0.5
         ELSE
            IF(T0.LT.100.) THEN
               BE=(1.-EXP(-T0))/T0
               DBEDTA=((T0+1.)*EXP(-T0)-1.)/T0**2
            ELSE
               BE=1./T0
               DBEDTA=-1./T0**2
            ENDIF
         ENDIF
      ENDIF
123   CONTINUE
      RETURN
      END


      SUBROUTINE MULTISIMPQ(iel,ion,Z,TE,XI,rltot)
c
c Obs! Keep total number of 'neutral' and 'ionized' atoms constant
c
c     IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT none
      INTEGER N
c      PARAMETER (MD=350,MDP1=MD+1)
      include 'param'
      integer nl,nlp1,nfel,nlp,ik,IAGN,ISTAT,ISOBOL,ionq,iel,ion,idep
      integer i,ierr,ilabcfx,im1,imiss,j,k,ioni,ipr
      integer NION,NLEV,NP1H,NMAX,NMIN,ipre,kmax,itime,ideb,kch
      real*8 C3,C33,dens,del,RIN,DRQ,R1Q,TDAY,EM,ESC,TTOT,TOP
      real*8 FRH,WEQ,SR,WOBS,e00,CS,CI,G,E
      real*8 WLN,DCDT,DCIDT,RECNET,RECNDT,RECCO,RTE
      real*8  PH,DRDT,RECT,PHET,RECO,XN,XN1,XN2,XN3
      real*8  xion,abn,time,dr_dv,c,be,xi,aa,x,xo,betot,tau
      real*8  pi,elch,amu,wlin,z,te,rltot,a,cinx,den,denel,qso
      real*8 taul,taulinex,tev,totopl,wl,wlix,x0,xplus,wtot
      real*8 au,besc,betadest,dbedta,dei,di,eps,err,errmax,fbeta,simul
      real*8 fhumm,opcon,opline,opline2,pd,refn,s,sig,t0,t0t,vterm,wlair
      integer ikk,iprint,l,ll,neq,nlhs
      PARAMETER (NL=340,NLP1=NL+1)
      PARAMETER (NFEL=3000)
      COMMON/NION/IONQ
      parameter (nlp=30000)
c      COMMON/COLD/COLTOT(14,27),COLTOTT(14,27),SRED(NFEL+nlp)
      COMMON/QSOM/qso,IAGN,ISTAT,ISOBOL
      COMMON/A7/C3,C33
      COMMON/PHY/DENS(MD)
      COMMON/ELEC/DEL(MD)
      COMMON/TPAR/RIN,DRQ,R1Q,TDAY
      COMMON/IND/IK
      COMMON/A19/EM(NL,NL),ESC(NL,NL),TTOT(NL,NL),TOP(NL,NL)
      COMMON/EQUIH/FRH(2,401),WEQ(401),SR(NFEL+nlp),WOBS(NL,NL)
      COMMON/NLEV/e00,NION,NLEV,NP1H,NMAX,NMIN
      COMMON/A14/CS(NL,NL),CI(NL),G(NL),E(NL),A(NL,NL),
     &     WLN(NL,NL),DCDT(NL,NL),DCIDT(NL)
      COMMON/A17/RECNET(NL),RECNDT(NL),RECCO(NL),RTE(NL),
     &     PH(NL),DRDT(NL),RECT(6,NL),PHET(NL)
      COMMON/RECAL/RECO(NL)
      common/w/wlin(401),wl(NL,NL),taul(401)
      common/line_em/cinx(14,26,401),taulinex(14,26,401),
     &     wlix(14,26,401),ilabcfx(14,26,401)
      common/lineopac/totopl(nlp)
      COMMON/HPOP/XN(6,NL),XN1,XN2,XN3
      common/ionx/xion(md,14,27)
      common/abl/abn(15)
      common/preion/ipre
      common/kmaxpop/kmax(14,27)
      common/timecheck/time,itime
      common/velgrad/dr_dv      
      common/debug/ideb
      DIMENSION BE(NL,NL),XI(NL),C(NL,NL)
      DIMENSION AA(nlp1,nlp1),X(nl),XO(NL),betot(nl,nl)
c      real*8 AA(nlp1,nlp1),X(nl),XO(NL),betot(nl,nl)
      dimension kch(nl,nl),tau(nl,nl)
      DATA PI/3.1415926E0/,ELCH/1.60219E-12/,AMU/1.660531E-24/
      ipr=0
      ideb=0
      if(iel.eq.4.and.ion.eq.2) then
         ipr=0
      elseif(iel.eq.5.and.ion.eq.2) then
         ipr=0
      endif
      if(iel==12.and.ion==5) then
         ipr=0
      endif
      if(iel==14.and.ion==2) then
         ipr=0
      endif
      if(ipr==1)  then

         do i=1,5
            do j=1,5
               WLin=DABS(12398.54d0/(E(I)-E(J)))
            enddo
         enddo
         
      endif

      DEN=DENS(IK)
      DENEL=DEN*DEL(IK)
      TEV=TE/1.1609E4      
      x0 = xion(ik,iel,ion)
      xplus = xion(ik,iel,ion+1)
      z=abn(iel)


      if(iel.eq.-14.and.ion==1) then
         write(6,*)istat,time
         write(6,9275)iel,ion,ik,z,xplus,
     &        x0,z*x0,den,del(ik),te
         write(6,2929)(xion(ik,iel,ioni),ioni=1,12)

 2929    format('Fe ion ',1pe12.3,20e12.3)
         do i=1,12
            write(6,9277)i,ph(i),reco(i)
 9277       format(i5,1pe12.3,10e12.3)
         enddo

 9275    format(' xplus coltot ',3i5,1pe12.3,10e12.3)
      endif
      if(ik.gt.2.and.(iel.eq.1.or.iel.eq.2)) then
         do i=1,nlev
            x(i)=xn(ionq,i)
         enddo
      endif
      do i=1,nl
         x(i)=1.e-10
      enddo
      x(1)=0.999
      do i=1,401
         weq(i)=0.
         sr(i)=0.
         wlin(i)=0.
         taul(i)=0.
      enddo

      DO I=1,NLP1
         DO J=1,NLP1
            aa(i,j)=0.d0
         ENDDO
      ENDDO
      DO I=1,NL
         DO J=1,NL
            c(I,J)=0.
c            be(I,J)=0.
c            betot(I,J)=0.
         ENDDO
      ENDDO
C     ADD ONE LEVEL FOR CONTINUUM
      N=NLEV+1
      E(N)=E00
      DO I=1,N
c!!!  Obs! skip all recomb an phoionizations since the ion equil is already determined
c     through x0 and xplus. This means no recombination cascade either! Questionable!!
         reco(i)=0.
         ph(i)=0.
         DO J=1,N            
            WL(I,J)=0.0
            IF(I.NE.J.AND.E(I).NE.E(J)) THEN
               WL(I,J)=DABS(12398.54d0/(E(I)-E(J)))
            ENDIF
         ENDDO
      ENDDO
      TEV=TE/1.1609D4
      TIME=8.64E4*TDAY
      x(1)=x0
      DO I=1,NLEV
        DO J=1,NLEV
          C(J,I)=DENEL*CS(J,I)
c          if(iel.eq.14.and.ion.le.14.and.ion.ge.10) then
c             if(j.le.2.and.i.le.2) then
c                c(j,i)=c(j,i)
c             else
c                c(j,i)=1.e-20*c(j,i)
c             endif
c          endif
        ENDDO
      ENDDO
      DO L=1,20
         DO I=1,Nlev
            DEI=X(I)*Z*DEN
            DO J=1,N
               AA(I,J)=0.
c               IF(I.EQ.1) THEN
c replace first level eqn. by number cons.
c                  AA(I,J)=1.D0
c               ELSE
                  BE(I,I)=0.
                  IF(I.NE.J) THEN
C     
C     CALCULATE THE ESCAPE PROBABILITY.
C
                     IF(ISTAT.EQ.1) THEN
c                        T0=C33*COLTOT(iel,ion)*X(I)*A(J,I)*G(J)*     
c     &                       WL(J,I)**3./G(I)
                        T0T=T0

                     ELSE
                        T0=1.e-24*WL(J,I)**3*A(J,I)*G(J)*DEI*TIME/
     &                       (8.*PI*G(I))

c vel grad

c                        T0=1.e-24*WL(J,I)**3*A(J,I)*G(J)*DEI*dr_dv/
c!!     &                       (8.*PI*G(I))
                     ENDIF
                     IF(L.EQ.1) T0=.0
C!!   ASSUME MEAN ATOMIC WEIGHT = 4. 
                     AU=4.
                     if(iel.eq.1) then
                        AU=1.
                     elseif(iel.eq.2) then
                        au=4.
                     elseif(iel.eq.3) then
                        au=12.
                     elseif(iel.eq.4) then
                        au=14.
                     elseif(iel.eq.5) then
                        au=16.
                     elseif(iel.eq.13) then
                        au=56.
                     elseif(iel.eq.12) then
                        au=40.
                     elseif(iel.eq.14) then
                        au=56.                        
                     endif
                     VTERM=1.285E6*SQRT(TE/(AU*1.E4))
                     CALL ESCAPE(T0,T0T,WL(J,I),A(J,I),VTERM,
     &                    BE(J,I),DBEDTA)
                     if(iel.eq.-14.and.ion.eq.1) 
     &                  then 
                        write(6,9231)j,i,istat,wl(j,i),
     &                       time,t0,
     &                       vterm,t0t,be(j,i)
                     endif

c                     ESC(J,I)=BE(J,I)
                     
c!!                     tau(j,i) = t0

                     if(j.gt.i) then 
                        be(i,j)=be(j,i)
                     endif
                     if(iel.eq.-14.and.ion==1) then
c                     if(iel.eq.1) then
                        write(6,9231)j,i,istat,wl(j,i),g(j),dei,time,t0,
     &                       a(j,i),be(j,i)
 9231                   format(' esc ',3i5,f10.2,1pe12.3,10e12.3)
                     endif
                     besc=be(j,i)
                     pd=0.
                     ll=0
                     if(ionq.eq.3.or.ionq.eq.5) then
c     first 584 A line
                        if(ionq.eq.3) then
                           if(i.eq.5.and.j.eq.1) then
                              ll=47
                           endif
                        endif
c     now Lyman alpha
                        if(ionq.eq.5) then
                           if(i.eq.2.and.j.eq.1) then
                              ll=63
                           endif
                        endif
                        if(ll.ne.0) then
c     line center opacity
c     line center opacity
                           opline=3.979E-26*wl(i,j)**3.*g(i)*a(i,j)*
     &                          dei/(g(j)*vterm) 
                           opline2=t0/(vterm*time*den)
c     continuum desruction probability
                           opcon=totopl(ll)
c     H-R expression
                           fhumm=8.5
                           if(opline.gt.0.) then
                              betadest=opcon/opline
                           else
                              betadest=1.d33
                           endif
                           fbeta=min(fhumm*betadest,1.d-3)
                           pd=fbeta/(fbeta+1.)
                        endif
                     endif
                     betot(j,i)=(1.-pd)*besc+pd
                     if(pd.ne.0..and.iel.eq.1) then
                        write(6,9371)ll,i,j,ionq,t0,opcon,opline,fbeta,
     &                       pd,besc,betot(j,i)
 9371                   format(' pdest ',4i5,1pe12.3,10e12.3)
                     endif
                  ENDIF
                  IF(I.EQ.J) THEN
                     S=0.
                     DO K=1,N
                        S=S+C(I,K)+BEtot(I,K)*A(I,K)
                        if(iel.eq.-14.and.ion==1) then
c                        if(iel.eq.1) then
                           write(6,928)i,k,c(i,k),betot(i,k),a(i,k)
                        endif
 928                    format(' s ',2i5,1pe12.4,10e12.4)
                     ENDDO
C     ADD PHOTOIONIZATION
                     AA(I,I)=-S-PH(I)
                  ELSEIF(J.LT.N) THEN
                     AA(I,J)=BEtot(J,I)*A(J,I)+C(J,I)
                  ELSEIF(J.EQ.N) THEN
C     RECOMBINATION CONTRIBUTION
                     AA(I,N)=RECO(I)*DENEL*xplus
                  ENDIF     
            ENDDO
         ENDDO



c replace i=1 eqn with NUMBER CONSERVATION
c replace first level eqn. by number cons.


         DO J=1,Nlev
            AA(1,J)=1.D0
         enddo

         AA(1,N)=x0

         EPS=1.D-30
c         DO I=1,nlev
         DO I=1,n
            XO(I)=X(I)
         ENDDO
c!!

         neq=n-1
         nlhs=n
         if(iel.eq.-14.and.ion.eq.1) then
c         if(iel.eq.1) then
            write(6,9222)(x(j),j=1,nlhs+1)
            do i=1,neq+1
               write(6,*)' i ',i
               write(6,9222)(aa(i,j),j=1,nlhs+1)
c               write(6,9222)(x(j)*aa(i,j),j=1,nlhs+1)
            enddo
 9222       format('aa ',1pe12.3,9e12.3)
         endif
         do i=1,n
            x(i)=0.
         enddo
c HIT
         
         DI=SIMUL(Neq,AA,X,EPS,1,Nlhs)
         ERRMAX=0.
         x(n)=xplus
         DO I=1,N
            IF(X(I).NE.0.) THEN
               ERR=ABS((XO(I)-X(I))/X(I))
            ELSE
               ERR=1.
            ENDIF
            IF(ERR.GT.ERRMAX) ERRMAX=ERR
         ENDDO

         IF(ERRMAX.LT.0.01.and.l.ge.2) GOTO 555

         if(iel.eq.-12.and.ion.eq.5) then
c         if(iel.eq.-1) then
            do i=1,n
               write(6,9292)iel,ion,i,x(i),te,ph(i),reco(i),denel
               do j=1,i-1
                  write(6,9292)iel,ion,j,x(j),be(i,j),c(i,j),c(j,i)
 9292             format(' multi ',3i5,1pe12.3,10e12.3)
               enddo
            enddo
         endif


      ENDDO
 555  CONTINUE
      DO I=1,N
         XI(I)=X(I)
      ENDDO
      K=0
c      write(6,*)' multismpq ',iel,nlev,nl,te,xplus,denel,x(1),x(1)*z


      if(iel.eq.1.and.ion.ge.1.and.x(1).lt.0.) then
         write(6,*)' x(H) < 0'
         write(6,*)istat,time
         write(6,9275)iel,ion,ik,z,xplus,
     &     x0,z*x0,den,del(ik),te
         do i=1,nlev
            write(6,9177)i,x(i),ph(i),reco(i)
 9177       format(i5,1pe12.3,10e12.3)
         enddo
      endif
      iprint=0
      if(iel.eq.4.and.ion.eq.2) then
         iprint=0
      elseif(iel.eq.4.and.(ion.eq.3.or.ion.eq.4)) then
         iprint=0
      elseif(iel.eq.3.and.ion.eq.3) then
         iprint=0
      elseif(iel.eq.5.and.ion.eq.7) then
         iprint=0
      elseif(iel.eq.8.and.ion.eq.2) then
         iprint=0
c         write(6,*)'Mg II x0,z,x(1) ',x0,z,x(1)
      elseif(iel.eq.12.and.ion.eq.5) then
         iprint=0
c         write(6,*)'Ar V x0,z,x(1) ',x0,z,x(1) 
      elseif(iel.eq.14.and.(ion.ge.9.and.ion.le.17)) then
c         write(6,*)' cooling of Fe ',ion,z,x0,x(1)
         iprint=0
      endif

      
      rltot=0.
      wtot=0.
      ierr=0
      DO I=2,NLEV
         IM1=I-1
         DO J=1,IM1
            kch(i,j)=0
            WOBS(I,J)=1.602E-12*Z*(E(I)-E(J))*X(I)*A(I,J)*BE(I,J)/DEN
            EM(I,J)=1.602E-12*Z*(E(I)-E(J))*
     &           (X(J)*C(J,I)-X(I)*C(I,J))/DENEL
            rltot=rltot+em(i,j)
            wtot=wtot+wobs(i,j)
c     if(iprint.eq.1.and.j.le.10.and.i.le.10) then
            if(iel==-14.and.ion==2.and.i<=16) then
               write(6,9228)i,j,k,wl(i,j),x(i),x(j)
     &            ,a(i,j),be(i,j),c(i,j),c(j,i),em(i,j),rltot,wobs(i,j)
            endif

            IF(K.LE.400.and.a(i,j).gt.0.) THEN
               if(iel.eq.1.and.(j.eq.2.or.j.eq.3)) then
                  ikk=1
               elseif(iel.eq.1.and.i.le.20) then
                  ikk=1
               elseif(iel.ne.1) then
                  ikk=1
               else
                  ikk=0
               endif
               if(ikk.eq.1) then
                  K=K+1
                  SR(K)=0.
                  WEQ(K)=WOBS(I,J)
                  wlin(k)=wl(i,j)
c                  taul(k)=tau(i,j)
                  sig=1.e4/wlin(k)
                  refn=1.+1.e-8*(6431.8+2949330./(146.-sig**2)+
     &                 25536./(41.-sig**2))
                  wlair=wl(i,j)/refn
                  if(wl(i,j).gt.2000.) then
                     wlin(k)=wlair
                  endif
                  wlix(iel,ion,k)=wlin(k)
c                  if(wlin(k).gt.3332422680.) then
c                     write(6,*)' zzzz ',iel,ion,i,j,k,wlin(k)
c                  endif
                  kch(i,j)=1
c     if(iel.eq.1.and.ion.eq.1.and.i.le.5) then
                  if(weq(k).lt.-1.d-30) then
                     ierr=1
                  endif

                  if(weq(k).lt.-1.d-30.or.weq(k).gt.1.d-10.or.
     &                 (iel.eq.-11.and.ion.ge.4)) then
                     write(6,9229)k,i,j,wl(i,j),x(i),e(i),e(j),te,
     &                    xplus,a(i,j),be(i,j),c(i,j),em(i,j),del(ik),weq(k)
                  endif
 9228             format('weq ',3i5,f15.2,1pe12.3,10e12.3)
 9229             format('weq<0 ',3i5,f12.2,1pe12.3,10e12.3)
               endif
            endif
         ENDDO
      ENDDO

c      write(6,*)' rltot,wtot ',iel,nlev,nl,te,rltot,wtot

      imiss = 0
      DO I=2,NLEV
         IM1=I-1
         DO J=1,IM1

            if(iel.eq.-2.and.ion.eq.1.and.i.le.16) then
c               write(6,9228)i,j,k,wl(i,j),x(i)
c     &              ,be(i,j)
c     &         ,wtot
c$$$     &              wobs(i,j),wobs(i,j)/wobs(8,4)
            endif

            if(wobs(i,j).gt.0.01*wtot.and.kch(i,j).eq.0) then
               imiss = 1
c               write(6,9264)iel,ion,i,j,wl(i,j),wobs(i,j),wtot
 9264          format('obs! missing line: iel,ion,i,j,wl(i,j),wobs(i,j),
     &wtot'
     &              ,4i5,f10.2,1pe12.3,10e12.3)
            endif
         enddo
      enddo

      if(imiss.eq.1) then
c         write(6,*)' there were missing lines in multiatom',iel,ion

      endif

      if(ierr.eq.1) then
         do i=1,n
c            write(6,929)i,x(i),te,ph(i),reco(i),denel
         enddo
      endif
      xi(n)=xplus
      kmax(iel,ion)= k-1
      if(iel.eq.1.and.ion.ge.1.and.x(1).lt.0.) stop
      RETURN
      END


      SUBROUTINE COLLEX(ION,TS)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NL=340,NLP1=NL+1)
      COMMON/NLEV/e00,NION,N,NP1H,NHMAX,NHMIN
      COMMON/A14/C(NL,NL),CI(NL),G(NL),E(NL),AS(NL,NL),
     &WL(NL,NL),DCDT(NL,NL),DCIDT(NL)
      DIMENSION CQ(NL,NL),CQP(NL,NL)
      NP1=N+1
      DT=0.0001
      TSPDT=(1.+DT)*TS
      DT=DT*TS
      CALL CEX(ION,TS,CQ)
      CALL CEX(ION,TSPDT,CQP)
      DO J=1,N
         IP1=J+1
         DO I=IP1,N
            EIJ=ABS(E(J)-E(I))
            TEV=TSPDT/1.1609E4
            ET=EIJ/TEV
            IF(EIJ/TEV.GT.700.) ET=700.
            CPIJ=G(J)*EXP(ET)*CQP(J,I)/G(I)
            C(J,I)=CQ(J,I)
            TEV=TS/1.1609E4
            ET=EIJ/TEV
            IF(EIJ/TEV.GT.700.) ET=700.
            C(I,J)=G(J)*EXP(ET)*C(J,I)/G(I)
            DCDT(J,I)=(CQP(J,I)-C(J,I))/DT
            DCDT(I,J)=(CPIJ-C(I,J))/DT
         enddo
      enddo
C
C     CALCULATE COLL. IONIZATION RATES
C
      DO I=1,N
         NLQ=I
         NMAX=NP1
         CALL CION(ION,NLQ,NMAX,TS,CI(I))
         CALL CION(ION,NLQ,NMAX,TSPDT,CIP)
         DCIDT(I)=(CIP-CI(I))/DT
      enddo
      RETURN
      END

