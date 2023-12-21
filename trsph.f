      SUBROUTINE SPHTRP(IMAX0,JF,R15,TEINT,RM,IMAXP)
C     ******************************************************************
C     *****
C
C     BINNING IN THIS ROUTINE COMPARED TO THE MAIN PROGRAM
C
C         RR        RR        RR                RR      RR      RR
C       TAUT      TAUT      TAUT              TAUT    TAUT    TAUT
C            DTAUT     DTAUT                      DTAUT   DTAUT(2)
C        X&S       X&S       X&S               X&S     X&S     X&S
C        I(+)       I                                    I     I(-)=0
C
C SPHTRP JTAU     JTAU-1    JTAU-2...............3       2       1
C        core                                                 surface
c        r -->
c
C INOUT = 1
C MAIN    1         2         3...............IMAX-2  IMAX-1   IMAX
C
C         R         R         R                  R       R        R
C            DE(2)     DE(3)                        DE     DE(IMAX)
C            TA(2)     TA(3)                        TA     TA(IMAX)
C         FD(1)                                               FD(IMAX)
C INOUT = 0
C            TA(IMAX)                               TA(3)  TA(2)
C         R(IMAX)                                        R(2)    R(1)
C         FD(IMAX)                                      FD(2)
c     r = rr = r(cgs)/(r15*1.e15)
c     dx = dr(cgs)
c     taupl = delta (tau)
c     tautot = tau(total)
c     ta = opacity in cm**-2
c     tauc = optical depth at jf=-3. used as optical depth scale
c     x = ta*r15*1.e15 ( x * (r(i)-r(i-1)) = opt. depth interval)
c
c     *****
C     ******************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
c      PARAMETER (MD=350,MDP1=MD+1,MPAR=1000,MPP1=1001)
      PARAMETER (MPAR=1500,MPP1=1501,m2pp1=2*mpp1)
c      include 'param'
      include 'parameters.h'
c      PARAMETER (NE1=-100,NE2=250,NE3=NE2+1)
      COMMON/SPECT/TEL,FD(MD,NE1:NE2),F0(NE1:NE2),IPARA
      COMMON/DIF/EM(MD,NE1:NE2),TAUPL(MD,NE1:NE2),TAUTOT(MD,NE1:NE2),
     &           EMC(MD,NE1:NE2)
      COMMON/RADIE/R(0:MD)
      COMMON/PHY/DEN(MD)
      COMMON/DXA/DX(MD)
      COMMON /CTRAN/X(MPAR),S(MPAR),XJ(MPAR),XH(MPAR),XK(MPAR)
     & ,FJ(MPAR),SOURCE(MPAR),JTAU0,JTAU1
      COMMON /TAUCQ/TAU(MPAR),TT(MPP1,MPAR),DTAULN(MPAR),JTAU
      COMMON/SPACE2/P(MPAR)
     & ,SP1(MPAR),SP2(MPAR),SP3(MPAR),AD(MPAR),BD(MPAR)
     & ,NIMPAC,KIMPAC(MPP1),PIMPAC(MPP1),MMU(MPAR)
     & ,TAUT(MPAR),DTAUT(MPAR)
     & ,PFEAU(MPP1,MPAR),XMU(MPP1,MPAR),FUN(MPP1),DER(M2PP1),DMU(MPP1)
      COMMON /CSPHER/NCORE,NATMOS,TAUM,RR(MPAR)
      COMMON/IOPAR/IOSHELL,IORAD1,IORAD2,IOFLUX
      COMMON/FRE/NINQ,JMIN,JJ
      COMMON/CINOUT/INOUT,IPULS
      common/contpoint/otscon(100),ncont,jcont(100),iioncont(100)
      COMMON/EDDFLUX/EDDFLUX(NE1:NE2),SURFJ(NE1:NE2)
      DIMENSION TEINT(MD),TA(MD),RM(MD),TAUC(MPAR)

      write(6,*)'inout,jtau',jf,inout,jtau
      write(0,*)'inout,jtau',jf,inout,jtau


      if(jf.eq.-14) then
         open(35,file='slasktrut')
         rewind 35
         write(35,*)IMAX0,JF,R15,IMAXP,MD,JMIN,JJ
         write(35,*)' rm'
         write(35,*)RM
         write(35,*)' te'
         write(35,*)TEINT
         do i=1,md
            write(35,*)'i = ',i
            write(35,*)(EM(i,j),j=JMIN,JJ)
            write(35,*)(EMC(i,j),j=JMIN,JJ)
            write(35,*)(TAUPL(i,j),j=JMIN,JJ)
            write(35,*)(TAUTOT(i,j),j=JMIN,JJ)
         enddo
         write(35,*)' r '
         write(35,*)r
         write(35,*)' den'
         write(35,*)den
         write(35,*)' dx'
         write(35,*)dx
         close (35)
         inittr=1
      endif
c add shells in the interior so that a reasonable coverage in radius is
c obtained
      if(inout.eq.0) then
         nadd=40
         dr=r(imax0)/nadd
         do i=1,nadd
            r(imax0+i)=(1.+1.e-5)*r(imax0)-i*dr
            den(imax0+i)=den(imax0)
            dx(imax0+i)=dr*r15*1.e15
            do j=jmin,jj
               EM(imax0+i,j)=em(imax0,j)*1.e-10
               EMc(imax0+i,j)=emc(imax0,j)*1.e-10
               taupl(imax0+i,j)=taupl(imax0,j)
            enddo
         enddo
         imax_old=imax0
         imax0=imax0+nadd
      endif

cfmod!!core      NCORE=2
c  this statement is probably overruled by ncore=2 below. this is not checked
      NCORE=50
cfmodhit
      NATMOS=100
c  determine depth to where J=S. 
      imax=imax0
      do i=2,imax0
         if(tautot(i,jf).gt.1.e8) then
            imax=i
            goto 133
         endif
      enddo
 133  continue
      DO I=2,IMAX0
         TA(I)=TAUPL(I,JF)/(DX(I)*DEN(I))
      enddo
      DO K=2,imax0
C SOURCE FCN IN SAME DIRECTION AS R (I.E. IN MAIN)
         XX=0.5*(TA(K+1)+TA(K))
         SS=0.5*(DEN(K)*EM(K,JF)+EM(K+1,JF)
     &        *DEN(K+1))/XX
         FD(K,JF)=SS
      ENDDO
      JTAU=IMAX
C     
C     OPACITY FOR ALL SHELLS
C
      DO I=2,IMAX
         TA(I)=TAUPL(I,JF)/(DX(I)*DEN(I))
      enddo
      RMIN=R(IMAX)
      IF(INOUT.EQ.1) RMIN=R(1)
C
C     CALCULATE OPACITY, SOURCE FUNCTION AND DISTANCE FROM THE SURFACE
C        FOR EACH SHELL
C
      IF(INOUT.EQ.0) then
         RR(1)=R(1)
         X(1)=TA(2)
         SOURCE(1)=DEN(2)*EM(2,JF)/X(1)
         TAU(1)=DEN(2)*(R(1)-R(2))/2.d0
         RR(JTAU)=R(IMAX)
         X(JTAU)=TA(IMAX)
         SOURCE(JTAU)=DEN(IMAX)*EM(IMAX,JF)/X(JTAU)
         nin=0
      ENDIF
      IF(INOUT.EQ.1) THEN
         RR(1)=R(IMAX)
         X(1)=TA(IMAX)
         SOURCE(1)=DEN(IMAX)*EM(IMAX,JF)/X(1)
         TAU(1)=DEN(IMAX)*(R(IMAX)-R(IMAX-1))/2.d0
         RR(JTAU)=R(1)
         X(JTAU)=TA(2)
         SOURCE(JTAU)=DEN(2)*EM(2,JF)/X(JTAU)
      ENDIF
      jfq=-99999
      if(jf.eq.jfq) then
         write(6,*)'inout,jtau',inout,jtau
      endif
      DO K=2,JTAU-1
         IF(INOUT.EQ.0) THEN
            RR(K)=R(K)
            X(K)=0.5*(TA(K+1)+TA(K))
            SOURCE(K)=0.5*(DEN(K)*EM(K,JF)+EM(K+1,JF)
     &           *DEN(K+1))/X(K)
            TAU(K)=TAU(K-1)+DEN(K)*(RR(K-1)-RR(K))
         ENDIF
         IF(INOUT.EQ.1) THEN
            RR(K)=R(IMAX-K+1)
            X(K)=0.5*(TA(IMAX-K+1)+TA(IMAX-K+2))
            SOURCE(K)=0.5*(DEN(IMAX-K+1)*EM(IMAX-K+1,JF)+EM(IMAX-K+2,JF)
     &           *DEN(IMAX-K+2))/X(K)
            TAU(K)=TAU(K-1)+DEN(IMAX-K+2)*(RR(K-1)-RR(K))
         ENDIF
         if(jf.eq.jfq) then
            PRINT 166,K,RR(K),DEN(k),TAU(K),x(k),em(k,jf),source(k)
 166        FORMAT(' tau ',2I4,10E11.4)
         endif
      ENDDO
c check if source fcn small 
      iso=0
      do k=1,jtau
         if(source(k).gt.1.d-30) then
            iso=1
         endif
      enddo
      if(iso.ne.1) then
         write(0,*)' small soure fcn at ',jf
         do k=1,jtau
            source(k)=1.e-30
         enddo
      endif
      IF(INOUT.EQ.1) THEN
         TAU(JTAU)=TAU(JTAU-1)+DEN(2)*(RR(JTAU-1)-RR(JTAU))
      ELSEIF(INOUT.EQ.0) THEN
         TAU(JTAU)=TAU(JTAU-1)+DEN(IMAX)*(RR(JTAU-1)-RR(JTAU))
      ENDIF
      if(r(jtau).gt.0..and.inout.eq.1) then
         nin=21
         dr=(rr(jtau)-0.00001)/nin
         do i=jtau+1,jtau+nin
            rr(i)=rr(jtau)-(i-jtau)*dr
            den(i)=den(jtau)
c!!            do j=jmin,jj
            source(i)=0.
            x(i)=1.e-5*x(jtau)
c!!            enddo
         enddo
         do k=jtau+1,jtau+nin
            TAU(K)=TAU(K-1)+DEN(jtau)*(RR(K-1)-RR(K))
         enddo
      endif
      jtau=jtau+nin
      DO  K=1,JTAU
         X(K)=R15*1.E15*X(K)
      enddo
C     CHECK IF SOURCE FUNCTION IS DIFFERENT FROM ZERO IN ANY SHELL
C     
      IQWE=0
      DO K=1,JTAU
         IF(SOURCE(K).GT.0.) IQWE=1
      ENDDO
      iof=1
      IF(IQWE.EQ.0.and.IOF.eq.1)PRINT*,' SOURCE FUNTION EQUAL TO ZERO'
     &     ,' FOR J= ',JF
C     CALCULATE THE TOTAL OPTICAL DEPTH FROM THE SURFACE
C     
      TAUC(1)=0.
      DO I=2,JTAU
         TAUC(I)=TAUPL(I,-3)+TAUC(I-1)
         IF(INOUT.EQ.1) TAUC(I)=TAUPL(IMAX-I+1,-3)+TAUC(I-1)
      ENDDO
c      DO I=1,JTAU
c     if(iqwe.eq.1) PRINT 221,I,RM(I),TAU(I),TAUC(I),RR(I),SOURCE(I)
c      ENDDO
c      DO  I=1,JTAU
c     if(iqwe.eq.1) PRINT 221,I,rm(i),em(I,jf),TA(I)       
c      ENDDO
 221  FORMAT(I5,8E11.4)
      IF(IQWE.EQ.0) GOTO 153
C     
C     ADD SHELLS SO THAT THEY ARE DISTRIBUTED UNIFORMILY IN OPTICAL
C     DEPTH AT ENERGY JF= -3.
C     
      I=2
 140  I=I+1
      IYX=0
      IYS=0
      IYTAUC=0
      IYTAU=0
 120  IF(IOF.EQ.0) GOTO 292
      IF(TAUC(I).LE.0..OR.TAUC(I-1).LE.0.) PRINT 828,I,TAUC(I),
     &     TAUC(I-1)
 828  FORMAT('TAUC LE ',I4,2E11.4)
      IF(X(I).LE.0..OR.X(I-1).LE.0.) PRINT 829,I,X(I),
     &     X(I-1)
 829  FORMAT('X LE ',I4,2E11.4)
      IF(TAU(I).LE.0..OR.TAU(I-1).LE.0.) PRINT 830,I,jf,TAU(I),
     &     TAU(I-1)
 830  FORMAT('TAU LE ',2I4,2E11.4)
      IF(SOURCE(I).LT.0..OR.SOURCE(I-1).LT.0.) PRINT 831,I,SOURCE(I),
     &     SOURCE(I-1)
 831  FORMAT('SOU LE ',I4,2E11.4)
 292  DTLOG=DLOG10(TAUC(I)/TAUC(I-1))
      IF(DTLOG.LT.0.50) GOTO 100
      IF(X(I).LE.0..OR.X(I-1).LE.0.)IYX=1
      IF(SOURCE(I).LE.0..OR.SOURCE(I-1).LE.0.) IYS=1
      IF(TAU(I).LE.0..OR.TAU(I-1).LE.0.) IYTAU=1
      IF(TAUC(I).LE.0..OR.TAUC(I-1).LE.0.) IYTAUC=1
      IF(IYS.EQ.1.AND.IOF.EQ.1)PRINT*,' NEG SOUR',I,SOURCE(I),
     &     SOURCE(I-1)
      IF(IYTAU.EQ.1.AND.IOF.EQ.1)PRINT*,' NEG TAU ',I,TAU(I),TAU(I-1)
      IF(IYTAUC.EQ.1.AND.IOF.EQ.1)PRINT*,' NEG TAUC ',I,TAUC(I),
     &     TAUC(I-1)
      IF(IYX.EQ.1.AND.IOF.EQ.1)PRINT*,' NEG X ',I,X(I),X(I-1)
      IF(IYTAUC.EQ.0) TAUM=0.5*(DLOG10(TAUC(I))+DLOG10(TAUC(I-1)))
      IF(IYX.EQ.0) XMEAN=0.5*(DLOG10(X(I))+DLOG10(X(I-1)))
      IF(IYTAU.EQ.0) TMEAN=0.5*(DLOG10(TAU(I))+DLOG10(TAU(I-1)))
      IF(IYS.EQ.0) SMEAN=0.5*(DLOG10(SOURCE(I))+DLOG10(SOURCE(I-1)))
      IF(IYTAUC.EQ.1) TAUM=0.5*(TAUC(I)+TAUC(I-1))
      IF(IYX.EQ.1) XMEAN=0.5*(X(I)+X(I-1))
      IF(IYTAU.EQ.1) TMEAN=0.5*(TAU(I)+TAU(I-1))
      IF(IYS.EQ.1) SMEAN=0.5*(SOURCE(I)+SOURCE(I-1))
      RRMEAN=0.5*(RR(I)+RR(I-1))
      IJM=JTAU-I+2
      IF(JF.EQ.jcont(1)-5) then
         RMEAN=0.5*(RM(IJM)+RM(IJM-1))
         TEMEAN=0.5*(TEINT(IJM)+TEINT(IJM-1))
      endif
      DO IJ=1,IJM
         II=JTAU-IJ+2
         X(II+1)=X(II)
         RR(II+1)=RR(II)
         SOURCE(II+1)=SOURCE(II)
         TAU(II+1)=TAU(II)
         TAUC(II+1)=TAUC(II)
      enddo
      IF(JF.eq.jcont(1)-5) then
         DO  IJ=JTAU,IJM,-1
            TEINT(IJ+1)=TEINT(IJ)
            RM(IJ+1)=RM(IJ)
         enddo
         TEINT(IJM)=TEMEAN
         RM(IJM)=RMEAN
      endif
C     PRINT 252,I,RMEAN,XMEAN,SMEAN,RM(I),RM(I-1)
 252  FORMAT(I4,8E11.4)
C     PRINT*,(RM(IK),IK=1,JTAU)
C     PRINT*,(RM(IK),IK=1,JTAU)
      IF(IYX.EQ.0) X(I)=10.**XMEAN
      IF(IYS.EQ.0) SOURCE(I)=10.**SMEAN
      IF(IYTAUC.EQ.0) TAUC(I)=10.**TAUM
      IF(IYTAU.EQ.0) TAU(I)=10.**TMEAN
      IF(IYX.EQ.1) X(I)=XMEAN
      IF(IYS.EQ.1) SOURCE(I)=SMEAN
      IF(IYTAUC.EQ.1) TAUC(I)=TAUM
      IF(IYTAU.EQ.1) TAU(I)=TMEAN
      RR(I)=RRMEAN
      JTAU=JTAU+1
C     PRINT*,JTAU,I
C     PRINT*,(RM(IK),IK=1,JTAU)
      GOTO 120
C100  PRINT 255,I,RM(I),TAUC(I),X(I),SOURCE(I)
 100  CONTINUE
 255  FORMAT(' I ',I4,5E11.4)
      IF(JTAU.GT.MD) GOTO 150
      IF(I.GE.JTAU) GOTO 150
      GOTO 140
 150  CONTINUE
      IMAXP=JTAU
C     PRINT 112,RR,TAU,SOURCE,X
 112  FORMAT(5E11.4)
C     IF(JF.EQ.-14) PRINT*,' RR'
C     IF(JF.EQ.-14) PRINT 901,(RR(K),K=1,JTAU)
C     IF(JF.EQ.-14) PRINT*,' TAU'
C     IF(JF.EQ.-14) PRINT 901,(TAU(K),K=1,JTAU)
C     PRINT *,' X'
C     PRINT 901,(X(K),K=1,JTAU)
C     PRINT *,' TAUC'
C     PRINT 901,(TAUC(K),K=1,JTAU)
C     PRINT*,' SOURCE'
C     PRINT 901,(SOURCE(K),K=1,JTAU)
 901  FORMAT(5E11.4)
C     PRINT 111,RR(JTAU),RR(1)
 111  FORMAT(' RMIN=',F8.3,3X,' RMAX=',F8.3)
C     
C     SOLV FOR CONTINUUM
C     
      IF(JF.EQ.22222) THEN
         WRITE(6,*)' k, rr, Source at J=',jf,inout
         DO K=1,JTAU
            WRITE(6,9288)K,RR(K),SOURCE(K),tautot(k,jf)
 9288       FORMAT(I5,1PE12.4,10E12.4)
         ENDDO
      ENDIF
      IF(IQWE.EQ.1) CALL TRANEQ
 153  CONTINUE
c      IF(IQWE.EQ.1.AND.IOF.EQ.1)PRINT*,' INTENSITIES AT ENERGY J=',JF
c      DO 151 I=1,JTAU
c         IF(IQWE.EQ.1.AND.IOF.EQ.1)PRINT 152,I,RM(I),RR(I),TT(1,I),
c     &        tauc(i),tautot(i,jf),SOURCE(I),XJ(I),XH(I)
c 151  CONTINUE
 152  FORMAT(I5,2E11.4,7E11.3)
      JTAU=JTAU-NIN
      DO  I=1,JTAU
         IF(IQWE.EQ.0) XJ(I)=0.
         IF(XJ(I).LT.0..AND.IOF.EQ.1)
     &        PRINT*,' NEG MEAN INT. J= ',i,JF,XJ(i),source(i),x(i),
     &  em(imax-i+1,jf)
         IF(XJ(I).LT.0.) XJ(I)=0.
         IF(XH(I).LT.0.) XH(I)=-XH(I)
         FD(I+1,JF)=XJ(I)
         IF(I.EQ.2) THEN
            EDDFLUX(JF)=4.*3.1416*XH(2)
            SURFJ(JF)=XJ(2)
         ENDIF
         IF(INOUT.EQ.1) FD(JTAU-I+1,JF)=XJ(I)
      enddo
c the shells for which J=S and H=0
      DO  I=JTAU+1,imax0
         XH(I)=0.
      enddo
c      IF(IQWE.EQ.1.AND.IOF.EQ.1) then
c      do i=1,imax0
c         write(6,9289)i,RM(I),RR(I),TT(1,I),
c     &        tauc(i),tautot(i,jf),SOURCE(I),XJ(I),XH(I),FD(I+1,JF)
c      enddo
c      endif
 9289 format(i5,1pe12.4,10e12.4)
      imax0=imax_old
      imaxp=jtau-nadd
      RETURN
      END

      SUBROUTINE TRANEQ
      IMPLICIT REAL*8(A-H,O-Z)
c      PARAMETER (MD=350,MDP1=MD+1,MPAR=1000,MPP1=1001)
      PARAMETER (MPAR=1500,MPP1=1501,m2pp1=2*mpp1)
c     include 'PARAM'
      include 'parameters.h'
C
C TRANEQ SOLVES THE TRANSFER EQUATION INCLUDING CONTINUUM SCATTERING.
C FEATURES:
C
C 1. CANNONS PERTURBATION TECHNIQUE IS USED ON THE ANGULAR QUADRATURE.
C    THE BASIC IDEA IN THIS TECHNIQUE IS TO REPLACE THE INVERSION OF
C    A COMPLICATED (MMU ORDER) OPERATOR WITH THE INVERSION OF A SIMPLE
C    OPERATOR (ONE POINT=EDDINGTON APPROXIMATION), PLUS ITERATION ON
C    THE ERROR.
C 2. A TRICK DUE TO ROBERT STEIN (PRIV. COMM., 1979) IS USED TO
C    ELIMINATE THE NEED FOR DOUBLE PRECISION STORAGE OF THE MATRIX
C    ELEMENTS. THE IDEA IS TO STORE THE (SMALL) SUM OF THE THREE
C    MATRIX ELEMENTS ON A ROW, INSTEAD OF THE (LARGE) DIAGONAL ELE-
C    MENT.
C 3. THE SOLUTION IS A CUBIC SPLINE, RATHER THAN A PIECE-WISE
C    QUADRATIC FUNCTION. THIS IS ACCOMPLISHED WITH THE CORRECTION
C    TERMS AD AND BD IN SUBROUTINE TRANFR.
C 4. A BOUNDARY CONDITION WHICH INCLUDES AN ESTIMATED INFALLING
C    RADIATION MAKES THE SOLUTION GOOD ALSO FOR VALUES OF X+S
C    LARGE COMPARED WITH 1./TAU(1). A LOGARITHMIC TAU-SCALE
C    SHOULD BE USED.
C
C THIS VERSION OF TRANEQ IS COMPATIBLE WITH PREVIOUS TRANEQ'S.
C 79.06.21 *NORD*
C
      COMMON /CTRAN/X(MPAR),S(MPAR),XJ(MPAR),XH(MPAR),XK(MPAR)
     & ,FJ(MPAR),SOURCE(MPAR),JTAU0,JTAU1
      COMMON /TAUCQ/TAU(MPAR),TT(MPP1,MPAR),DTAULN(MPAR),JTAU
      COMMON/SPACE2/P(MPAR)
     & ,SP1(MPAR),SP2(MPAR),SP3(MPAR),AD(MPAR),BD(MPAR)
     & ,NIMPAC,KIMPAC(MPP1),PIMPAC(MPP1),MMU(MPAR)
     & ,TAUT(MPAR),DTAUT(MPAR)
     & ,PFEAU(MPP1,MPAR),XMU(MPP1,MPAR),FUN(MPP1),DER(M2PP1),DMU(MPP1)
      COMMON /CSPHER/NCORE,NATMOS,TAUM,RR(MPAR)
      COMMON /TRDBUG/IDEBUG
      LOGICAL DEBUG
      DIMENSION A(9)
      DATA ITMAX/7/,DEBUG/.FALSE./
c      write(6,*)' traneq nimpac ',nimpac,jtau
C
109   IDEBUG=0
      IF (DEBUG) PRINT 103,X,SOURCE,XJ,XH
103   FORMAT(' X,S,XJ,D2B='/(/4(1X,10E10.4/)))
C
C     PRINT 999,RR
999   FORMAT(' TRR ',5E11.4)
C CALCULATE THE MATRIX ELEMENTS
      CALL TRRAYS
c!!      CALL TRANGS
c      write(6,*)'jtau b ',jtau
      CALL TRANFR
C     CALL FORMAL
      NIMP1=NIMPAC+1
      IF (IDEBUG.GT.1) GO TO 150
C
C
C CALCULATE MOMENTS, AND CHECK DEBUG CONTROL
c      CALL TRMO
C150   IF (DEBUG.AND.IDEBUG.GT.1) STOP
 150   CONTINUE
      IF (DEBUG.AND.IDEBUG.EQ.1) IDEBUG=0
      DEBUG=IDEBUG.GT.1
C!!!! IF (DEBUG) GO TO 109
C
      RETURN
      END

      SUBROUTINE TRANFR
       IMPLICIT REAL*8(A-H,O-Z)
c      PARAMETER (MD=350,MDP1=MD+1,MPAR=1000,MPP1=1001)
      PARAMETER (MPAR=1500,MPP1=1501,m2pp1=2*mpp1)
c     include 'PARAM'
      include 'parameters.h'
C
C FORMAL SOLVES THE TRANSFER EQUATION WITH GIVEN SOURCE FUNCTION 'SOURCE
C 'ERROR' IS THE RESULTING ERROR IN THE DEFINITION OF THE CONTINUUM
C SCATTERING SOURCE FUNCTION. TRANSFR CALCULATES THE MATRIX ELEMENTS
C OF THE PROBLEM. INTENSITIES AT TAU=0 ARE RETURNED IN /CSURF/.
C 80.08.05 *NORD*
C SPLINE APPROXIMATION USED FOR DIFF. EQNS.
C A*P(K-1)+B*P(K)+C*P(K+1)=L(K)
C SP1=-A
C SP3=-C
C SP2=-A-B-C (=1., EXCEPT AT BOUNDARIES) BOB STEINS TRICK
C P=-L
C
      COMMON /CTRAN/X(MPAR),S(MPAR),XJ(MPAR),XH(MPAR),XK(MPAR)
     & ,FJ(MPAR),SOURCE(MPAR),JTAU0,JTAU1
      COMMON /TAUCQ/TAU(MPAR),TT(MPP1,MPAR),DTAULN(MPAR),JTAU
      COMMON/SPACE2/P(MPAR)
     & ,SP1(MPAR),SP2(MPAR),SP3(MPAR),AD(MPAR),BD(MPAR)
     & ,NIMPAC,KIMPAC(MPP1),PIMPAC(MPP1),MMU(MPAR)
     & ,TAUT(MPAR),DTAUT(MPAR)
     & ,PFEAU(MPP1,MPAR),XMU(MPP1,MPAR),FUN(MPP1),DER(M2PP1),DMU(MPP1)
      COMMON /CSPHER/NCORE,NATMOS,TAUM,RR(MPAR)
      COMMON /TRDBUG/IDEBUG
      common/frek/jf
      DIMENSION DTCNT(MPP1),VFE(MPP1,MPAR)
c      write(0,*)' tranfr ',jf,nimpac
c      write(6,*)' tranfr ',jf,nimpac
c      if(jf.eq.-249.or.jf.eq.-10) then
c         write(6,*)'rr traneq a ',rr(2),rr(100),rr(500)
c      endif
C

      xjmax=0.

      jfq=-99999

C MU LOOP
      DO  I=1,NIMPAC
c         write(0,*)' mu ',i
C INNERMOST SHELL FOR IMPACT RAY I (FROM SURFACE)
         NTAU=KIMPAC(I)
         NTAU1=NTAU-1
C
C CALCULATE DTAUT ALONG THE RAY
         ZOLD=0.0
         DO K=1,NTAU
c            write(0,*)'k ',k,ntau-k+1,i
            if((RR(NTAU-K+1)-PIMPAC(I)).le.0.d0) then
               z=0.
               iq = 1
            else
               Z=DSQRT(RR(NTAU-K+1)**2-PIMPAC(I)**2)
               iq = 2
            endif
            if(abs(rr(ntau-k+1)).le.1.e-33) 
     &           write(6,*)' overflow i trans!!! ',k,ntau,z,rr(ntau-k+2)
            XMU(I,NTAU-K+1)=-Z/RR(NTAU-K+1)
            IF (K.EQ.1) GO TO 100
            DZ=Z-ZOLD
            DZDR=DZ/(RR(NTAU-K+1)-RR(NTAU-K+2))
            DTAUT(NTAU-K+2)=DZDR*0.5d0*(X(NTAU-K+1)+S(NTAU-K+1)
     &           +X(NTAU-K+2)+S(NTAU-K+2))*(TAU(NTAU-K+2)-TAU(NTAU-K+1))
 100        ZOLD=Z
            if(jf.eq.jfq) then
               write(6,9227)k,ntau-k+2,RR(NTAU-K+1)-RR(NTAU-K+2),
     &              RR(NTAU-K+1),RR(NTAU-K+2),pimpac(i),z,dz,
     &              RR(NTAU-K+1)-PIMPAC(I),TAU(NTAU-K+2),TAU(NTAU-K+1)
               write(6,9224)k,ntau-k+2,dzdr,x(ntau-k+1),s(ntau-k+1),
     &              x(ntau-k+2),s(ntau-k+2),dtaut(ntau-k+2)
 9227          format('rr    ',2i5,1pe13.5,20e13.5)
 9224          format('dtaut ',2i5,1pe13.5,20e13.5)
            endif
         ENDDO
         TAUT(1)=DZDR*(X(1)+S(1))*TAU(1)
         DO K=2,NTAU
            TAUT(K)=TAUT(K-1)+DTAUT(K)
            TT(I,K)=TAUT(K)
         ENDDO
C
C K=1 = SURFACE CONDITION
         A=1./DTAUT(2)
         B=A**2
         SP2(1)=1.+2.d0*A
         SP3(1)=-2.d0*B
         EX=TAUT(1)*(1.-0.5d0*TAUT(1)*(1.-0.3333333333*TAUT(1)))
         IF (TAUT(1).GT.0.1) EX=1.-DEXP(-TAUT(1))
         SP2(1)=SP2(1)/(1.+2.d0*A*EX)
         SP3(1)=SP3(1)/(1.+2.d0*A*EX)
C     
C     K=2,NTAU-1
         DO  K=2,NTAU1
            DTAUC=0.5d0*(DTAUT(K)+DTAUT(K+1))
            AD(K)=0.166666666667*DTAUT(K)/DTAUC
            BD(K)=0.166666666667*DTAUT(K+1)/DTAUC
C     C..  AD(K)=0.
C...  BD(K)=0.
            SP1(K)=-1./(DTAUT(K)*DTAUC)+AD(K)
            SP2(K)=1.
            SP3(K)=-1./(DTAUT(K+1)*DTAUC)+BD(K)

            if((i.le.10.or.i.ge.ntau-10).and.jf.eq.jfq) then
               write(6,9223)k,ad(k),bd(k),sp1(k),dtaut(k),dtaut(k+1)
 9223          format('dtau ',i5,1pe13.5,10e13.5)
            endif

         ENDDO
C     
C     K=NTAU
C     FIRST THE RAYS HITTING THE CORE

         AD(NTAU)=0.0d0
         BD(NTAU)=0.0d0
         SP1(NTAU)=-1.d0
         SP2(NTAU)=DTAUC+0.5d0*DTAUC**2
         SP3(NTAU)=0.0d0
         DTCNT(I)=0.5d0*(DTAUT(NTAU-1)+DTAUT(NTAU))
C!!!  If Cloudy!! 
c     fmod
         IF (I.GT.NCORE) THEN
c     fmodhit
c     IF (I.LE.NCORE) GO TO 120
C     NOW THE RAYS MISSING THE CORE
            AD(NTAU)=0.3333333
            SP1(NTAU)=0.3333333-2.d0/DTAUT(NTAU)**2
            SP2(NTAU)=1.
         ENDIF
c         write(0,*)'tr b',dtaut(ntau)
         if(jf.eq.jfq) write(6,*)'tr b',jf,dtaut(ntau)
C ELIMINATE SUBDIAGONAL
         DO  K=1,NTAU1
            SP1(K)=-SP1(K+1)/(SP2(K)-SP3(K))
            SP2(K+1)=SP2(K+1)+SP1(K)*SP2(K)
            SP2(K)=SP2(K)-SP3(K)
            if(jf.eq. jfq) then
               write(6,9225)k,sp1(k),sp1(k+1),sp2(k),sp3(k),dtaut(k)
 9225          format(' sp ',i5,1pe12.3,10e12.3)
            endif
         ENDDO
C
C INITIATE
         P(1)=SOURCE(1)
         DO  K=2,NTAU
            P(K)=(1.-AD(K)-BD(K))*SOURCE(K)+AD(K)*SOURCE(K-1)+
     &           BD(K)*SOURCE(K+1)
            if(k.ge.ntau-3.and.jf.eq.jfq) then
               write(6,9222)k,ad(k),bd(k),source(k),source(k+1),p(k)
 9222          format('p(k) ',i5,1pe13.5,10e13.5)
            endif
         ENDDO
C INNER BOUNDARY VALUE FOR RAYS HITTING THE CORE
C OUTGOING RADIATION.
C IF-DIFFUSION APPROX. USE I(+)=S(N)+(S(N)-S(N-1))/DTAUC
         IF(I.LE.NCORE) DIPLUS=0.
         DTAUC=DTCNT(I)
cfmod
         if(i.le.ncore) diplus=SOURCE(jtau1)+(SOURCE(jtau1)-
     &        SOURCE(jtau1-1))/DTAUC
cfmodhit

c!!! mod 070819
         diplus = 0.

         IF(I.LE.NCORE) THEN
            P(JTAU1)=0.5d0*DTAUC**2*SOURCE(JTAU1)+DTAUC*DIPLUS
            if(jf.eq.jfq) then
               write(6,9232)jtau1,i,dtauc,diplus,source(jtau1),p(jtau1)
            endif
c??? only for diff apprx.??
            diplus=SOURCE(ntau)+(SOURCE(ntau)-
     &        SOURCE(ntau-1))/DTAUC

c!!! mod 070819
            diplus = 0.

            P(NTAU)=0.5d0*DTAUC**2*SOURCE(NTAU)+DTAUC*DIPLUS
            if(jf.eq.jfq) then
               write(6,9232)ntau,i,dtauc,diplus,source(ntau),p(ntau)
 9232          format('p(k)nc ',2i5,1pe13.5,10e13.5)
            endif
         ENDIF

C
C ACCUMULATE RIGHT HAND SIDE
         DO  K=1,NTAU1
            P(K+1)=P(K+1)+SP1(K)*P(K)
            if(jf.eq.jfq.and.k.ge.ntau1-3) then
               write(6,9236)i,k,p(ntau),sp2(ntau),
     &              PFEAU(i,ntau)
 9236          format(' pfeau bsq ',2i5,1pe12.3,10e12.3)
            endif
         ENDDO
C
C BACKSUBSTITUTE
         PFEAU(I,NTAU)=P(NTAU)/SP2(NTAU)
         if(jf.eq.jfq) then
            write(6,9226)i,ntau,p(ntau),sp2(ntau),
     &           PFEAU(i,ntau)
 9226       format(' pfeau bs1 ',2i5,1pe12.3,10e12.3)
         endif
         PFMAX=0.
         DO K=1,NTAU1
            PFEAU(I,NTAU-K)=(P(NTAU-K)-
     &           SP3(NTAU-K)*PFEAU(I,NTAU-K+1))/SP2(NTAU-K)
            if(jf.eq.jfq.and.k.le.ntau-3) then
               write(6,9229)i,k,p(ntau-k),sp2(ntau-k),sp3(ntau-k),
     &              PFEAU(I,ntau-k+1),PFEAU(I,ntau-k)
 9229          format(' pfeau bs ',2i5,1pe12.3,10e12.3)
            endif
c     IF (PFEAU(I,NTAU-K).LE.0.0) GO TO 230
            PFMAX=DMAX1(PFEAU(I,NTAU-K),PFMAX)
         ENDDO
         DO  K=1,NTAU1

            IF (PFEAU(I,NTAU-K).LE.0.) THEN
c            IF (PFEAU(I,NTAU-K).LE.-1.E-5*ABS(PFMAX)) THEN
c               WRITE(43,935)I,NTAU-K,PFEAU(I,NTAU-K),PFMAX
 935           FORMAT(' PF LE 0 ',2I4,4E12.4)
               PFEAU(I,NTAU-K)=1.d-30
            ENDIF
         ENDDO

         IF(JF.EQ.jfq) THEN
            WRITE(6,*)' I=',I,'PIMPAC=',PIMPAC(I),'NCORE',NCORE,
     &           'ntau',ntau
            write(6,9279)diplus,SOURCE(jtau1),SOURCE(jtau1),
     &        SOURCE(jtau1-1),DTAUC
 9279       format(' diplus ',1pe12.4,10e12.4)
c            DO K=1,JTAU1
            DO K=1,jTAU
               if(k.le.20.or.k.ge.jtau-10) then
                  WRITE(6,9288)i,K,RR(K),SOURCE(K),DTAUT(K),P(K),
     &                 pfeau(i,k),sp1(k),sp2(k)
 9288             FORMAT(' pf ',2I5,1PE12.4,10E12.4)
               endif
            ENDDO
         ENDIF

C
C END MU LOOP

         DO K=1,NTAU1
            VFE(I,K)=(PFEAU(I,K+1)-PFEAU(I,K))/
     &           (TT(I,K+1)-TT(I,K))
            IF(VFE(I,K).LE.0.) VFE(I,K)=1.E-30
         ENDDO
C     PRINT*,I,PIMPAC(I)
         DO  K=1,NTAU1,2
            FMIN=PFEAU(I,K)-VFE(I,K)
            FPLU=PFEAU(I,K)+VFE(I,K)
C     PRINT 923,K,TT(I,K),PFEAU(I,K),VFE(I,K),FMIN,FPLU
         ENDDO
 923  FORMAT(I4,5E11.4)
      ENDDO
C
C INTERPOLATE TO PFEAU AT MU=0 FOR THOSE K THAT HAVE NO RAY
      DO  K=1,JTAU1
         II=MMU(K)
         IF (KIMPAC(II).ne.K) then
            PX=-XMU(II-2,K)/(XMU(II-1,K)-XMU(II-2,K))
            QX=1.-PX
            IF(VFE(II-2,K).LE.0..OR.VFE(II-1,K).LE.0.) PRINT*,
     &           II,K,VFE(II-2,K),VFE(II-1,K)
            IF(PFEAU(II-2,K).LE.0..OR.PFEAU(II-1,K).LE.0.) PRINT*,
     &           II,K,PFEAU(II-2,K),PFEAU(II-1,K)
            PFEAU(II,K)=DEXP(DLOG(PFEAU(II-2,K))*QX+DLOG(PFEAU(II-1,K))
     &           *PX)
            if(PFEAU(II-2,K).gt.1.d-30.and.PFEAU(II-1,K).gt.1.d-30) then
               PFEAU(II,K)=DEXP(DLOG(PFEAU(II-2,K))*QX+
     &              DLOG(PFEAU(II-1,K))*PX)
            else
               PFEAU(II,K)=1.d-30
            endif
            VFE(II,K)=DEXP(DLOG(VFE(II-2,K))*QX+DLOG(VFE(II-1,K))*PX)
            if(VFE(II-2,K).gt.1.d-30.and.VFE(II-1,K).gt.1.d-30) then
               VFE(II,K)=DEXP(DLOG(VFE(II-2,K))*QX+
     &              DLOG(VFE(II-1,K))*PX)
            else
               VFE(II,K)=1.d-30
            endif
            if(jf.eq. jfq) then
c               write(6,9228)ii,k,PFEAU(II-2,K),PFEAU(II-1,K),px
 9228          format(' pfeau inter ',2i5,1pe12.3,10e12.3)
            endif
         endif
      enddo

C
C CALCULATE MEAN INTENSITY
      DO K=1,JTAU1-1
         XJ(K)=0.
         NMU=MMU(K)
         DO  I=2,NMU
            DMU(I)=XMU(I,K)-XMU(I-1,K)
            DER(I)=(PFEAU(I,K)-PFEAU(I-1,K))/DMU(I)
            XJ(K)=XJ(K)+DMU(I)*(PFEAU(I,K)+PFEAU(I-1,K))
            if(k.eq.2.or.k.eq.100.or.k.eq.500) then
c               write(6,9277)k,i,xmu(i,k),dmu(i),pfeau(i,k)
 9277          format(' pfaa ',2i5,1pe12.3,10e12.3)
            endif
         enddo
         XJ(K)=XJ(K)*6.
         NMU1=NMU-1
         DO  I=2,NMU1
            XJ(K)=XJ(K)+(DMU(I+1)-DMU(I))*(DMU(I)*DER(I+1)+
     &           DMU(I+1)*DER(I))
         enddo
         XJ(K)=XJ(K)+DMU(2)**2*DER(2)-DMU(NMU)**2*DER(NMU)
         XJ(K)=XJ(K)*0.083333333
      enddo

C CALCULATE FLUX
      DO K=1,JTAU1-1
         XH(K)=0.
         NMU=MMU(K)
         DO I=2,NMU
            DMU(I)=XMU(I,K)-XMU(I-1,K)
            DER(I)=(XMU(I,K)*VFE(I,K)-XMU(I-1,K)*VFE(I-1,K))/DMU(I)
            XH(K)=XH(K)+DMU(I)*(XMU(I,K)*VFE(I,K)+XMU(I-1,K)*VFE(I-1,K))
         enddo
         XH(K)=XH(K)*6.
         NMU1=NMU-1
         DO I=2,NMU1
            XH(K)=XH(K)+(DMU(I+1)-DMU(I))*(DMU(I)*DER(I+1)+DMU(I+1)*
     &           DER(I))
         ENDDO
         XH(K)=XH(K)+DMU(2)**2*DER(2)-DMU(NMU)**2*DER(NMU)
         XH(K)=XH(K)*0.083333333
      enddo

      NCORP1=NCORE+1
      DO K=1,JTAU1-1
         XH(K)=0.
         XJ(K)=0.
         NMU=MMU(K)
         DO  I=2,NMU
            DMU(I)=XMU(I,K)-XMU(I-1,K)
            IF(I.NE.NCORP1)XJ(K)=XJ(K)+DMU(I)*(PFEAU(I,K)+PFEAU(I-1,K))
            IF(I.NE.NCORP1)XH(K)=XH(K)+DMU(I)*(XMU(I,K)*VFE(I,K)+
     &           XMU(I-1,K)*VFE(I-1,K))
C     AVOID DISCONTINUITY AT THE LIMB FOR THE MOMENTS
            IF(I.EQ.NCORP1)XJ(K)=XJ(K)+DMU(I)*(PFEAU(I,K)+PFEAU(I,K))
            IF(I.EQ.NCORP1)XH(K)=XH(K)+DMU(I)*(XMU(I,K)*VFE(I,K)+
     &           XMU(I-1,K)*VFE(I,K))
C     IF(K.GE.37) PRINT 828,I,DMU(I),VFE(I,K),PFEAU(I,K),XH(K),XJ(K)
 828        FORMAT(I4,5E12.5)
         enddo
         XH(K)=XH(K)*0.5d0
         XJ(K)=XJ(K)*0.5d0
      enddo

C--------------------------------------------------------------------
C
C
C FLUX AT TAU(1)
      XH(1)=TRQUAD(MMU(1),XMU,FUN,DER)
C
C CALCULATE SECOND MOMENT XK
      DO  K=1,JTAU1-1
         NMU=MMU(K)
         DO I=1,NMU
c            write(0,*)i,k,nmu
            FUN(I)=PFEAU(I,K)*XMU(I,K)**2
         enddo
c         write(0,*)' to trq',k,nmu
         XK(K)=TRQUAD(NMU,XMU(1,K),FUN,DER)
      enddo
c      write(0,*)'tr j'
C
C CALCULATE FIRST MOMENT, XH, FROM MOMENT RELATION
      DO  K=2,JTAU
         XH(K)=(XK(K)-XK(K-1)+(XJ(K)+XJ(K-1)-3.*(XK(K)+XK(K-1)))
     &        *(RR(K-1)-RR(K))/(RR(K)+RR(K-1)))
     &        *2.d0/((TAU(K)-TAU(K-1))*(X(K)+S(K)+X(K-1)+S(K-1)))
      enddo
C
C CALCULATE FIRST MOMENT XH BY QUADRATURE
      DO K=1,JTAU1-1
         NMU=MMU(K)
         DO I=1,NMU
            FUN(I)=PFEAU(I,K)*XMU(I,K)
         enddo
         XH(K)=TRQUAD(NMU,XMU(1,K),FUN,DER)
      enddo
c      write(0,*)'tr k'
c      write(6,*)' out from aa',jtau
      do k=1,jtau
         r2 = 1.0
         r1 = 0.9
         r = rr(k)
         xja = 0.5 * r * (r2/r + (r2**2-r**2)/(2.*r**2)*log((r+r2)/
     &       abs(r-r2))  - r1/r - 
     &        (r1**2-r**2)/(2.*r)*log((r+r1)/abs(r-r1)))

         if(abs(xj(k)).gt.xjmax) then
            xjmax=abs(xj(k))
            xjabs=xj(k)
            kmax=k
         endif
         if(jf.eq.jfq) then
            write(6,913)k,rr(k),taut(k),xj(k),xh(k)
 913        format(' aan',i5,1pe19.10,10e12.3)
         endif
c         write(66,914)k,rr(k),taut(k),xj(k),xh(k)
 914     format(i5,1pe19.10,10e12.3)
      enddo

      if(xjmax.ge.1.e10) then
         write(6,*)' xjmax, jmax ',jf,kmax,xjmax,xjabs
         write(0,*)' xjmax, jmax ',jf,kmax,xjmax,xjabs
         do k=1,jtau
            write(6,913)k,jf,rr(k),taut(k),xj(k),xh(k)
         enddo
      endif

C

c      if(jf.eq.-249.or.jf.eq.-10) then
c         write(6,*)'rr traneq b ',rr(2),rr(100),rr(500)
c      endif

c         write(0,*)'tr l'

      RETURN
C------------------------------------------------------------------
C
C EMERGENCY EXIT
230   KK=NTAU-K
      IDEBUG=2
      PRINT 232,I,KK,JTAU0,JTAU1,NTAU
232   FORMAT('0NON-POSITIVE RESULT AT I,K,J0,J1,N=',5I3)
      GO TO 233
231   KK=0
      IDEBUG=3
C     PRINT 232,I,KK,JTAU0,JTAU1,NTAU
233   PRINT 237,NCORE,EX
237   FORMAT(' NCORE,EX(I)=',I3,3e12.4)
      PRINT 238
238   FORMAT(' DEBUG INFORMATION WRITTEN ON UNIT 13')
      RETURN
      END

      DOUBLE PRECISION FUNCTION TRQUAD(N,X,F,W)
       IMPLICIT REAL*8(A-H,O-Z)
c      PARAMETER (MD=350,MDP1=MD+1,MPAR=1500,MPP1=1501,M2PP1=2*MPP1)			
      PARAMETER (MPAR=1500,MPP1=1501,M2PP1=2*MPP1)
c     include 'PARAM'
      include 'parameters.h'
      DIMENSION X(N),F(N),W(M2PP1)
c skall det verkligen vara mpp1? ursporungligen m2pp1. 
c dock inte konsistent med der(mpp1)
c      DIMENSION X(N),F(N),W(MPP1)
C
C TRAPEZOIDAL QUADRATURE PLUS NEXT ORDER CORRECTION FOR NON-
C -EQUIDISTANT GRID.
c      write(0,*)'trq ',n
      N1=N-1
      Q=0.
      DO 100 K=2,N
      W(K)=X(K)-X(K-1)
      W(N+K)=(F(K)-F(K-1))/W(K)
100   Q=Q+W(K)*(F(K-1)+F(K))
      Q=Q*6.
      DO 101 K=2,N1
101   Q=Q+(W(K+1)-W(K))*(W(K)*W(N+K+1)+W(K+1)*W(N+K))
      W1=((W(2)+0.5d0*W(3))*W(N+2)-0.5d0*W(2)*W(N+3))*2.0d0/(W(2)+W(3))
      WN=((W(N)+0.5d0*W(N1))*W(N+N)-0.5d0*W(N)*W(N+N1))*2.0d0/(W(N)+
     &     W(N1))
      Q=0.083333333*(Q+W(2)**2*W1-W(N)**2*WN)
      TRQUAD=Q
c      write(0,*)'trq ut'
      RETURN
      END

      SUBROUTINE TRRAYS
       IMPLICIT REAL*8(A-H,O-Z)
c      PARAMETER (MD=350,MDP1=MD+1,MPAR=1500,MPP1=1501)
      PARAMETER (MPAR=1500,MPP1=1501,m2pp1=2*mpp1)
c     include 'PARAM'
      include 'parameters.h'
      COMMON /CTRAN/X(MPAR),S(MPAR),XJ(MPAR),XH(MPAR),XK(MPAR)
     & ,FJ(MPAR),SOURCE(MPAR),JTAU0,JTAU1
      COMMON /TAUCQ/TAU(MPAR),TT(MPP1,MPAR),DTAULN(MPAR),JTAU
      COMMON/SPACE2/P(MPAR)
     & ,SP1(MPAR),SP2(MPAR),SP3(MPAR),AD(MPAR),BD(MPAR)
     & ,NIMPAC,KIMPAC(MPP1),PIMPAC(MPP1),MMU(MPAR)
     & ,TAUT(MPAR),DTAUT(MPAR)
     & ,PFEAU(MPP1,MPAR),XMU(MPP1,MPAR),FUN(MPP1),DER(M2PP1),DMU(MPP1)
      COMMON /CSPHER/NCORE,NATMOS,TAUM,RR(MPAR)
      common/frek/jf
C
C FIND JTAU0 (THERMALIZATION DEPTH) AND JTAU1 ( TAUM TIMES LARGER)
C TAUM=DEPTH WHERE DIFFUSION STARTS
c      write(6,*)' trrays ',jtau
      TAUT(1)=TAU(1)*(X(1)+S(1))
      JTAU0=JTAU
      taum=1.e4
      DO K=2,JTAU
         TAUT(K)=TAUT(K-1)+0.5d0*(X(K)+S(K)+X(K-1)+
     &        S(K-1))*(TAU(K)-TAU(K-1))
         if(TAUT(K).gt.taum) goto 221
      enddo
 221  JTAU1=MAX0(k,3)
      ncore=jtau-jtau1+1
c!!!
      ncore=2
c!! mod 070819
      ncore = 0
      jtau1=jtau
C
C DISTRIBUTE RAYS IN CORE  WITH IMPACT INSIDE JTAU1
C I=1  P=0.
C I=NCORE P=R(JTAU1)
C KIMPAC(I)= LABEL OF INNERMOST SHELL WHICH RAY I HITS
C            INSIDE NCORE = JTAU1
C            OUTSIDE RADIAL SHELL AT Z=0.
      RRK=RR(JTAU0)
      NCORE1=NCORE-1
C DISTRIBUTE RAYS UNIFORMILY IN Z (AND THUS MU) FOR P=0 TO P=R(JTAU1)
C
      DR=RR(JTAU1)/NCORE1
      DO 110 I=1,NCORE1
C RRK = Z
      PIMPAC(I)=DSQRT(RR(JTAU1)**2-RRK**2)
      RRK=RRK-DR
110   KIMPAC(I)=JTAU1
c use radial shells also for core rays
      do i=1,ncore1
         pimpac(i)=rr(i)
         kimpac(i)=i
      enddo
C
C RAYS IN ATMOSPHERE WITH IMPACT OUTSIDE JTAU1
C MEAN TAU FACTOR IN ATMOSPHERE
      TAUFCT=0.
C     RADFCT=10.**(-DLOG10(RR(JTAU1)/RR(1))/(NATMOS-1.))
CVD$  NOCONCUR
      I=NCORE
cf mod
      i=1
c      PRINT*,NCORE,NATMOS,JTAU1,TAUFCT
      KI=JTAU1
120   KIMPAC(I)=KI
      PIMPAC(I)=RR(KI)
      KIP=KI
121   KI=KI-1
C IS RATIO OF SHELL DEPTHS LARGER THAN THE CHOOSEN INTERVAL?
CVD$  NOCONCUR
      IF (TAU(KIP)/TAU(KI).LT.TAUFCT.AND.KI.GT.1) GO TO 121
C     IF (RR(KI)/RR(KIP).LT.RADFCT.AND.KI.GT.1) GO TO 121
      I=I+1
      IF (KI.GE.3) GO TO 120
      KIMPAC(I)=0
      PIMPAC(I)=RR(1)
      NIMPAC=I-2
cf mod
      nimpac=i-3
      if(jf.eq.jfq) then
         write(6,*)'nimpac',nimpac,ncore,jtau,jtau1
         do ik=1,jtau
            PRINT 955,ik,kimpac(ik),PIMPAC(Ik),rr(ik)
         enddo
 955     FORMAT(' PI ',2i5,1pe19.10,5E19.10)
      endif
      IF (NIMPAC.LT.MPP1) GO TO 131
      PRINT 122,NIMPAC,NCORE,NATMOS
122   FORMAT('0**** NMBR OF RAYS TOO LARGE =',I3,'  NCORE,NATMOS =',
     & 2I5)
      STOP
C
C FIND THE NUMBER OF MU-PNTS FOR EACH K, PLUS ONE EXTRA FOR MU=0.0
131   II=NIMPAC+1
      DO 130 K=1,JTAU1
         MMU(K)=II
         XMU(II,K)=0.0
         if(ii.ge.2) then
            IF (K+1.EQ.KIMPAC(II-1)) II=II-1
         endif
130   CONTINUE
C     PRINT 999,RR
999   FORMAT(' TRR ',5E11.4)
      write(0,*)' exit trrays '
      RETURN
      END

      SUBROUTINE TRANGS
       IMPLICIT REAL*8(A-H,O-Z)
c      PARAMETER (MD=350,MDP1=MD+1,MPAR=1000,MPP1=1001)
      PARAMETER (MPAR=1500,MPP1=1501,m2pp1=2*mpp1)
c     include 'PARAM'
      include 'parameters.h'
C
C TRANGS SOLVES THE TRANSFER EQUATION WITH GOERANS METHOD
C  FORMAL SOLVER FOR THE TRANSFER EQUATION USING DIRECT INTEGRATION
C  FOR CUBIC SPLINE SOURCES
C
C  ITRAN DETERMINES THE MODE OF FORMAL SOLUTION
C  ITRAN  = 3  INTEGRAL CUBIC SPLINE METHOD ACCORDING TO SCHARMER
C           4  INTEGRAL CUBIC SPLINE METHOD WITH LOCALLY DETERMINED
C              FIRST DERIVATIVE
C
C  REF G.B.SCHARMER
C
      REAL*8 SPRIM(MPAR),SBISS(MPAR),STRISS(MPAR),
     & EXPD(MPAR),WW(MPAR)
      INTEGER KP(MPAR)
      REAL*8 IP(0:MPAR),IM(0:MPAR),IPLUS(0:MPAR),IMINUS(0:MPAR)
      COMMON /CTRAN/X(MPAR),SI(MPAR),XJ(MPAR),XH(MPAR),XK(MPAR)
     & ,FJ(MPAR),S(MPAR),JTAU0,JTAU1
      COMMON /TAUCQ/TAU(MPAR),TT(MPP1,MPAR),DTAULN(MPAR),JTAU
      COMMON /SPACE2/P(MPAR)
     & ,SP1(MPAR),SP2(MPAR),SP3(MPAR),AD(MPAR),BD(MPAR)
     & ,NIMPAC,KIMPAC(MPP1),PIMPAC(MPP1),MMU(MPAR)
     & ,TAUT(MPAR),DTAUT(MPAR)
     & ,PFEAU(MPP1,MPAR),XMU(MPP1,MPAR),FUN(MPP1),DER(M2PP1),DMU(MPP1)
      COMMON /CSPHER/NCORE,NATMOS,TAUM,RR(MPAR)
      COMMON /TRDBUG/IDEBUG
      DIMENSION DTCNT(MPP1),VFE(MPP1,MPAR)
      DATA DT1/1.E-1/,DT2/1.E2/
C
C  THE PARAMETER DT1 DETERMINES WHEN I HAS TO BE CALCULATED IN A
C  THIN WAY DUE TO NUMERICAL REASONS. THE PARAMETER SHOULD BE SET
C  TO ABOUT THE THIRD ROOT OF THE MACHINE ACCURACY.
C
C
C MU LOOP
C     PRINT*,' I TRANGS ',NTAU
      DO 131 I=1,NIMPAC
C INNERMOST SHELL FOR IMPACT RAY I (FROM SURFACE)
      NTAU=KIMPAC(I)
      NTAU1=NTAU-1
C
C CALCULATE DTAUT ALONG THE RAY
      ZOLD=0.0
      DO 100 K=1,NTAU
      Z=DSQRT(RR(NTAU-K+1)**2-PIMPAC(I)**2)
      if(abs(rr(ntau-k+1)).le.1.e-33) 
     &         write(6,*)' overflow i trangs!!! ',k,ntau,z,rr(ntau-k+1)
      XMU(I,NTAU-K+1)=-Z/RR(NTAU-K+1)
      IF (K.EQ.1) GO TO 100
      DZ=Z-ZOLD
      DZDR=DZ/(RR(NTAU-K+1)-RR(NTAU-K+2))
      DTAUT(NTAU-K+2)=DZDR*0.5d0*(X(NTAU-K+1)+SI(NTAU-K+1)
     & +X(NTAU-K+2)+SI(NTAU-K+2))*(TAU(NTAU-K+2)-TAU(NTAU-K+1))
2626  FORMAT(I4,7E11.4)
100   ZOLD=Z
      TAUT(1)=DZDR*(X(1)+SI(1))*TAU(1)
      DO 101 K=2,NTAU
      TAUT(K)=TAUT(K-1)+DTAUT(K)
      TT(I,K)=TAUT(K)
101   CONTINUE
c      DO 888 K=1,NTAU
C     IF(I.EQ.1) PRINT 2626,K,RR(K),X(K),S(K),TAU(K),DTAUT(K)
c888   CONTINUE
C
C
C
      IMINUS(0)=0.0
      T=TAUT(1)
C
      ITRAN=4
      IF(ITRAN.EQ.3) THEN
        CALL SPLIN0(NTAU,DTAUT,S,SPRIM,SBISS,STRISS,WW)
      ELSE
        CALL SPLIN1(NTAU,DTAUT,S,SPRIM,SBISS,STRISS)
      ENDIF
C
C  CALCULATE ONLY NECESSARY EXPONENTIALS, AND DO IT IN A SEPARATE
C  LOOP TO ALLOW THE CRAY TO VECTORIZE THE CALLS
C
      NM=0
      DO 400 K=1,NTAU-1
        IF(DTAUT(K+1).GT.DT1 .AND. DTAUT(K+1).LE.DT2) THEN
          NM=NM+1
          WW(NM)=DTAUT(K+1)
          KP(NM)=K
        ENDIF
  400 CONTINUE
      DO 410 M=1,NM
        WW(M)=DEXP(-WW(M))
  410 CONTINUE
      DO 420 M=1,NM
        EXPD(KP(M))=WW(M)
  420 CONTINUE
C
C  CALCULATE INCOMING INTENSITY
C
      IF(T.LT.0.01) THEN
        IM(1)=IMINUS(0)+T*(1.-T*(.5-T*(1./6.-T/24.)))*
     &          (S(1)-IMINUS(0))
      ELSE IF(T.LT.20.) THEN
        IM(1)=IMINUS(0)+(1.-DEXP(-T))*(S(1)-IMINUS(0))
      ELSE
        IM(1)=S(1)
      ENDIF
      IMINUS(1)=IM(1)
      DO 600 K=1,NTAU-1
        IF(DTAUT(K+1).LE.DT1) THEN
          IM(K+1)=IM(K)-DTAUT(K+1)*(IM(K)-S(K)-
     &     .5*DTAUT(K+1)*(IM(K)-S(K)+SPRIM(K)-.333333333*
     &     DTAUT(K+1)*(IM(K)-S(K)+SPRIM(K)-SBISS(K)-.25*
     &     DTAUT(K+1)*(IM(K)+SPRIM(K)-SBISS(K)+STRISS(K)))))
        ELSE IF(DTAUT(K+1).LE.DT2) THEN
          SBISS1=SBISS(K)+DTAUT(K+1)*STRISS(K)
          IM(K+1)=S(K+1)-SPRIM(K+1)+SBISS1-STRISS(K)+EXPD(K)*
     &     (IM(K)-S(K)+SPRIM(K)-SBISS(K)+STRISS(K))
        ELSE
          SBISS1=SBISS(K)+DTAUT(K+1)*STRISS(K)
          IM(K+1)=S(K+1)-SPRIM(K+1)+SBISS1-STRISS(K)
        ENDIF
        IMINUS(K)=IM(K)
  600 CONTINUE
C
C     INNER BOUNDARY CONDITION I+=S(NTAU)+...
C
      IPLUS(NTAU)=S(NTAU)+SPRIM(NTAU)+SBISS(NTAU)+STRISS(NTAU)
C
C     INNER BOUNDARY CONDITION I+=I-
C
      IPLUS(NTAU)=IMINUS(NTAU)
C
C  CALCULATE OUTGOING INTENSITY IPLUS(K)
C
      IP(NTAU)=IPLUS(NTAU)
      DO 500 K=NTAU-1,1,-1
        IF(DTAUT(K+1).GT.DT2) THEN
          IP(K)=S(K)+SPRIM(K)+SBISS(K)+STRISS(K)
        ELSE IF (DTAUT(K+1).GT.DT1) THEN
          SBISS1=SBISS(K)+DTAUT(K+1)*STRISS(K)
          IP(K)=S(K)+SPRIM(K)+SBISS(K)+STRISS(K)+EXPD(K)*
     &     (IP(K+1)-S(K+1)-SPRIM(K+1)-SBISS1-STRISS(K))
        ELSE
          SBISS1=SBISS(K)+DTAUT(K+1)*STRISS(K)
          IP(K)=IP(K+1)-DTAUT(K+1)*(IP(K+1)-S(K+1)-
     &     0.5d0*DTAUT(K+1)*(IP(K+1)-S(K+1)-SPRIM(K+1)-
     &     0.333333333*DTAUT(K+1)*(IP(K+1)-S(K+1)-SPRIM(K+1)-
     &     SBISS1-0.25*DTAUT(K+1)*(IPLUS(K+1)-SPRIM(K+1)-SBISS1-
     &     STRISS(K)))))
        ENDIF
        IPLUS(K)=IP(K)
  500 CONTINUE
C
C  SURFACE INTENSITY
C
      IPLUS(0)=IPLUS(1)+T*(1.-T*(.5-T*(1./6.-T/24.)))*(S(1)-IPLUS(1))
C
C  FEAUTRIERS P
C
      DO 700 K=1,NTAU
        PFEAU(I,K)=0.5d0*(IPLUS(K)+IMINUS(K))
        VFE(I,K)=0.5d0*(IPLUS(K)-IMINUS(K))
      IF(VFE(I,K).LE.0.) VFE(I,K)=1.E-30
      IF(PFEAU(I,K).LE.0.) PFEAU(I,K)=1.E-30
C     PRINT 955,K,IPLUS(K),IMINUS(K),PFEAU(I,K),S(K),DTAUT(K)
      IF(PFEAU(I,K).LE.0.) PFEAU(I,K)=1.E-30
      IF(VFE(I,K).LE.0.) VFE(I,K)=1.E-30
955    FORMAT(I5,9E11.4)
  700 CONTINUE
C
C
C END MU LOOP
      DO 977 K=1,NTAU1
C     VFE(I,K)=(PFEAU(I,K+1)-PFEAU(I,K))/(TT(I,K+1)-TT(I,K))
      IF(VFE(I,K).LE.0.) VFE(I,K)=1.E-30
 977  CONTINUE
      IF(ICOU.LE.4) GOTO 398
      ICOU=0
C     PRINT*,I,PIMPAC(I)
      DO 999 K=1,NTAU1,2
      FMIN=PFEAU(I,K)-VFE(I,K)
      FPLU=PFEAU(I,K)+VFE(I,K)
C      PRINT 923,K,TT(I,K),PFEAU(I,K),VFE(I,K),FMIN,FPLU
 999  CONTINUE
 923  FORMAT(I4,5E11.4)
 398  ICOU=ICOU+1
131   CONTINUE
C
C INTERPOLATE TO PFEAU AT MU=0 FOR THOSE K THAT HAVE NO RAY
      DO 181 K=1,JTAU1
      II=MMU(K)
      IF (KIMPAC(II).EQ.K) GO TO 181
      PX=-XMU(II-2,K)/(XMU(II-1,K)-XMU(II-2,K))
      QX=1.-PX
      IF(VFE(II-2,K).LE.0.) VFE(II-2,K)=1.E-30
      IF(VFE(II-1,K).LE.0.) VFE(II-1,K)=1.E-30
      IF(PFEAU(II-2,K).LE.0.) PFEAU(II-2,K)=1.E-30
      IF(PFEAU(II-1,K).LE.0.) PFEAU(II-1,K)=1.E-30
      PFEAU(II,K)=DEXP(DLOG(PFEAU(II-2,K))*QX+DLOG(PFEAU(II-1,K))*PX)
      VFE(II,K)=DEXP(DLOG(VFE(II-2,K))*QX+DLOG(VFE(II-1,K))*PX)
      IF(VFE(II-2,K).LE.0..OR.VFE(II-1,K).LE.0.) PRINT*,
     &II,K,VFE(II-2,K),VFE(II-1,K)
      IF(PFEAU(II-2,K).LE.0..OR.PFEAU(II-1,K).LE.0.) PRINT*,
     &II,K,PFEAU(II-2,K),PFEAU(II-1,K)
181   CONTINUE
C
C CALCULATE MEAN INTENSITY
      DO 190 K=1,JTAU1
      XJ(K)=0.
      NMU=MMU(K)
      DO 191 I=2,NMU
      DMU(I)=XMU(I,K)-XMU(I-1,K)
      DER(I)=(PFEAU(I,K)-PFEAU(I-1,K))/DMU(I)
191   XJ(K)=XJ(K)+DMU(I)*(PFEAU(I,K)+PFEAU(I-1,K))
      XJ(K)=XJ(K)*6.
      NMU1=NMU-1
      DO 192 I=2,NMU1
192   XJ(K)=XJ(K)+(DMU(I+1)-DMU(I))*(DMU(I)*DER(I+1)+DMU(I+1)*DER(I))
      XJ(K)=XJ(K)+DMU(2)**2*DER(2)-DMU(NMU)**2*DER(NMU)
      XJ(K)=XJ(K)*0.083333333
190   CONTINUE
C CALCULATE FLUX
      DO 195 K=1,JTAU1
      XH(K)=0.
      NMU=MMU(K)
      DO 196 I=2,NMU
      DMU(I)=XMU(I,K)-XMU(I-1,K)
      DER(I)=(XMU(I,K)*VFE(I,K)-XMU(I-1,K)*VFE(I-1,K))/DMU(I)
196   XH(K)=XH(K)+DMU(I)*(XMU(I,K)*VFE(I,K)+XMU(I-1,K)*VFE(I-1,K))
      XH(K)=XH(K)*6.
      NMU1=NMU-1
      DO 197 I=2,NMU1
197   XH(K)=XH(K)+(DMU(I+1)-DMU(I))*(DMU(I)*DER(I+1)+DMU(I+1)*DER(I))
      XH(K)=XH(K)+DMU(2)**2*DER(2)-DMU(NMU)**2*DER(NMU)
      XH(K)=XH(K)*0.083333333
195   CONTINUE
      NCORP1=NCORE+1
      DO 295 K=1,JTAU1
      XH(K)=0.
      XJ(K)=0.
      NMU=MMU(K)
      DO 296 I=2,NMU
      DMU(I)=XMU(I,K)-XMU(I-1,K)
      IF(I.NE.NCORP1)XJ(K)=XJ(K)+DMU(I)*(PFEAU(I,K)+PFEAU(I-1,K))
      IF(I.NE.NCORP1)XH(K)=XH(K)+DMU(I)*(XMU(I,K)*VFE(I,K)+
     &XMU(I-1,K)*VFE(I-1,K))
C AVOID DISCONTINUITY AT THE LIMB FOR THE MOMENTS
      IF(I.EQ.NCORP1)XJ(K)=XJ(K)+DMU(I)*(PFEAU(I,K)+PFEAU(I,K))
      IF(I.EQ.NCORP1)XH(K)=XH(K)+DMU(I)*(XMU(I,K)*VFE(I,K)+
     &XMU(I-1,K)*VFE(I,K))
C     IF(K.GE.37) PRINT 828,I,DMU(I),VFE(I,K),PFEAU(I,K),XH(K),XJ(K)
828   FORMAT(I4,5E12.5)
296   CONTINUE
      XH(K)=XH(K)*0.5d0
      XJ(K)=XJ(K)*0.5d0
295   CONTINUE
C--------------------------------------------------------------------
C
C
C FLUX AT TAU(1)
      XH(1)=TRQUAD(MMU(1),XMU,FUN,DER)
C
C CALCULATE SECOND MOMENT XK
      DO 201 K=1,JTAU1
      NMU=MMU(K)
      DO 200 I=1,NMU
200   FUN(I)=PFEAU(I,K)*XMU(I,K)**2
201   XK(K)=TRQUAD(NMU,XMU(1,K),FUN,DER)
C
C CALCULATE FIRST MOMENT, XH, FROM MOMENT RELATION
      DO 211 K=2,JTAU
211   XH(K)=(XK(K)-XK(K-1)+(XJ(K)+XJ(K-1)-3.*(XK(K)+XK(K-1)))
     & *(RR(K-1)-RR(K))/(RR(K)+RR(K-1)))
     & *2.d0/((TAU(K)-TAU(K-1))*(X(K)+S(K)+X(K-1)+S(K-1)))
C
C CALCULATE FIRST MOMENT XH BY QUADRATURE
      DO 251 K=1,JTAU1
      NMU=MMU(K)
      DO 250 I=1,NMU
250   FUN(I)=PFEAU(I,K)*XMU(I,K)
251   XH(K)=TRQUAD(NMU,XMU(1,K),FUN,DER)
C
      RETURN
C------------------------------------------------------------------
C
C EMERGENCY EXIT
230   KK=NTAU-K
      IDEBUG=2
      PRINT 232,I,KK,JTAU0,JTAU1,NTAU
232   FORMAT('0NON-POSITIVE RESULT AT I,K,J0,J1,N=',5I3)
      GO TO 233
231   KK=0
      IDEBUG=3
      PRINT 232,I,KK,JTAU0,JTAU1,NTAU
233   PRINT 237,NCORE,EX
237   FORMAT(' NCORE,EX(I)=',I3,3G12.4)
      PRINT 238
238   FORMAT(' DEBUG INFORMATION WRITTEN ON UNIT 13')
      RETURN
      END

      SUBROUTINE SPLIN0(N,DX,F,D,D2,D3,WW)
C
C  CALCULATES STANDARD CUBIC SPLINE WITH CONTINOUS SECOND DERIVATIVE
C
C  CODED BY G.SCHARMER, 1982
C
C  DX(K)=X(K)-X(K-1)
C  WW IS A WORKING ARRAY
C
       IMPLICIT REAL*8(A-H,O-Z)
      REAL*8  DX(N),F(N)
      REAL*8 D(N),D2(N),D3(N),WW(N)
C
      FAC=-DX(2)/DX(3)
      D(2)=(2.d0-FAC)*(DX(2)+DX(3))
      C2=DX(3)+FAC*DX(2)
      WW(2)=(F(2)-F(1))/DX(2)
      WW(3)=(F(3)-F(2))/DX(3)
      D3(2)=6.*(WW(3)-WW(2))
      DO 100 K=3,N-1
        D(K)=2.d0*(DX(K)+DX(K+1))
        WW(K+1)=(F(K+1)-F(K))/DX(K+1)
        FAC=-DX(K)/D(K-1)
        D(K)=D(K)+FAC*DX(K)
        IF(K.EQ.3) D(K)=D(K)+FAC*(C2-DX(K))
        D3(K)=6.*(WW(K+1)-WW(K))+FAC*D3(K-1)
  100 CONTINUE
      FAC=-DX(N)/D(N-2)
      AN=-DX(N-1)-DX(N)+FAC*DX(N-1)
      D3(N)=FAC*D3(N-2)
      FAC=-AN/D(N-1)
      D2(N)=(D3(N)+FAC*D3(N-1))/(DX(N-1)+FAC*DX(N))
      DO 150 K=N-1,3,-1
        D2(K)=(D3(K)-DX(K+1)*D2(K+1))/D(K)
  150 CONTINUE
      D2(2)=(D3(2)-C2*D2(3))/D(2)
      D2(1)=((DX(2)+DX(3))*D2(2)-DX(2)*D2(3))/DX(3)
      DO 180 K=1,N-1
        D(K)=WW(K+1)-(D2(K+1)+D2(K)+D2(K))*DX(K+1)/6.
        D3(K)=(D2(K+1)-D2(K))/DX(K+1)
  180 CONTINUE
      D3(N)=D3(N-1)
      D(N)=D(N-1)+DX(N)*(D2(N-1)+0.5d0*DX(N)*D3(N-1))
C
      RETURN
      END

      SUBROUTINE SPLIN1(N,DX,F,D,D2,D3)
C**********************************************************************
C
C  CALCULATES CUBIC SPLINES WITH LOCALLY DETERMINED FIRST DERIVATIVE
C  STABILITY BETTER THAN FOR STANDARD CUBIC SPLINE, AT COST OF
C  DISCONTINOUS SECOND DERIVATIVE.
C
C  CODED BY AA.NORDLUND, MAR-83
C
C  DX(K)=X(K)-X(K-1)
C
       IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 DX(N),F(N)
      REAL*8 D(N),D2(N),D3(N)
C
C  FIRST DERIVATIVE BY CENTERED DIFFERENCE
C
      DO 100 K=2,N-1
        D(K)=(F(K+1)-F(K-1))/(DX(K+1)+DX(K))
  100 CONTINUE
      D(1)=(F(2)-F(1))/DX(2)
      D(N)=(F(N)-F(N-1))/DX(N)
C
C  SECOND AND THIRD DERIVATIVE FROM SPLINE CONDITIONS
C
      DO 110 K=1,N-1
        CX=1.0/DX(K+1)
        DFDX=(F(K+1)-F(K))*CX
        D2(K)=(6.*DFDX-4.*D(K)-2.d0*D(K+1))*CX
        D3(K)=6.*(D(K)+D(K+1)-2.d0*DFDX)*CX*CX
  110 CONTINUE
      D2(N)=(4.*D(N)+2.d0*D(N-1)-6.*DFDX)*CX
      D3(N)=D3(N-1)
C
      RETURN
      END
