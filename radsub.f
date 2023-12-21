      SUBROUTINE RLOSS(iel,ion,WL,OM,G1,G2,A21,R,XEL,Z,TE,RL,WEQ
     &     ,CINT,RFL,L,WLI)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 MTOT,MNI
      PARAMETER (NFEL=3000)
C     **************************************************************
C     *****
C     THIS ROUTINE CALCULATES THE HEAT LOSS AND OBSERVED FLUX FROM
C     A LINE WITH THE ESCAPE PROBABILITY METHOD.
C     WL = WAVELEGTH (A)
C     OM = COLLISION STRENGTH
C     G1,G2 = STAT.WEIGHTS OF THE LOWER AND UPPER LEVELS
C     A21 = TRANSITION PROBABILITY
C     R = RADIUS (IN TERMS OF SHOCK RADIUS)
C     XEL = ELECTRON FRACTION
C     Z = TOTAL ABUNDANCE OF THE ION
C     TE = TEMPERATURE
C     RL = LOSS RATE IN ERG CM**3 /S
C     WEQ = FLUX EMITTED TO INFINITY (ERG/S)
C     CINT = OBSERVED FLUX
C     RFL = RADIATION FORCE (NOT USED HERE)
C     SRED = (1-EXP(-TAU))*SOURCE FCN FOR LINE PROFILE
C     REVISED 87-05-15
C
C     *****
C     ****************************************************************
      COMMON/NBA/NBACK
      COMMON/A2/SIg,T0,TAUL
      parameter (nlp=30000)
      COMMON/COLD/COLTOT(14,27),COLTOTT(14,27),SRED(NFEL+nlp)
      COMMON/QSOM/qso,IAGN,ISTAT,ISOBOL
      COMMON/A7/C3,C33
      COMMON/A10/AV
      COMMON/SPH/ISPH
      COMMON/TPAR/RIN,DRQ,R1Q,TDAY
      COMMON/A11/R11,CV,FLUX
      COMMON/A5/TAU,ALFA,ENE,TSN,XL40,TXEV,RMAX
      COMMON/MPAR/MTOT,MNI,VEXP
      COMMON/HYDROS/HSCALE,DEN1,XQL,AMEAN
      common/lineopac/totopl(nlp)
      common/linephoto/phline(nlp)
      common/escapemax/be00
      common/pescpd/pdbe(nlp,2)
      common/contdest/idest
      common/velgrad/dr_dv
      dimension ami(100)
      data ami/1.,2*4.,4*16.,6*12.,
     &         4*16.,7*14.,14*28.,
     &         3*24.,4*56.,4*26.,3*40.,3*23.,16*32.,10*20.,
     &         2*36.,11*56.,6*0/
      dimension atw(15)
      data atw/1.,4.,12.,14.,16.,20.,22.,24.,26.,
     &     28.,32.,36.,40.,56.,58./
      DATA PI/3.1415926E0/,ELCH/1.60219E-12/,AMU/1.660531E-24/
      WLI=WL
      SIg=AV-1.
      EN=1.4413E8/WL
      ENCGS=1.989E-8/WL
      EEV=ENCGS/ELCH
      dtot=de(r)
      DENEL=XEL*dtot
      C21=8.63E-6*OM/(G2*SQRT(TE))
      E=C21*DENEL/A21

C
C     THIS PART APPLIES FOR NO VELOCITY GRADIENT AND NO BACKGROUND
C     CONTINUUM.
C
      IF(ISTAT.EQ.1) THEN
        T0=C33*COLTOT(iel,ion)*A21*G2*WL**3./G1
        t00=t0
        T0T=(COLTOTT(iel,ion)-COLTOT(iel,ion))*T0/COLTOTT(iel,ion)
        IF(T0T.LE.0.) T0T = 0.
      ENDIF
      W=0.E0
      nback=0
      IF(NBACK.EQ.1) THEN 
        FIN=FMEAN(3,EEV)
      ELSEIF(NBACK.EQ.0) THEN
        FIN=0.
      ENDIF
      IF(ISTAT.EQ.1) GOTO 200
C
C     THIS PART APPLIES FOR A VELOCITY GRADIENT AND A BACKGROUND
C     CONTINUUM, FC.
C
      ENX=EN/TSN
      IF(ISPH.EQ.1) W=.5*(1.-SQRT(1.-1./(R*R)))
      w=wd(r)
      TIME=8.64E4*TDAY
C     OPT. DEPTH. COEFF = 1E-24/8 PI
      T00=3.979E-26*TIME*WL**3.*G2*A21*Z*dtot/G1

c vel grad

      T00=3.979E-26*dr_dv*WL**3.*G2*A21*Z*dtot/G1

      T0T=0.
 200  CONTINUE
C
C     CALCULATE THE ESCAPE PROBABILITY.
C
      MM=0
      T0=T00
232   BE=1.
      TOP=T0
      ENTE=EN/TE
      AU=atw(iel)
      VTERM=1.285E6*SQRT(TE/(AU*1.E4))
      CALL ESCAPE(TOP,T0T,WL,A21,VTERM,BE,DBEDTA)
      IF(EN/TE.GT.100.) ENTE=100.
C
C     N1*G2/N2*G1
C
c     line center opacity
      opline=3.979E-26*wl**3.*g2*a21*Z/(g1*vterm)
c      opline=top/(vterm*time*dtot)
c     continuum desruction probability
      opcon=totopl(l)
c     factor 4.9 is cecilias estimate. may be wrong by factor 2 from exact
c     H-R expression
      fhumm=8.5
      betadest=opcon/opline
      fbeta=min(fhumm*betadest,1.d-3)
      pd=fbeta/(fbeta+1.)
c  if idest = 0 supress continuum destruction
c!!!! neglect cont. dest.
      idest=0

      if(idest.eq.0) then
         pd=0.
      endif
c
c      pd=fhumm*opcon/(fhumm*opcon+opline)
      RN1G2N2G1=(BE*(1.-pd)+pd*(be00-be)+E)/(E*EXP(-ENTE))
      RN1N2=G1*RN1G2N2G1/G2
      rn2=1./(1.+rn1n2)
C
C     CORRECT OPTICAL DEPTH FOR STIMULATED EMISSION
C
      IF(RN1N2.LT.10.) T0NEW=T00*(1.-1./RN1G2N2G1)/(1.+1./RN1N2)
      IF(RN1N2.LT.10.) T0=T0NEW
      IF(RN1N2.LT.10.) MM=MM+1
      IF(RN1N2.LT.10..AND.MM.LT.4) GOTO 232
C
C     COOLING RATE OF ELECTRONS TAKING THERMALIZATION INTO ACCOUNT
C
      RL=ENCGS*G2*C21*Z*EXP(-ENTE)*(1.-EXP(ENTE)/RN1G2N2G1)
     &/(G1*(1.+1./RN1N2))
      XX=G1*(A21*BE+C21*DENEL)/(G2*C21*DENEL*EXP(-ENTE))
      WEQ=Z*A21*ENCGS*BE*dtot*((1.-W))/(1.+XX)

      if(iel.eq.14.and.ion.ge.10) then

c         write(6,9286)om,g2,z,be,a21*be/(c21*denel),xx,w,dtot,weq
 9286    format(' rloss a ',1pe12.3,10e12.3)
      endif

c     photoionization due to line l
c     the actual rate to ion k is obtained from phline(l)*sk(k,l)
c     and the total phline(l)*opcon
c
c      phline(l)=fhumm*rn2*z*a21*(be00-be)/(opline+fhumm*opcon)
      if(opcon.gt.0.) then
         phline(l)=(pd/opcon)*rn2*z*a21*(be00-be)
         if(phline(l).lt.0.) then
            write(6,*)' phline lt 0 i radsub',l,
     &           pd,rn2,z,be00,be,phline(l)
            write(6,*)TOP,T0T,WL,A21,VTERM,BE,DBEDTA
         endif
      else
         phline(l)=0.
      endif
      if(l.ge.85.and.l.le.90) then
c         write(6,922)l,opline,opcon,rn2,t0,be,pd,phline(l)
      endif
 922  format(' opdest ',i5,1pe12.4,10e12.4)
      IF(EN/TE.LT.100.) BPL=2.031E-04*EEV**3/(EXP(EN/TE)-1.)
C     SOURCE FUNCTION INCL. STIM. EMISSION
      EP=E
      IF(EN/TE.LT.100.) EP=E*(1.-EXP(-EN/TE))
      IF(EN/TE.LT..001) EP=E*EN/TE
      SOURCE=(EP*BPL)/(EP+BE)
C     OBSERVED INTENSITY
      EXPM1=1.-EXP(-T0)
      IF(T0.LT.1.E-3) EXPM1=T0-0.5*T0**2
      IF(ISTAT.NE.1) THEN
        CINT=4.*PI*1.E8*EXPM1*SOURCE/(TIME*WL*DENEL*dtot)

c vel grad

        CINT=4.*PI*1.E8*EXPM1*SOURCE/(dr_dv*WL*DENEL*dtot)

      ELSE
        CINT=WEQ/(DENEL*dtot)
      ENDIF

      if(iel.eq.14.and.ion.ge.10) then

c         write(6,9287)istat,ion,wl,c21,t0,expm1,cint
 9287    format(' rloss ',2i5,f12.2,1pe12.3,10e12.3)
      endif

C     CONVERT TO INTENSITY/ ANGSTROM
      SOURCE=3.E18*SOURCE/WL**2
      SRED(L)=EXPM1*SOURCE
c      write(6,98)K,WL,L,WLI
 98   format(' rlqq ',i7,f10.1,i6,f10.1)
      pdbe(l,1)=be
      pdbe(l,2)=pd
      RETURN
      END




      SUBROUTINE DIEL(A,T0,B,T1,TE,ALD)                                 
C
      IMPLICIT REAL*8(A-H,O-Z)                                          
C     *******************************************************
C     DIELECTRIC RECOMBINATION RATE ACCORDING TO FORMULA IN 
C     ALD. AND PEQ.
C     *******************************************************
      EX=T0/TE                                                          
      ALD=0.                                                            
      IF(EX.LT.100)ALD=A*DEXP(-T0/TE)*(1.+B*DEXP(-T1/TE))/              
     &(TE*DSQRT(TE))                                                    
C
      RETURN                                                            
      END                                                               
C************************************************************           
C************************************************************           
      SUBROUTINE DIELB(A,B,C,D,F,T4,ALDB)                               
C
      IMPLICIT REAL*8(A-H,O-Z)                                          
C     *******************************************************
C     DIELECTRIC RECOMBINATION RATE ACCORDING TO FORMULA IN 
C     NUSSBAUMER AND STOREY (A&A 126,75 (1983))            
C     *******************************************************
      
      EX=F/T4                                                           
      ALDB=0.                                                           
      IF(T4.GT.6.) GOTO 100                                             
      IF(EX.GT.100.) GOTO 100                                           
      ALDB=1.D-12*(A/T4+B+C*T4+D*T4*T4)*DEXP(-F/T4)/T4**1.5             
 100  CONTINUE                                                          
      RETURN                                                            
      END                                                               


      real*8 function dielsupp(iz,te,eden)
c     iz=charge
      IMPLICIT REAL*8 (A-H,O-Z)
      common/isupp/isupp
*     dielectronic burgess recombination suppresion
*     default is this to be true
      if(isupp.eq.1 ) then
*         following is effective density for scaling from Davidson's plot
*         first do temperature scaling, term in () is sqrt(te/15,000)
          effden =  eden / (sqrt(te) / 122.47)
*         this is rough charge dependence, z^7 from Davidson
          effden = effden / ((float(iz)/3.)**7)
          dielsupp = min(1.d0,  (1.-0.092*log10(effden) )  )
          dielsupp = max( 0.08d0, dielsupp )
      else
          dielsupp = 1.
      endif
      return
      end



      SUBROUTINE FIVELEV_dp(KI,iel,ion,N,E1,E2,E3,E4,E5,G1,G2,
     &G3,G4,G5,A21,A31,A41,A51,A32,A42,A52,A43,A53,
     &     A54,C21,C31,C41,C51,C32,C42,C52,C43,C53,C54,TE,Z,XI)
      IMPLICIT REAL*8 (A-H,O-Z)
c      REAL*4 A21,A31,A41,A51,A32,A42,A52,A43,A53,A54
c     &,C21,C31,C41,C51,C32,C42,C52,C43,C53,C54
c     &,E1,E2,E3,E4,E5,G1,G2,G3,G4,G5
      INTEGER N
      include "parameters.h"
      parameter (nlp=30000)
      COMMON/COLD/COLTOT(14,27),COLTOTT(14,27),SRED(NFEL+nlp)
c      COMMON/COLD/OPT(nio),COLTOT(nio),COLTOTT(nio),SRED(NFEL+100)
      COMMON/QSOM/qso,IAGN,ISTAT,ISOBOL
      COMMON/A7/C3,C33
      COMMON/PHY/DENS(MD)
      COMMON/ELEC/DEL(MD)
      COMMON/TPAR/RIN,DRQ,R1Q,TDAY
      COMMON/IND/IK
      COMMON/A19/EM(NL,NL),ESC(NL,NL),TTOT(NL,NL),TOP(NL,NL)
      COMMON/EQUIH/FRH(2,401),WEQ(401),SR(NFEL+nlp),WOBS(NL,NL)
      common/w/wlin(401),wl(NL,NL),taul(401)
      common/kmaxpop/kmaxp(14,27)
      DIMENSION  C(5,5),G(5),E(5),A(5,5),BE(5,5),XI(15)
      DIMENSION AA(nlp1,nlp1),X(nl),RL(5,5),xold(nl),topt(5,5)
      DATA PI/3.1415926E0/,ELCH/1.60219E-12/,AMU/1.660531E-24/
      DEN=DENS(IK)
      DENEL=DEN*DEL(IK)
      NM1=N-1
      NP1=N+1
      DO I=1,NP1
         DO J=1,NP1
            AA(I,J)=0.
         enddo
      enddo
      E(1)=DBLE(E1)
      E(2)=DBLE(E2)
      E(3)=DBLE(E3)
      E(4)=DBLE(E4)
      E(5)=DBLE(E5)
      G(1)=DBLE(G1)
      G(2)=DBLE(G2)
      G(3)=DBLE(G3)
      G(4)=DBLE(G4)
      G(5)=DBLE(G5)
      DO I=1,N
         E(I)=E(I)/8065.46
      enddo
      DO 5395 I=1,N
      DO 5394 J=1,N
      WL(I,J)=0.0
      IF(I.EQ.J) GOTO 5394
      IF(E(I).EQ.E(J)) GOTO 5394
      WL(I,J)=DABS(12398.54/(E(I)-E(J)))
5394  CONTINUE
5395  CONTINUE
      DO I=1,N
         DO  J=1,N
            C(I,J)=0.d0
            A(I,J)=0.d0
            be(I,J)=0.d0
         enddo
      enddo
      A(2,1)=DBLE(A21)
      A(3,1)=DBLE(A31)
      A(4,1)=DBLE(A41)
      A(5,1)=DBLE(A51)
      A(3,2)=DBLE(A32)
      A(4,2)=DBLE(A42)
      A(5,2)=DBLE(A52)
      A(4,3)=DBLE(A43)
      A(5,3)=DBLE(A53)
      A(5,4)=DBLE(A54)
      C(2,1)=DBLE(C21)
      C(3,1)=DBLE(C31)
      C(4,1)=DBLE(C41)
      C(5,1)=DBLE(C51)
      C(3,2)=DBLE(C32)
      C(4,2)=DBLE(C42)
      C(5,2)=DBLE(C52)
      C(4,3)=DBLE(C43)
      C(5,3)=DBLE(C53)
      C(5,4)=DBLE(C54)
      TEV=TE/1.1609D4
      TIME=8.64E4*TDAY
      DO I=1,NM1
         IP1=I+1
         DO J=IP1,N
            C(J,I)=8.63D-6*DENEL*C(J,I)/(DSQRT(TE)*G(J))
            C(I,J)=C(J,I)*DEXP(-(E(J)-E(I))/TEV)*G(J)/G(I)
         enddo
      enddo
      do k=1,5
         if(k==1) then
            x(1)=1.
         else
            x(k)=0.
         endif
         xold(k)=x(k)
      enddo
      DO L=1,20
      DO 100 I=1,N
      DO 100 J=1,N
      AA(I,J)=0.
      IF(I.EQ.N) AA(I,J)=1.D0
      IF(I.EQ.N) GOTO 100
      DEI=X(I)*Z*DEN
C
C     CALCULATE THE ESCAPE PROBABILITY.
C
      IF(ISTAT.EQ.1) THEN
        T0=C33*COLTOT(iel,ion)*A(J,I)*G(J)*WL(J,I)**3./G(I)
        T0T=(COLTOTT(iel,ion)-COLTOT(iel,ion))*T0/COLTOTT(iel,ion)
        IF(T0T.LE.0.) T0T = 0.
      ELSE
        T0=1.e-24*WL(J,I)**3*A(J,I)*G(J)*DEI*TIME/(8.*PI*G(I))
        T0T=0.
      ENDIF
      IF(L.EQ.1) T0=.0
C!!   ASSUME MEAN ATOMIC WEIGHT = 20. 
      AU=20.
      VTERM=1.285E6*SQRT(TE/(AU*1.E4))
      CALL ESCAPE(T0,T0T,WL(J,I),A(J,I),VTERM,BE(J,I),DBEDTA)
      topt(j,i)=t0
      IF(I.NE.J) GOTO 2
      S=0.d0
      DO K=1,N
         if(k.ne.i) then
            S=S+C(I,K)+BE(I,K)*A(I,K)
         endif
      enddo
      AA(I,I)=-S
      GOTO 100
2     AA(I,J)=BE(J,I)*A(J,I)+C(J,I)
 100  CONTINUE
      NRC=N+1
      AA(N,NRC)=1.
      EPS=1.D-30
      NA=N
      DI=SIMUL(NA,AA,X,EPS,1,NRC)
      errmax=0.
      do k=1,n
         deltax=abs(x(k)-xold(k))/x(k)
         if(deltax>errmax) then
            errmax=deltax
         endif
      enddo
      if(errmax < 1.e-5) goto 44
      do k=1,5
         xold(k)=x(k)
      enddo
      ENDDO
 44   continue
      DO I=1,N
        XI(I)=X(I)
      ENDDO
      K=0
      DO I=2,N         
         DO J=1,I-1
            WOBS(I,J)=1.602E-12*Z*(E(I)-E(J))*X(I)*A(I,J)*BE(I,J)/DEN
            EM(I,J)=WOBS(I,J)
            IF(K.le.400.and.a(i,j).gt.0.) then
               K=K+1
               SR(K)=0.
               WEQ(K)=WOBS(I,J)
               wlin(k)=wl(i,j)
               taul(k)=topt(i,j)
               kmaxp(iel,ion)=k
            endif
         enddo
      
         kmaxp(iel,ion)= k
      enddo
      RETURN
      END


      SUBROUTINE FORB3(IJ,iel,ion,E21,E32,A21,A31,A32,O21,O31,O32,S1,
     &     S2,S3,F21,F31,F32,TE,Z,F)

C     COOLING/VOL= F*X(EL)*AB(I)*X(I)*DEN**2. -> cin = weq 

      IMPLICIT REAL*8(A-H,O-Z)
c      PARAMETER (MD=350,MDP1=MD+1)
      include 'PARAM'
      COMMON/TAUFBO/TAFB(30,3)
      COMMON/PHY/DENS(MD)
      COMMON/ELEC/DEL(MD)
      COMMON/TPAR/RIN,DRQ,R1Q,TDAY
      COMMON/IND/I
      common/wlforb3/wl21,wl31,wl32
      COMMON/QSOM/qso,IAGN,ISTAT,ISOBOL
      PARAMETER (NL=340,NLP1=NL+1)
      common/w/wlin(401),wl(NL,NL),taul(401)
      parameter (nlp=30000)
      PARAMETER (NFEL=3000)
      COMMON/EQUIH/FRH(2,401),WEq(401),SR(NFEL+nlp),WOBS(NL,NL)
      DATA PI/3.1415926E0/,ELCH/1.60219E-12/,AMU/1.660531E-24/
      DEN=DENS(I)
      XEL=DEL(I)
      E31=E21+E32
      C21=8.63E-6*DEN*XEL*O21/(S2*SQRT(TE))
      C31=8.63E-6*DEN*XEL*O31/(S3*SQRT(TE))
      C32=8.63E-6*DEN*XEL*O32/(S3*SQRT(TE))
      TEV=TE/1.1609E4
      C12=0.
      C13=0.
      C23=0.
      IF(E21/TEV.LT.100.) C12=S2*C21*EXP(-E21*1.1609E4/TE)/S1
      IF(E31/TEV.LT.100.) C13=S3*C31*EXP(-E31*1.1609E4/TE)/S1
      IF(E32/TEV.LT.100.) C23=S3*C32*EXP(-E32*1.1609E4/TE)/S2
      WL21=12398.*1.E-8/E21
      WL31=12398.*1.E-8/E31
      WL32=12398.*1.E-8/E32
      wlin(1)=wl21*1.e8
      wlin(2)=wl31*1.e8
      wlin(3)=wl32*1.e8
      DENI=Z*DEns(i)
      DN1=1.
      DN2=1.e-10
      DN3=1.e-10
      mw=0
98    DE1=DN1*DENI
      DE2=DN2*DENI
      DE3=DN3*DENI
      mw=mw+1
      DN1O=DN1
      DN2O=DN2
      DN3O=DN3
      TIME=8.64E4*TDAY
      CALL ESCFOR(iel,ion,WL21,A21,S1,S2,DE1,TIME,T0,BE21)
      t21=t0
      CALL ESCFOR(iel,ion,WL31,A31,S1,S3,DE1,TIME,T0,BE31)
      t31=t0
      CALL ESCFOR(iel,ion,WL32,A32,S2,S3,DE2,TIME,T0,BE32)
      t32=t0
c!!! cloudy!! escape both inward and outwards for forb. lines if static!
      if(istat.eq.1) then
         be21=2.*be21
         be32=2.*be32
         be31=2.*be31
      endif
      Q1=(C12+C13)*(A32*BE32+C32)+C12*(A31*BE31+C31)
      Q2=(A21*BE21+C21)*(A32*BE32+C32)+(A21*BE21+C21+C23)*(A31*BE31+C31)
      DN21=Q1/Q2
      DN31=(C12+C13-DN21*(A21*BE21+C21))/(A31*BE31+C31)
      DN1=1./(1.+DN21+DN31)
      DN2=DN1*DN21
      DN3=DN1*DN31
      if(dn1o.le.0.and.mw.gt.20)WRITE(6,92)mw,DN1,DN2,DN3
      if(dn2o.le.0..and.mw.gt.20)WRITE(6,92)mw,DN1,DN2,DN3
      if(dn3o.le.0..and.mw.gt.20)WRITE(6,92)mw,wl31,DN1,DN2,DN3
      dr1=0.
      dr2=0.
      dr3=0.
      if(dn1o.gt.0.)DR1=ABS(DN1O-DN1)/DN1O
      if(dn2o.gt.0.)DR2=ABS(DN2O-DN2)/DN2O
      if(dn3o.gt.0.)DR3=ABS(DN3O-DN3)/DN3O
      DNMAX=MAX(DR1,DR2,DR3)
92    FORMAT(' DN ',i5,6E11.4)
      if(mw.gt.20)write(6,*)' no conv in forb'
      if(mw.gt.20) goto 231
      IF(DNMAX.GT.0.1) GOTO 98
231   IF(E21/TEV.LT.100.) B21=S1*EXP(E21/TEV)*DN21/S2
      IF(E31/TEV.LT.100.) B31=S1*EXP(E31/TEV)*DN31/S3
      F21=ELCH*DN1*DN21*A21*BE21*E21/(XEL*DEN)
      F31=ELCH*DN1*DN31*A31*BE31*E31/(XEL*DEN)
      F32=ELCH*DN1*DN31*A32*BE32*E32/(XEL*DEN)
      weq(1)=z*f21
      weq(2)=z*f31
      weq(3)=z*f32
      F=F21+F31+F32
      tafb(ij,1)=t21
      tafb(ij,2)=t31
      tafb(ij,3)=t32
      taul(1)=t21
      taul(2)=t31
      taul(3)=t32
      if(iel.eq.14.and.ion.eq.12) then
         write(6,981)o21,o31,o32
         write(6,981)a21,a31,a32,be21,be31,be32
         write(6,981)c21,c31,c32,den,xel,z
         write(6,981)dn1,dn21*dn1,dn31*dn1
         write(6,981)(weq(k),k=1,3)
981     format(' Fe rloss F ',1pe12.3,10e12.3)
      endif
      RETURN
      END

      SUBROUTINE ESCFOR(iel,ion,WL,A21,G1,G2,DE1,TIME,T0,BE)
C
C     CALCULATE THE ESCAPE PROBABILITY.
C
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NFEL=3000)
      parameter (nlp=30000)
      COMMON/COLD/COLTOT(14,27),COLTOTT(14,27),SRED(NFEL+nlp)
      COMMON/QSOM/qso,IAGN,ISTAT,ISOBOL
      COMMON/A7/C3,C33
      common/velgrad/dr_dv
      DATA PI/3.1415926E0/,ELCH/1.60219E-12/,AMU/1.660531E-24/
      BE=1.
      IF(ISTAT.EQ.1) THEN
        WL21=WL/1.E-8
        T0=C33*COLTOT(iel,ion)*A21*G2*WL21**3./G1
        T0T=(COLTOTT(iel,ion)-COLTOT(iel,ion))*T0/COLTOTT(iel,ion)
        IF(T0T.LE.0.) T0T = 0.
      ELSE
        T0=WL**3*A21*G2*DE1*TIME/(8.*PI*G1)

c vel grad

        T0=WL**3*A21*G2*DE1*dr_dv/(8.*PI*G1)

        T0T=0.
      ENDIF
C!!   ASSUME MEAN ATOMIC WEIGHT = 20. 
      AU=20.
c!!
      te=1.e4
      VTERM=1.285E6*SQRT(TE/(AU*1.E4))
      CALL ESCAPE(T0,T0T,WL,A21,VTERM,BE,DBEDTA)
      RETURN
      END



      SUBROUTINE PSI(IZ,IA,N,ZS,PS)
      IMPLICIT REAL*8(A-H,O-Z)
C     IZ=LOWER INDEX Z
C     IA=UPPER INDEX (1...5)
C     N= C OF EL. LEFT
C     ZS=IONIZATION RATES
C     PS=PSI
      DIMENSION ZS(18,4)
C     ZS(IZ,I), I=1 = K-SHELL, 2= 2S, 3= 2P, 4= 3S.
      PS=0.
      IF(IA.NE.1) GOTO 100
      PS=ZS(IZ,4)
      IF(N.GT.10) GOTO 110
      PS=PS+ZS(IZ,2)+ZS(IZ,3)
  110 IF(N.NE.11) GOTO 120
      PS=PS+ZS(IZ,3)
  120 GOTO 200
  100 IF(IA.NE.2) GOTO 210
      IF(N.LT.12) GOTO 130
      PS=ZS(IZ,3)
  130 IF(N.NE.11) GOTO 140
      PS=PS+ZS(IZ,2)
  140 IF(N.NE.12) GOTO 150
      PS=PS+ZS(IZ,2)
  150 IF(N.GT.10) GOTO 160
      PS=PS+ZS(IZ,1)
  160 IF(N.NE.11) GOTO 170
      PS=PS+2.*ZS(IZ,1)/3.
  170 GOTO 200
  210 IF(IA.NE.3) GOTO 220
      IF(N.LT.13) GOTO 310
      PS=PS+ZS(IZ,2)
  310 IF(N.NE.13) GOTO 320
      PS=PS+.666667*ZS(IZ,1)
  320 IF(N.NE.12) GOTO 330
      PS=PS+ZS(IZ,1)
  330 IF(N.NE.11) GOTO 340
      PS=PS+.3333333*ZS(IZ,1)
  340   GOTO 200
  220 IF(IA.NE.4) GOTO 230
      IF(N.NE.13) GOTO 350
      PS=PS+ZS(IZ,1)/3.
  350 IF(N.NE.14) GOTO 360
      PS=PS+ZS(IZ,1)
  360 IF(N.LT.15) GOTO 370
      PS=PS+2.*ZS(IZ,1)/3.
  370 GOTO 200
  230 IF(IA.NE.5) GOTO 200
      IF(N.LT.15) GOTO 200
      PS=ZS(IZ,1)/3.
  200 CONTINUE
      RETURN
      END

      SUBROUTINE COSI
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/COSIL/COS(18)
C     COLLISIONAL IONIZATION RATES FOR SILICON FROM MEWE ET AL
      CALL COLLI(2.E0,0.E0,8.25E0,C1)
      CALL COLLI(2.E0,0.E0,14.7E0,C2)
      COS(1)=C1+C2
      CALL COLLI(1.E0,1.E0,16.3E0,C1)
      CALL COLLI(2.E0,1.E0,22.8E0,C2)
      COS(2)=C1+C2
      CALL COLLI(2.E0,2.E0,33.5E0,COS(3))
      CALL COLLI(1.E0,3.E0,45.1E0,COS(4))
      CALL COLLI(6.E0,4.E0,167.E0,C1)
      CALL COLLI(2.E0,4.E0,220.E0,C2)
      COS(5)=C1+C2
      CALL COLLI(5.E0,5.E0,205.E0,C1)
      CALL COLLI(2.E0,5.E0,256.E0,C2)
      COS(6)=C1+C2
      CALL COLLI(4.E0,6.E0,246.E0,C1)
      CALL COLLI(2.E0,6.E0,294.E0,C2)
      COS(7)=C1+C2
      CALL COLLI(3.E0,7.E0,303.E0,C1)
      CALL COLLI(2.E0,7.E0,335.E0,C2)
      COS(8)=C1+C2
      CALL COLLI(2.E0,8.E0,351.E0,C1)
      CALL COLLI(2.E0,8.E0,378.E0,C2)
      COS(9)=C1+C2
      CALL COLLI(1.E0,9.E0,401.E0,C1)
      CALL COLLI(2.E0,9.E0,423.E0,C2)
      COS(10)=C1+C2
      CALL COLLI(2.E0,10.E0,476.E0,COS(11))
      CALL COLLI(1.E0,11.E0,523.E0,COS(12))
      CALL COLLI(2.E0,12.E0,2430.E0,COS(13))
      CALL COLLI(1.E0,13.E0,2670.E0,COS(14))
      RETURN
      END

      SUBROUTINE COLLI(S,Z,E0,CO)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*4 S,Z,E0
      COMMON/TEQQ/TE
      X=1.578E5*(E0/13.6)/TE
      CO=0.
      F=3.1-1.2/(Z+1.)-.9/(Z+1)**2
      IF(X.LT.100.)CO=2.223E-8*S*F*SQRT(TE)*EXP(-X)*(1.-EXP(-X))
     &/E0**2
      RETURN
      END



      SUBROUTINE COMG
      IMPLICIT REAL*8(A-H,O-Z)
C      ***********************************************************
C      ******
C      COLLISIONAL IONIZATION FROM SHULL & STEENBERG (APJ SUPP 48:95)
C      *******
C      ************************************************************
      COMMON/COM/CMG(3)
      COMMON/TEQQ/TE
      T4=TE/1.E4
      CMG(1)=8.90E-11*SQRT(TE)*EXP(-8.87/T4)/(1.+0.1*T4/8.87)
      CMG(2)=5.90E-11*SQRT(TE)*EXP(-17.4/T4)/(1.+0.1*T4/17.4)
      CMG(3)=0.
      IF(T4.LT.2) GOTO 10
      CMG(3)=1.10E-11*SQRT(TE)*EXP(-93.0/T4)/(1.+0.1*T4/93.0)
10    CONTINUE
      RETURN
      END



      SUBROUTINE COAL
      IMPLICIT REAL*8(A-H,O-Z)
C      ***********************************************************
C      ******
C      COLLISIONAL IONIZATION
C      *******
C      ************************************************************
      COMMON/COLAL/CAL(4)
      COMMON/TEQQ/TE
      T4=TE/1.E4
C
C     RATES ESTIMATED BY FORMULA (5-79) IN MIHALAS AND PHOTO-
C     IONIZATION CROSSECTIONS FROM REILMAN AND MANISON
C
      CAL(1)=2.89E-10*SQRT(TE)*EXP(-6.95/T4)
      CAL(2)=3.14E-12*SQRT(TE)*EXP(-21.83/T4)
      CAL(3)=3.30E-12*SQRT(TE)*EXP(-32.97/T4)
      CAL(4)=0.
      IF(T4.GT.1.) CAL(4)=1.61E-11*SQRT(TE)*EXP(-139.3/T4)
      RETURN
      END

      SUBROUTINE COSU
      IMPLICIT REAL*8(A-H,O-Z)
C  ***************************************************************
C     SULPHUR COLL.ION. RATES.
C     (Shull and van Steenberg).
C
C      ************************************************************
      COMMON/COSUL/CSU(18)
      COMMON/TEQQ/TE
      T4=TE/1.E4
C
      do i=1,10
         csu(i)=0.
      enddo
      if(te.lt.5.d4) goto 7406
      CSU(10)=4.75D-13*DSQRT(TE)*DEXP(-519./T4)/(1.+0.1*T4/519.)
      CSU(9)=6.97D-13*DSQRT(TE)*DEXP(-440./T4)/(1.+0.1*T4/440.)
      CSU(8)=1.01D-12*DSQRT(TE)*DEXP(-381./T4)/(1.+0.1*T4/381.)
      CSU(7)=1.71D-12*DSQRT(TE)*DEXP(-326./T4)/(1.+0.1*T4/326.)
 7406 continue
      if(te.lt.2.d4) goto 7407
      CSU(6)=1.85D-12*DSQRT(TE)*DEXP(-102./T4)/(1.+0.1*T4/102.)
      CSU(5)=6.43D-12*DSQRT(TE)*DEXP(-84.7/T4)/(1.+0.1*T4/84.7)
      CSU(4)=6.39D-12*DSQRT(TE)*DEXP(-54.9/T4)/(1.+0.1*T4/54.9)
 7407 continue
      if(te.lt.2.d3) goto 7408
      CSU(3)=2.12D-11*DSQRT(TE)*DEXP(-40.6/T4)/(1.+0.1*T4/40.6)
      CSU(2)=7.11D-11*DSQRT(TE)*DEXP(-27.1/T4)/(1.+0.1*T4/27.1)
      CSU(1)=1.45D-10*DSQRT(TE)*DEXP(-12.0/T4)/(1.+0.1*T4/12.0)
 7408 continue
C
      RETURN
      END



