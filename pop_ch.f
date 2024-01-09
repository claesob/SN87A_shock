c     ion = ionel: 1 = Ca II, 2 = O I, 3 = He I, 4 = Fe II, 5 = H, 6 = Fe I,
c     7 = He II ??               
      SUBROUTINE POPchianti(iel,ion,te,XQ,rltot)
      IMPLICIT REAL*8(A-H,O-Z)
c      PARAMETER (MD=350,MDP1=MD+1)
      include 'param'
      PARAMETER (NL=340,NLP1=NL+1)
      PARAMETER(NFEL=3000)
      parameter (nlp=30000)
      common/abl/abn(15)
      common/ionx/xion(md,14,27)
      common/ichia/ipopch
      COMMON/INUT/IUY
      COMMON/NION/IONQ
      COMMON/INI/INITH,INIT,INITCA,INITHE,INITFE,initfei
      COMMON/IND/IK
      COMMON/ELEC/DEL(MD)
      COMMON/NLEV/e00,NION,N,NP1H,NMAX,NMIN
      COMMON/BOLD/TOLD,TOLDH,TOLDCA,TOLDO,TOLDHE,TOLDFE,toldfei,BB(NL),
     &       BOLH(NL),BOL(NL),BOLCA(NL),BOLHE(NL),BOLFE(NL),bolfei(nl)
      COMMON/RECAL/RECO(NL)
      COMMON/HPOP/XN(6,NL),XN1,XN2,XN3
      COMMON/A14/C(NL,NL),CI(NL),G(NL),E(NL),AS(NL,NL),
     &     WL(NL,NL),DCDT(NL,NL),DCIDT(NL)
      COMMON/A17/RECNET(NL),RECNDT(NL),RECCO(NL),RTE(NL),
     &     PH(NL),DRDT(NL),RECT(6,NL),PHET(NL)
      COMMON/EQUIH/FRH(2,401),WEH(401),SR(NFEL+nlp),WOBS(NL,NL)
      common/timecheck/time,itime
      DIMENSION XINC(NL),XI(NL),
     &     A(NLP1,NLP1),BS(NL)

      DATA CSAHA/2.0708d-16/

      IONQ=ion

      abxa = abn(iel)*xion(ik,iel,ion)
      ipopch=0
      if(itime.le.3.or.abxa.gt.1.e-8) then
         ipopch=1
         if(iel.eq.3.and.ion.eq.5) then
            N=49
         elseif(iel.eq.3.and.ion.eq.6) then
            N=25
         elseif(iel.eq.4.and.ion.eq.6) then
            N=49
         elseif(iel.eq.4.and.ion.eq.7) then
            N=25
         elseif(iel.eq.5.and.ion.eq.7) then
            N=49
         elseif(iel.eq.5.and.ion.eq.8) then
            N=25
         elseif(iel.eq.6.and.ion.eq.6) then
            N=180
         elseif(iel.eq.6.and.ion.eq.7) then
            N=46
         elseif(iel.eq.6.and.ion.eq.8) then
            N=40
         elseif(iel.eq.6.and.ion.eq.9) then
            N=49
         elseif(iel.eq.6.and.ion.eq.10) then
            N=36
         elseif(iel.eq.8.and.ion.eq.8) then
            N=125
         elseif(iel.eq.8.and.ion.eq.9) then
            N=46
         elseif(iel.eq.8.and.ion.eq.10) then
            N=40
         elseif(iel.eq.8.and.ion.eq.11) then
            N=49
         elseif(iel.eq.8.and.ion.eq.12) then
            N=25
         elseif(iel.eq.10.and.ion.eq.8) then
            N=72
         elseif(iel.eq.10.and.ion.eq.9) then
            N=46
         elseif(iel.eq.10.and.ion.eq.10) then
            N=125
         elseif(iel.eq.10.and.ion.eq.11) then
            N=46
         elseif(iel.eq.10.and.ion.eq.12) then
            N=40
         elseif(iel.eq.10.and.ion.eq.13) then
            N=49
         elseif(iel.eq.10.and.ion.eq.14) then
            N=25
         elseif(iel.eq.11.and.ion.eq.2) then
            N=43
         elseif(iel.eq.11.and.ion.eq.3) then
            N=49
         elseif(iel.eq.11.and.ion.eq.13) then
            N=46
         elseif(iel.eq.11.and.ion.eq.14) then
            N=40
         elseif(iel.eq.11.and.ion.eq.15) then
            N=49
         elseif(iel.eq.11.and.ion.eq.16) then
            N=25
         elseif(iel.eq.13.and.ion.eq.9) then
            N=16
         elseif(iel.eq.13.and.ion.eq.10) then
            N=21
         elseif(iel.eq.13.and.ion.eq.11) then
            N=89
         elseif(iel.eq.13.and.ion.eq.12) then
            N=3
         elseif(iel.eq.13.and.ion.eq.13) then
            N=86
         elseif(iel.eq.13.and.ion.eq.14) then
            N=91
         elseif(iel.eq.14.and.ion.eq.10) then
            N=172
         elseif(iel.eq.14.and.ion.eq.11) then
            N=47
         elseif(iel.eq.14.and.ion.eq.12) then
            N=143
         elseif(iel.eq.14.and.ion.eq.13) then
            N=27
         elseif(iel.eq.14.and.ion.eq.14) then
            N=40
         elseif(iel.eq.14.and.ion.eq.15) then
            N=283
         elseif(iel.eq.14.and.ion.eq.16) then
            N=49
         elseif(iel.eq.14.and.ion.eq.17) then
            N=267
         elseif(iel.eq.14.and.ion.eq.18) then
            N=337
         elseif(iel.eq.14.and.ion.eq.19) then
            N=200
         elseif(iel.eq.14.and.ion.eq.20) then
            N=302
         elseif(iel.eq.14.and.ion.eq.21) then
            N=200
         elseif(iel.eq.14.and.ion.eq.22) then
            N=200
         elseif(iel.eq.14.and.ion.eq.23) then
            N=58
         elseif(iel.eq.14.and.ion.eq.24) then
            N=40
         elseif(iel.eq.14.and.ion.eq.25) then
            N=49
         elseif(iel.eq.14.and.ion.eq.26) then
            N=25
         endif
         CALL ATDATchianti(iel,ion,te)
         NMIN=n
         NQ=n
         NP1=N+1
         NP1H=NP1
         NVAR=NP1
         NRC=NP1+1
         nmax=nmin
         xelec=del(ik)

c!! no recombination or photoionization to these

         do k=1,nl
            reco(k)=0.
            ph(k)=0.
         enddo


C
C     CALC. THE FRACTIONAL POPULATIONS (RELATIVE TO THE TOTAL H
C     DENSITYXN1,XN2,XN3, OF LEVELS 1S,2S AND 2P.
C
         C1=CSAHA
         DEN=DE(RS)
         TEV=Te/1.1609E4
         

         iout=0
         if(ik.gt.222222) then
            WRITE(iout,6236)itera,te,den,del(ik),(bb(kh),kh=1,np1)
 6236       format(' o7-pop',i5,1pe12.4,30e12.4)
            write(iout,6237)(ph(kh),kh=1,np1)
 6237       format(' o7 ',1pe12.4,30e12.4)
         endif
         CALL MULTISIMPq(iel,ion,ZI,Te,XI,rltot)
         BB(NP1)=XI(NP1)
         DENEL=DEL(IK)*DE(RS)
         DO J=1,N
            BB(J)=XI(J)/(G(J)*BB(NP1)*DENEL*C1*EXP((E00-E(J))/TEV)
     &           /TS**1.5)
         ENDDO
         ITCON=1
         
      else


         do k=1,401
            weh(k)=1.d-40
         enddo

         rltot = 1.d-40

      endif

      RETURN
      END


      SUBROUTINE POPsimp(iel,ion,Te,XQ,rltot)
      IMPLICIT REAL*8(A-H,O-Z)
c      PARAMETER (MD=350,MDP1=MD+1)
      include 'param'
      PARAMETER (NL=340,NLP1=NL+1)
      COMMON/INUT/IUY
      COMMON/NION/IONQ
      COMMON/INI/INITH,INIT,INITCA,INITHE,INITFE,initfei
      COMMON/IND/IK
      COMMON/ELEC/DEL(MD)
      COMMON/NLEV/e00,NION,N,NP1H,NMAX,NMIN
      COMMON/BOLD/TOLD,TOLDH,TOLDCA,TOLDO,TOLDHE,TOLDFE,toldfei,BB(NL),
     &       BOLH(NL),BOL(NL),BOLCA(NL),BOLHE(NL),BOLFE(NL),bolfei(nl)
      COMMON/RECAL/RECO(NL)
      COMMON/HPOP/XN(6,NL),XN1,XN2,XN3
      COMMON/A14/C(NL,NL),CI(NL),G(NL),E(NL),AS(NL,NL),
     &     WL(NL,NL),DCDT(NL,NL),DCIDT(NL)
      COMMON/A17/RECNET(NL),RECNDT(NL),RECCO(NL),RTE(NL),
     &     PH(NL),DRDT(NL),RECT(6,NL),PHET(NL)
      DIMENSION XINC(NL),XI(NL),
     &     A(NLP1,NLP1),BS(NL)

      DATA CSAHA/2.0708d-16/

      IONQ=ion
      if(iel.eq.4.and.ion.eq.2) then
         N=12
         CALL ATDATn2(te)
      endif
      if(iel.eq.12.and.ion.eq.5) then
         N=5
         CALL ATDATar5(te)
      endif
      if(iel.eq.14.and.ion.eq.7) then
         N=9
         CALL ATDATfe7(te)
      endif
      NMIN=n
      NQ=n
      NP1=N+1
      NP1H=NP1
      NVAR=NP1
      NRC=NP1+1
      nmax=nmin

      xelec=del(ik)

C
C     CALC. THE FRACTIONAL POPULATIONS (RELATIVE TO THE TOTAL H
C     DENSITYXN1,XN2,XN3, OF LEVELS 1S,2S AND 2P.
C
      C1=CSAHA
      DEN=DE(RS)
      TEV=Te/1.1609E4


      iout=0
      if(ik.gt.22222) then
         WRITE(iout,6236)itera,te,den,del(ik),(bb(kh),kh=1,np1)
 6236    format(' he-pop',i5,1pe12.4,30e12.4)
         write(iout,6237)(ph(kh),kh=1,np1)
 6237    format(' ph ',1pe12.4,30e12.4)
         write(iout,6238)(rte(kh),kh=1,np1)
 6238    format(' rte ',1pe12.4,30e12.4)
         write(iout,6239)(reco(kh),kh=1,np1)
 6239    format(' rec ',1pe12.4,30e12.4)
      endif
      CALL MULTISIMPq(iel,ion,ZI,Te,XI,rltot)
      BB(NP1)=XI(NP1)
      DENEL=DEL(IK)*DE(RS)
      DO J=1,N
c         XN(ION,J)=XI(J)
         BB(J)=XI(J)/(G(J)*BB(NP1)*DENEL*C1*EXP((E00-E(J))/TEV)
     &        /TS**1.5)
      ENDDO
      ITCON=1

      TEV=TS/1.1609E4

      RETURN
      END


      SUBROUTINE POPO(RS,TSIN,ZS,XQ,PHO,IFPOP)
      IMPLICIT REAL*8(A-H,O-Z)
c      PARAMETER (MD=350,MDP1=MD+1)
      include 'param'
      PARAMETER (NL=340,NLP1=NL+1)
      COMMON/A19/EM(NL,NL),ESC(NL,NL),TTOT(NL,NL),TOP(NL,NL)
      common/abl/abn(15)
      COMMON/INUT/IUY
      COMMON/HRAT/RHYD,ZMETAL,XRQ,HYR,HEAT,COOL
      COMMON/NION/IONQ
      COMMON/A25/DHDT,DCODT,DZDT,DPDT
      COMMON/A17/RECNET(NL),RECNDT(NL),RECCO(NL),RTE(NL),
     &PH(NL),DRDT(NL),RECT(6,NL),PHET(NL)
      COMMON/MULTI/IHMUL,IHEMUL,IOMUL,ICAMUL,IFEMUL
      COMMON/INI/INITH,INIT,INITCA,INITHE,INITFE,initfei
      COMMON/ITER/ITE
      COMMON/SECX/CSEC(20),DCSDX(20),CISEC(20),DCISDX
      COMMON/NBA/NBACK
      COMMON/IND/IK
      COMMON/ELEC/DEL(MD)
      COMMON/NLEV/e00,NION,N,NP1H,NMAX,NMIN
      COMMON/BOLD/TOLD,TOLDH,TOLDCA,TOLDO,TOLDHE,TOLDFE,toldfei,BB(NL),
     &       BOLH(NL),BOL(NL),BOLCA(NL),BOLHE(NL),BOLFE(NL),bolfei(nl)
      COMMON/ELDEN/IDENS
      COMMON/PHEAT/PHE(NL)
      COMMON/BIN/B1IN,ELD
      COMMON/HPOP/XN(6,NL),XN1,XN2,XN3
      COMMON/A14/C(NL,NL),CI(NL),G(NL),E(NL),AS(NL,NL),
     &     WL(NL,NL),DCDT(NL,NL),DCIDT(NL)
      COMMON/A12/ITMAX
      COMMON/LYBFLUOR/OIFLUOR,OIFLUORH,OIFLUORO,OILYB,OIREC
      common/multisol/imultih,imultihe,imultio,imultica,imultifei,
     &     imultifeii
      DIMENSION XINC(NL),XI(NL),bbo(nl),
     &                                    A(NLP1,NLP1),BS(NL)
      DATA CSAHA/2.0708E-16/
      write(6,*)'POPO called '
      
      ifpop=0
      ITE=-1
      ION=2
      IONQ=2
      N=9
      NMIN=N
      NMAX=N
      CALL ATDATO
      DO 15 NQ=NMIN,NMAX
      NP1=N+1
      NP1H=NP1
      NP2=N+2
      NVAR=NP2
      IIT=0
1398  ICO=0
      IW=1
c      IF(IK.EQ.2) INIT=1
      TS=TSIN
      IF(IIT.GT.10) STOP
      RCM=1.E15*RS
      IF(ITE.EQ.-1) DS=DE(RS)
      IF(N.GT.NMIN) GOTO 600
      DO 1963 KK=1,N
      IF(INIT.EQ.1) GOTO 1963
      BB(KK)=BOL(KK)*(TS/TOLDO)**1.5*EXP((E00-E(KK))*
     &1.1609E4*(1./TOLDO-1./TS))
1963  BS(KK)=BB(KK)
      BB(NP1)=BOL(NP1)
      BS(NP1)=BB(NP1)
      IF(IK.NE.2) GOTO 600
 600  NRC=NP1+1
1399  CONTINUE
      EPS1=1.D-300
      EPS2=1.E-3
      IPRINT=2
      J=1
      IQQ=0
      IF(IK.EQ.2.AND.IUY.EQ.1) XELEC=ZS
      IF(IK.EQ.2.AND.IUY.EQ.1) XN1=0.
      IF(IK.NE.2) XELEC=DEL(IK-1)
      DO 9 ITERA=1,ITMAX
      IF(ITERA.EQ.ITMAX) WRITE(6,*)' NO CONVERGENCE FOR O!!'
C     IF(ITE.EQ.1) TS=BB(NP2)
      RCM=1.E15*RS
      IF(ITE.EQ.-1) DS=DE(RS)
C
C     CALCULATE COLL EXCITATION RATES
C
        IQW=0
C
C     CALC. THE FRACTIONAL POPULATIONS (RELATIVE TO THE TOTAL H
C     DENSITY XN1,XN2,XN3, OF LEVELS 1S,2S AND 2P.
C
      C1=CSAHA
      DEN=DE(RS)
      G1=9.
      G2=5.
      G3=1.
      TEV=TS/1.1609E4
      IF(13.6/TEV.GT.700.) GOTO 777
      XN1=G1*BB(NP1)*DEL(IK)*C1*DEN*EXP(13.60/TEV)*BB(1)/TS**1.5
      PHN1=XN1*PHE(1)
      XN2=G2*BB(NP1)*DEL(IK)*C1*DEN*EXP(11.657/TEV)*BB(2)/TS**1.5
      PHN2=XN2*PHE(2)
C     XN3=G3*BB(NP1)*DEL(IK)*C1*DEN*EXP(9.446/TEV)*BB(3)/TS**1.5
      PHN3=XN3*PHE(3)
777   PHH=PHN1+PHN2+PHN3
      IF(ITERA.GT.1) then
         XELEC=ABn(5)*BB(NP1)*XQ+ZMETAL
      else
         xelec = 0.5
      endif
      CALL COLLEX(ION,TS)
      CALL HORATE(ION,RS,TS)
      icc=1
11    XELOLD=XELEC
      icc=icc+1
      if(icc.gt.500) stop
      CALL SECEX(XELEC)
      IF(N.GT.NMIN) GOTO 6382
      IF(IUY.NE.1) GOTO 6382
      IF(ITERA.NE.1) GOTO 6382
      do nk=1,n
      BB(nk)=nk*0.1
      enddo
      T4=TS/1.E4
C!!   OTS-APP. DOes not work for external rad field!!
cnots      PH(1)=0.
      ALO=0.
      DO J=1,N
         CALL RECOMB(ION,J,TS,AL,RI)
      ALO=ALO+AL
      ENDDO
      RI=ALO*EXP((E(1)-E00)/TEV)*TS**1.5/(CSAHA*G(1))
      BB(1)=RI/(PH(1)+DEN*XELEC*CI(1)+CISEC(2)/DEN)
      QQ=(PH(1)/DEN+CI(1)*XELEC+CISEC(2)/DEN**2)/(ALO*ABn(5)*XQ)
      ZQ=ZS/(ABn(5)*XQ)
      BB(NP1)=SQRT(QQ+(QQ+ZQ)**2/4.)-(ZQ+QQ)/2.
      XELECn=ZS+BB(NP1)*ABn(5)*XQ
      xelec=sqrt(xelecn*xelold)
      IF(ABS(XELEC-XELOLD)/XELOLD.GT.0.01) GOTO 11
      INIT=0
      BB(2)=BB(1)/2.
      BB(3)=BB(1)/3.
      ZI=ABn(5)*XQ
      IDEP=0
      CALL MULTISIMPq(5,1,ZI,TS,XI,rltot)
      BB(NP1)=XI(NP1)
      DENEL=DEL(IK)*DE(RS)
      IF(IDEP.EQ.1) THEN
         DO J=1,N
            BB(J)=XI(J)
            XN(ION,J)=XI(J)*G(J)*BB(NP1)*DENEL*C1*EXP((E00-E(J))/TEV)
     &           /TS**1.5
         ENDDO
      ELSE
         DO J=1,N
            XN(ION,J)=XI(J)
            EXT=DMIN1(700.D0,(E00-E(J))/TEV)
            bbo(j)=bb(j)
            BB(J)=XI(J)/(G(J)*BB(NP1)*DENEL*C1*EXP(EXT)/TS**1.5)
            xinc(j)=bb(j)-bbo(j)
         ENDDO
      ENDIF
      do NK=1,NP1
         BS(NK)=BB(NK)
      enddo
      IHLO=0
      ITCON=1
      IF(TS.LT.-2000..OR.IOMUL.EQ.0) GOTO 111
c!! changed for pulsar calc.
c!!      IF(TS.LT.2800..OR.IOMUL.EQ.0) GOTO 111
6382  CONTINUE
      IHLO=1
      IF(N.GT.NMIN.AND.ITERA.EQ.1) BB(NP1)=BB(NP1-1)
      IF(N.GT.NMIN.AND.ITERA.EQ.1) BB(N)=BB(N-1)
      IF(N.GT.NMIN.AND.ITERA.EQ.1) BS(N)=BB(N)
      IF(N.GT.NMIN.AND.ITERA.EQ.1) BS(NP1)=BB(NP1)
      if(imultio.eq.1) then
         goto 111
      endif

9     CONTINUE
      IF(ITCON.EQ.0) STOP
111   CONTINUE
C
C     CALC. THE FRACTIONAL POPULATIONS (RELATIVE TO THE TOTAL H
C     DENSITY) XN1,XN2,XN3, OF LEVELS 1S,2S AND 2P.
C
      C1=CSAHA
      DEN=DE(RS)
c      DEL(IK)=ABn(5)*XQ*BB(NP1)+ZMETAL
      TEV=TS/1.1609E4
      DO 1928 J=1,N
      EXT=DMIN1(700.D0,(E00-E(J))/TEV)
1928  XN(2,J)=G(J)*BB(NP1)*DEL(IK)*C1*DEN*EXP(EXT)*BB(J)/
     &TS**1.5
 5555 continue
      XN1=0.
      XN2=0.
      XN3=0.
      IF(13.6/TEV.GT.700.) GOTO 778
      XN1=G1*BB(NP1)*DEL(IK)*C1*DEN*EXP(13.60/TEV)*BB(1)/TS**1.5
      XN2=G2*BB(NP1)*DEL(IK)*C1*DEN*EXP(11.657/TEV)*BB(2)/TS**1.5
      XN3=G3*BB(NP1)*DEL(IK)*C1*DEN*EXP(9.446/TEV)*BB(3)/TS**1.5
778   PHN1=XN1*PHE(1)
      PHN2=XN2*PHE(2)
      PHN3=XN3*PHE(3)
      PHH=PHN1+PHN2+PHN3
      TOLDO=TS
      DO 3767 IS=1,NVAR
3767  BOL(IS)=BB(IS)
15    CONTINUE
      IF(ITCON.EQ.0) IFPOP=1
      IF(ITCON.EQ.0) write(6,204)
      XELEC=ABn(5)*BB(NP1)*XQ+ZMETAL
c      IF(IHLO.EQ.1) CALL HLOSS(2,BB,RS,TS,XELEC,XQ,HCA)
      IF(IFPOP.EQ.1.AND.IK.EQ.1) INIT=1
C     LYMAN BETA FLUORESENCE A LA KWAN AND KROLIK
      IF(XN(5,1).GT.0.) THEN
        IF(ABn(1).GT.0.01) THEN
          OIFLUOR=1.41E6*6.65*ESC(7,8)/(5.55+1.1*ESC(7,8))
        ELSE
          OIFLUOR=0.
        ENDIF
        OIFLUORH=OIFLUOR*ABn(5)*XQ*XN(2,1)/(ABn(1)*XN(5,1))
        OIFLUORO=OIFLUOR*XN(5,3)/XN(5,1)
      ENDIF
      RETURN
201   FORMAT(1X,'MATRIX ILL-CONDITIONED OR SINGULAR',E11.4)
202   FORMAT(1X,I8,E10.4,I8,8E12.5)
204   FORMAT(1X,'NO CONVERGENCE FOR O!')
      END


      SUBROUTINE POPO_new_old(RS,TSIN,ZS,XQ,PHO,IFPOP)
      IMPLICIT REAL*8(A-H,O-Z)
c      PARAMETER (MD=350,MDP1=MD+1)
      include 'param'
      PARAMETER (NL=340,NLP1=NL+1)
      COMMON/A19/EM(NL,NL),ESC(NL,NL),TTOT(NL,NL),TOP(NL,NL)
      common/abl/abn(15)
      COMMON/INUT/IUY
      COMMON/HRAT/RHYD,ZMETAL,XRQ,HYR,HEAT,COOL
      COMMON/NION/IONQ
      COMMON/A25/DHDT,DCODT,DZDT,DPDT
      COMMON/A17/RECNET(NL),RECNDT(NL),RECCO(NL),RTE(NL),
     &PH(NL),DRDT(NL),RECT(6,NL),PHET(NL)
      COMMON/MULTI/IHMUL,IHEMUL,IOMUL,ICAMUL,IFEMUL
      COMMON/INI/INITH,INIT,INITCA,INITHE,INITFE,initfei
      COMMON/ITER/ITE
      COMMON/SECX/CSEC(20),DCSDX(20),CISEC(20),DCISDX
      COMMON/NBA/NBACK
      COMMON/IND/IK
      COMMON/ELEC/DEL(MD)
      COMMON/NLEV/e00,NION,N,NP1H,NMAX,NMIN
      COMMON/BOLD/TOLD,TOLDH,TOLDCA,TOLDO,TOLDHE,TOLDFE,toldfei,BB(NL),
     &       BOLH(NL),BOL(NL),BOLCA(NL),BOLHE(NL),BOLFE(NL),bolfei(nl)
      COMMON/ELDEN/IDENS
      COMMON/PHEAT/PHE(NL)
      COMMON/BIN/B1IN,ELD
      COMMON/HPOP/XN(6,NL),XN1,XN2,XN3
      COMMON/A14/C(NL,NL),CI(NL),G(NL),E(NL),AS(NL,NL),
     &     WL(NL,NL),DCDT(NL,NL),DCIDT(NL)
      COMMON/A12/ITMAX
      COMMON/LYBFLUOR/OIFLUOR,OIFLUORH,OIFLUORO,OILYB,OIREC
      common/multisol/imultih,imultihe,imultio,imultica,imultifei,
     &     imultifeii
      DIMENSION XINC(NL),XI(NL),bbo(nl),
     &                                    A(NLP1,NLP1),BS(NL)
      DATA CSAHA/2.0708E-16/

      ifpop=0
      ITE=-1
      ION=2
      IONQ=2
      N=13
      NMIN=N
      NMAX=N
      do i=1,nl
         bol(i)=0.
      enddo
      CALL ATDATO_new

      DO  NQ=NMIN,NMAX
         NP1=N+1
         NP1H=NP1
         NP2=N+2
         NVAR=NP2
         IIT=0
         ICO=0
         IW=1
c     IF(IK.EQ.2) INIT=1
         TS=TSIN
         IF(IIT.GT.10) STOP
         RCM=1.E15*RS
         IF(ITE.EQ.-1) DS=DE(RS)
         NRC=NP1+1
         EPS1=1.D-300
         EPS2=1.E-3
         IPRINT=2
         J=1
         IQQ=0

         IF(IK.EQ.2.AND.IUY.EQ.1) XELEC=ZS
         IF(IK.EQ.2.AND.IUY.EQ.1) XN1=0.
         IF(IK.ge.2) XELEC=DEL(IK-1)

         RCM=1.E15*RS
         IF(ITE.EQ.-1) DS=DE(RS)
C     
C     CALCULATE COLL EXCITATION RATES
C     
         IQW=0
C     
C     CALC. THE FRACTIONAL POPULATIONS (RELATIVE TO THE TOTAL H
C     DENSITY XN1,XN2,XN3, OF LEVELS 1S,2S AND 2P.
C     
         C1=CSAHA
         DEN=DE(RS)
         G1=9.
         G2=5.
         G3=1.
         TEV=TS/1.1609E4
         IF(13.6/TEV.GT.700.) GOTO 777
         XN1=G1*BB(NP1)*DEL(IK)*C1*DEN*EXP(13.60/TEV)*BB(1)/TS**1.5
         PHN1=XN1*PHE(1)
         XN2=G2*BB(NP1)*DEL(IK)*C1*DEN*EXP(11.657/TEV)*BB(2)/TS**1.5
         PHN2=XN2*PHE(2)
C     XN3=G3*BB(NP1)*DEL(IK)*C1*DEN*EXP(9.446/TEV)*BB(3)/TS**1.5
         PHN3=XN3*PHE(3)
 777     PHH=PHN1+PHN2+PHN3
         IF(ITERA.GT.1) then
            XELEC=ABn(5)*BB(NP1)*XQ+ZMETAL
         else
            xelec = 0.5
         endif

         CALL COLLEX(ION,TS)

         CALL HORATE(ION,RS,TS)

 11      XELOLD=XELEC

         CALL SECEX(XELEC)

         T4=TS/1.E4
C!!   OTS-APP. DOes not work for external rad field!!
c     nots      PH(1)=0.
         ALO=0.


         DO J=1,N
            CALL RECOMB(ION,J,TS,AL,RI)
            ALO=ALO+AL
         ENDDO

         RI=ALO*EXP((E(1)-E00)/TEV)*TS**1.5/(CSAHA*G(1))
         BB(1)=RI/(PH(1)+DEN*XELEC*CI(1)+CISEC(2)/DEN)
         QQ=(PH(1)/DEN+CI(1)*XELEC+CISEC(2)/DEN**2)/(ALO*ABn(5)*XQ)
         ZQ=ZS/(ABn(5)*XQ)
         BB(NP1)=SQRT(QQ+(QQ+ZQ)**2/4.)-(ZQ+QQ)/2.
         XELECn=ZS+BB(NP1)*ABn(5)*XQ
         xelec=sqrt(xelecn*xelold)

         IF(ABS(XELEC-XELOLD)/XELOLD.GT.0.01) GOTO 11
         INIT=0

         ZI=ABn(5)*XQ
         IDEP=0

         CALL MULTISIMPq(5,1,ZI,TS,XI,rltot)

         BB(NP1)=XI(NP1)
         DENEL=DEL(IK)*DE(RS)

         IF(IDEP.EQ.1) THEN
            DO J=1,N
               BB(J)=XI(J)
               XN(ION,J)=XI(J)*G(J)*BB(NP1)*DENEL*C1*EXP((E00-E(J))/TEV)
     &              /TS**1.5
            ENDDO
         ELSE
            DO J=1,N
               XN(ION,J)=XI(J)
               EXT=DMIN1(700.D0,(E00-E(J))/TEV)
               bbo(j)=bb(j)
               BB(J)=XI(J)/(G(J)*BB(NP1)*DENEL*C1*EXP(EXT)/TS**1.5)
               xinc(j)=bb(j)-bbo(j)
            ENDDO
         ENDIF
         DO 1827 NK=1,NP1
 1827       BS(NK)=BB(NK)
            IHLO=0
            ITCON=1


C     
C     CALC. THE FRACTIONAL POPULATIONS (RELATIVE TO THE TOTAL H
C     DENSITY) XN1,XN2,XN3, OF LEVELS 1S,2S AND 2P.
C     
            C1=CSAHA
            DEN=DE(RS)

            TEV=TS/1.1609E4
            DO  J=1,N
               EXT=DMIN1(700.D0,(E00-E(J))/TEV)
               XN(2,J)=G(J)*BB(NP1)*DEL(IK)*C1*DEN*EXP(EXT)*BB(J)/TS**1.5
            enddo
 5555       continue
            XN1=0.
            XN2=0.
            XN3=0.
            IF(13.6/TEV.GT.700.) GOTO 778
            XN1=G1*BB(NP1)*DEL(IK)*C1*DEN*EXP(13.60/TEV)*BB(1)/TS**1.5
            XN2=G2*BB(NP1)*DEL(IK)*C1*DEN*EXP(11.657/TEV)*BB(2)/TS**1.5
            XN3=G3*BB(NP1)*DEL(IK)*C1*DEN*EXP(9.446/TEV)*BB(3)/TS**1.5
 778        PHN1=XN1*PHE(1)
            PHN2=XN2*PHE(2)
            PHN3=XN3*PHE(3)
            PHH=PHN1+PHN2+PHN3
            TOLDO=TS
            DO  IS=1,NVAR
               BOL(IS)=BB(IS)
            enddo
         enddo
         IF(ITCON.EQ.0) IFPOP=1
      IF(ITCON.EQ.0) write(6,204)
      XELEC=ABn(5)*BB(NP1)*XQ+ZMETAL
c      IF(IHLO.EQ.1) CALL HLOSS(2,BB,RS,TS,XELEC,XQ,HCA)
      IF(IFPOP.EQ.1.AND.IK.EQ.1) INIT=1
C     LYMAN BETA FLUORESENCE A LA KWAN AND KROLIK
      IF(XN(5,1).GT.0.) THEN
        IF(ABn(1).GT.0.01) THEN
          OIFLUOR=1.41E6*6.65*ESC(7,8)/(5.55+1.1*ESC(7,8))
        ELSE
          OIFLUOR=0.
        ENDIF
        OIFLUORH=OIFLUOR*ABn(5)*XQ*XN(2,1)/(ABn(1)*XN(5,1))
        OIFLUORO=OIFLUOR*XN(5,3)/XN(5,1)
      ENDIF     
      RETURN
201   FORMAT(1X,'MATRIX ILL-CONDITIONED OR SINGULAR',E11.4)
202   FORMAT(1X,I8,E10.4,I8,8E12.5)
204   FORMAT(1X,'NO CONVERGENCE FOR O!')
      END


      SUBROUTINE POPO_new(RS,TSIN,ZS,XQ,PHO,IFPOP)
      IMPLICIT REAL*8(A-H,O-Z)
c      PARAMETER (MD=350,MDP1=MD+1)
c      PARAMETER (NL=340,NLP1=NL+1)
      include "parameters.h"
      COMMON/A19/EM(NL,NL),ESC(NL,NL),TTOT(NL,NL),TOP(NL,NL)
      common/abl/abn(15)
      COMMON/INUT/IUY
      COMMON/HRAT/RHYD,ZMETAL,XRQ,HYR,HEAT,COOL
      COMMON/NION/IONQ
      COMMON/A25/DHDT,DCODT,DZDT,DPDT
      COMMON/A17/RECNET(NL),RECNDT(NL),RECCO(NL),RTE(NL),
     &PH(NL),DRDT(NL),RECT(6,NL),PHET(NL)
      COMMON/MULTI/IHMUL,IHEMUL,IOMUL,ICAMUL,IFEMUL
      COMMON/INI/INITH,INIT,INITCA,INITHE,INITFE,initfei
      COMMON/ITER/ITE
      COMMON/SECX/CSEC(20),DCSDX(20),CISEC(20),DCISDX
      COMMON/NBA/NBACK
      COMMON/IND/IK
      COMMON/ELEC/DEL(MD)
      COMMON/NLEV/e00,NION,N,NP1H,NMAX,NMIN
      COMMON/BOLD/TOLD,TOLDH,TOLDCA,TOLDO,TOLDHE,TOLDFE,toldfei,BB(NL),
     &       BOLH(NL),BOL(NL),BOLCA(NL),BOLHE(NL),BOLFE(NL),bolfei(nl)
      COMMON/ELDEN/IDENS
      COMMON/PHEAT/PHE(NL)
      COMMON/BIN/B1IN,ELD
      COMMON/HPOP/XN(6,NL),XN1,XN2,XN3
      COMMON/A14/C(NL,NL),CI(NL),G(NL),E(NL),AS(NL,NL),
     &     WL(NL,NL),DCDT(NL,NL),DCIDT(NL)
      COMMON/A12/ITMAX
      COMMON/LYBFLUOR/OIFLUOR,OIFLUORH,OIFLUORO,OILYB,OIREC
      common/multisol/imultih,imultihe,imultio,imultica,imultifei,
     &     imultifeii
      DIMENSION XINC(NL),XI(NL),bbo(nl),
     &                                    A(NLP1,NLP1),BS(NL)
      DATA CSAHA/2.0708E-16/
      iel=5
      ion=1
      ifpop=0
      ITE=-1
      IONQ=2
      N=13
      NMIN=N
      NMAX=N
      do i=1,nl
         bol(i)=0.
      enddo
      CALL ATDATO_new

      DO  NQ=NMIN,NMAX
         NP1=N+1
         NP1H=NP1
         NP2=N+2
         NVAR=NP2
         IIT=0
         ICO=0
         IW=1
c     IF(IK.EQ.2) INIT=1
         TS=TSIN
         IF(IIT.GT.10) STOP
         RCM=1.E15*RS
         IF(ITE.EQ.-1) DS=DE(RS)
         NRC=NP1+1
         EPS1=1.D-300
         EPS2=1.E-3
         IPRINT=2
         J=1
         IQQ=0

         IF(IK.EQ.2.AND.IUY.EQ.1) XELEC=ZS
         IF(IK.EQ.2.AND.IUY.EQ.1) XN1=0.
         IF(IK.ge.2) XELEC=DEL(IK-1)

         RCM=1.E15*RS
         IF(ITE.EQ.-1) DS=DE(RS)
C     
C     CALCULATE COLL EXCITATION RATES
C     
         IQW=0
C     
C     CALC. THE FRACTIONAL POPULATIONS (RELATIVE TO THE TOTAL H
C     DENSITY XN1,XN2,XN3, OF LEVELS 1S,2S AND 2P.
C     
         C1=CSAHA
         DEN=DE(RS)
         G1=9.
         G2=5.
         G3=1.
         TEV=TS/1.1609E4

         CALL COLLEX(IONQ,TS)

         CALL HORATE(IONQ,RS,TS)

 11      XELOLD=XELEC

         CALL SECEX(XELEC)

         T4=TS/1.E4
C!!   OTS-APP. DOes not work for external rad field!!
c     nots      PH(1)=0.
         ALO=0.


         DO J=1,N
            CALL RECOMB(IONQ,J,TS,AL,RI)
            ALO=ALO+AL
         ENDDO

         RI=ALO*EXP((E(1)-E00)/TEV)*TS**1.5/(CSAHA*G(1))
         BB(1)=RI/(PH(1)+DEN*XELEC*CI(1)+CISEC(2)/DEN)
         QQ=(PH(1)/DEN+CI(1)*XELEC+CISEC(2)/DEN**2)/(ALO*ABn(5)*XQ)
         ZQ=ZS/(ABn(5)*XQ)
         BB(NP1)=SQRT(QQ+(QQ+ZQ)**2/4.)-(ZQ+QQ)/2.
         XELECn=ZS+BB(NP1)*ABn(5)*XQ
         xelec=sqrt(xelecn*xelold)

         IF(ABS(XELEC-XELOLD)/XELOLD.GT.0.01) GOTO 11
         INIT=0

         ZI=ABn(5)*XQ
         IDEP=0

         CALL MULTISIMPq(iel,ion,ZI,TS,XI,rltot)

         BB(NP1)=XI(NP1)
         DENEL=DEL(IK)*DE(RS)

         IF(IDEP.EQ.1) THEN
            DO J=1,N
               BB(J)=XI(J)
               XN(IONQ,J)=XI(J)*G(J)*BB(NP1)*DENEL*C1*EXP((E00-E(J))/TEV)
     &              /TS**1.5
            ENDDO
         ELSE
            DO J=1,N
               XN(IONQ,J)=XI(J)
               EXT=DMIN1(700.D0,(E00-E(J))/TEV)
               bbo(j)=bb(j)
               BB(J)=XI(J)/(G(J)*BB(NP1)*DENEL*C1*EXP(EXT)/TS**1.5)
               xinc(j)=bb(j)-bbo(j)
            ENDDO
         ENDIF
         DO  NK=1,NP1
            BS(NK)=BB(NK)
         enddo
         IHLO=0
         ITCON=1


C     
C     CALC. THE FRACTIONAL POPULATIONS (RELATIVE TO THE TOTAL H
C     DENSITY) XN1,XN2,XN3, OF LEVELS 1S,2S AND 2P.
C     
         C1=CSAHA
         DEN=DE(RS)

         TEV=TS/1.1609E4
         DO  J=1,N
            EXT=DMIN1(700.D0,(E00-E(J))/TEV)
            XN(2,J)=G(J)*BB(NP1)*DEL(IK)*C1*DEN*EXP(EXT)*BB(J)/TS**1.5
         enddo
 5555    continue
         XN1=0.
         XN2=0.
         XN3=0.
         IF(13.6/TEV.GT.700.) GOTO 778
         XN1=G1*BB(NP1)*DEL(IK)*C1*DEN*EXP(13.60/TEV)*BB(1)/TS**1.5
         XN2=G2*BB(NP1)*DEL(IK)*C1*DEN*EXP(11.657/TEV)*BB(2)/TS**1.5
         XN3=G3*BB(NP1)*DEL(IK)*C1*DEN*EXP(9.446/TEV)*BB(3)/TS**1.5
 778     PHN1=XN1*PHE(1)
         PHN2=XN2*PHE(2)
         PHN3=XN3*PHE(3)
         PHH=PHN1+PHN2+PHN3
         TOLDO=TS
         DO  IS=1,NVAR
            BOL(IS)=BB(IS)
         enddo
      enddo
      IF(ITCON.EQ.0) IFPOP=1
      IF(ITCON.EQ.0) write(6,204)
      XELEC=ABn(5)*BB(NP1)*XQ+ZMETAL
c     IF(IHLO.EQ.1) CALL HLOSS(2,BB,RS,TS,XELEC,XQ,HCA)
      IF(IFPOP.EQ.1.AND.IK.EQ.1) INIT=1
C     LYMAN BETA FLUORESENCE A LA KWAN AND KROLIK
      IF(XN(5,1).GT.0.) THEN
         IF(ABn(1).GT.0.01) THEN
            OIFLUOR=1.41E6*6.65*ESC(7,8)/(5.55+1.1*ESC(7,8))
         ELSE
            OIFLUOR=0.
         ENDIF
         OIFLUORH=OIFLUOR*ABn(5)*XQ*XN(2,1)/(ABn(1)*XN(5,1))
         OIFLUORO=OIFLUOR*XN(5,3)/XN(5,1)
      ENDIF     
      RETURN
 201  FORMAT(1X,'MATRIX ILL-CONDITIONED OR SINGULAR',E11.4)
 202  FORMAT(1X,I8,E10.4,I8,8E12.5)
 204  FORMAT(1X,'NO CONVERGENCE FOR O!')
      END


      SUBROUTINE       POPSI_I(RS,TSIN,ZS,XQ,PHH,IFPOP)
      IMPLICIT REAL*8(A-H,O-Z)
c     PARAMETER (MD=350,MDP1=MD+1)
      include "parameters.h"
      common/abl/abn(15)
      common/timecheck/time,itime
      common/ionx/xion(md,14,27)
c      PARAMETER(NFEL=3000)
      parameter (nlp=30000)
      COMMON/EQUIH/FRH(2,401),WEH(401),SR(NFEL+nlp),WOBS(NL,NL)
      COMMON/INUT/IUY
      COMMON/HRAT/RHYD,ZMETAL,XRQ,HYR,HEAT,COOL
      COMMON/NION/IONQ
      COMMON/A25/DHDT,DCODT,DZDT,DPDT
      COMMON/IND/IK
      COMMON/ELEC/DEL(MD)
      COMMON/NLEV/e00,NION,N,NP1H,NMAX,NMIN
      COMMON/BOLD/TOLD,TOLDH,TOLDCA,TOLDO,TOLDHE,TOLDFE,toldfei,BB(NL),
     &       BOLH(NL),BOL(NL),BOLCA(NL),BOLHE(NL),BOLFE(NL),bolfei(nl)
      COMMON/ELDEN/IDENS
      COMMON/PHEAT/PHE(NL)
      COMMON/BIN/B1IN,ELD
      COMMON/HPOP/XN(6,NL),XN1,XN2,XN3
      COMMON/FELEV/NFEII,nfei
      common/multisol/imultih,imultihe,imultio,imultica,imultifei,
     &     imultifeii
      COMMON/SILEV/NSIII,nsi1
      DIMENSION XINC(NL),XI(NL),A(NLP1,NLP1),BS(NL)

      DATA CSAHA/2.0708E-16/
      ION=7
      IONQ=7

      nsi1=65

      nsi1=56

      n=nsi1
      NMIN=nsi1
      NMAX=NMIN

      CALL ATDATSI_I
      call si_i_recomb(tsin)
      ifpop=0
      ITE=-1
      IIT=0

      TS=TSIN

C
C     CALCULATE COLL EXCITATION RATES
C
      IQW=0
C
C     CALC. THE FRACTIONAL POPULATIONS (RELATIVE TO THE TOTAL H
C     DENSITYXN1,XN2,XN3, OF LEVELS 1S,2S AND 2P.
C

      call collsi_i(ts)

      ZI=1.

      ZI=ABn(10)*XQ
      IDEP=0

      CALL MULTISIMPq(10,1,zi,TS,XI,rltot)
      
      return

      END

      subroutine collsi_i(te)
      implicit real*8 (a-h,o-z)
      include "parameters.h"
      COMMON/NLEV/e00,NION,N,NP1H,NHMAX,NHMIN
      COMMON/A14/C(NL,NL),CI(NL),G(NL),E(NL),AS(NL,NL),WL(NL,NL),DCDT(NL,NL),DCIDT(NL)
      common/abl/abn(15)
      tev=te/1.1609d4
      te5000=te/5000
      t4=te/1.e4
      tesqrt=dsqrt(te)
c *** Si I ***
      ie=iesi
      ion=1
      do i=1,nl
        ci(i)=0.
        do j=1,nl
          c(i,j)=0.
        enddo
      enddo
c
      if (abn(10).gt.0.) then
c     Si I 129.68 68.474 my omegas guessed!! from o i
         o32=9.76e-6*(te-228)+3.46e-11*(te-228.)**2
         o31=3.39e-6*(te-326.)-2.90e-11*(te-326.)**2
         o21=1.89e-6*(te-326.)+8.00e-11*(te-326.)**2
         if(o21.lt.0.) o21=1.e-10
         if(o31.lt.0.) o31=1.e-10
         if(o32.lt.0.) o32=1.e-10
c
         et=(e(3)-e(2))/tev
         if (et.gt.650.) et=650.
         etexp=dexp(-et)
         c(2,3)=8.63e-6*o32*etexp/(g(2)*tesqrt)
         c(3,2)=g(2)*c(2,3)/(g(3)*etexp)
c
         et=(e(3)-e(1))/tev
         if (et.gt.650.) et=650.
         etexp=dexp(-et)
         c(1,3)=8.63e-6*o31*etexp/(g(1)*tesqrt)
         c(3,1)=g(1)*c(1,3)/(g(3)*etexp)
c
         et=(e(2)-e(1))/tev
         if (et.gt.650.) et=650.
         etexp=dexp(-et)
         c(1,2)=8.63e-6*o21*etexp/(g(1)*tesqrt)
         c(2,1)=g(1)*c(1,2)/(g(2)*etexp)
c
c     (Si I) 6591, 10995, 16360  omega for 3p - 1d from
c     m.s. pindzola, a.k. bhatia, a.temkin, phys.rev. 15, 35 (1977)
c     rest scaled from c i with scaling = 3
c      o21=3.97*t4**0.9
         o4_123=3.97*t4**0.9
         o41=o4_123/9.d0
         o42=o4_123/3.d0
         o43=o4_123*5.d0/9.d0
c
         et=(e(4)-e(1))/tev
         if (et.gt.650.) et=650.
         etexp=dexp(-et)
         c(1,4)=8.63e-6*o41*etexp/(g(1)*tesqrt)
         c(4,1)=g(1)*c(1,4)/(g(4)*etexp)
c
         et=(e(4)-e(2))/tev
         if (et.gt.650.) et=650.
         etexp=dexp(-et)
         c(2,4)=8.63e-6*o42*etexp/(g(2)*tesqrt)
         c(4,2)=g(2)*c(2,4)/(g(4)*etexp)
c
         et=(e(4)-e(3))/tev
         if (et.gt.650.) et=650.
         etexp=dexp(-et)
         c(3,4)=8.63e-6*o43*etexp/(g(3)*tesqrt)
         c(4,3)=g(3)*c(3,4)/(g(4)*etexp)
c
c      o31=3.*0.149*(te/5000.)**0.871
         o5_123=3.*0.149*(te5000)**0.871
         o51=o5_123/9.d0
         o52=o5_123/3.d0
         o53=o5_123*5.d0/9.d0
c
         et=(e(5)-e(1))/tev
         if (et.gt.650.) et=650.
         etexp=dexp(-et)
         c(1,5)=8.63e-6*o51*etexp/(g(1)*tesqrt)
         c(5,1)=g(1)*c(1,5)/(g(5)*etexp)
c
         et=(e(5)-e(2))/tev
         if (et.gt.650.) et=650.
         etexp=dexp(-et)
         c(2,5)=8.63e-6*o52*etexp/(g(2)*tesqrt)
         c(5,2)=g(2)*c(2,5)/(g(5)*etexp)
c
         et=(e(5)-e(3))/tev
         if (et.gt.650.) et=650.
         etexp=dexp(-et)
         c(3,5)=8.63e-6*o53*etexp/(g(3)*tesqrt)
         c(5,3)=g(3)*c(3,5)/(g(5)*etexp)
c      o32=3.*0.196*(te/5000.)**0.499
         o54=3.*0.196*(te5000)**0.499
c
         et=(e(5)-e(4))/tev
         if (et.gt.650.) et=650.
         etexp=dexp(-et)
         c(4,5)=8.63e-6*o54*etexp/(g(4)*tesqrt)
         c(5,4)=g(4)*c(4,5)/(g(5)*etexp)

         do j=6,n
            do i=1,j-1
               gi=g(i)
               gj=g(j)
               eij=abs(e(j)-e(i))
               et=eij/tev
               if(et.gt.650.) et=650.
               etexp=dexp(-et)
c     omij=omsi1(i,j)*(te5000**alfsi1(i,j))
c     For the remaining transitions omega=0.1. Just a wild guess!
c     Better for the lower transitions from  Vernazza, Avrett, & Loeser 1976, ApJS, 30, 1
c     see Cissis file
               omij=0.1
               c(i,j)=8.63e-6*omij*etexp/(gi*tesqrt)
               c(j,i)=gi*c(i,j)/(gj*etexp)
            enddo
         enddo

c Collisional ionization rate (cm**3 s**-1) from 
c Vernazza, Avrett, & Loeser 1976, ApJS, 30, 1
         et=(e00-e(1))/tev
         if(et.gt.650.) et=650.
         ci(1)=2.41d-8*(te5000**0.649)*dexp(-et)
         ci(2)=ci(1)
         ci(3)=ci(1)
c
         et=(e00-e(4))/tev
         if(et.gt.650.) et=650.
         ci(4)=0.27d-8*(te5000**0.688)*dexp(-et)
c
         et=(e00-e(5))/tev
         if(et.gt.650.) et=650.
         ci(5)=5.39d-8*(te5000**0.353)*dexp(-et)
c
         et=(e00-e(6))/tev
         if(et.gt.650.) et=650.
c change from CK to CF  numbering         
c     ci(6)=0.35d-8*(te5000**0.582)*dexp(-et)
         ci(7)=0.35d-8*(te5000**0.582)*dexp(-et)
         ci(8)=ci(7)
         ci(9)=ci(7)
c
         et=(e00-e(7))/tev
         if(et.gt.650.) et=650.
c     ci(7)=1.23d-8*(te5000**0.547)*dexp(-et)
         ci(10)=1.23d-8*(te5000**0.547)*dexp(-et)
c
         et=(e00-e(8))/tev
         if(et.gt.650.) et=650.
c     ci(8)=6.63d-8*(te5000**0.525)*dexp(-et)
         ci(16)=6.63d-8*(te5000**0.525)*dexp(-et)
         ci(17)=ci(16)
         ci(18)=ci(16)
      endif

      return
      end

      SUBROUTINE ATDATsi_i
      IMPLICIT REAL*8(A-H,O-Z)
      character*80 dum
      character*2 ch2
      character*4 ch4
      character*5 ch5
      SAVE
      include "parameters.h"
      COMMON/A14/C(NL,NL),CI(NL),G(NL),E(NL),A(NL,NL),WL(NL,NL),DCDT(NL,NL),DCIDT(NL)
      COMMON/A27/OM(NL,NL)
      COMMON/NBA/NBACK
      COMMON/A16/SIG(NL),GA(NL),BP(NL)
      COMMON/NLEV/e00,NION,NH,NP1H,NHMAX,NHMIN
      COMMON/SILEV/NSIII,nsi1
      common/initat/initfeii,initfei,initsi1 
      integer nlsi1
      parameter(nlsi1=66)
      common/csi_idat/gsi1(nlsi1)
      common/si1_reco/tesi1r(81),alsi1r(65,81),totrecsi1(81)
      DIMENSION GS(NL),ES(NL),AS(NL,NL),WLS(NL,NL),OMS(NL,NL)
      ntot=355
      nsi1=65
      N=Nsi1
      E00=8.15168
      IF(Initsi1.EQ.0) THEN
         OPEN(23,FILE='./ATDAT/Si_I_levels_Aij_clean.dat',status='old')
         rewind 23
         do i=1,ntot
            read(23,*)nn,wn,ji,gsi
            if(i<=nsi1) then
               gsi1(i)=gsi
               gs(i)=gsi
               es(i)=wn/8065.46d0
            endif
         enddo
         DO 5395 I=1,N
            DO 5394 J=1,N
               WLS(I,J)=0.0
               IF(I.EQ.J) GOTO 5394
               IF(ES(I).EQ.ES(J)) GOTO 5394
               WLS(I,J)=ABS(12398.54/(ES(I)-ES(J)))
 5394       CONTINUE
 5395    CONTINUE
         DO I=1,N
            DO  J=1,N
               OMS(I,J)=0.
               AS(I,J)=0.
            enddo
         enddo
c     read(23,987)dum
 987     format(a)
         do i=1,1000
            read(23,*,end=88)iu,il,enl,enu,wlq,aij
            if(iu < n) then
               as(iu,il)=aij+as(iu,il)
            endif
         enddo
 88      continue
         close (23)
         open(23,file='./ATDAT/SiI_recomb_nahar_sorted.dat',status='old')
         read(23,*)(tesi1r(k),k=1,81)
         do i=1,nsi1
            read(23,*)il,(alsi1r(i,k),k=1,81)
         enddo
         read(23,987)dum
         read(23,*)(totrecsi1(k),k=1,81)
         do iu=2,n
            do il=1,iu-1
c               write(62,955)iu,il,wls(iu,il),as(iu,il)
 955           format(2i4,f10.1,1pe12.4,e12.4)
            enddo
         enddo
         Initsi1=1
      ENDIF
      DO I=1,N
         E(I)=ES(I)
         G(I)=GS(I)
         DO J=1,N
            A(I,J)=AS(I,J)
            OM(I,J)=OMS(I,J)
            WL(I,J)=WLS(I,J)
         ENDDO
      ENDDO
C
      SIG(1)=0.
      SIG(2)=0.
      DO 8457 I=3,N
 8457 SIG(I)=0.
      GA(1)=2.99
      GA(2)=2.5
      DO 1009 I=3,N
1009  GA(I)=3.
      GA(8)=3.56
      GA(9)=3.56
      DO 165 I=1,N
165   C(I,I)=0.
      RETURN
      END


      subroutine si_i_recomb(te)
      implicit real*8 (a-h,o-z)
c
      integer mie,mion

      parameter(mie=15,mion=5)
c
      data ieh/1/,iehe/2/,iec/3/,ien/4/,ieo/5/,iene/6/,iena/7/,
     &     iemg/8/,iesi/9/,ies/10/,iear/11/,ieca/12/,iefe/13/,
     &     ieco/14/,ieni/15/,iemax/15/
cqqqq
      COMMON/NLEV/e00,NION,N,NP1H,NMAX,NMIN
      common/cnlev/nlev(mie,mion)
      COMMON/SILEV/NSIII,nsi1
      integer nlsi1
      parameter(nlsi1=66)
      common/csi_idat/gsi1(nlsi1)
      common/sirecombi/recosi1(nlsi1)
      common/si1_reco/tesi1r(81),alsi1r(65,81),totrecsi1(81)
c     
c     Si I
      do il=1,nsi1
         recosi1(il)=0.d0
      enddo
      if (te.lt.tesi1r(1)) then
         tconst=0.d0
         it1=1
         it2=1
      elseif (te.ge.tesi1r(81)) then
         tconst=0.d0
         it1=81
         it2=81
      else
         if(te.lt.tesi1r(1)) then
            it1=1
            it2=2
         elseif(te.gt.tesi1r(80)) then
            it1=79
            it2=80
         else
            do i=1,80
               if (te.ge.tesi1r(i).and.te.lt.tesi1r(i+1)) then
                  it1=i
                  it2=i+1
                  goto 9
               endif
            enddo
 9          continue
         endif
         tconst=(te-tesi1r(it1))/(tesi1r(it2)-tesi1r(it1))
      endif

      do il=1,nsi1
         recosi1(il)=alsi1r(il,it1)+tconst*(alsi1r(il,it2)-alsi1r(il,it1))
      enddo
c     
      totrec=totrecsi1(it1)+tconst*(totrecsi1(it2)-totrecsi1(it1))
c     
c     Add the remaining recombinations to the 8 highest levels 
c     (114 - 121; y5P, y3F, y3D) weighted
c     with their statistical weights so that the total recombination 
c     coefficient is correct ; gsum=42 (7,5,7,3,5,7,5,3)
c     
      recsum=0.d0
      do il=1,n
         recsum=recsum+recosi1(il)
      enddo
      recadd=totrec-recsum
      if (recadd.lt.0.) then
         write(50,*) 'Error in calc. rec. coeff. for Fe I !!'
      endif
c     add the remaining recomb. contribution to
c     all the levels, not only the highest
      gsum=0.d0
      do il=1,nsi1
         gsum=gsum+gsi1(il)
      enddo
      do il=1,nsi1
         recosi1(il)=recosi1(il)+recadd*gsi1(il)/gsum
      enddo
c     
      return
      end

      



