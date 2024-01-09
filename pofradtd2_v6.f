      subroutine enum_cf(ik)
      IMPLICIT REAL*8(A-H,O-Z)

      include 'param'
      common/ionx/xion(md,14,27)
      COMMON/ABUN/AB(20)
      common/abl/abn(15)
      dimension iabcf(14)

c
c conversion from logical to cf  abundance enumeration      
c
      data iabcf/1,2,4,5,3,13,10,7,9,6,12,14,11,8/

c  CF  normal

c  1 = 1,

c  2 = 2,

c  3 = 5,   O

c  4 = 3,   C

c  5 = 4,   N

c  6 = 10,  Si

c  7 = 8,   Mg

c  8 = 14,  Fe

c  9 = 9,   Al

c 10 = 7,   Na

c 11 = 13,  Ca

c 12 = 11,  S

c 13 = 6,   Ne

c 14 0 12   Ar


      xe = 0.

      return
      end

      subroutine cfenum
c
c translate from normal system to cf enumeration
c note that first index is ionization stage and second is element
c
      integer inum
      common/cfum/inum(27,14)
      data inum/1,26*0,2,3,25*0,12,13,8,9,10,11,21*0,
     &         18,19,20,21,22,23,24,20*0,14,15,16,17,7,4,5,6,19*0,
     &         72,73,74,75,76,77,78,79,80,81,17*0,
     &         53,54,55,24*0,39,40,41,24*0,46,47,48,49,23*0,
     &         25,26,27,28,29,30,31,32,33,34,35,36,37,38,13*0,
     &         56,57,58,59,60,61,62,63,64,65,66,67,68,69,71,12*0,
     &         82,83,96,97,98,99,100,20*0,50,51,52,24*0,
     &         42,43,44,45,84,85,86,87,88,89,90,91,92,93,94,12*0/
      return
      end      


      DOUBLE PRECISION FUNCTION RAD(TE,XELEC,den,dr,IFPOP)
      IMPLICIT REAL*8(A-H,O-Z)
c      real*4 om21,om31,om41,om51,om32,om42,om52,om43,om53,om54
      real*8 jlah,j2gh,jrech,jemh
      character*12 lev
      character*80 lab
       include 'param'
c      PARAMETER (MD=350,MDP1=MD+1)
      PARAMETER (NL=340,NLP1=NL+1)
      PARAMETER(NFEL=3000)
c      PARAMETER (NE1=-100,NE2=250,NE3=NE2+1)
C     ****************************************************************
C     *****
C     THIS SUBROUTINE CALCULATES ALL RADIATIVE AND COLLISIONAL RATES
C     AS WELL AS THE IONIZATION STRUCTURE (EXCEPT FOR SILICON AND
C     MAGNESIUM ), FOR A GIVEN TEMPERATURE TE.
C

      COMMON/EMHY/RECEM(10,NL),TWOPH,TWOPHHEI,COHI,PHI,PHIT,PO,POT
      common/oiiilines/em4959,em5007
      COMMON/DIF/EM(MD,NE1:NE2),TAU(MD,NE1:NE2),TAUTOT(MD,NE1:NE2),
     &           EMC(MD,NE1:NE2)
      common/heiires/x2he,he2coll12,he2jla,he2j2g,he2jrec(4),
     &                                                he2jem(7,4)
      common/hres/x2h,coll12h,jlah,j2gh,jrech(4),jemh(7,4)
      COMMON/RECMET/recappc,recappo,recappsi,recapps,recappfe
      COMMON/GSREC/ALGS(100),ALEX(100),ALTOT(100),RECEX(100),RECGS(100)
      COMMON/CSEX/CSLYA
      COMMON/HYPOP/XNH(NL)
      COMMON/INT/FL(2,NE1:NE2),SI(14,27,NE1:NE2),E1(NE1:NE3),E(NE1:NE3)
      COMMON/FRE/NINT,JMIN,JJ
      COMMON/SECX/CSEC(20),DCSDX(20),CISEC(20),DCISDX
      COMMON/SPOT/OTSP(16)
      COMMON/NLEV/e00,NION,NHY,NP1H,NMA,NMI
      COMMON/GAMPAR/GAMLUM,GAELH,GAHE
      COMMON/COL/RTE(4),FF,HY,HE,C131,C142,C151,C161,C231,C241,COH,
     &COHE,C351,C321,C331,C341,C222,C232,C221,C332,C441
      parameter (nlp=30000)
      COMMON/EQUIV/RADF,FR(2,nlp),W(nlp),CIN(nlp),WLI(nlp+NFEL),FB(300),tauline(nlp+nfel)
      common/lineioncf/ilabcf(nlp+nfel)
      COMMON/EQUIH/FRH(2,401),WEH(401),SR(NFEL+nlp),WOBS(NL,NL)
      COMMON/COLD/COLTOT(14,27),COLTOTT(14,27),SRED(NFEL+nlp)
      COMMON/PHOTOHEAT/PHEAT(10,NL),PHEATT(10,NL)
      COMMON/HEA/CO,PHE0,PHEI,POX(8),PC(6),PN(7),PMG,PSI,PFE,
     &     pnei(10)
      COMMON/TEQQ/TEZ
      common/phycon/te_h_he
      COMMON/REHEL/REHE21,REHE22,REC31,REN41
      COMMON/IND/IK
      COMMON/mdIND/kradius
      COMMON/PHQ/ZEA(14,27),GEA(14,27),ZKA(14,27)
      COMMON/COSUL/COSUL(18)
      COMMON/REC/AL2(16)
      parameter (mh=2,mhe=3,mc=7,mn=8,mo=9,mne=11,mna=12,mmg=13,
     &     mal=14,msi=15,ms=17,mar=19,mca=21,mfe=27,mni=29,
     &     mmat=mh+mhe+mc+mn+mo+mne+mna+mmg+mal+msi+ms+mar+mca+mfe+mni,
     &     ispecy=15)
      common/abc1/alc(mc),aln(mn),alo(mo),alne(mne),alna(mna),
     &     almg(mmg),alal(mal),alsi(msi),alsu(ms),alar(mar),alca(mca),
     &     alfe(mfe),alni(mni)

      COMMON/ABC/AL(16)

      common/augalu/alaug(27),coaug(26)

      COMMON/ELEC/DEL(MD)
      common/abl/abn(15)
      COMMON/RADIE/R(0:MD)
      COMMON/NHY/NH
      COMMON/HPOP/XNQ(6,NL),XN1,XN2,XN3
      COMMON/BOLD/TOLD,TOLDH,TOLDCA,TOLDO,TOLDHE,TOLDFE,toldfei,BB(NL),
     &       BOLH(NL),BOL(NL),BOLCA(NL),BOLHE(NL),BOLFE(NL),bolfei(nl)
      COMMON/HRAT/RHYD,ZEL,XEL,HYR,HEAT,COOL
      COMMON/REC1/RE1(7),RHE2B
      COMMON/A14/CQ(NL,NL),CIQ(NL),GQ(NL),EI(NL),AQFE(NL,NL),WL(NL,NL)
     &,DCDT(NL,NL),DCIDT(NL)
      COMMON/A17/RECNET(NL),RECNDT(NL),RECCO(NL),RTEQ(NL),
     &PHN(NL),DRDT(NL),RECT(6,NL),PHET(NL)
      COMMON/RECAL/RECO(NL)
      COMMON/A19/EMCA(NL,NL),ESCCA(NL,NL),TTOT(NL,NL),TOP(NL,NL)
      COMMON/A31/TOPI(10,10,10),EMHY(NL,NL)
      COMMON/IOINF/IOSPEC,IOLEVEL,IOEMISS,IOTERM,IOION,IOELITER
      COMMON/NION/IONQ
      COMMON/IOPAR/IOSHELL,IORAD1,IORAD2,IOFLUX
      COMMON/MULTI/IHMUL,IHEMUL,IOMUL,ICAMUL,IFEMUL
      COMMON/LOWION/ILOWION
      COMMON/CINOUT/INOUT,IPULS
      COMMON/CICS/RSHOCK,ics,ics87a
      COMMON/FELEV/NFE,nfei
      COMMON/WFE/WOBFE(NL,NL),wobfei(nl,nl),IQPOPH
      COMMON/FESAVE/INIFE(4),inifei(4),BOLFES(4,NL),TOLDFES(4),
     &     bolfeis(4,nl),toldfeis(4)
      common/ctionfe/phfe(2)
      common /initcoll/initcoll
      common/nperm/nlinetot
      common/ichia/ipopch
      common/coolants/cl(nlp)
      common/contdest/idest
      common/linepoint/nlines,jline(nlp),iion(nlp)
      common/bowem/ emb644, emb703, emb3210, emb3287, em374e,
     &              emb4642, emb691, emb990, emb4100, emb452
      common/w/wlin(401),wlij(NL,NL),taul(401)
      parameter(nspec=91,neq=nspec-3)
      common/ionabun/iab(100)
      COMMON/SILEV/NSIII,nsi1
      DIMENSION YI(7),DELA(MD),FD(MD),YEX1(7),CEX(500)
      DIMENSION PF(27),ZE(14,27),ZK(14,27),GE(14,27)
      DIMENSION rech(10,3),g(100),emfe4(20,20)
      DIMENSION XN(30),ZC(7),XO(30),XR(30),ZA(8),XNA(30),XSU(30),XAR(30)
      DIMENSION OXR(10),FBO(30)
      DIMENSION XC(30),ZB(6),XF(30),XSI(30),XM(30),XAL(30),XCA(30),
     &XNE(30),LHEI(16),WLHEI(16),LFEII(16),WLFEII(16),ZS(16),
     &WLHI(10),X(100),CHI(100)
      DIMENSION ZNEO(10),XI(15),ZSUL(16),A21HNU(5),gai(14),xico(nl)
      dimension supr_high(26)
      dimension fextot(15)      
      common/numbow/nbow
      common/ionx/xion(md,14,27)
      common/raug/zion(30,27,7),geion(30,27)
      common/ca4/teca4(50),coll_ca4(3,50)
      real*8 line_cool
      common/line_c/line_cool(14,27)
      common/fbapp/fbapp6300,fbapp5007
      common/kmaxpop/kmaxp(14,27)
      common/timecheck/time,itime
      common/position/rcomm
      common/debug/ideb
      common/recoxtest/RECCA,RECOX,RECSI,RECS,RECFE
      common/intcosio/initco_sio
      common/number/ielec,igam,icr
      common/cophoto/cophr
      COMMON/A2/SIg,T0,TAUrL
      common/line_em/cinx(14,26,401),taulinex(14,26,401),wlix(14,26,401)
     &     ,ilabcfx(14,26,401)
      dimension omar(5,5),omca4(3)
      common/rec_coll/alrec(30,30),collion(30,30),ction(14,27)
      common/collion/coolx(14,27)
      dimension photoion(14,26),dheat(14,27)
      dimension nline_min(14,26)
      character*10 ch10
      character*8 hatom      
      real*8 wlmgs(1000),glmgs(1000),gumgs(1000),a21mgs(1000)
     &           ,cs1e4(1000),alfmgs(1000)
      integer ielmgs(1000),ionmgs(1000),stop_tmin
      common/mi_te/stop_tmin
      save wlmgs,glmgs,gumgs,a21mgs,cs1e4,alfmgs
      save ielmgs,ionmgs,nline,nline_min
      parameter(nchi22=20)
      PARAMETER (nups=65)
      common/chianti_2022_cs/upsil(nchi22,nups,nups)
      DATA PI/3.1415926E0/,ELCH/1.60219E-12/,AMU/1.660531E-24/
      DATA CSAHA/2.0708E-16/
      DATA A21HNU/1.,7.663E-3,1.301E-4,3.390E-5,1.144E-5/
      DATA LHEI/7,5,17,38,9,14,25,32,49,21,42,1,35,52,4,63/
      DATA WLHEI/584.4,10830.,3889.,3188.,20581.,7065.,5876.,4713.,
     &4471.,42944.,12547.,625.6,21120.,17002.,591.4,18648./
      DATA WLHI/1216.,1025.,6563.,972.,4861.,18751.,949.,4340.,
     &12818.,40512./
      DATA LFEII/7,5,17,38,9,14,25,32,49,21,42,28,35,52,44,63/
      DATA WLFEII/584.,10830.,3889.,3188.,20581.,7065.,5876.,4713.,
     &4471.,42944.,12547.,186233.,21120.,17002.,19543.,18648./
      DATA IAB/1,2,2,5,5,5,5,3,3,3,
     &         3,3,3,5,5,5,5,4,4,4,
     &         4,4,4,4,10,10,10,10,10,10,
     &         10,10,10,10,10,10,10,10,8,8,
     &         8,14,14,14,14,9,9,9,9,13,
     &         13,13,7,7,7,11,11,11,11,11,
     &         11,11,11,11,11,11,11,11,11,11,
     &         11,6,6,6,6,6,6,6,6,6,
     &         6,12,12,14,14,14,14,14,14,14,
     &         14,14,14,14,5,5*12/

      dimension nionel(14)
      data nionel/1,2,6,7,8,10,11,12,13,14,16,18,20,26/
      real*8 ar2_te(11),ar2_coll(11)
      data ar2_te/1.58E+03,2.51E+03,3.98E+03,6.31E+03,1.00E+04,1.58E+04,
     &     2.51E+04,3.98E+04,6.31E+04,1.00E+05,1.58E+05/
      data ar2_coll/2.48,2.54,2.63,2.77,2.93,3.09,3.19,3.2,3.13,2.97,2.71/
      real*8 ar6_te(43),ar6_coll(43)
      data ar6_te/2.00E+01,5.00E+01,1.00E+02,2.00E+02,3.00E+02,4.00E+02,
     &     5.00E+02,6.00E+02,7.00E+02,8.00E+02,9.00E+02,1.00E+03,1.50E+03,
     &     2.00E+03,2.50E+03,3.00E+03,3.50E+03,4.00E+03,4.50E+03,5.00E+03,
     &     5.50E+03,6.00E+03,6.50E+03,7.00E+03,7.50E+03,8.00E+03,8.50E+03,
     &     9.00E+03,9.50E+03,1.00E+04,1.10E+04,1.20E+04,1.30E+04,1.40E+04,
     &     1.50E+04,1.60E+04,1.80E+04,2.00E+04,2.40E+04,2.80E+04,3.20E+04,
     &     3.60E+04,4.00E+04/
      data ar6_coll/3.04,3.02,3.0,2.97,2.95,2.96,2.97,3.0,3.03,3.06,3.08,
     &     3.11,3.23,3.4,3.62,3.86,4.11,4.35,4.57,4.78,4.96,5.13,5.27,5.4,5.51,
     &     5.61,5.7,5.77,5.84,5.9,6.0,6.08,6.14,6.19,6.23,6.26,6.3,6.33,
     &     6.36,6.36,6.36,6.35,6.33/
      real*8 cpu_start, cpu_tsec
      common/cpu_t/cpu_start,cpu_tsec      
      save feifs_c,fecon
      write(6,*)' te,xelec',te,xelec
      te_h_he=te
      do ii=1,nlp
         cin(ii)=0.
      enddo
      if(itime.eq.1) then
         do ii=1,nlp+NFEL
            WLI(ii)=0.
            tauline(ii)=0.
         enddo
         do ielem=1,14
            do ionq=1,26
               nline_min(ielem,ionq)=0
            enddo
         enddo
      endif

      ielec = 88

      if(initphd.eq.0) then
         cophr=0.
      endif


      TLIMO=100.
c!!   use popca to very low T!
      TLIMCA=100.
      TEZ=TE
      T4=TE/1.E4
      T3=TE/1.E3
      XEL=XELEC
      DEL(IK)=XEL
c      write(6,9262)ik,xel
 9262 format('xel in rad at 9262 ',i5,1pe12.3)
      eden=xel*den
      MIMAX=50
C
C     DIVIDE ALL IONIZATION AND HEATING RATES BY THE DENSITY
C


      do iel=3,14
         do iz=1,nionel(iel)

c     ge(iel,iz) = gea(iel,iz)/den
            ge(iel,iz) = geion(iel,iz)/den
            ze(iel,iz) = zea(iel,iz)/den
            zk(iel,iz) = zka(iel,iz)/den
         enddo            
      enddo            

      do i=1,7
         xc(i)=xion(ik,3,i)
      enddo

      do i=1,8
         xn(i)=xion(ik,4,i)
      enddo

      do i=1,9
         xo(i)=xion(ik,5,i)
      enddo

      do i=1,11
         xne(i)=xion(ik,6,i)
      enddo

      xna(1)=xion(ik,7,1)
      xna(2)=xion(ik,7,2)

      xm(1)=xion(ik,8,1)
      xm(2)=xion(ik,8,2)

      do i=1,4
         xal(i)=xion(ik,9,i)
      enddo

      do i=1,15
         xsi(i)=xion(ik,10,i)
      enddo

      do i=1,17
         xsu(i)=xion(ik,11,i)
      enddo

      do i=1,8
         xar(i)=xion(ik,12,i)
      enddo

      do i=1,10
         xca(i)=xion(ik,13,i)
      enddo

      do i=1,27
         xf(i)=xion(ik,14,i)         
      enddo

      do iel=3,14
         do ion=1,27
            line_cool(iel,ion)=0.
         enddo
      enddo

      MQ=0
      MI=1
  776 CONTINUE
      DELA(1)=xelec
      DEEL=DELA(MI)
      ionq=3
c      CALL SECEX(DELA(MI))
      IF(DELA(MI).LE.0.)write(6,*)'NEG EL ',IK,MI,DELA(MI),DEL(IK)

C     ***************************************************************
C     *****
C     RECOMBINATION COOLING (CASE A) FOR HYDROGEN AND HELIUM
C     RATES ARE ADJUSTED TO AGREE WITH OSTERBROCKS VALUES FOR
C     T<2.E4*Z**2
C     *****
C     ***************************************************************
      YI(4)=64.*YI(1)
      YI(5)=36.*YI(1)
      YI(6)=49.*YI(1)
      YI(7)=196.*YI(1)
C
C     CARBON
C
      RECCA=0.
      DO K=1,6
         RECCA=RECCA+ABn(3)*XC(K+1)*1.5*1.38E-16*TE*ALC(K+1)
      ENDDO
C
C     OXYGEN
C
      RECOX=0.
c!!! neglect charge transfer
      alcto2=0.
      alcto3=0.
      DO K=1,8
         ALFO=ALO(K+1)
         IF(K.EQ.1) ALFO=ALO(2)-ALCTO2
         IF(K.EQ.2) ALFO=ALO(3)-ALCTO3
         RECOX=RECOX+ABn(5)*XO(K+1)*1.5*1.38E-16*TE*ALFO
         if(del(ik).lt.0.6) then
c            write(6,9198)abn(5),xo(k+1),alfo,recox
 9198       format('recox ',1pe12.3,10e12.3)
         endif
      ENDDO
C
C     SILICON
C
      RECSI=0.
      DO K=1,14
         CSI=ABn(10)*XSI(K+1)*1.5*1.38E-16*TE*ALSi(K+1)
         RECSI=RECSI+CSI
      enddo
C     SULPHUR
C
C
      RECS=0.
      DO K=1,10
         CSU=ABn(11)*XSU(K+1)*1.5*1.38E-16*TE*ALSU(K+1)
         RECS=RECS+CSU
      enddo
C     IRON
C
C
      RECFE=0.
      DO K=1,15
         CFE=ABn(14)*XF(K+1)*1.5*1.38E-16*TE*ALFE(K+1)
         RECFE=RECFE+CFE
      enddo
C     ***************************************************************
C     *****
C     FREE-FREE COOLING,(COX&TUCKER AP.J. 157:1157) FOR T>1+7
C     AND SPITZER TABLE 3.3 FOR GAUNT FACTOR
C     *****
C     ***************************************************************
      FF=0.
      DO KK=4,7
         ZZ=1.
         IF(KK.EQ.3) ZZ=4.
         IF(KK.GE.4) ZZ=(KK-3.)**2
         GAUNT=-1.08+0.925*LOG10(TE/ZZ)-0.085*(LOG10(TE/ZZ))**2.
         IF(KK.EQ.7) ABX=ABn(14)*(1.-XF(1)-XF(2)-XF(3)-XF(4))
         IF(KK.GE.4.AND.KK.LE.6) ABX=ABn(14)*XF(KK-2)
 8437    FF=FF+1.426E-27*ZZ*ABX*GAUNT*SQRT(TE)
      enddo

C     ***************************************************************
C     *****
C     LINE COOLING RATES AND EMITTED FLUX TO INFINITY ( W ).
C     *****
C     ***************************************************************
C     **************************************************************
C     *****
C     FORBIDDEN LINE COOLING
C     *****
C     **************************************************************
C
C     COLL. EXIT. OF CARBON
C
C     (C I) 4619,8729,9812 A. NUSSBAUMER & RUSCA AA 72, 129.
      O21=0.603*(TE/5000.)**0.96
      O31=0.149*(TE/5000.)**0.871
      O32=0.196*(TE/5000.)**0.499
      Z=ABN(3)*XC(1)
       CALL FORB3(1,3,1,1.262D0,1.420D0,3.26D-4,2.73D-3,0.528D0,O21,
     &     O31,O32,9.D0,5.D0,1.D0,FB(20),FB(21),FB(22),TE,Z,F)
      if(te.gt.1.e3.and.te.lt.6.e4) then
         CALL DIELB(-2.02d-2,.3799d0,0.0890d0,-0.57d-2,.9237d0,T4,ALDB)
      else
         aldb=0.
      endif
      DIFB=ABN(3)*XC(2)*ALDB*1.2633*1.602E-12
c add dielectronic contrib.
      wlin(1)=wlin(1)+ABN(3)*XC(2)*ALDB*1.2633*1.602E-12
      C211=ABN(3)*XC(1)*F
      ll=0
      if(itime.eq.1) nline_min(3,1)=ll
      ll=nline_min(3,1)
      do k=1,3
         ll=ll+1
         cin(ll)=weh(k)
         tauline(ll)=taul(k)
         if(itime.le.888888) wli(ll)=wlin(k)
         ilabcf(ll)=301
         cinx(3,1,k)=weh(k)/del(ik)
         taulinex(3,1,k)=taul(k)
         wlix(3,1,k)=wlin(k)
         ilabcfx(3,1,k)=301
         kx=k
      enddo

      DO NM=20,22
         FB(NM)=ABN(3)*XC(1)*FB(NM)
      enddo
      FB(20)=FB(20)+ABN(3)*XC(2)*ALDB*1.2633*1.602E-12


C     C I 2966-68 A  OMEGA FROM HAYES & NUSSBAUMER AA 134:193
C
      Z=ABN(3)*xion(ik,3,1)
      OM=0.475*(TE/5000.)**0.5
      CALL RLOSS(3,1,2967.D0,OM,9.D0,5.D0,6.D8,RS,XEL,Z,TE,cl(27)
     &     ,W(27),CIN(27),FR(2,27),27,WLI(27))
      tauline(27)=t0
      ilabcf(27)=301
      cinx(3,1,kx+1)=cin(27)
      taulinex(3,1,kx+1)=t0
      wlix(3,1,kx+1)=wli(27)
      ilabcfx(3,1,kx+1)=301
      C212=cl(27)
C     C I 370.4 609.1 MY MAZDA COLLISION STRENGTHS FROM O I!!
      O32=9.76E-6*(TE-228)+3.46E-11*(TE-228.)**2
      O31=3.39E-6*(TE-326.)-2.90E-11*(TE-326.)**2
      O21=1.89E-6*(TE-326.)+8.00E-11*(TE-326.)**2
      IF(O21.LT.0.) O21=1.e-10
      IF(O31.LT.0.) O31=1.e-10
      IF(O32.LT.0.) O32=1.e-10
      Z=ABN(3)*XC(1)
      CALL FORB3(13,3,1,3.348D-3,2.03D-3,7.93D-8,1.71D-14,2.65D-7,O21,
     &     O31,O32,1.D0,3.D0,5.D0,FB(32),FBQQ,FB(33),TE,Z,F)

      do k=1,3
         ll=ll+1
         cin(ll)=weh(k)
         tauline(ll)=taul(k)
         if(itime.le.888888) wli(ll)=wlin(k)
         ilabcf(ll)=301
         cinx(3,1,kx+2+k)=weh(k)
         taulinex(3,1,kx+2+k)=taul(k)
         wlix(3,1,kx+2+k)=wlin(k)
         ilabcfx(3,1,kx+2+k)=301
         C212=cl(27)
      enddo


      DO NM=32,33
         FB(NM)=ABN(3)*XC(1)*FB(NM)
      enddo
      C212=ABN(3)*XC(1)*F
      line_cool(3,1)=c211+c212
      

C     C II 157.73 MY HAYES & NUSSBAUMER  A A 134, 193. Multi 
      Z=ABN(3)*XC(2)
      IF(TE.LT.1.E3) OM=1.82
      IF(TE.Ge.1.E3) OM=1.87*T3**0.21
      CALL RLOSS(3,1,1.58D6,OM,2.D0,4.D0,2.29D-6,RS,XEL,Z,TE,cl(26)
     &     ,W(26),CIN(26),FRQQ,26,WLI(26))
      tauline(26)=t0
      ilabcf(26)=302
      C223=cl(26)
c      cinx(3,2,3)=cin(26)
c      taulinex(3,2,3)=t0
c      wlix(3,2,3)=wli(26)
c      ilabcfx(3,2,3)=302
C
C     C II 1334 A  OMEGA FROM HAYES & NUSSBAUMER AA 134:193 . Multi (?)
C
      Z=ABN(3)*xion(ik,3,2)
      OM=5.79*T4**0.05
      CALL RLOSS(3,2,1334.D0,OM,6.D0,10.D0,6.D8,RS,XEL,Z,TE,cl(1)
     &     ,W(1),CIN(1),FR(2,1),1,WLI(1))
      tauline(1)=t0
      ilabcf(1)=302
      C221=cl(1)
c      cinx(3,2,2)=cin(1)
c      taulinex(3,2,2)=t0
c      wlix(3,2,2)=wli(1)
c      ilabcfx(3,2,2)=302
C
C     C II] 2326 OMEGA FROM HAYES & NUSSBAUMER AA 134:193 A=44.1 MULTILEV!
C
      CALL RLOSS(3,2,2326.D0,2.8D0,6.D0,12.D0,4.41D1,RS,XEL,Z,TE,cl(2)
     &     ,W(2),CIN(2),FR(2,2),2,WLI(2))
      tauline(2)=t0
      ilabcf(2)=302
      C222=cl(2)
      kmaxp(3,2)=3
      line_cool(3,2)=c221+c222+c223
c      cinx(3,2,1)=cin(2)
c      taulinex(3,2,1)=t0
c      wlix(3,2,1)=wli(2)
c      ilabcfx(3,2,1)=302
C
C     C III 977 A Pradhan comp. MULTILEV!
C
      Z=ABN(3)*xion(ik,3,3)
      OM=4.34*T4**0.112
C III 977, cs from Berrington '85, J Phys B 18, L395.
      om = min( 7.0d0 , 1.556d0 * te**0.1 )
      CALL RLOSS(3,3,977.D0,OM,1.D0,3.D0,1.79D9,RS,XEL,Z,TE,cl(3)
     &     ,W(3),CIN(3),FR(2,3),3,WLI(3))
      tauline(3)=t0
      ilabcf(3)=303
      C232=cl(3)

C     CIII] 1909 (MENDOZA) MULTILEV!
C
      om=1.01/t4**0.1
      CALL RLOSS(3,3,1909.D0,om,1.D0,9.D0,4.03D1,RS,XEL,Z,TE,cl(4)
     &     ,W(4),CIN(4),FR(2,4),4,WLI(4))
      tauline(4)=t0
      ilabcf(4)=303
      C231=cl(4)

C
C     C IV 1546 A,(Pradhan compil.) 
C
      Z=ABN(3)*xion(ik,3,4)
      om=8.88d0*t4**0.0113
      CALL RLOSS(3,4,1550.D0,om,2.D0,6.D0,2.65D8,RS,XEL,Z,TE,cl(5)
     &     ,W(5),CIN(5),FR(2,5),5,WLI(5))
      tauline(5)=t0
      ilabcf(5)=304
      C241=cl(5)
      cinx(3,4,1)=cin(5)
      taulinex(3,4,1)=t0
      wlix(3,4,1)=wli(5)
      ilabcfx(3,4,1)=304
      kmaxp(3,4)=1
      line_cool(3,4)=cl(5)
C
C     NITROGEN
C
C
C     N II 2143 .OMEGA=1.29 (JACKSON 73)  ONLY GUESS OF A21 FROM OIII 16 MULTILEV!
C
      Z=ABN(4)*XN(2)
      CALL RLOSS(4,2,2143.D0,1.29D0,9.D0,5.D0,2.D2,RS,XEL,Z,TE,cl(10)
     &     ,W(10),CIN(10),FR(2,10),10,WLI(10))
      tauline(10)=t0
      ilabcf(10)=402
      C321=cl(10)
c!! done with 12 level atom
      c321=0.
      Z=ABN(4)*XN(3)
      om=3.657*t4**0.08
      CALL RLOSS(4,3,1750.D0,om,6.D0,12.D0,5.89D2,RS,XEL,Z,TE,cl(11)
     &     ,W(11),CIN(11),FR(2,11),11,WLI(11))
 9275 format('N III 1750 2level',f12.2,1pe12.3,10e12.3)
C

C
C     N III 990 A (MENDOZA)C ALPINE) OMEGA=5.22 G1=6 MULTILEV!
C
      Z=ABN(4)*XN(3)
      CALL RLOSS(4,3,990.D0,7.12D0,6.D0,10.D0,7.3D8,RS,XEL,Z,TE,cl(12)
     &     ,W(12),CIN(12),FR(2,12),12,WLI(12))
      tauline(12)=t0
      ilabcf(12)=403
      C332=cl(12)
c     Add diel. rec. from NS for 0.1 < T4 < 6
      if(t4.gt.0.1.and.t4.lt.6) then
         adeff=1.d-12*(2.8315/t4+12.9695+16.8995*t4-0.5167*t4**2)*
     &        exp(-0.8162/t4)/t4**1.5
      else
         adeff=0.
      endif
      rec990=abn(4)*xn(4)*adeff*2.01e-11
C
C     N V 1238 A and OMEGA from Pradhan compil
C
      Z=ABN(4)*XN(5)
      om=6.65*t4**0.015
      CALL RLOSS(4,5,1238.D0,om,2.D0,6.D0,3.40D8,RS,XEL,Z,TE,cl(14)
     &     ,W(14),CIN(14),FR(2,14),14,WLI(14))
      tauline(14)=t0
      ilabcf(14)=405
      C351=cl(14)
      cinx(4,5,1)=cin(14)
      taulinex(4,5,1)=t0
      wlix(4,5,1)=wli(14)
      ilabcfx(4,5,1)=405
      kmaxp(4,5)=1
      line_cool(4,5)=cl(14)
C
C     COLL. EXIT. OF OXYGEN
C

c     O I 7949. DIEL. REC. 
      if(te.gt.1.e3.and.te.lt.6.e4) then
         CALL DIELB(0.d0,0.0400d0,0.d0,0.d0,0.5587d0,T4,ALDB)
      else
         aldb=0.
      endif
      CIN(87)=ABN(5)*xion(ik,5,2)*ALDB*1.9864E-8/7949.
      WLI(87)=7949.
      ilabcf(87)=501
c     O I 6319. REC. 
      if(te.gt.1.e3.and.te.lt.6.e4) then
         CALL DIELB(0.d0,0.0280d0,0.d0,0.d0,1.0257d0,T4,ALDB)
      else
         aldb=0.
      endif
      CIN_6319=ABN(5)*xion(ik,5,2)*ALDB*1.9864E-8/6319.
      WLI_6319=6319.

c     O II 4651 REC. (BORKOWSKI AND SHULL 1990 AP J 348,169)
      CIN(84)=ABN(5)*xion(ik,5,3)*ALEX(15)*1.9864E-8/4651.*0.16
      WLI(84)=4651.
      ilabcf(84)=502
C
C     O III 1664  Pradhan compil. MULTILEV!
C
      Z=ABN(5)*XO(3)
      om=1.21*t4**0.058
      CALL RLOSS(5,3,1664.D0,om,9.D0,5.D0,4.06D2,RS,XEL,Z,TE,cl(9)
     &     ,W(9),CIN(9),FR(2,9),9,WLI(9))
      tauline(9)=t0
      ilabcf(9)=503
      C132=cl(9)
c     O III 3762. REC. (BORKOWSKI AND SHULL 1990 AP J 348,169)
      if(te.gt.1.e3.and.te.lt.6.e4) then
         CALL DIELB(5.30d-3,-9.830d-2,.4693d0,5.40d-3,0.5684d0,T4,ALDB)
      else
         aldb=0.
      endif
      CIN(85)=ABN(5)*xion(ik,5,4)*(ALEX(16)*0.17+ALDB)*1.9864E-8/3762.
      WLI(85)=3762.
      ilabcf(85)=503
c     O III 3266. DIEL. REC. 
      if(te.gt.1.e3.and.te.lt.6.e4) then
         CALL DIELB(1.d-4,.9493d0,.1777d0,.0623d0,1.2049d0,T4,ALDB)
      else
         aldb=0.
      endif
      CIN(89)=ABN(5)*xion(ik,5,4)*ALDB*1.9864E-8/3266.
      WLI(89)=3266.
      ilabcf(89)=503
C
C     O IV 1400 OMEGA=1.37 (O&W -77)
C
      Z=ABN(5)*XO(4)
      CALL RLOSS(5,4,1400.D0,1.37D0,6.D0,12.D0,1.6D3,RS,XEL,Z,TE,cl(8)
     &     ,W(8),CIN(8),FR(2,8),8,WLI(8))
      tauline(8)=t0
      ilabcf(8)=504
      C142=cl(8)
      cinx(5,4,1)=cin(8)
      taulinex(5,4,1)=t0
      wlix(5,4,1)=wli(8)
      ilabcfx(5,4,1)=504
C
C     O IV 789 A from Ferland Cl
C
      CALL RLOSS(5,4,789.d0,0.66D0,6.D0,10.D0,7.07d8,RS,XEL,Z,TE,cl(74)
     &     ,W(74),CIN(74),FRQQ,74,WLI(74))
      tauline(74)=t0
      ilabcf(74)=504
      C141=cl(74)
      cinx(5,4,2)=cin(74)
      taulinex(5,4,2)=t0
      wlix(5,4,2)=wli(74)
      ilabcfx(5,4,2)=504
C
C     O IV 25.88 MU MAZDA
C
      Z=ABN(5)*XO(4)
      CALL RLOSS(5,4,2.588D5,2.33D0,2.D0,4.D0,5.20D-4,RS,XEL,Z,TE,cl(83)
     &     ,W(83),CIN(83),FRQQ,83,WLI(83))
c      tauline(83)=t0
c      ilabcf(83)=504
c      C143=cl(83)
      cinx(5,4,3)=cin(83)
      taulinex(5,4,3)=t0
      wlix(5,4,3)=wli(83)
      ilabcfx(5,4,3)=504
c     O IV 3066. REC. (BORKOWSKI AND SHULL 1990 AP J 348,169)
      CIN(86)=ABN(5)*xion(ik,5,5)*ALEX(17)*1.9864E-8/3066.*0.16
      WLI(86)=3066.
      ilabcf(86)=504
      cinx(5,4,4)=cin(86)
      taulinex(5,4,4)=t0
      wlix(5,4,4)=wli(86)
      ilabcfx(5,4,4)=504
c     O IV 3027. DIEL. REC. 
      if(te.gt.1.e3.and.te.lt.6.e4) then
         CALL DIELB(0.0d0,1.5217d0,0.d0,0.d0,0.7085d0,T4,ALDB)
      else
         aldb=0.
      endif
      CIN(90)=ABN(5)*xion(ik,5,5)*ALDB*1.9864E-8/3027.
      WLI(90)=3027.
      ilabcf(90)=504
      cinx(5,4,5)=cin(90)
      taulinex(5,4,5)=t0
      wlix(5,4,5)=wli(90)
      ilabcfx(5,4,5)=504
      kmaxp(5,4)=5
      line_cool(5,4)=cl(8)+cl(74)+cl(83)
C
C     EXIT. OF O V (MENDOZA)
C
      Z=ABN(5)*XO(5)
      CALL RLOSS(5,5,1216.D0,0.721D0,1.D0,9.D0,2.25D3,RS,XEL,Z,TE,cl(7)
     &     ,W(7),CIN(7),FR(2,7),7,WLI(7))
      tauline(7)=t0
      ilabcf(7)=505
      C151=cl(7)
      cinx(5,5,1)=cin(7)
      taulinex(5,5,1)=t0
      wlix(5,5,1)=wli(7)
      ilabcfx(5,5,1)=505
C     O V 629.7 A Pradhan compil 95
      om=2.76*t4**0.046
      CALL RLOSS(5,5,629.7D0,om,1.D0,3.D0,2.80D9,RS,XEL,Z,TE,cl(18)
     &     ,W(18),CIN(18),FRqq,18,WLI(18))
      tauline(18)=t0
      ilabcf(18)=505
      C152=cl(18)
      cinx(5,5,2)=cin(18)
      taulinex(5,5,2)=t0
      wlix(5,5,2)=wli(18)
      ilabcfx(5,5,2)=505
c     O V 6488. DIEL. REC. 
      if(te.gt.1.e3.and.te.lt.6.e4) then
         CALL DIELB(-.0765d0,4.168d0,1.2966d0,-.1261d0,2.8234d0,T4,ALDB)
      else
         aldb=0.
      endif
      CIN(91)=ABN(5)*XO(6)*ALDB*1.9864E-8/6488.
      WLI(91)=6488.      
      ilabcf(91)=505
      cinx(5,5,3)=cin(91)
      taulinex(5,5,3)=t0
      wlix(5,5,3)=wli(91)
      ilabcfx(5,5,3)=505
      kmaxp(5,5)=3
      line_cool(5,5)=cl(7) + cl(18)
         
C     EXIT. OF O VI, A & Omega from Pradhan
C
      Z=ABN(5)*XO(6)
      om=5.00*t4**0.014
      CALL RLOSS(5,6,1034.D0,om,2.D0,6.D0,4.12D8,RS,XEL,Z,TE,cl(6)
     &     ,W(6),CIN(6),FR(2,6),6,WLI(6))
      tauline(6)=t0
      ilabcf(6)=506
      C161=cl(6)
      cinx(5,6,1)=cin(6)
      taulinex(5,6,1)=t0
      wlix(5,6,1)=wli(6)
      ilabcfx(5,6,1)=506
      kmaxp(5,6)=1
      line_cool(5,6)=cl(6)
C
C     NEON
C
C     NE II 12.814 MY Pradhan. Updated coll. str. Wang+ 2017
      Z=ABN(6)*XNE(2)
      OM=0.303*t4**0.065
      CALL RLOSS(6,2,12.814D4,OM,4.D0,2.D0,8.55D-3,RS,XEL,Z,TE,cl(46)
     &     ,W(46),CIN(46),FRQQ,46,WLI(46))
      tauline(46)=t0
      ilabcf(46)=602
      C1322=cl(46)
      cinx(6,2,1)=cin(46)
      taulinex(6,2,1)=t0
      wlix(6,2,1)=wli(46)
      ilabcfx(6,2,1)=602
      kmaxp(6,2)=1
      line_cool(6,2)=cl(46)

      
C
C     NE VII 895.1 Pradhan 95 compil. DONE with multilev
C
      Z=ABN(6)*XNE(7)
      if(t4.lt.4.) then
         om1=0.0352*te**0.17
      else
         om1=0.736/te**0.09
      endif
      OM=0.172*T4**0.41
      CALL RLOSS(6,7,895.1D0,OM,1.D0,3.D0,1.92D4,RS,XEL,Z,TE,cl(75)
     &     ,W(75),CIN(75),FR(2,75),75,WLI(75))
      tauline(75)=t0
      ilabcf(75)=607
c      line_cool(6,7)=cl(75)
      
c      C1371=cl(75)
C
C     SODIUM
C
C
C     NA I 5889 OM FROM KENNEDY ET AL J. PHYS. B 10, 3759
C
      Z=ABN(7)*XNA(1)
      OM=10.45*T4**0.837
      CALL RLOSS(7,1,5889.D0,OM,2.D0,6.D0,1.07D8,RS,XEL,Z,TE,cl(22)
     &     ,W(22),CIN(22),FR(2,22),22,WLI(22))
      tauline(22)=t0
      ilabcf(22)=701
      C1011=cl(22)
      cinx(7,1,1)=cin(22)
      taulinex(7,1,1)=t0
      wlix(7,1,1)=wli(22)
      ilabcfx(7,1,1)=701
C
C     NA I 5889 RECOMBINATION (ALL REC. COMES OUT IN 5889)
C
      CIN(41)=ABN(7)*XNA(2)*2.52E-13*3.37E-12/T4**.682
      wli(41)=5889.
      ilabcf(41)=701
      cinx(7,1,2)=cin(41)
      taulinex(7,1,2)=t0
      wlix(7,1,2)=5889
      ilabcfx(7,1,2)=701
      kmaxp(7,1)=2
      line_cool(7,1)=cl(22)
C
C     MAGNESIUM
C
C     MG I 5177 DIELECTRONIC REC. (NUSSBAUMERY&.STOREY 1985)
C
      if(te.gt.1.e3.and.te.lt.6.e4) then
         CALL DIELB(0.0d0,.1658d0,0.5114d0,-5.69d-2,2.1817d0,T4,ALDB)
      else
         aldb=0.
      endif

      CIN(28)=ABN(8)*XM(2)*ALDB*3.87E-12
      WLI(28)=5177.
      ilabcf(28)=801
      cinx(8,1,1)=cin(28)
      taulinex(8,1,1)=t0
      wlix(8,1,1)=wli(28)
      ilabcfx(8,1,1)=801
C
C     MG I 8806 DIELECTRONIC REC. (NUSSBAUMERY&.STOREY 1985)
C
      if(te.gt.1.e3.and.te.lt.6.e4) then
         CALL DIELB(.1259d0,-.6903d0,1.1785d0,-9.60d-2,2.0651d0,T4,ALDB)
      else
         aldb=0.
      endif
      CIN(29)=ABN(8)*XM(2)*ALDB*2.27E-12
      WLI(29)=8806.
      ilabcf(29)=801
      cinx(8,1,2)=cin(29)
      taulinex(8,1,2)=t0
      wlix(8,1,2)=wli(29)
      ilabcfx(8,1,2)=801
C
C     MG I 3834 DIELECTRONIC REC. (NUSSBAUMERY&.STOREY 1985)
C
      if(te.gt.1.e3.and.te.lt.6.e4) then
         CALL DIELB(.1283d0,-1.2733d0,4.6334d0,-4.73d-1,2.7691d0,
     &        T4,ALDB)
      else
         aldb=0.
      endif
      CIN(30)=ABN(8)*XM(2)*ALDB*5.22E-12
      WLI(30)=3834.
      ilabcf(30)=801
      cinx(8,1,3)=cin(30)
      taulinex(8,1,3)=t0
      wlix(8,1,3)=wli(30)
      ilabcfx(8,1,3)=801
C
C     MG I 4571 OM FROM FABRIKANT (J. PHYS. B 7:91 -74), AS GIVEN
C          BY OSTERBROCK (-74) P 245.
C               A FROM
C
      Z=ABN(8)*XM(1)
      OM=1.60*T4**0.56
      CALL RLOSS(8,1,4572.D0,OM,1.D0,9.D0,7.20D1,RS,XEL,Z,TE,cl(20)
     &     ,W(20),CIN(20),FR(2,20),20,WLI(20))
      tauline(20)=t0
      ilabcf(20)=801
      C511=cl(20)
      cinx(8,1,4)=cin(20)
      taulinex(8,1,4)=t0
      wlix(8,1,4)=wli(20)
      ilabcfx(8,1,4)=801

C
C     MG I 4572 RECOMBINATION (SEE NOTES). DIEL. FROM NUSSBAUMER & STORE
C
      if(te.gt.1.e3.and.te.lt.6.e4) then
         CALL DIELB(.5116d0,-2.8906d0,7.445d0,-.72340d0,2.414d0,T4,ALDB)
      else
         aldb=0.
      endif
      CIN(40)=ABN(8)*XM(2)*(1.96E-13/T4**0.855+ALDB)*4.34E-12
      wli(40)=4571.
      ilabcf(40)=801
      cinx(8,1,1)=cin(40)
      taulinex(8,1,5)=t0
      wlix(8,1,5)=wli(40)
      ilabcfx(8,1,5)=801
C
C     MG I 2852 OM FROM FABRIKANT (J. PHYS. B 7:91)
C
      Z=ABN(8)*XM(1)
      OM=2.11*T4**1.187
      CALL RLOSS(8,1,2852.D0,OM,1.D0,3.D0,4.93D8,RS,XEL,Z,TE,cl(21)
     &     ,W(21),CIN(21),FR(2,21),21,WLI(21))
      tauline(21)=t0
      ilabcf(21)=801
      C512=cl(21)
      cinx(8,1,6)=cin(21)
      taulinex(8,1,6)=t0
      wlix(8,1,6)=wli(21)
      ilabcfx(8,1,6)=801
      kmaxp(8,1)=6
      line_cool(8,1)=cl(20) + cl(21)
            
C
C     MG II 2800 OM (MENDOZA 83). Done in multlev
C
      Z=ABN(8)*XM(2)
      om=16.9*t4**0.1
      CALL RLOSS(8,2,2795.D0,om,2.D0,6.D0,2.55D8,RS,XEL,Z,TE,cl(15)
     &     ,W(15),CIN(15),FR(2,15),15,WLI(15))
      tauline(15)=t0
      ilabcf(15)=802
      C521=cl(15)
c      cinx(8,2,1)=cin(15)
c      taulinex(8,2,1)=t0
c      wlix(8,2,1)=wli(15)
c      ilabcfx(8,2,1)=802


C
C     Al II 1671 Pradhan compli.
C
      Z=ABN(9)*XAL(2)
      OM=3.251*T4**0.112
      CALL RLOSS(9,2,1671.D0,OM,1.D0,3.D0,1.46D9,RS,XEL,Z,TE,cl(42)
     &     ,W(42),CIN(42),FR(2,42),42,WLI(42))
c      tauline(42)=t0
c      ilabcf(42)=902
      C921=cl(42)
      cinx(9,2,1)=cin(42)
      taulinex(9,2,1)=t0
      wlix(9,2,1)=wli(42)
      ilabcfx(9,2,1)=902
C
C     AL II 2670. Pradhan comp.
C
      Z=ABN(9)*XAL(2)
      om=3.251*t4**0.25
      CALL RLOSS(9,2,2670.D0,om,1.D0,9.D0,1.11d3,RS,XEL,Z,TE,cl(81)
     &     ,W(81),CIN(81),FR(2,81),81,WLI(81))
c      tauline(81)=t0
c      ilabcf(81)=902
      C922=cl(81)
      cinx(9,2,2)=cin(81)
      taulinex(9,2,2)=t0
      wlix(9,2,2)=wli(81)
      ilabcfx(9,2,2)=902
      line_cool(9,2)=cl(42) + cl(81)

C
C     AL III 1854.7 1862.8 OM and A Kingdon 
C
      Z=ABN(9)*XAL(3)
      om=16.0*t4**0.1
      CALL RLOSS(9,3,1858.D0,om,2.D0,6.D0,5.55D8,RS,XEL,Z,TE,cl(96)
     &     ,W(96),CIN(96),FR(2,96),96,WLI(96))
c      tauline(96)=t0
c      ilabcf(96)=903
      C931=cl(96)
      cinx(9,3,1)=cin(96)
      taulinex(9,3,1)=t0
      wlix(9,3,1)=wli(96)
      ilabcfx(9,3,1)=903
      line_cool(9,3)=cl(96)

C
C     SILICON
C

C
C     SI II] 2335 (MENDOZA 83)
C
      Z=ABN(10)*xion(ik,10,2)
      CALL RLOSS(10,2,2335.D0,5.14D0,6.D0,12.D0,4.14D1,RS,XEL,Z,TE,
     &     cl(17),W(17),CIN(17),FR(2,17),17,WLI(17))
      tauline(17)=t0
      ilabcf(17)=1002
      C421=cl(17)
      cinx(10,2,1)=cin(17)
      taulinex(10,2,1)=t0
      wlix(10,2,1)=wli(17)
      ilabcfx(10,2,1)=1002
C     SI II 34.81 MY A FROM(MENDOZA)
C     omega from Keenan et al MNRAS 214, 37p, 1985.
      Z=ABN(10)*XSI(2)
      OM=5.6
      CALL RLOSS(10,2,3.481D5,OM,2.D0,4.D0,2.17D-4,RS,XEL,Z,TE,cl(44)
     &     ,W(44),CIN(44),FRQQ,44,WLI(44))
      tauline(44)=t0
      ilabcf(44)=1002
      C422=cl(44)
      cinx(10,2,2)=cin(44)
      taulinex(10,2,2)=t0
      wlix(10,2,2)=wli(44)
      ilabcfx(10,2,2)=1002
      kmaxp(10,2)=2
      line_cool(10,2)=cl(17)+cl(44)
C
C     SI III 1892 OMEGA FROM MENDOZA 83
C
      Z=ABN(10)*xion(ik,10,3)
      OM=5.43/T4**0.34
      CALL RLOSS(10,3,1892.D0,OM,1.D0,9.D0,5.5D3,RS,XEL,Z,TE,cl(73)
     &     ,W(73),CIN(73),FR(2,73),73,WLI(73))
      tauline(73)=t0
      ilabcf(73)=1003
      C431=cl(73)
      cinx(10,3,1)=cin(73)
      taulinex(10,3,1)=t0
      wlix(10,3,1)=wli(73)
      ilabcfx(10,3,1)=1003
      kmaxp(10,3)=1
      line_cool(10,3)=cl(73)
C
C     SI IV 1400 OMEGA from Pradhan comp
C
      Z=ABN(10)*xion(ik,10,4)
      CALL RLOSS(10,4,1400.D0,16.D0,2.D0,6.D0,7.7D8,RS,XEL,Z,TE,cl(16)
     &     ,W(16),CIN(16),FR(2,16),16,WLI(16))
      tauline(16)=t0
      ilabcf(16)=1004
      C441=cl(16)
      cinx(10,4,1)=cin(16)
      taulinex(10,4,1)=t0
      wlix(10,4,1)=wli(16)
      ilabcfx(10,4,1)=1004
      kmaxp(10,4)=1
      line_cool(10,4)=cl(16)
C
C     SULPHUR
C
C     S IV 1393.4    (G & S) in multi
C
      Z=ABN(11)*xion(ik,11,4)
      OM=6.99D0*T4**(-0.129)
      CALL RLOSS(11,4,1393.4D0,OM,6.D0,12.D0,2.06D4,RS,XEL,Z,TE,cl(80)
     &     ,W(80),CIN(80),FR(2,80),80,WLI(80))
      tauline(80)=t0
      ilabcf(80)=1104
      C1242=cl(80)
      
C     S IV 1069.6    (G & S) in multi
C
      OM=5.24D0*T4**0.028
      CALL RLOSS(11,4,1069.6D0,OM,6.D0,10.D0,1.31D8,RS,XEL,Z,TE,cl(78)
     &     ,W(78),CIN(78),FR(2,78),78,WLI(78))
      tauline(78)=t0
      ilabcf(78)=1104
      C1241=cl(78)
      line_cool(11,4)=cl(78)+cl(80)
C
C     S V 786.5 Pradhan 95 compil. In multi
C
      Z=ABN(11)*xion(ik,11,5)
      CALL RLOSS(11,5,786.5D0,7.30d0,1.d0,3.d0,5.25d9,RS,XEL,Z,TE,cl(76)
     &     ,W(76),CIN(76),FR(2,76),76,WLI(76))
      tauline(76)=t0
      ilabcf(76)=1105
      C1251=cl(76)
      line_cool(11,5)=cl(76)
C                                            
C     S VI 937.1  Pradhan
C
      Z=ABN(11)*xion(ik,11,6)
      OM=11.9d0
      CALL RLOSS(11,6,937.1D0,OM,2.D0,6.D0,1.61D9,RS,XEL,Z,TE,cl(77)
     &     ,W(77),CIN(77),FR(2,77),77,WLI(77))
      tauline(77)=t0
      ilabcf(77)=1106
      C1261=cl(77)
      line_cool(11,6)=cl(77)

C     AR II 6.985 MY from Pelan & Berrington 1995
      telog=log10(te)
      if(te < ar2_te(1)) then
         al_ar2=log10(ar2_coll(2)/ar2_coll(1))/log10(ar2_te(2)/ar2_te(1))
         om=ar2_coll(1)*(te/ar2_te(1))**al_ar2
         ii=1
      elseif(te > ar2_te(11)) then
         al_ar2=log10(ar2_coll(11)/ar2_coll(10))/log10(ar2_te(11)/ar2_te(10))
         om=ar2_coll(11)*(te/ar2_te(11))**al_ar2
         ii=10
      else
         do i=1,10
            if(te > ar2_te(i) .and. te <= ar2_te(i+1)) then
               al_ar2=log10(ar2_coll(i+1)/ar2_coll(i))/log10(ar2_te(i+1)/ar2_te(i))
               om=ar2_coll(i)*(te/ar2_te(i))**al_ar2
               ii=i
            endif
         enddo
      endif
      Z=ABN(12)*XAR(2)
      CALL RLOSS(12,2,6.985D4,OM,4.D0,2.D0,5.27D-2,RS,XEL,Z,TE,cl(45)
     &     ,W(45),CIN(45),FRQQ,45,WLI(45))
      tauline(45)=t0
      ilabcf(45)=1202
      C1222=cl(45)
      cinx(12,2,1)=cin(45)
      taulinex(12,2,1)=t0
      wlix(12,2,1)=wli(45)
      ilabcfx(12,2,1)=1202
      kmaxp(12,2)=1
      line_cool(12,2)=cl(45)
C     ************************************************************
C     *****
c   Ca IV 3 level atom up to 3p6 2S     
C     *****
C     *************************************************************

c     interpolate coll strength
      if(te.lt.teca4(1)) then
         o21=coll_ca4(1,1)
         o31=coll_ca4(2,1)
         o32=coll_ca4(3,1)
      elseif(te.ge.teca4(50)) then
         o21=coll_ca4(1,50)
         o31=coll_ca4(2,50)
         o32=coll_ca4(3,50)
      else
         do kc=1,49
            if(te>=teca4(50)) then
               omca4(k)=coll_ca4(k,50)
            elseif(te<=teca4(1)) then
               omca4(k)=coll_ca4(k,1)
            elseif(te>=teca4(kc) .and. te<=teca4(kc+1)) then
               do k=1,3
                  omca4(k)=coll_ca4(k,kc) + (te-teca4(kc))*
     &                 (coll_ca4(k,kc+1)-coll_ca4(k,kc))/(teca4(kc+1)-teca4(kc))
               enddo
               om21=omca4(1)
               om31=omca4(2)
               om32=omca4(3)
            endif
         enddo
      endif
      Z=ABn(13)*xca(4)
      e21=0.38661d0
      e32=18.90010d0-e21
      CALL FORB3(21,12,4,e21,e32,0.545d0,5.11d8,2.466d8,om21,om31,om32, 
     &     4.D0,2.D0,2.D0,FB21,FB31,FB32,TE,Z,F)      
      DO NM=64,66                                                     
         FB(NM)=Z*FB(NM)
      enddo
      do k4=1,3         
         if(k4==1) then
            cinx(13,4,k4)=z*fb21
            wlix(13,4,k4)=wlin(1)
         elseif(k4==2) then
            cinx(13,4,k4)=z*fb31                     
            wlix(13,4,k4)=wlin(2)
         elseif(k4==3) then
            cinx(13,4,k4)=z*fb32            
            wlix(13,4,k4)=wlin(3)
         endif
         ilabcfx(13,4,k4)=1304
      enddo
      kmaxp(13,4)=3
      Cocaiv=Z*F
      cocal=cocaiv
      line_cool(13,4)=cocaiv

c Ca II
      ilcf=1302
      z=abn(13)*xca(2)
      if(te <= 1.e5) then
         call popchianti_new(13,2,te,coca2)
         do k=1,kmaxp(13,2)
            if(ipopch.eq.1) then
               cinx(13,2,k)=weh(k)/del(ik)
               taulinex(13,2,k)=taul(k)
               wlix(13,2,k)=wlin(k)
            endif
            ilabcfx(13,2,k)=1302
         enddo         
         line_cool(13,2)=coca2
      else
         do k=1,kmaxp(13,2)
            if(ipopch.eq.1) then
               cinx(13,2,k)=0.
               taulinex(13,2,k)=taul(k)
               wlix(13,2,k)=wlin(k)
            endif
            ilabcfx(13,2,k)=1302
         enddo         
         line_cool(13,2)=0.
      endif
      
C     
C     Ca V
C
      om21=2.8e0
      om31=0.66e0
      om32=0.996e0
      om1d3p=3.1e0
      om43=1.*om1d3p/9.
      om42=3.*om1d3p/9.
      om41=5.*om1d3p/9.
      om1s3p=0.147e0
      om53=1.*om1s3p/9.
      om52=3.*om1s3p/9.
      om51=5.*om1s3p/9.
      om54=1.09e0
      call upsilon_chianti(0,17,13,5,5,te)
      om21=upsil(17,1,2)
      om31=upsil(17,1,3)
      om41=upsil(17,1,4)
      om51=upsil(17,1,5)
      om32=upsil(17,2,3)
      om42=upsil(17,2,4)
      om52=upsil(17,2,5)
      om43=upsil(17,3,4)
      om53=upsil(17,3,5)
      om54=upsil(17,4,5)
      z=abn(13)*xca(5)

      CALL FIVELEV_dp(113,13,5,5,0.d0,2404.7d0,3275.6d0,18830.3d0,
     &     43836.5d0,5.d0,3.d0,1.d0,5.d0,1.d0,
     &     3.079d-01,3.960d-05,1.955d+00,1.224d-01,3.681d-02,4.356d-01,
     &     2.345d+01,7.986d-05,0.d0, 3.697d+00,
     &     om21,om31,om41,om51,om32,om42,om52,
     &     om43,om53,om54,TE,Z,XI)
      DO I=1,5
        DO J=1,5
         EMCA(I,J)=EMCA(I,J)/DEL(IK)
        ENDDO
      ENDDO
      coca5=0.
      DO  I=1,5
        DO  J=I+1,5
           COca5=COca5+EMCA(J,I)
        ENDDO
      ENDDO
      do k=1,kmaxp(13,5)
         wlix(13,5,k)=wlin(k)
         cinx(13,5,k)=weh(k)/(del(ik))
         ilabcfx(13,5,k)=1305
      enddo
      line_cool(13,5)=coca5
      cocal=cocal+coca5
           
C
C     FE I FIT TO AXELROD  AT N=1.E7
C
      W18=ABN(14)*XF(1)*5.1E-19*EXP(-2.83/T4)/SQRT(T4)
C     FINE STRUCTURE
      FEIFS=1.15E-23*EXP(-0.6499/T3)/(T3**0.358*(DEN*XEL/1.E7)**
     &     0.96)

c normalize to the last convergent calc. of Fe I in multi

c      FEIFS=EXP(-0.6499/T3)/(T3**0.358*(DEN*XEL/1.E7)**
c     &     0.96)

      FEIFS=FEIFS*ABN(14)*XF(1)
c      WLI(18)=-1.
c      ilabcf(18)=1401
      CIN18=fecon*(W18+FEIFS)
      feian=(W18+FEIFS)
c      cl(18)=cin(18)
C
C     FE II FIT TO AXELROD AT N=1.E7
C
      W19=ABN(14)*XF(2)*3.7E-19*EXP(-2.91/T4)/SQRT(T4)
C     FINE STRUCTURE
      FEIIFS=1.03E-23*EXP(-0.5077/T3)/(T3**0.162*(DEN*XEL/1.E7)**
     &     0.96)
      FEIIFS=FEIIFS*ABN(14)*XF(2)
c      WLI(19)=-1.
c      ilabcf(19)=1402
      CIN19=W19+FEIIFS

c     cl(19)=cin(19)

C     FE III Done as multi!
      feiii=0.
      cin42=0.
      IF(ABN(14)*XF(3).GT.1.E-6) THEN
         denel=den*xel
         z=abn(14)*xf(3)
         if(ilowion.eq.0) call feiiiatom(te,denel,z,rtot)
         if(ilowion.eq.1) then
C     
C     FE III FIT TO AXELROD AT N=1.E7
C     
            CIN42=ABN(14)*XF(3)*3.7E-19*EXP(-31.82/T3)/T3**0.634
            FEIIIFS=2.92E-21*EXP(-0.6313/T3)/T3**0.0313
            FEIIIFS=FEIIIFS/(1.+(DEN*XEL)/4.04E4)
            FEIIIFS=FEIIIFS*ABN(14)*XF(3)
            CIN42=CIN42+FEIIIFS
         endif
         feiii=rtot*abn(14)*xf(3)
      endif
c      WLI(42)=-1.
c      ilabcf(42)=1403
c      cl(42)=feiii
C     FE IV
      denel=den*xel
cfqqq
      rtot=1.80e-19*exp(-4.626e4/te)/(te**0.079*(1.+denel/5.675e7)**.6)
cfqqq
      z=abn(14)*xf(4)
c!! done as multi!      
      IF(ABN(14)*XF(4).GT.1.E-6.AND.ilowion.eq.0) 
     &                  call feivatom(te,denel,z,rtot)
      feiv=rtot*abn(14)*xf(4)

c!!!  use same rate for Fe V
      fev=rtot*abn(14)*xf(5)
      goto 4565
      REWIND 27
      read(27,1921)lab
      read(27,1921)lab
1921   format(a)
      iout=12
      do i=1,iout
         read(27,912)ix,ex,g(i),lev
 912     format(i3,f12.2,f7.1,a11)
      enddo
      fe4=0.
      read(27,1921)lab
      read(27,1921)lab
      do i=2,iout
         im1=i-1
         do j=1,im1
            read(27,*)iq,jq,dlair,aij,omega
            hny=12398.54/dlair
            hnykt=hny*1.1604e4/te
            emfe4(i,j)=8.63e-6*omega*1.602e-12*hny*exp(-hnykt)/
     &           (sqrt(te)*g(j))
            emfe4(i,j)=abn(14)*xf(4)*emfe4(i,j)
            fe4=fe4+emfe4(i,j)
         enddo
      enddo
4565  continue
      DO  IW=1,14
         CEX(IW)=0.
      enddo
c      IF(TE.LT.0.1E4) GOTO 5005

      TEV=TE/1.1609E4

      idest=1

      nline=96

      ll=nline

C
C     O I 3P2, 3P1, 3P0, 1D, 1S Done as multi in popo
C
      om21=0.0987*t4**1.063
      om31=0.0292*t4**0.928
      om32=0.0265*t4**1.315
      om1d3p=0.266*t4**1.007
      om41=5.*om1d3p/9.
      om42=3.*om1d3p/9.
      om43=1.*om1d3p/9.
      om1s3p=0.0324*t4**0.994
      om51=5.*om1s3p/9.
      om52=3.*om1s3p/9.
      om53=1.*om1s3p/9.
      om54=0.105*t4**0.507
      Z=ABN(5)*XO(1)

      exen=3.16e-12
      texi=2.27e4
      fbapp6300=8.83e-6*abn(5)*xo(1)*exen*om1d3p*exp(-texi/te)/
     &     (9.*te**0.5)
      CALL FIVELEV_dp(14,5,1,5,0.d0,158.265d0,226.977d0,15867.86d0,33792.58d0,5.d0,3.d0,1.d0,
     &     5.d0,1.d0,8.92d-5,1.0d-10,6.34d-3,2.88d-4,1.74d-5,6.74d-3,
     &     7.60d-2,8.92d-7,0.0d-0,1.22d0,om21,om31,om41,om51,om32,
     &     om42,om52,om43,om53,om54,TE,Z,XI)      

      DO I=1,5
        DO J=1,5
         EMCA(I,J)=EMCA(I,J)/DEL(IK)
        ENDDO
      ENDDO
      C111b=0.
      DO  I=1,5
        DO  J=I+1,5
          C111b=C111b+EMCA(J,I)
        ENDDO
      ENDDO

c 6300-64
      FB(10)=EMCA(4,1)+EMCA(4,2)
c 2959
      FB(11)=EMCA(5,2)
c 5577
      FB(12)=EMCA(5,4)
c 44.1 mu
      FB(17)=EMCA(3,1)
c 63.2 mu
      FB(18)=EMCA(2,1)
c 146 mu
      FB(19)=EMCA(3,2)

C     (O III) 2321,4363,4959 & 5007 MENDOZA
      Z=ABN(5)*XO(3)
C
C     O III 3P,1D,1S OMEGAS from Burke et al MN 236,353 and Mazda (FS)
C
      om21=0.545*t4**0.047
      om31=0.271*t4**0.093
      om32=1.29*t4**0.0655
      om1d3p=2.29*t4**0.132
      om41=1.*om1d3p/9.
      om42=3.*om1d3p/9.
      om43=5.*om1d3p/9.
      om1s3p=0.293*t4**0.167
      om51=1.*om1s3p/9.
      om52=3.*om1s3p/9.
      om53=5.*om1s3p/9.
      om54=0.582*t4**0.0677

      exen=4.0e-12
      texi=2.895e4
      fbapp5007=8.83e-6*abn(5)*xo(3)*exen*om1d3p*exp(-texi/te)/
     &     (9.*te**0.5)

      CALL FIVELEV_dp(16,5,3,5,0.d0,113.2d0,306.2d0,20273.3d0,43185.7d0,1.d0,3.d0,5.d0,5.d0,1.d0,
     &     2.62d-5,3.02d-11,2.74d-6,0.d0,9.69d-5,6.74d-3,0.327d0,
     &     1.96d-2,7.85d-4,2.65d0,om21,om31,om41,om51,om32,om42,om52,
     &     om43,om53,om54,TE,Z,XI)
c     POPULATION OF EXCITED STATES
      XO31D=XI(4)
C     TOTAL O III 1D POPULATION
      xionoiii1d=XI(4)*xion(ik,5,3)
      DO I=1,5
        DO J=1,5
         EMCA(I,J)=EMCA(I,J)/DEL(IK)
        ENDDO
      ENDDO
      COOIII=0.
      DO  I=1,5
        DO  J=I+1,5
          COOIII=COOIII+EMCA(J,I)
        ENDDO
      ENDDO
      C131=COOIII
      line_cool(5,3)=cooiii
      FB(1)=EMCA(4,2)+EMCA(4,3)
      FB(2)=EMCA(5,2)
      FB(3)=EMCA(5,4)
      FB(49)=EMCA(2,1)
      FB(50)=EMCA(3,2)
      em4959=EMCA(4,2)
      em5007=EMCA(4,3)
      if(itime.eq.1) nline_min(5,3)=ll
      ll=nline_min(5,3)
      do k=1,kmaxp(5,3)
         Ll=ll+1
         cin(ll)=weh(k)/del(ik)
         tauline(ll)=taul(k)
         cinx(5,3,k)=weh(k)/del(ik)
         taulinex(5,3,k)=taul(k)
         wlix(5,3,k)=wlin(k)
         ilabcfx(5,3,k)=503         
         if(itime.le.888888) wli(ll)=wlin(k)
         ilabcf(ll)=503
      enddo

C                                                             
C       Ne III       (1814.6, 3342, 3869-3968) (Pradhan & Peng) 
C                                                           
      Z=ABN(6)*XNE(3)               

c
c     Coll. strengths from Butler & Zeippen 1994, A:s from Pradhan comp.
c    updated 230601 Mao et al 2021,McLaughlin+ 2011      
c
      call FIVELEV_dp(74,6,3,5,
     &0.0d0,642.9d0,920.4d0,25840.8d0,55750.6d0,
     &5.d0,3.d0,1.d0,1.d0,1.d0,
     &     5.97d-3,2.18d-8,1.71d-1,3.94d-3,1.15d-3,5.42d-2,
     &     2.06d0,8.51d-6,0.0d0,0.271d0,
     &     0.60d0,0.170d0,0.55d0,0.084d0,0.167d0,0.32d0,
     &     0.050d0,0.151d0,0.017d0,0.55d0,
     &     TE,Z,XI)

      DO I=1,5
        DO J=1,5
         EMCA(I,J)=EMCA(I,J)/DEL(IK)
        ENDDO
      ENDDO
      CONEIII=0.
      DO  I=1,5
        DO  J=I+1,5
          CONEIII=CONEIII+EMCA(J,I)
        ENDDO
      ENDDO
      line_cool(6,3)=coneiii
      if(itime.eq.1) nline_min(6,3)=ll
      ll=nline_min(6,3)
      do k=1,kmaxp(6,3)
         ll=ll+1
         cin(ll)=weh(k)/del(ik)
         tauline(ll)=taul(k)
         if(itime.le.888888) wli(ll)=wlin(k)
         ilabcf(ll)=603
         cinx(6,3,k)=weh(k)/del(ik)
         taulinex(6,3,k)=taul(k)
         wlix(6,3,k)=wlin(k)
         ilabcfx(6,3,k)=603
      enddo
   
      FB(37)=EMCA(4,2)+EMCA(4,1)
      FB(38)=EMCA(5,2)
      FB(39)=EMCA(5,4)
      FB(46)=EMCA(2,1)
      FB(47)=EMCA(3,1)
      FB(48)=EMCA(3,2)
      C1332=CONEIII

C                                                                   
C       Ne V    (1575, 2975, 3346-3426)                                
C                                                                  
      Z=ABN(6)*XNE(5)                                     

      DO NM=43,45                       
         FB(NM)=Z*FB(NM)
      enddo
      CALL FIVELEV_dp(76,6,5,5,0.d0,411.227d0,1109.5d0,30290.7d0,63915.4d0,1.d0,3.d0,5.d0,5.d0,1.d0,
     &     1.28d-3,5.08d-9,2.37d-5,0.d00,4.59d-3,1.31d-1,4.21d0,
     &     0.365d0,6.69d-3,2.85d0,
     &     1.41d0,1.81d0,0.232d0,0.027d0,5.83d0,0.695d0,0.082d0,1.159d0,
     &     0.137d0,0.58d0,
     &     TE,Z,XI)
      DO I=1,5
        DO J=1,5
         EMCA(I,J)=EMCA(I,J)/DEL(IK)
        ENDDO
      ENDDO
      CONEV=0.
      DO  I=1,5
        DO  J=I+1,5
          CONEV=CONEV+EMCA(J,I)
        ENDDO
      ENDDO
      FB(45)=EMCA(4,2)+EMCA(4,3)
      FB(43)=EMCA(5,2)
      FB(44)=EMCA(5,4)
      C1351=CONEV
      line_cool(6,5)=conev
      if(itime.eq.1) nline_min(6,5)=ll
      ll=nline_min(6,5)
      do k=1,kmaxp(6,5)
         ll=ll+1
         cin(ll)=weh(k)/del(ik)
         tauline(ll)=taul(k)
         if(itime.le.888888) wli(ll)=wlin(k)
         ilabcf(ll)=605
         cinx(6,5,k)=weh(k)/del(ik)
         taulinex(6,5,k)=taul(k)
         wlix(6,5,k)=wlin(k)
         ilabcfx(6,5,k)=605
      enddo
      
c 3 level atoms!  should all be replaced with 5 level!

      
C     (N I ) 3468, 5201, 10406 (MENDOZA 83)
      O21=(0.48*T4)
      O31=(0.17*T4)
      O32=(0.62*T4)
      Z=ABN(4)*XN(1)
      
      CALL FORB3(7,4,1,2.38396D0,1.19164D0,1.24D-5,5.28D-3,8.85D-2,O21,
     &     O31,O32,4.D0,10.D0,6.D0,FB(13),FB(14),FB(15),TE,Z,F)
      C311=ABN(4)*XN(1)*F

      if(itime.eq.1) nline_min(4,1)=ll
      ll=nline_min(4,1)
      do k=1,3
         ll=ll+1
         cin(ll)=weh(k)
         tauline(ll)=taul(k)
         if(itime.le.888888) wli(ll)=wlin(k)
         ilabcf(ll)=401
         cinx(4,1,k)=weh(k)
         taulinex(4,1,k)=taul(k)
         wlix(4,1,k)=wlin(k)
         ilabcfx(4,1,k)=401
         C212=cl(27)
      enddo
      line_cool(4,1)=c212
      DO NM=13,15
         FB(NM)=ABN(4)*XN(1)*FB(NM)
      enddo

      c322=0.

      DO NM=7,9
         FB(NM)=ABN(4)*XN(2)*FB(NM)
      enddo
C     (O I ) 2964, 5581, 6302-64 (MENDOZA 83)
      IN=N
      TEV=TE/1.1609E4
      IF(T4.GT.0.5) GOTO 10
      O21=.0151*T3**1.31
      O31=.00184*T3**1.32
      O32=.031*T3**.534
      GOTO 50
 10   IF(T4.GT.1.) GOTO 20
      O21=.266*T4**1.10
      O31=.034*T4**1.08
      O32=.105*T4**.52
      GOTO 50
 20   O21=.266*T4**.91
      O31=.0324*T4**.91
      O32=.105*T4**.50
 50   CONTINUE
      Z=ABN(5)*XO(1)



C     O I 63.18, 145.5 MY MAZDA
      O21=(.0018*(TE/1000.)**1.169)
      IF(TE.LE.1.E3) O31=(.0022*(TE/1000.)**1.874)
      IF(TE.GT.1.E3) O31=(.0292*(TE/10000.)**1.123)
      IF(TE.GT.1.E3) O32=(.0987*(TE/10000.)**1.113)
      IF(TE.LE.1.E3) O32=(.0076*(TE/1000.)**1.493)
      Z=ABN(5)*XO(1)

      DO NM=17,19         
         FB(NM)=ABN(5)*XO(1)*FB(NM)
      enddo

C     (O II) 2471,3727,7327 (MENDOZA 83) Done as multi!
      Z=ABN(5)*XO(2)
      CALL FORB3(3,5,2,3.32D0,1.69D0,8.89D-5,4.53D-2,1.73D-1,1.335D0,
     &0.405D0,1.30D0,4.D0,10.D0,6.D0,FB(4),FB(5),FB(6),TE,Z,F)
      C121=ABN(5)*XO(2)*F
      if(itime.eq.1) nline_min(5,2)=ll
      ll=nline_min(5,2)
      do k=1,3
         ll=ll+1
         cin(ll)=weh(k)
         tauline(ll)=taul(k)
         if(itime.le.888888) wli(ll)=wlin(k)
         ilabcf(ll)=502
      enddo


      DO NM=4,6
         FB(NM)=ABN(5)*XO(2)*FB(NM)
      enddo

C                                                                       
C       Ne IV    (1602, 2423-25, 4714-24)                                 
C                                                                       
      Z=ABN(6)*XNE(4) 
      CALL FORB3(11,6,4,5.12D0,2.63D0,1.42D-3,1.12D0,1.04D0,1.402D0
     &,0.469D0,0.98D0,4.D0,10.D0,6.D0,FB(40),FB(41),FB(42),TE,Z,F)  
      if(itime.eq.1) nline_min(6,4)=ll
      ll=nline_min(6,4)
      do k=1,3
         ll=ll+1
         cin(ll)=weh(k)
         tauline(ll)=taul(k)
         if(itime.le.888888) wli(ll)=wlin(k)
         ilabcf(ll)=604
         cinx(6,4,k)=weh(k)
         taulinex(6,4,k)=taul(k)
         wlix(6,4,k)=wlin(k)
         ilabcfx(6,4,k)=604
      enddo

      DO NM=40,42        
         FB(NM)=Z*FB(NM)
      enddo
      C1341=Z*F
      line_cool(6,4)=c1341      


C     (SI I) 6591, 10995, 16360 A. MENDOZA. OMEGA FOR 3P - 1D (1.64 mu) FROM
C     M.S. PINDZOLA, A.K. BHATIA, A.TEMKIN, PHYS.REV. 15, 35 (1977)
C     REST SCALED FROM C I WITH SCALING = 3
      O21=3.97*T4**0.9
      O31=3.*0.149*(TE/5000.)**0.871
      O32=3.*0.196*(TE/5000.)**0.499
      Z=ABN(10)*XSI(1)
      CALL FORB3(5,10,1,0.7580D0,1.1277D0,3.043D-3,3.22D-2,1.14D0,O21
     &     ,O31,O32,9.D0,5.D0,1.D0,FB(30),FBQQ,FB(31),TE,Z,F)
      C1011=ABN(10)*XSI(1)*F

      if(itime.eq.1) nline_min(10,1)=ll
      ll=nline_min(10,1)
      do k=1,3
         ll=ll+1
         cin(ll)=weh(k)
         tauline(ll)=taul(k)
         if(itime.le.888888) wli(ll)=wlin(k)
         ilabcf(ll)=1001
         cinx(10,1,k)=weh(k)
         taulinex(10,1,k)=taul(k)
         wlix(10,1,k)=wlin(k)
         ilabcfx(10,1,k)=1001
         kx=k
      enddo

      DO NM=30,31
         FB(NM)=ABN(10)*XSI(1)*FB(NM)
      enddo
C     SI I 129.68 68.474 MY A:s FROM MAZDA OMEGAS GUESSED!! FROM O I
      O32=9.76E-6*(TE-228)+3.46E-11*(TE-228.)**2
      O31=3.39E-6*(TE-326.)-2.90E-11*(TE-326.)**2
      O21=1.89E-6*(TE-326.)+8.00E-11*(TE-326.)**2
      IF(O21.LT.0.) O21=1.e-10
      IF(O31.LT.0.) O31=1.e-10
      IF(O32.LT.0.) O32=1.e-10
      Z=ABN(10)*XSI(1)
      CALL FORB3(15,10,1,9.559D-3,1.810D-2,8.25D-6,0.D0,4.21D-5,O21,O31,
     &     O32,1.D0,3.D0,5.D0,FB(28),FBQQ,FB(29),TE,Z,F)

      do k=1,3
         ll=ll+1
         cin(ll)=weh(k)
         tauline(ll)=taul(k)
         if(itime.le.888888) wli(ll)=wlin(k)
         ilabcf(ll)=1001
         cinx(10,1,kx+1+k)=weh(k)
         taulinex(10,1,kx+1+k)=taul(k)
         wlix(10,1,kx+1+k)=wlin(k)
         ilabcfx(10,1,kx+1+k)=1001
      enddo

      
      DO NM=28,29
         FB(NM)=ABN(10)*XSI(1)*FB(NM)
      enddo
      C1012=ABN(10)*XSI(1)*F
      line_cool(10,1)=c1011+c1012       

C     S I 25.25 , 56.31 MY A:s FROM MAZDA OMEGAS GUESSED!! FROM O I
      O32=9.76E-6*(TE-228)+3.46E-11*(TE-228.)**2
      O31=3.39E-6*(TE-326.)-2.90E-11*(TE-326.)**2
      O21=1.89E-6*(TE-326.)+8.00E-11*(TE-326.)**2
      IF(O21.LT.0.) O21=1.e-10
      IF(O31.LT.0.) O31=1.e-10
      IF(O32.LT.0.) O32=1.e-10
      Z=ABN(11)*XSU(1)
      CALL FORB3(16,11,1,4.912D-2,2.202D-2,1.39D-3,6.71D-8,3.02D-4,O21,
     &O31,O32,5.D0,3.D0,1.D0,FB(26),FBQQ,FB(27),TE,Z,F)
      if(itime.eq.1) nline_min(11,1)=ll
      ll=nline_min(11,1)
      do k=1,3
         ll=ll+1
         cin(ll)=weh(k)
         tauline(ll)=taul(k)
         if(itime.le.888888) wli(ll)=wlin(k)
         ilabcf(ll)=1101
         cinx(11,1,k)=weh(k)
         taulinex(11,1,k)=taul(k)
         wlix(11,1,k)=wlin(k)
         ilabcfx(11,1,k)=1101
      enddo


      DO NM=26,27
         FB(NM)=ABN(11)*XSU(1)*FB(NM)
      enddo
      C111=ABN(11)*XSU(1)*F
      line_cool(11,1)=c111

C     (S II) 6628,4071,10320 A. (MENDOZA 83) Done as multi!
      O21=6.98/T4**0.065
      O31=2.28/T4**0.076
      O32=14.63/T4**0.015
      Z=ABN(11)*XSU(2)
      CALL FORB3(6,11,2,1.8427D0,1.201D0,5.10D-4,1.80D-1,0.288D0,O21
     &,O31,O32,4.D0,10.D0,6.D0,FB(23),FB(24),FB(25),TE,Z,F)
      do k=1,3
         ilabcf(ll)=1102
      enddo



      DO NM=23,25
         FB(NM)=ABN(11)*XSU(2)*FB(NM)
      enddo
      C1221=ABN(11)*XSU(2)*F

C     S II 314.5 MY (MENDOZA) Done as multi!
      Z=ABN(11)*XSU(2)
      OM=7.50/T4**0.103
      ilabcf(43)=1102

C
C      S III   (3722, 6312, 9069 - 9532) MENDOZA done by multilevel atom!
C                                                                       
      Z=ABN(11)*XSU(3)                                                  
      O21=DMAX1(8.20D0,8.39D0*T4**(-.073))
      O31=1.20D0*t4**0.048
      O32=DMIN1(2.10D0,1.88D0*T4**.275)
      CALL FORB3(9,11,3,1.365D0,1.965D0,11.02D-2,6.915D-1,3.22D0
     &,O21,O31,O32,9.D0,5.D0,1.D0,FB(34),FB(35),FB(36),TE,Z,F)      

      DO NM=34,36            
         FB(NM)=Z*FB(NM)
      enddo
C
C     Ar III 3P,1D,1S OMEGAS and A:s from Pradhan
C
      om21=2.24e0
      om31=0.531e0
      om32=1.18e0
      om1d3p=4.74e0
      om43=1.*om1d3p/9.
      om42=3.*om1d3p/9.
      om41=5.*om1d3p/9.
      om1s3p=0.680e0
      om53=1.*om1s3p/9.
      om52=3.*om1s3p/9.
      om51=5.*om1s3p/9.
      om54=0.823e0
      z=abn(12)*xar(3)
      z3=z
      CALL FIVELEV_dp(96,12,3,5,0.d0,1112.175d0,1570.229d0,14010.004d0,33265.724d0,
     &     5.d0,3.d0,1.d0,5.d0,1.d0,
     &     3.08d-2,2.37d-6,3.14d-1,4.17d-2,5.17d-3,8.23d-2,3.91d0,
     &     2.21d-5,0.d0,2.59d0,
     &     om21,om31,om41,om51,om32,om42,om52,
     &     om43,om53,om54,TE,Z,XI)
      DO I=1,5
        DO J=1,5
         EMCA(I,J)=EMCA(I,J)/DEL(IK)
        ENDDO
      ENDDO
      coarIII=0.
      DO  I=1,5
        DO  J=I+1,5
          COariii=COariii+EMCA(J,I)
        ENDDO
      ENDDO
      line_cool(12,3)=coariii
      if(itime.eq.1) nline_min(12,3)=ll
      ll=nline_min(12,3)
      do k=1,kmaxp(12,3)
         ll=ll+1
         cin(ll)=weh(k)/del(ik)
         tauline(ll)=taul(k)
         if(itime.le.888888) wli(ll)=wlin(k)
         ilabcf(ll)=1203
         cinx(12,3,k)=weh(k)/del(ik)
         taulinex(12,3,k)=taul(k)
         wlix(12,3,k)=wlin(k)
         ilabcfx(12,3,k)=1203
      enddo

      fb(55)=emca(2,1)
      fb(56)=emca(3,2)
      fb(57)=emca(4,1)
      fb(58)=emca(4,2)
      fb(59)=emca(5,1)
      fb(60)=emca(5,2)
      fb(61)=emca(5,4)
C
C     Ar IV 3P,1D,1S OMEGAS from Ramsbottom et al '97 and A:s from Pradhan
C
      om21=1.144*t4**0.027
      om31=0.762*t4**0.028
      om32=7.055*t4**0.086
      om1d3p=2.29*t4**0.132
      om41=0.785*t4**0.140
      om42=3.939*t4**0.025
      om43=2.139*t4**0.027
      om1s3p=0.293*t4**0.167
      om51=0.393*t4**0.140
      om52=1.533*t4**0.028
      om53=1.507*t4**0.018
      om54=2.065*t4**0.383
      z=abn(12)*xar(4)
      z4=z
 9388 format(' Ar IV z ',1pe12.3,10e12.3)
      CALL FIVELEV_dp(97,12,4,5,0.d0,21090.4d0,21219.3d0,34855.5d0,35032.6d0,
     &     4.d0,4.d0,6.d0,2.d0,4.d0,
     &     1.77d-3,2.23d-2,2.11d0,0.862d0,2.30d-5,0.598d0,
     &     0.119d0,0.789d0,0.603d0,4.94d-5,
     &     om21,om31,om41,om51,om32,om42,om52,
     &     om43,om53,om54,TE,Z,XI)
      DO I=1,5
        DO J=1,5
         EMCA(I,J)=EMCA(I,J)/DEL(IK)
        ENDDO
      ENDDO
      coariv=0.
      DO  I=1,5
        DO  J=I+1,5
          coariv=coariv+EMCA(J,I)
        ENDDO
      ENDDO
      line_cool(12,4)=coariv       
      if(itime.eq.1) nline_min(12,4)=ll
      ll=nline_min(12,4)
      do k=1,kmaxp(12,4)
         ll=ll+1
         cin(ll)=weh(k)/del(ik)         
         tauline(ll)=taul(k)
         if(itime.le.888888) wli(ll)=wlin(k)         
         ilabcf(ll)=1204
         cinx(12,4,k)=weh(k)/del(ik)
         taulinex(12,4,k)=taul(k)
         wlix(12,4,k)=wlin(k)
         ilabcfx(12,4,k)=1204
      enddo
      fb(62)=emca(2,1)
      fb(63)=emca(3,1)
      fb(64)=emca(4,1)
      fb(65)=emca(5,1)
      fb(66)=emca(3,2)
      fb(67)=emca(5,4)
      fb(68)=emca(4,2)
      fb(69)=emca(4,3)
      fb(70)=emca(5,2)
      fb(71)=emca(5,3)
C     *************************************************************
C     *****
C     CALCULATE CORONAL AND FS EMISSION FROM IRON
C     *****
C     *************************************************************

c fit to S XII FS from Chianti
      aq=0.15
      t1=1.3e5
      t2=1.e6
      om7612=aq+(te/t1)**1.7/(1.+25.*(te/t2)**1.7)
      

      z=abn(11)*xion(ik,11,12)
c a21 from chianti
      a21=19.60
      if(itime.eq.1) nline_min(11,12)=ll+1
      ll=nline_min(11,12)
      call rloss(11,12,7612.68d0,om7612,2.d0,6.d0,a21,rs,xel,z,
     &     te,cl(ll),w(ll),cin(ll),frqq,ll,wli(ll))
      tauline(ll)=t0
      ilabcf(ll)=1112
      cinx(11,12,1)=cin(ll)
      taulinex(11,12,1)=t0
      wlix(11,12,1)=wli(ll)
      ilabcfx(11,12,1)=1112

c fit to Ar XIV FS 

c om from Chianti
      aq=2.3
      t1=1.e5
      t2=2.e8
      om4414=aq-(te/t1)**1.05/(1.+0.73*(te/t1)**1.05)

c a21 from Träbert 2000
      a21=103.1

      z=abn(14)*xf(12)

c Ne V 1134 Ne V 5S - 3P, OM Lennon and Burke, 1994, A&A Supp 103, 273.
c     A from Mendoza 
      if(itime.eq.1) nline_min(6,5)=ll
      ll=nline_min(6,5)
      ll=ll+1
      om1134=11.9/te**0.23
      z=abn(6)*xne(5)
      call rloss(6,5,1134.d0,om1134,9.d0,5.d0,4.67d3,rs,xel,z,
     &     te,cl(ll),w(ll),cin(ll),frqq,ll,wli(ll))
      tauline(ll)=t0
      ilabcf(ll)=605
      kx=kmaxp(6,5)
      cinx(6,5,kx+1)=cin(ll)
      taulinex(6,5,kx+1)=t0
      wlix(6,5,kx+1)=wli(ll)
      ilabcfx(6,5,kx+1)=605
      line_cool(6,5)=cl(ll)+line_cool(6,5)       
c
c Bowen lines
c

c O III Bowen 374 A
      if(itime.eq.1) then
         nline_min(12,4)=ll
      endif
      ll=nline_min(12,4)
      ll=ll+1
      nbow=ll
      cin(ll)=em374e
      wli(ll)=374.
      ilabcf(ll)=503
      kx=kmaxp(5,3)
      cinx(5,3,kx+1)=em374e
      taulinex(5,3,kx+1)=0.
      wlix(5,3,kx+1)=wli(ll)
      ilabcfx(5,3,kx+1)=503
c O III Bowen 644 A
      ll=ll+1
      cin(ll)=emb644
      wli(ll)=644.
      ilabcf(ll)=503
      cinx(5,3,kx+2)=emb644
      taulinex(5,3,kx+2)=0.
      wlix(5,3,kx+2)=wli(ll)
      ilabcfx(5,3,kx+2)=503
c O III Bowen 703 A
      ll=ll+1
      cin(ll)=emb703
      wli(ll)=703.
      ilabcf(ll)=503
      cinx(5,3,kx+3)=emb703
      taulinex(5,3,kx+3)=0.
      wlix(5,3,kx+3)=wli(ll)
      ilabcfx(5,3,kx+3)=503
c O III Bowen 3210
      ll=ll+1
      cin(ll)=emb3210
      wli(ll)=3210.
      ilabcf(ll)=503
      cinx(5,3,kx+4)=emb3210
      taulinex(5,3,kx+4)=0.
      wlix(5,3,kx+4)=wli(ll)
      ilabcfx(5,3,kx+4)=503
c O III Bowen 3287
      ll=ll+1
      cin(ll)=emb3287
      wli(ll)=3287.
      ilabcf(ll)=503
      cinx(5,3,kx+5)=emb3287
      taulinex(5,3,kx+5)=0.
      wlix(5,3,kx+5)=wli(ll)
      ilabcfx(5,3,kx+5)=503
 9341 format(' O-Bow ',7(0pf7.1,1pe11.3))
c Replace the highest multilevel lines  with the Bowen lines
c     N III Bowen 4642
      ll=ll+1
      nbown=ll
      cin(ll)=emb4642
      wli(ll)=4642.
      ilabcf(ll)=403
      cinx(4,3,401)=emb4642
      taulinex(4,3,401)=0.
      wlix(4,3,401)=4642.
      ilabcfx(4,3,401)=403
c N III Bowen 4100
      ll=ll+1
      cin(ll)=emb4100
      wli(ll)=4100.
      ilabcf(ll)=403
      cinx(4,3,400)=emb4100
      taulinex(4,3,400)=0.
      wlix(4,3,400)=4100.
      ilabcfx(4,3,400)=403
c N III Bowen 691
      ll=ll+1
      wli(ll)=691.
      ilabcf(ll)=403
      cinx(4,3,6)=emb403
      taulinex(4,3,6)=0.
      wlix(4,3,6)=691.
      ilabcfx(4,3,399)=403
c N III Bowen 990
      ll=ll+1
      cin(ll)=emb990
      wli(ll)=990.
      ilabcf(ll)=403
      cinx(4,3,399)=emb990
      taulinex(4,3,399)=0.
      wlix(4,3,399)=990.
      ilabcfx(4,3,399)=403
c N III Bowen 452
      ll=ll+1
      cin(ll)=emb452
      wli(ll)=452.
      ilabcf(ll)=403
      cinx(4,3,398)=emb452
      taulinex(4,3,398)=0.
      wlix(4,3,398)=452.
      ilabcfx(4,3,398)=403
c  total number of lines excluding forb. and Fe lines
      nlinetotb=ll
c  total cooling due to these lines      
      coolcor=0.
      do l1=nline+1,nlinetotb
         coolcor=coolcor+cl(l1)
      enddo

C     *************************************************************
C     *****
C     CALCULATE POPULATIONS AND LINE EMISSION FROM H I
C     *****
C     *************************************************************
      line_cool(1,1)=0.   
C     *************************************************************
C     *****
C     CALCULATE POPULATIONS AND LINE EMISSION FROM HE I
C     *****
C     *************************************************************
C     ************************************************************
C     *****
c   Ar V 5 level atom      
C     *****
C     *************************************************************
      ilcf=1205
c     Omegas from  Galavis, M.~E., Mendoza, C., \& Zeippen, C.~J.\ 1995, \aaps, 111, 347
      call popsimp(12,5,te,XQ,coarv)
      if(itime.eq.1) nline_min(12,5)=ll
      ll=nline_min(12,5)
      do k=1,kmaxp(12,5)
         ll=ll+1
         cin(ll)=weh(k)/del(ik)
         tauline(ll)=taul(k)
         if(itime.le.888888) wli(ll)=wlin(k)
         ilabcf(ll)=1205
         cinx(12,5,k)=weh(k)/del(ik)
         taulinex(12,5,k)=taul(k)
c!!         wlix(12,5,k)=wlin(k)
         ilabcfx(12,5,k)=1205
      enddo
      line_cool(12,5)=coarv       
C     ************************************************************
C     *****
c   Ar VI 2 level atom      
C     *****
C     *************************************************************
      Z=ABn(12)*xar(6)
      OM=0.303*t4**0.065
      telog=log10(te)
c Omegas from Saraph, H. E., & Storey, P. J. 1996, A&AS, 115, 151      
      if(te < ar6_te(1)) then
         al_ar6=log10(ar6_coll(2)/ar6_coll(1))/log10(ar6_te(2)/ar6_te(1))
         om=ar6_coll(1)*(te/ar6_te(1))**al_ar6
         ii=1
      elseif(te > ar6_te(43)) then
         al_ar6=log10(ar6_coll(43)/ar6_coll(42))/log10(ar6_te(43)/ar6_te(42))
         om=ar6_coll(43)*(te/ar6_te(43))**al_ar6
      else
         do i=1,42
            if(te > ar6_te(i) .and. te <= ar6_te(i+1)) then
               al_ar6=log10(ar6_coll(i+1)/ar6_coll(i))/log10(ar6_te(i+1)/ar6_te(i))
               om=ar6_coll(i)*(te/ar6_te(i))**al_ar6
               ii=i
            endif
         enddo
      endif
      CALL RLOSS(12,6,4.52922D4,OM,2.D0,4.D0,9.7D-2,RS,XEL,Z,TE,cl(45)
     &     ,W(45),CIN(45),FRQQ,45,WLI(45))
      tauline(45)=t0
      ilabcf(45)=1206
      Coarvi=cl(45)
      cinx(12,6,1)=cin(45)
      taulinex(12,6,1)=t0
      wlix(12,6,1)=wli(45)
      ilabcfx(12,6,1)=1206
      kmaxp(12,6)=1
      line_cool(12,6)=cl(45)       
C     *************************************************************
C     *****
c   C II 53 level atom      
C     *****
C     *************************************************************
      ilcf=302
      call cpu_time(cpu_tsec)
      call popchianti_new(3,2,te,coc2)
      if(itime.eq.1) nline_min(3,2)=ll
      ll=nline_min(3,2)
      do k=1,kmaxp(3,2)
         ll=ll+1
         if(ipopch.eq.1) then
            wli(ll)=wlin(k)
            cin(ll)=weh(k)/del(ik)
            tauline(ll)=taul(k)
            cinx(3,2,k)=weh(k)/del(ik)
            taulinex(3,2,k)=taul(k)
c!!            wlix(5,7,k)=wlin(k)
         endif
         ilabcf(ll)=ilcf
         ilabcfx(3,2,k)=302
      enddo
      line_cool(3,2)=coc2
      
C     *************************************************************
C     *****
c   C V  49 level atom      
C     *****
C     *************************************************************
      ilcf=305
      call cpu_time(cpu_tsec)
      call popchianti(3,5,te,XQ,coc5)
      if(itime.eq.1) nline_min(3,5)=ll
      ll=nline_min(3,5)
      do k=1,kmaxp(3,5)
         ll=ll+1
         if(ipopch.eq.1) then
            wli(ll)=wlin(k)
            cin(ll)=weh(k)/del(ik)
            tauline(ll)=taul(k)
            cinx(3,5,k)=weh(k)/del(ik)
            taulinex(3,5,k)=taul(k)
c!!            wlix(3,5,k)=wlin(k)
            ilabcfx(3,5,k)=305
         endif
         ilabcf(ll)=ilcf
      enddo
      line_cool(3,5)=coc5
C     *************************************************************
C     *****
c   C VI 25  level atom      
C     *****
C     *************************************************************
      ilcf=306
      call cpu_time(cpu_tsec)
      call popchianti(3,6,te,XQ,coc6)
      if(itime.eq.1) nline_min(3,6)=ll
      ll=nline_min(3,6)
      do k=1,kmaxp(3,6)
         ll=ll+1
         if(ipopch.eq.1) then
            wli(ll)=wlin(k)
            cin(ll)=weh(k)/del(ik)
            tauline(ll)=taul(k)
            cinx(3,6,k)=weh(k)/del(ik)
            taulinex(3,6,k)=taul(k)
            wlix(3,6,k)=wlin(k)
            ilabcfx(3,6,k)=306
         endif
         ilabcf(ll)=ilcf
      enddo
      line_cool(3,6)=coc6
C     *************************************************************
C     *****
c   N VI  49 level atom      
C     *****
C     *************************************************************
      ilcf=406
      call cpu_time(cpu_tsec)
      call popchianti(4,6,te,XQ,con6)
      if(itime.eq.1) nline_min(4,6)=ll
      ll=nline_min(4,6)
      do k=1,kmaxp(4,6)
         ll=ll+1
         if(ipopch.eq.1) then
            wli(ll)=wlin(k)
            cin(ll)=weh(k)/del(ik)
            tauline(ll)=taul(k)
         endif
         ilabcf(ll)=ilcf
         cinx(4,6,k)=weh(k)/del(ik)
         taulinex(4,6,k)=taul(k)
         wlix(4,6,k)=wlin(k)
         ilabcfx(4,6,k)=406
      enddo
      line_cool(4,6)=con6
C     *************************************************************
C     *****
c   N VII  25 level atom      
C     *****
C     *************************************************************
      ilcf=407
      call popchianti(4,7,te,XQ,con7)
      if(itime.eq.1) nline_min(4,7)=ll
      ll=nline_min(4,7)
      do k=1,kmaxp(4,7)
         ll=ll+1
         if(ipopch.eq.1) then
            wli(ll)=wlin(k)
            cin(ll)=weh(k)/del(ik)
            tauline(ll)=taul(k)
         endif
         ilabcf(ll)=ilcf
         cinx(4,7,k)=weh(k)/del(ik)
         taulinex(4,7,k)=taul(k)
         wlix(4,7,k)=wlin(k)
         ilabcfx(4,7,k)=407
      enddo
      line_cool(4,7)=con7
C     *************************************************************
C     *****
c   O II 35 level atom      
C     *****
C     *************************************************************
      ilcf=502
      call popchianti_new(5,2,te,coo2)
      if(itime.eq.1) nline_min(5,2)=ll
      ll=nline_min(5,2)
      do k=1,kmaxp(5,2)
         ll=ll+1
         if(ipopch.eq.1) then
            wli(ll)=wlin(k)
            cin(ll)=weh(k)/del(ik)
            tauline(ll)=taul(k)
            cinx(5,2,k)=weh(k)/del(ik)
            taulinex(5,2,k)=taul(k)
            wlix(5,7,k)=wlin(k)
         endif
         ilabcf(ll)=ilcf
         ilabcfx(5,2,k)=502
      enddo
      line_cool(5,2)=coo2
C     *************************************************************
C     *****
c   O VII 47 level atom      
C     *****
C     *************************************************************
      ilcf=507
      call popchianti(5,7,te,XQ,coo7)
      if(itime.eq.1) nline_min(5,7)=ll
      ll=nline_min(5,7)
      do k=1,kmaxp(5,7)
         ll=ll+1
         if(ipopch.eq.1) then
            wli(ll)=wlin(k)
            cin(ll)=weh(k)/del(ik)
            tauline(ll)=taul(k)
            cinx(5,7,k)=weh(k)/del(ik)
            taulinex(5,7,k)=taul(k)
            wlix(5,7,k)=wlin(k)
         endif
         ilabcf(ll)=ilcf
         ilabcfx(5,7,k)=507
      enddo
      line_cool(5,7)=coo7
C     *************************************************************

C     *************************************************************
C     *****
c   O VIII 25 level atom      
C     *****
C     *************************************************************
      ilcf=508
      call popchianti(5,8,te,XQ,coo8)
      if(itime.eq.1) nline_min(5,8)=ll
      ll=nline_min(5,8)
      do k=1,kmaxp(5,8)
         ll=ll+1
         if(ipopch.eq.1) then
            wli(ll)=wlin(k)
            cin(ll)=weh(k)/del(ik)
            tauline(ll)=taul(k)
            cinx(5,8,k)=weh(k)/del(ik)
            taulinex(5,8,k)=taul(k)
            wlix(5,8,k)=wlin(k)
         endif
         ilabcfx(5,8,k)=508
         ilabcf(ll)=ilcf
      enddo
      line_cool(5,8)=coo8
C     *************************************************************
C     *****
c   Ne VI 180 level atom      
C     *****
C     *************************************************************
      ilcf=606
      call popchianti(6,6,te,XQ,cone6)
      if(itime.eq.1) nline_min(6,6)=ll
      ll=nline_min(6,6)
      do k=1,kmaxp(6,6)
         ll=ll+1
         if(ipopch.eq.1) then
            wli(ll)=wlin(k)
            cin(ll)=weh(k)/del(ik)
            tauline(ll)=taul(k)
            cinx(6,6,k)=weh(k)/del(ik)
            taulinex(6,6,k)=taul(k)
            wlix(6,6,k)=wlin(k)
         endif
         ilabcfx(6,6,k)=606
         ilabcf(ll)=ilcf
      enddo
      line_cool(6,6)=cone6
C     *************************************************************
C     *****
c   Ne VII 46 level atom      
C     *****
C     *************************************************************
      ilcf=607
      call popchianti(6,7,te,XQ,cone7)
      if(itime.eq.1) nline_min(6,7)=ll
      ll=nline_min(6,7)
      do k=1,kmaxp(6,7)
         ll=ll+1
         if(ipopch.eq.1) then
            wli(ll)=wlin(k)
            cin(ll)=weh(k)/del(ik)
            tauline(ll)=taul(k)
            cinx(6,7,k)=weh(k)/del(ik)
            taulinex(6,7,k)=taul(k)
            wlix(6,7,k)=wlin(k)
         endif
         ilabcfx(6,7,k)=607
         ilabcf(ll)=ilcf
      enddo
      line_cool(6,7)=cone7
C     *************************************************************
C     *****
c   Ne VIII 40 level atom      
C     *****
C     *************************************************************
      ilcf=608
      call popchianti(6,8,te,XQ,cone8)
      if(itime.eq.1) nline_min(6,8)=ll
      ll=nline_min(6,8)
      do k=1,kmaxp(6,8)
         ll=ll+1
         if(ipopch.eq.1) then
            wli(ll)=wlin(k)
            cin(ll)=weh(k)/del(ik)
            tauline(ll)=taul(k)
            cinx(6,8,k)=weh(k)/del(ik)
            taulinex(6,8,k)=taul(k)
            wlix(6,8,k)=wlin(k)
         endif
         ilabcfx(6,8,k)=608
         ilabcf(ll)=ilcf
      enddo
      line_cool(6,8)=cone8
C     *************************************************************
C     *****
c   Ne IX 49 level atom      
C     *****
C     *************************************************************
      ilcf=609
      call popchianti(6,9,te,XQ,cone9)
      if(itime.eq.1) nline_min(6,9)=ll
      ll=nline_min(6,9)
      do k=1,kmaxp(6,9)
         ll=ll+1
         if(ipopch.eq.1) then
            wli(ll)=wlin(k)
            cin(ll)=weh(k)/del(ik)
            tauline(ll)=taul(k)
            cinx(6,9,k)=weh(k)/del(ik)
            taulinex(6,9,k)=taul(k)
            wlix(6,9,k)=wlin(k)
         endif
         ilabcfx(6,9,k)=609
         ilabcf(ll)=ilcf
      enddo
      line_cool(6,9)=cone9
C     *************************************************************
C     *****
c   Ne X 36 level atom      
C     *****
C     *************************************************************
      ilcf=610
      call popchianti(6,10,te,XQ,cone10)
      if(itime.eq.1) nline_min(6,10)=ll
      ll=nline_min(6,10)
      do k=1,kmaxp(6,10)
         ll=ll+1
         if(ipopch.eq.1) then
            wli(ll)=wlin(k)
            cin(ll)=weh(k)/del(ik)
            tauline(ll)=taul(k)
            cinx(6,10,k)=weh(k)/del(ik)
            taulinex(6,10,k)=taul(k)
            wlix(6,10,k)=wlin(k)
         endif
         ilabcfx(6,10,k)=610         
         ilabcf(ll)=ilcf
      enddo
      line_cool(6,10)=cone10
C     *************************************************************
C     *****
c   Mg II 41 level atom      
C     *****
C     *************************************************************
      ilcf=802
      call popchianti_new(8,2,te,comg2)
      if(itime.eq.1) nline_min(8,2)=ll
      ll=nline_min(8,2)
      do k=1,kmaxp(8,2)
         ll=ll+1
         if(ipopch.eq.1) then
            wli(ll)=wlin(k)
            cin(ll)=weh(k)/del(ik)
            tauline(ll)=taul(k)
            cinx(8,2,k)=weh(k)/del(ik)
            taulinex(8,2,k)=taul(k)
            wlix(8,2,k)=wlin(k)
         endif
         ilabcf(ll)=ilcf
         ilabcfx(8,2,k)=802
      enddo
      line_cool(8,2)=comg2
      
C     *************************************************************
C     *****
c   Mg VIII 125 level atom      
C     *****
C     *************************************************************
      ilcf=808
      call popchianti(8,8,te,XQ,comg8)
      if(itime.eq.1) nline_min(8,8)=ll
      ll=nline_min(8,8)
      do k=1,kmaxp(8,8)
         ll=ll+1
         if(ipopch.eq.1) then
            wli(ll)=wlin(k)
            cin(ll)=weh(k)/del(ik)
            tauline(ll)=taul(k)
            cinx(8,8,k)=weh(k)/del(ik)
            taulinex(8,8,k)=taul(k)
            wlix(8,8,k)=wlin(k)
         endif
         ilabcfx(8,8,k)=808
         ilabcf(ll)=ilcf
      enddo
      line_cool(8,8)=comg8
C     *************************************************************
C     *****
c   Mg IX 46 level atom      
C     *****
C     *************************************************************
      ilcf=809
      call popchianti(8,9,te,XQ,comg9)
      if(itime.eq.1) nline_min(8,9)=ll
      ll=nline_min(8,9)
      do k=1,kmaxp(8,9)
         ll=ll+1
         if(ipopch.eq.1) then
            wli(ll)=wlin(k)
            cin(ll)=weh(k)/del(ik)
            tauline(ll)=taul(k)
            cinx(8,9,k)=weh(k)/del(ik)
            taulinex(8,9,k)=taul(k)
            wlix(8,9,k)=wlin(k)
         endif
         ilabcfx(8,9,k)=809
         ilabcf(ll)=ilcf
      enddo
      line_cool(8,9)=comg9
C     *************************************************************
C     *****
c   Mg X 40 level atom      
C     *****
C     *************************************************************
      ilcf=810
      call popchianti(8,10,te,XQ,comg10)
      if(itime.eq.1) nline_min(8,10)=ll
      ll=nline_min(8,10)
      do k=1,kmaxp(8,10)
         ll=ll+1
         if(ipopch.eq.1) then
            wli(ll)=wlin(k)
            cin(ll)=weh(k)/del(ik)
            tauline(ll)=taul(k)
            cinx(8,10,k)=weh(k)/del(ik)
            taulinex(8,10,k)=taul(k)
            wlix(8,10,k)=wlin(k)
         endif
         ilabcfx(8,10,k)=810
         ilabcf(ll)=ilcf
      enddo
      line_cool(8,10)=comg10
C     *************************************************************
C     *****
c   Mg XI 49 level atom      
C     *****
C     *************************************************************
      ilcf=811
      call popchianti(8,11,te,XQ,comg11)
      if(itime.eq.1) nline_min(8,11)=ll
      ll=nline_min(8,11)
      do k=1,kmaxp(8,11)
         ll=ll+1
         if(ipopch.eq.1) then
            wli(ll)=wlin(k)
            cin(ll)=weh(k)/del(ik)
            tauline(ll)=taul(k)
            cinx(8,11,k)=weh(k)/del(ik)
            taulinex(8,11,k)=taul(k)
            wlix(8,11,k)=wlin(k)
         endif
         ilabcfx(8,11,k)=811
         ilabcf(ll)=ilcf
      enddo
      line_cool(8,11)=comg11
C     *************************************************************
C     *****
c   Mg XII 25 level atom      
C     *****
C     *************************************************************
      ilcf=812
      call popchianti(8,12,te,XQ,comg12)
      if(itime.eq.1) nline_min(8,12)=ll
      ll=nline_min(8,12)
      do k=1,kmaxp(8,12)
         ll=ll+1
         if(ipopch.eq.1) then
            wli(ll)=wlin(k)
            cin(ll)=weh(k)/del(ik)
            tauline(ll)=taul(k)
            cinx(8,12,k)=weh(k)/del(ik)
            taulinex(8,12,k)=taul(k)
            wlix(8,12,k)=wlin(k)
         endif
         ilabcfx(8,12,k)=812            
         ilabcf(ll)=ilcf
      enddo
      line_cool(8,12)=comg12
C     *************************************************************
C     *****
c   Si VIII 72 level atom      
C     *****
C     *************************************************************
      ilcf=1008
      call popchianti(10,8,te,XQ,cosi8)
      if(itime.eq.1) nline_min(10,8)=ll
      ll=nline_min(10,8)
      do k=1,kmaxp(10,8)
         ll=ll+1
         if(ipopch.eq.1) then
            wli(ll)=wlin(k)
            cin(ll)=weh(k)/del(ik)
            tauline(ll)=taul(k)
            cinx(10,8,k)=weh(k)/del(ik)
            taulinex(10,8,k)=taul(k)
            wlix(10,8,k)=wlin(k)
         endif
         ilabcfx(10,8,k)=1008
         ilabcf(ll)=ilcf
      enddo
      line_cool(10,8)=cosi8
C     *************************************************************
C     *****
c   Si IX 46 level atom      
C     *****
C     *************************************************************
      ilcf=1009
      call popchianti(10,9,te,XQ,cosi9)
      if(itime.eq.1) nline_min(10,9)=ll
      ll=nline_min(10,9)
      do k=1,kmaxp(10,9)
         ll=ll+1
         if(ipopch.eq.1) then
            wli(ll)=wlin(k)
            cin(ll)=weh(k)/del(ik)
            tauline(ll)=taul(k)
            cinx(10,9,k)=weh(k)/del(ik)
            taulinex(10,9,k)=taul(k)
            wlix(10,9,k)=wlin(k)
         endif
         ilabcfx(10,9,k)=1009
         ilabcf(ll)=ilcf
      enddo
      line_cool(10,9)=cosi9
C     *************************************************************
C     *****
c   Si X 125 level atom      
C     *****
C     *************************************************************
      ilcf=1010
      call popchianti(10,10,te,XQ,cosi10)
      if(itime.eq.1) nline_min(10,10)=ll
      ll=nline_min(10,10)
      do k=1,kmaxp(10,10)
         ll=ll+1
         if(ipopch.eq.1) then
            wli(ll)=wlin(k)
            cin(ll)=weh(k)/del(ik)
            tauline(ll)=taul(k)
            cinx(10,10,k)=weh(k)/del(ik)
            taulinex(10,10,k)=taul(k)
            wlix(10,10,k)=wlin(k)
         endif
         ilabcfx(10,10,k)=1010
         ilabcf(ll)=ilcf
      enddo
      line_cool(10,10)=cosi10
C     *************************************************************
C     *****
c   Si XI 46 level atom      
C     *****
C     *************************************************************
      ilcf=1011
      call popchianti(10,11,te,XQ,cosi11)
      if(itime.eq.1) nline_min(10,11)=ll
      ll=nline_min(10,11)
      do k=1,kmaxp(10,11)
         ll=ll+1
         if(ipopch.eq.1) then
            wli(ll)=wlin(k)
            cin(ll)=weh(k)/del(ik)
            tauline(ll)=taul(k)
            cinx(10,11,k)=weh(k)/del(ik)
            taulinex(10,11,k)=taul(k)
            wlix(10,11,k)=wlin(k)
         endif
         ilabcfx(10,11,k)=1011
         ilabcf(ll)=ilcf
      enddo
      line_cool(10,11)=cosi11
C     *************************************************************
C     *****
c   Si XII 40 level atom      
C     *****
C     *************************************************************
      ilcf=1012
      call popchianti(10,12,te,XQ,cosi12)
      if(itime.eq.1) nline_min(10,12)=ll
      ll=nline_min(10,12)
      do k=1,kmaxp(10,12)
         ll=ll+1
         if(ipopch.eq.1) then
            wli(ll)=wlin(k)
            cin(ll)=weh(k)/del(ik)
            tauline(ll)=taul(k)
            cinx(10,12,k)=weh(k)/del(ik)
            taulinex(10,12,k)=taul(k)
            wlix(10,12,k)=wlin(k)
         endif
         ilabcfx(10,12,k)=1012
         ilabcf(ll)=ilcf
      enddo
      line_cool(10,12)=cosi12
C     *************************************************************
C     *****
c   Si XIII 49 level atom      
C     *****
C     *************************************************************
      ilcf=1013
      call popchianti(10,13,te,XQ,cosi13)
      if(itime.eq.1) nline_min(10,13)=ll
      ll=nline_min(10,13)
      do k=1,kmaxp(10,13)
         ll=ll+1
         if(ipopch.eq.1) then
            wli(ll)=wlin(k)
            cin(ll)=weh(k)/del(ik)
            tauline(ll)=taul(k)
            cinx(10,13,k)=weh(k)/del(ik)
            taulinex(10,13,k)=taul(k)
            wlix(10,13,k)=wlin(k)
         endif
         ilabcfx(10,13,k)=1013
         ilabcf(ll)=ilcf
      enddo
      line_cool(10,13)=cosi13
C     *************************************************************
C     *****
c   Si XIV  25 level atom      
C     *****
C     *************************************************************
      ilcf=1014
      call popchianti(10,14,te,XQ,cosi14)
      do k=1,kmaxp(10,14)
         ll=ll+1
         if(ipopch.eq.1) then
            wli(ll)=wlin(k)
            cin(ll)=weh(k)/del(ik)
            tauline(ll)=taul(k)
            cinx(10,14,k)=weh(k)/del(ik)
            taulinex(10,14,k)=taul(k)
            wlix(10,14,k)=wlin(k)
         endif
         ilabcfx(10,14,k)=1014
         ilabcf(ll)=ilcf
      enddo
      line_cool(10,14)=cosi14
C     *************************************************************
C     *****
c   S II 43 level atom      
C     *****
C     *************************************************************
      ilcf=1102
      call popchianti(11,2,te,XQ,cos2)
      do k=1,kmaxp(11,2)
         ll=ll+1
         if(ipopch.eq.1) then
            wli(ll)=wlin(k)
            cin(ll)=weh(k)/del(ik)
            tauline(ll)=taul(k)
            cinx(11,2,k)=weh(k)/del(ik)
            taulinex(11,2,k)=taul(k)
            wlix(11,2,k)=wlin(k)
         endif
         ilabcfx(11,2,k)=1102
         ilabcf(ll)=ilcf
      enddo
      line_cool(11,2)=cos2
      siiterm=1.8*1.6e-12*1.2e-3*exp(-2.086e4/te)*abn(11)*xion(ik,11,2)
     &     /denel
  
C     *************************************************************
C     *****
c   S III 49 level atom      
C     *****
C     *************************************************************
      ilcf=1103
      call popchianti(11,3,te,XQ,cos3)
      if(itime.eq.1) nline_min(11,3)=ll
      ll=nline_min(11,3)
      do k=1,kmaxp(11,3)
         ll=ll+1
         if(ipopch.eq.1) then
            wli(ll)=wlin(k)
            cin(ll)=weh(k)/del(ik)
            tauline(ll)=taul(k)
            cinx(11,3,k)=weh(k)/del(ik)
            taulinex(11,3,k)=taul(k)
            wlix(11,3,k)=wlin(k)
         endif
         ilabcfx(11,3,k)=1103
         ilabcf(ll)=ilcf
      enddo
      line_cool(11,3)=cos3
      
c S IV 
      ilcf=1104
      z=abn(11)*xsu(4)
      call popchianti_new(11,4,te,cos4)

      do k=1,kmaxp(11,4)
         if(ipopch.eq.1) then
            cinx(11,4,k)=weh(k)/del(ik)
            taulinex(11,4,k)=taul(k)
            wlix(11,4,k)=wlin(k)
         endif
         ilabcfx(11,4,k)=1104
      enddo
      line_cool(11,4)=cos4

c S V 
      ilcf=1105
      z=abn(11)*xsu(5)
      call popchianti_new(11,5,te,cos5)
      do k=1,kmaxp(11,5)
         if(ipopch.eq.1) then
            cinx(11,5,k)=weh(k)/del(ik)
            taulinex(11,5,k)=taul(k)
            wlix(11,5,k)=wlin(k)
         endif
         ilabcfx(11,5,k)=1105
      enddo
      line_cool(11,5)=cos5
C     *************************************************************
C     *****
c   S XIII 46 level atom      
C     *****
C     *************************************************************
      ilcf=1113
      call popchianti(11,13,te,XQ,cos13)
      if(itime.eq.1) nline_min(11,13)=ll
      ll=nline_min(11,13)
      do k=1,kmaxp(11,13)
         ll=ll+1
         if(ipopch.eq.1) then
            wli(ll)=wlin(k)
            cin(ll)=weh(k)/del(ik)
            tauline(ll)=taul(k)
            cinx(11,13,k)=weh(k)/del(ik)
            taulinex(11,13,k)=taul(k)
            wlix(11,13,k)=wlin(k)
         endif
         ilabcfx(11,13,k)=1113
         ilabcf(ll)=ilcf
      enddo
      line_cool(11,13)=cos13
C     *************************************************************
C     *****
c   S XIV 40 level atom      
C     *****
C     *************************************************************
      ilcf=1114
      call popchianti(11,14,te,XQ,cos14)
      if(itime.eq.1) nline_min(11,14)=ll
      ll=nline_min(11,14)
      do k=1,kmaxp(11,14)
         ll=ll+1
         if(ipopch.eq.1) then
            wli(ll)=wlin(k)
            cin(ll)=weh(k)/del(ik)
            tauline(ll)=taul(k)
            cinx(11,14,k)=weh(k)/del(ik)
            taulinex(11,14,k)=taul(k)
            wlix(11,14,k)=wlin(k)
         endif
         ilabcfx(11,14,k)=1114
         ilabcf(ll)=ilcf
      enddo
      line_cool(11,14)=cos14
C     *************************************************************
C     *****
c   S XV 49 level atom      
C     *****
C     *************************************************************
      ilcf=1115
      call popchianti(11,15,te,XQ,cos15)
      if(itime.eq.1) nline_min(11,15)=ll
      ll=nline_min(11,15)
      do k=1,kmaxp(11,15)
         ll=ll+1
         if(ipopch.eq.1) then
            wli(ll)=wlin(k)
            cin(ll)=weh(k)/del(ik)
            tauline(ll)=taul(k)
            cinx(11,15,k)=weh(k)/del(ik)
            taulinex(11,15,k)=taul(k)
            wlix(11,15,k)=wlin(k)
         endif
         ilabcfx(11,15,k)=1115
         ilabcf(ll)=ilcf
      enddo
      line_cool(11,15)=cos15
C     *************************************************************
C     *****
c   S XVI 25 level atom      
C     *****
C     *************************************************************
      ilcf=1116
      call popchianti(11,16,te,XQ,cos16)
      if(itime.eq.1) nline_min(11,16)=ll
      ll=nline_min(11,16)
      do k=1,kmaxp(11,16)
         ll=ll+1
         if(ipopch.eq.1) then
            wli(ll)=wlin(k)
            cin(ll)=weh(k)/del(ik)
            tauline(ll)=taul(k)
            cinx(11,16,k)=weh(k)/del(ik)
            taulinex(11,16,k)=taul(k)
            wlix(11,16,k)=wlin(k)
         endif
         ilabcfx(11,16,k)=1116
         ilabcf(ll)=ilcf
      enddo
      line_cool(11,16)=cos16
C     *************************************************************
C     *****
c   Ca IX 16 level atom      
C     *****
C     *************************************************************
      ilcf=1309
      call popchianti(13,9,te,XQ,coca9)
      if(itime.eq.1) nline_min(13,9)=ll
      ll=nline_min(13,9)
      do k=1,kmaxp(13,9)
         ll=ll+1
         if(ipopch.eq.1) then
            wli(ll)=wlin(k)
            cin(ll)=weh(k)/del(ik)
            tauline(ll)=taul(k)
            cinx(13,9,k)=weh(k)/del(ik)
            taulinex(13,9,k)=taul(k)
            wlix(13,9,k)=wlin(k)
         endif
         ilabcfx(13,9,k)=1309
         ilabcf(ll)=ilcf
      enddo
      line_cool(13,9)=coca9
C     *************************************************************
C     *****
c   Ca X 21 level atom      
C     *****
C     *************************************************************
      ilcf=1310
      call popchianti(13,10,te,XQ,coca10)
      if(itime.eq.1) nline_min(13,10)=ll
      ll=nline_min(13,10)
      do k=1,kmaxp(13,10)
         ll=ll+1
         if(ipopch.eq.1) then
            wli(ll)=wlin(k)
            cin(ll)=weh(k)/del(ik)
            tauline(ll)=taul(k)
            cinx(13,10,k)=weh(k)/del(ik)
            taulinex(13,10,k)=taul(k)
            wlix(13,10,k)=wlin(k)
         endif
         ilabcfx(13,10,k)=1310
         ilabcf(ll)=ilcf
      enddo
      line_cool(13,10)=coca10
C     *************************************************************
C     *****
c   Ca XI 89 level atom      
C     *****
C     *************************************************************
      ilcf=1311
      call popchianti(13,11,te,XQ,coca11)
      if(itime.eq.1) nline_min(13,11)=ll
      ll=nline_min(13,11)
      do k=1,kmaxp(13,11)
         ll=ll+1
         if(ipopch.eq.1) then
            wli(ll)=wlin(k)
            cin(ll)=weh(k)/del(ik)
            tauline(ll)=taul(k)
            cinx(13,11,k)=weh(k)/del(ik)
            taulinex(13,11,k)=taul(k)
            wlix(13,11,k)=wlin(k)
         endif
         ilabcfx(13,11,k)=1311
         ilabcf(ll)=ilcf
      enddo
      line_cool(13,11)=coca11
C     *************************************************************
C     *****
c   Ca XII 3 level atom      
C     *****
C     *************************************************************
      ilcf=1312
      call popchianti(13,12,te,XQ,coca12)
      if(itime.eq.1) nline_min(13,12)=ll
      ll=nline_min(13,12)
      do k=1,kmaxp(13,12)
         ll=ll+1
         if(ipopch.eq.1) then
            wli(ll)=wlin(k)
            cin(ll)=weh(k)/del(ik)
            tauline(ll)=taul(k)
            cinx(13,12,k)=weh(k)/del(ik)
            taulinex(13,12,k)=taul(k)
            wlix(13,12,k)=wlin(k)
         endif
         ilabcfx(13,12,k)=1312
         ilabcf(ll)=ilcf
      enddo
      line_cool(13,12)=coca12
C     *************************************************************
C     *****
c   Ca XIII 86 level atom      
C     *****
C     *************************************************************
      ilcf=1313
      call popchianti(13,13,te,XQ,coca13)
      if(itime.eq.1) nline_min(13,13)=ll
      ll=nline_min(13,13)
      do k=1,kmaxp(13,13)
         ll=ll+1
         if(ipopch.eq.1) then
            wli(ll)=wlin(k)
            cin(ll)=weh(k)/del(ik)
            tauline(ll)=taul(k)
            cinx(13,13,k)=weh(k)/del(ik)
            taulinex(13,13,k)=taul(k)
            wlix(13,13,k)=wlin(k)
         endif
         ilabcfx(13,13,k)=1313
         ilabcf(ll)=ilcf
      enddo
      line_cool(13,13)=coca13
C     *************************************************************
C     *****
c   Ca XIV 91 level atom      
C     *****
C     *************************************************************
      ilcf=1314
      call popchianti(13,14,te,XQ,coca14)
      if(itime.eq.1) nline_min(13,14)=ll
      ll=nline_min(13,14)
      do k=1,kmaxp(13,14)
         ll=ll+1
         if(ipopch.eq.1) then
            wli(ll)=wlin(k)
            cin(ll)=weh(k)/del(ik)
            tauline(ll)=taul(k)
            cinx(13,14,k)=weh(k)/del(ik)
            taulinex(13,14,k)=taul(k)
            wlix(13,14,k)=wlin(k)
         endif
         ilabcfx(13,14,k)=1314
         ilabcf(ll)=ilcf
      enddo
      line_cool(13,14)=coca14

C     *************************************************************
C     *****
c   Fe III 55 level atom (should be larger!)
C     *****
C     *************************************************************
      ilcf=1403
      call popchianti_new(14,3,te,cofe3)
      if(itime.eq.1) nline_min(14,3)=ll
      ll=nline_min(14,3)
      do k=1,kmaxp(14,3)
         ll=ll+1
         if(ipopch.eq.1) then
            wli(ll)=wlin(k)
            cin(ll)=weh(k)/del(ik)
            tauline(ll)=taul(k)
            cinx(14,3,k)=weh(k)/del(ik)
            taulinex(14,3,k)=taul(k)
            wlix(14,3,k)=wlin(k)
         endif
         ilabcf(ll)=ilcf
         ilabcfx(14,3,k)=1403
      enddo


      do io=3,9
         if(io==3) then
            ilcf=1403
         elseif(io==4) then
            ilcf=1404
         elseif(io==5) then
            ilcf=1405
         elseif(io==6) then
            ilcf=1406
         elseif(io==7) then
            ilcf=1407
         elseif(io==8) then
            ilcf=1408
         elseif(io==9) then
            ilcf=1409
         endif
         call popchianti_new(14,io,te,cofeio)
         if(itime.eq.1) nline_min(14,io)=ll
         ll=nline_min(14,io)
         do k=1,kmaxp(14,io)
            ll=ll+1
            if(ipopch.eq.1) then
               wli(ll)=wlin(k)
               cin(ll)=weh(k)/del(ik)
               tauline(ll)=taul(k)
               cinx(14,io,k)=weh(k)/del(ik)
               taulinex(14,io,k)=taul(k)
               wlix(14,io,k)=wlin(k)
            endif
            ilabcf(ll)=ilcf
            ilabcfx(14,io,k)=ilcf
         enddo
         if(io==3) then
            cofe3=cofeio
         elseif(io==4) then
            cofe4=cofeio
         elseif(io==5) then
            cofe5=cofeio
         elseif(io==6) then
            cofe6=cofeio
         elseif(io==7) then
            cofe7=cofeio
         elseif(io==8) then
            cofe8=cofeio
         elseif(io==9) then
            cofe9=cofeio
         endif
         line_cool(14,io)=cofeio
      enddo
         

C     *************************************************************
C     *****
c   Fe VII 9 level atom      
C     *****
C     *************************************************************
      call popsimp(14,7,te,XQ,cofe7)
      if(itime.eq.1) nline_min(14,7)=ll
      ll=nline_min(14,7)
      do k=1,kmaxp(14,7)
         ll=ll+1
         cin(ll)=weh(k)/del(ik)
         tauline(ll)=taul(k)
         if(itime.le.888888) wli(ll)=wlin(k)
         ilabcf(ll)=1407
         cinx(14,7,k)=weh(k)/del(ik)
         taulinex(14,7,k)=taul(k)
         wlix(14,7,k)=wlin(k)
         ilabcfx(14,7,k)=1407
      enddo
      line_cool(14,7)=cofe7
C     *************************************************************
C     *****
c   Fe X 172 level atom      
C     *****
C     *************************************************************
      ilcf=1410
      call popchianti(14,10,te,XQ,cofe10)
      if(itime.eq.1) nline_min(14,10)=ll
      ll=nline_min(14,10)
      do k=1,kmaxp(14,10)
         ll=ll+1
         if(ipopch.eq.1) then
            wli(ll)=wlin(k)
            cin(ll)=weh(k)/del(ik)
            tauline(ll)=taul(k)
            cinx(14,10,k)=weh(k)/del(ik)
            taulinex(14,10,k)=taul(k)
            wlix(14,10,k)=wlin(k)
         endif
         ilabcfx(14,10,k)=1410
         ilabcf(ll)=ilcf
      enddo
      line_cool(14,10)=cofe10
C     *************************************************************
C     *****
c   Fe XI 47 level atom      
C     *****
C     *************************************************************
      ilcf=1411
      call popchianti(14,11,te,XQ,cofe11)
      if(itime.eq.1) nline_min(14,11)=ll
      ll=nline_min(14,11)
      do k=1,kmaxp(14,11)
         ll=ll+1
         if(ipopch.eq.1) then
            wli(ll)=wlin(k)
            cin(ll)=weh(k)/del(ik)
            tauline(ll)=taul(k)
            cinx(14,11,k)=weh(k)/del(ik)
            taulinex(14,11,k)=taul(k)
            wlix(14,11,k)=wlin(k)
         endif
         ilabcfx(14,11,k)=1411
         ilabcf(ll)=ilcf
      enddo
      line_cool(14,11)=cofe11
C     *************************************************************
C     *****
c   Fe XII 143 level atom      
C     *****
C     *************************************************************
      ilcf=1412
      call popchianti(14,12,te,XQ,cofe12)
      if(itime.eq.1) nline_min(14,12)=ll
      ll=nline_min(14,12)
      do k=1,kmaxp(14,12)
         ll=ll+1
         if(ipopch.eq.1) then
            wli(ll)=wlin(k)
            cin(ll)=weh(k)/del(ik)
            tauline(ll)=taul(k)
            cinx(14,12,k)=weh(k)/del(ik)
            taulinex(14,12,k)=taul(k)
            wlix(14,12,k)=wlin(k)
         endif
         ilabcfx(14,12,k)=1412
         ilabcf(ll)=ilcf
      enddo
      line_cool(14,12)=cofe12
C     *************************************************************
C     *****
c   Fe XIII 27 level atom      
C     *****
C     *************************************************************
      ilcf=1413
      call popchianti(14,13,te,XQ,cofe13)
      if(itime.eq.1) nline_min(14,13)=ll
      ll=nline_min(14,13)
      do k=1,kmaxp(14,13)
         ll=ll+1
         if(ipopch.eq.1) then
            wli(ll)=wlin(k)
            cin(ll)=weh(k)/del(ik)
            tauline(ll)=taul(k)
            cinx(14,13,k)=weh(k)/del(ik)
            taulinex(14,13,k)=taul(k)
            wlix(14,13,k)=wlin(k)
         endif
         ilabcfx(14,13,k)=1413         
         ilabcf(ll)=ilcf
      enddo
      line_cool(14,13)=cofe13
C     *************************************************************
C     *****
c   Fe XIV 40 level atom      
C     *****
C     *************************************************************
      ilcf=1414
      call popchianti(14,14,te,XQ,cofe14)
      if(itime.eq.1) nline_min(14,14)=ll
      ll=nline_min(14,14)
      do k=1,kmaxp(14,14)
         ll=ll+1
         if(ipopch.eq.1) then
            wli(ll)=wlin(k)
            cin(ll)=weh(k)/del(ik)
            tauline(ll)=taul(k)
            cinx(14,14,k)=weh(k)/del(ik)
            taulinex(14,14,k)=taul(k)
            wlix(14,14,k)=wlin(k)
         endif
         ilabcfx(14,14,k)=1414
         ilabcf(ll)=ilcf
      enddo
      line_cool(14,14)=cofe14
C     *************************************************************
C     *****
c   Fe XV 283 level atom      
C     *****
C     *************************************************************
      ilcf=1415
      call popchianti(14,15,te,XQ,cofe15)
      if(itime.eq.1) nline_min(14,15)=ll
      ll=nline_min(14,15)
      do k=1,kmaxp(14,15)
         ll=ll+1
         if(ipopch.eq.1) then
            wli(ll)=wlin(k)
            cin(ll)=weh(k)/del(ik)
            tauline(ll)=taul(k)
            cinx(14,15,k)=weh(k)/del(ik)
            taulinex(14,15,k)=taul(k)
            wlix(14,15,k)=wlin(k)
         endif
         ilabcfx(14,15,k)=1415
         ilabcf(ll)=ilcf
      enddo
      line_cool(14,15)=cofe15
C     *************************************************************
C     *****
c   Fe XVI 267 level atom      
C     *****
C     *************************************************************
      ilcf=1416
      call popchianti(14,16,te,XQ,cofe16)
      if(itime.eq.1) nline_min(14,16)=ll
      ll=nline_min(14,16)
      do k=1,kmaxp(14,16)
         ll=ll+1
         if(ipopch.eq.1) then
            wli(ll)=wlin(k)
            cin(ll)=weh(k)/del(ik)
            tauline(ll)=taul(k)
            cinx(14,16,k)=weh(k)/del(ik)
            taulinex(14,16,k)=taul(k)
            wlix(14,16,k)=wlin(k)
         endif
         ilabcfx(14,16,k)=1416
         ilabcf(ll)=ilcf
      enddo
      line_cool(14,16)=cofe16
C     *************************************************************
C     *****
c   Fe XVII 267 level atom      
C     *****
C     *************************************************************
      ilcf=1417
      call popchianti(14,17,te,XQ,cofe17)
      if(itime.eq.1) nline_min(14,17)=ll
      ll=nline_min(14,17)
      do k=1,kmaxp(14,17)
         ll=ll+1
         if(ipopch.eq.1) then
            wli(ll)=wlin(k)
            cin(ll)=weh(k)/del(ik)
            tauline(ll)=taul(k)
            ilabcf(ll)=ilcf
            cinx(14,17,k)=weh(k)/del(ik)
            taulinex(14,17,k)=taul(k)
            wlix(14,17,k)=wlin(k)
         endif
         ilabcfx(14,17,k)=1417
      enddo
      line_cool(14,17)=cofe17
C     *************************************************************
C     *****
c   Fe XVIII 337 level atom      
C     *****
C     *************************************************************
      ilcf=1418
      call popchianti(14,18,te,XQ,cofe18)
      if(itime.eq.1) nline_min(14,18)=ll
      ll=nline_min(14,18)
      do k=1,kmaxp(14,18)
         ll=ll+1
         if(ipopch.eq.1) then
            wli(ll)=wlin(k)
            cin(ll)=weh(k)/del(ik)
            tauline(ll)=taul(k)
            cinx(14,18,k)=weh(k)/del(ik)
            taulinex(14,18,k)=taul(k)
            wlix(14,18,k)=wlin(k)
         endif
         ilabcfx(14,18,k)=1418
         ilabcf(ll)=ilcf
      enddo

      line_cool(14,18)=cofe18
C     *************************************************************
C     *****
c   Fe XIX truncated 200 level atom (origial Chianti has 637 levels)
C     *****
C     *************************************************************
      ilcf=1419
      call popchianti(14,19,te,XQ,cofe19)
      if(itime.eq.1) nline_min(14,19)=ll
      ll=nline_min(14,19)
      do k=1,kmaxp(14,19)
         ll=ll+1
         if(ipopch.eq.1) then
            wli(ll)=wlin(k)
            cin(ll)=weh(k)/del(ik)
            tauline(ll)=taul(k)
            cinx(14,19,k)=weh(k)/del(ik)
            taulinex(14,19,k)=taul(k)
            wlix(14,19,k)=wlin(k)
         endif
         ilabcfx(14,19,k)=1419
         ilabcf(ll)=ilcf
      enddo
      line_cool(14,19)=cofe19
C     *************************************************************
C     *****
c   Fe XX truncated 200 level atom (origial Chianti has 637 levels)
C     *****
C     *************************************************************
      ilcf=1420
      call popchianti(14,20,te,XQ,cofe20)
      if(itime.eq.1) nline_min(14,20)=ll
      ll=nline_min(14,20)
      do k=1,kmaxp(14,20)
         ll=ll+1
         if(ipopch.eq.1) then
            wli(ll)=wlin(k)
            cin(ll)=weh(k)/del(ik)
            tauline(ll)=taul(k)
            cinx(14,20,k)=weh(k)/del(ik)
            taulinex(14,20,k)=taul(k)
            wlix(14,20,k)=wlin(k)
         endif
         ilabcfx(14,20,k)=1420         
         ilabcf(ll)=ilcf
      enddo
      line_cool(14,20)=cofe20
C     *************************************************************
C     *****
c   Fe XXI truncated 200 level atom (origial Chianti has 637 levels)
C     *****
C     *************************************************************
      ilcf=1421
      call popchianti(14,21,te,XQ,cofe21)
      if(itime.eq.1) nline_min(14,21)=ll
      ll=nline_min(14,21)
      do k=1,kmaxp(14,21)
         ll=ll+1
         if(ipopch.eq.1) then
            wli(ll)=wlin(k)
            cin(ll)=weh(k)/del(ik)
            tauline(ll)=taul(k)
            cinx(14,21,k)=weh(k)/del(ik)
            taulinex(14,21,k)=taul(k)
            wlix(14,21,k)=wlin(k)
         endif
         ilabcfx(14,21,k)=1421
         ilabcf(ll)=ilcf
      enddo
      line_cool(14,21)=cofe21
C     *************************************************************
C     *****
c   Fe XXII truncated 200 level atom (origial Chianti has 637 levels)
C     *****
C     *************************************************************
      ilcf=1422
      call popchianti(14,22,te,XQ,cofe22)
      if(itime.eq.1) nline_min(14,22)=ll
      ll=nline_min(14,22)
      do k=1,kmaxp(14,22)
         ll=ll+1
         if(ipopch.eq.1) then
            wli(ll)=wlin(k)
c            cin(ll)=weh(k)/del(ik)
            tauline(ll)=taul(k)
            cinx(14,22,k)=weh(k)/del(ik)
            taulinex(14,22,k)=taul(k)
            wlix(14,22,k)=wlin(k)
         endif
         ilabcfx(14,22,k)=1422
         ilabcf(ll)=ilcf
      enddo
      line_cool(14,22)=cofe22
C     *************************************************************
C     *****
c   Fe XXIII truncated 200 level atom (origial Chianti has 637 levels)
C     *****
C     *************************************************************
      ilcf=1423
      call popchianti(14,23,te,XQ,cofe23)
      if(itime.eq.1) nline_min(14,23)=ll
      ll=nline_min(14,23)
      do k=1,kmaxp(14,23)
         ll=ll+1
         if(ipopch.eq.1) then
            wli(ll)=wlin(k)
c            cin(ll)=weh(k)/del(ik)
            tauline(ll)=taul(k)
            cinx(14,23,k)=weh(k)/del(ik)
            taulinex(14,23,k)=taul(k)
            wlix(14,23,k)=wlin(k)
         endif
         ilabcfx(14,23,k)=1423         
         ilabcf(ll)=ilcf
      enddo
      line_cool(14,23)=cofe23
C     *************************************************************
C     *****
c   Fe XXIV
C     *****
C     *************************************************************
      ilcf=1424
      call popchianti(14,24,te,XQ,cofe24)
      if(itime.eq.1) nline_min(14,24)=ll
      ll=nline_min(14,24)
      do k=1,kmaxp(14,24)
         ll=ll+1
         if(ipopch.eq.1) then
            wli(ll)=wlin(k)
            cin(ll)=weh(k)/del(ik)
            tauline(ll)=taul(k)
            cinx(14,24,k)=weh(k)/del(ik)
            taulinex(14,24,k)=taul(k)
            wlix(14,24,k)=wlin(k)
         endif
         ilabcfx(14,24,k)=1424         
         ilabcf(ll)=ilcf
      enddo
      line_cool(14,24)=cofe24
C     *************************************************************
C     *****
c   Fe XXV
C     *****
C     *************************************************************
      ilcf=1425
      call popchianti(14,25,te,XQ,cofe25)
      if(itime.eq.1) nline_min(14,25)=ll
      ll=nline_min(14,25)
      do k=1,kmaxp(14,25)
         ll=ll+1
         if(ipopch.eq.1) then
            wli(ll)=wlin(k)
            cin(ll)=weh(k)/del(ik)
            tauline(ll)=taul(k)
            cinx(14,25,k)=weh(k)/del(ik)
            taulinex(14,25,k)=taul(k)
            wlix(14,25,k)=wlin(k)
         endif
         ilabcfx(14,25,k)=1425         
         ilabcf(ll)=ilcf
      enddo
      line_cool(14,25)=cofe25
C     *************************************************************
C     *****
c   Fe XXVI 
C     *****
C     *************************************************************
      ilcf=1426
      call popchianti(14,26,te,XQ,cofe26)
      if(itime.eq.1) nline_min(14,26)=ll
      ll=nline_min(14,26)
      do k=1,kmaxp(14,26)
         ll=ll+1
         if(ipopch.eq.1) then
            wli(ll)=wlin(k)
c            cin(ll)=weh(k)/del(ik)
            tauline(ll)=taul(k)
            cinx(14,26,k)=weh(k)/del(ik)
            taulinex(14,26,k)=taul(k)
            wlix(14,26,k)=wlin(k)
         endif
         ilabcfx(14,26,k)=1426         
         ilabcf(ll)=ilcf
      enddo
      line_cool(14,26)=cofe26
C     *************************************************************
C     *****
c   N II 12 level atom based on V 10 Chianti
C     *****
C     *************************************************************
      ilcf=402
      call popchianti_new(4,2,te,con2)
      if(itime.eq.1) nline_min(4,2)=ll
      ll=nline_min(4,2)
      do k=1,kmaxp(4,2)
         ll=ll+1
         if(ipopch.eq.1) then
            wli(ll)=wlin(k)
            cin(ll)=weh(k)/del(ik)
            tauline(ll)=taul(k)
            cinx(4,2,k)=weh(k)/del(ik)
            taulinex(4,2,k)=taul(k)
            wlix(4,2,k)=wlin(k)
         endif
         ilabcf(ll)=ilcf
         ilabcfx(4,2,k)=402
      enddo
      c321=con2
      line_cool(4,2)=con2
      
C     *************************************************************
C     *****
c   N II 12 level atom      
C     *****
C     *************************************************************
C     *************************************************************
C     *****
c   C III 53 level atom based on V 10 Chianti
C     *****
C     *************************************************************
      ilcf=303
      call popchianti_new(3,3,te,coc3)
      if(itime.eq.1) nline_min(3,3)=ll
      ll=nline_min(3,3)
      do k=1,kmaxp(3,3)
         ll=ll+1
         if(ipopch.eq.1) then
            wli(ll)=wlin(k)
            cin(ll)=weh(k)/del(ik)
            tauline(ll)=taul(k)
            cinx(3,3,k)=weh(k)/del(ik)
            taulinex(3,3,k)=taul(k)
            wlix(3,3,k)=wlin(k)
         endif
         ilabcf(ll)=ilcf
         ilabcfx(3,3,k)=303
      enddo
      line_cool(3,3)=coc3
C     *************************************************************
C     *****
c   N III 53 level atom based on V 10 Chianti
C     *****
C     *************************************************************
      ilcf=403
      call popchianti_new(4,3,te,con3)
      if(itime.eq.1) nline_min(4,3)=ll
      ll=nline_min(4,3)
      tot1750=0.
      do k=1,kmaxp(4,3)
         ll=ll+1
         if(ipopch.eq.1) then
            wli(ll)=wlin(k)
            cin(ll)=weh(k)/del(ik)
            tauline(ll)=taul(k)
            cinx(4,3,k)=weh(k)/del(ik)
            taulinex(4,3,k)=taul(k)
            wlix(4,3,k)=wlin(k)
            if(wlin(k) > 1740. .and. wlin(k) < 1760.) then
               tot1750=tot1750 + cinx(4,3,k)
            endif
         endif
         ilabcf(ll)=ilcf
         ilabcfx(4,3,k)=403
      enddo
      line_cool(4,3)=con3
c  c331=con3
C     *************************************************************
C     *****
c   N IV 53 level atom based on V 10 Chianti
C     *****
C     *************************************************************
      ilcf=404
      call popchianti_new(4,4,te,con4)
      if(itime.eq.1) nline_min(4,4)=ll
      ll=nline_min(4,4)
      do k=1,kmaxp(4,4)
         ll=ll+1
         if(ipopch.eq.1) then
            wli(ll)=wlin(k)
            cin(ll)=weh(k)/del(ik)
            tauline(ll)=taul(k)
            cinx(4,4,k)=weh(k)/del(ik)
            taulinex(4,4,k)=taul(k)
            wlix(4,4,k)=wlin(k)
         endif
         ilabcf(ll)=ilcf
         ilabcfx(4,4,k)=404
      enddo
c      c331=con3
      line_cool(4,4)=con4
      
C     *************************************************************
C     *****
c  Collex from Mewe and from Gaetz and Salpeter
C     *****
C     *************************************************************

      cxtot = 0.
      do k=1,14
         cxtot = cxtot + fextot(k)
      enddo

C     *************************************************************
C     *****
C     CALCULATE POPULATIONS AND LINE EMISSION FROM O I
C     *****
C     *************************************************************

      ZEL=DEL(IK)-ABN(5)*XO(2)
      IFPOP=0
      XQ=XO(1)+XO(2)
      DO I=1,65
         SR(I)=0.
         WEH(I)=0.
      enddo      
      DO  I=1,12
         DO  J=I+1,13
            EMCA(J,I)=0.
         enddo
      enddo
      TOI=TE
      IF(TE.LT.TLIMO) TOI=TLIMO
C      IF(ABN(5).GT.1.E-5) THEN
      ipopoi=0
      IF((ABN(5)*XO(1).GT.1.E-10.and.te.lt.1.e5).or.itime.eq.1) then
         ipopoi=1
         CALL POPO_new(RS,TOI,ZEL,XQ,PHO,IFPOP)

      else
         do j=1,nl
            xnq(2,j)=0.
         enddo
      endif
      if(xo(2).gt.0.1) then
         recto=0.
         do k=1,9
            recto=recto+reco(k)
         enddo
      endif
      IF(IFPOP.EQ.1) GOTO 8278
C     31 = 6300-64, 32 = 2972, 33 = 5572, 34 = 1356, 35 = 1302
C     36 = 7774   , 37 = 8446, 38 = 1027, 39 =  990
      DO IJ=1,4
         SRED(IJ+30)=SR(IJ)
         CIN(IJ+30)=WEH(IJ+2)/DEL(IK)
      enddo
      CIN(35)=WEH(5)/DEL(IK)
      tauline(35)=taul(5)
      CIN(36)=WEH(8)/DEL(IK)
      tauline(36)=taul(8)
      CIN(37)=WEH(11)/DEL(IK)
      tauline(37)=taul(11)
      CIN(38)=WEH(12)/DEL(IK)
      tauline(38)=taul(12)
      CIN(39)=WEH(16)/DEL(IK)
      tauline(39)=taul(16)
      SRED(35)=SR(5)
      SRED(36)=SR(8)
      SRED(37)=SR(11)
      SRED(38)=SR(12)
      SRED(39)=SR(16)

      kmaxp(5,1) = 39
      if(itime.eq.1) nline_min(5,1)=ll
      ll=nline_min(5,1)
      do k=1,kmaxp(5,1)
         if(k.le.9) then
c            cin(30+k)=weh(k)/del(ik)
            tauline(30+k)=taul(k)
c            if(itime.le.888888) wli(30+k)=wlin(k)
            if(ipopoi.eq.1) wli(30+k)=wlin(k)
            ilabcf(30+k)=501
         else
            ll=ll+1
c            cin(ll)=weh(k)/del(ik)
            tauline(ll)=taul(k)
            if(ipopoi.eq.1) wli(ll)=wlin(k)
            ilabcf(ll)=501
         endif
         cinx(5,1,k)=weh(k)/del(ik)
         taulinex(5,1,k)=taul(k)
c!!         wlix(5,1,k)=wlin(k)         
         ilabcfx(5,1,k)=501         
      enddo

      COOI=0.
      DO I=1,12
         DO J=I+1,13
            COOI=COOI+EMCA(J,I)
         ENDDO
      ENDDO
      C111=COOI
      line_cool(5,1)=cooi

      IF(TE.LT.TLIMO) THEN
C     LOW TEMPERATURE APPROX.
         T3OI=TOI/1.E3
         O21=.0151*T3OI**1.31
         O31=.00184*T3OI**1.32
         O32=.031*T3OI**.534
         Z=ABN(5)*XO(1)
         CALL FORB3(4,5,1,1.957D0,2.223D0,8.45D-3,7.35D-2,1.22D0,O21,
     &        O31,O32,9.D0,5.D0,1.D0,FBO(10),FBO(11),FBO(12),TOI,Z,F)
         DO NM=10,12
            FBO(NM)=ABN(5)*XO(1)*FBO(NM)
         ENDDO
         CALL OXREC(OXR,TE)
         DO IJ=1,8
            OXR(IJ)=ABN(5)*xion(ik,5,2)*OXR(IJ)
         enddo
         CIN(31)=FB(10)*CIN(31)/FBO(10)
         CIN(32)=FB(11)*CIN(32)/FBO(11)
         CIN(33)=FB(12)*CIN(33)/FBO(12)
         CIN(34)=OXR(5)
         CIN(35)=OXR(1)
         CIN(36)=OXR(6)
         CIN(37)=OXR(2)
         CIN(38)=0.
         CIN(39)=0.
         if(itime.eq.1) nline_min(5,1)=ll
         ll=nline_min(5,1)
         do k=1,kmaxp(5,1)
            if(k.le.9) then
               ilabcf(30+k)=501
            else
               ll=ll+1
               cin(ll)=0.
               if(itime.le.888888) wli(ll)=wlin(k)
               ilabcf(ll)=501
            endif
         enddo


         DO IJ=1,8
            SRED(IJ+30)=0.
         enddo
         COOI=COOI*(FB(10)+FB(11)+FB(12))/(FBO(10)+FBO(11)+FBO(12))
      ENDIF

CF1016
C     RECOMB. EMISSION
      RECEMO=0.
      NP1H=10
      DO J=1,9
         CALL RECOMB(2,J,TE,ALOXI,RI)      
         RECEM(2,J)=ABN(5)*BOL(NP1H)*ALOXI*(E00-EI(J))*1.602e-12
         RECEMO=RECEMO+RECEM(2,J)
      ENDDO
      C111=COOI

c!!!! 333200!!
      CALL OXREC(OXR,TE)
      DO IJ=1,8
         OXR(IJ)=ABN(5)*xion(ik,5,2)*OXR(IJ)
      enddo

7327  IF(ABN(5).LT.1.E-5.OR.TE.LT.333200.) goto 7127
      C111=0.
      DO IJ=31,39
         C111=C111+CIN(IJ)
      enddo
 7127 continue
C     *************************************************************
C     *****
C     CALCULATE POPULATIONS AND LINE EMISSION FROM Si I
C     *****
C     *************************************************************
      ZEL=DEL(IK)-ABn(10)*(Xsi(1)+2.*Xsi(2))
      IFPOP=0
      DO I=2,nsi1
        DO J=1,I-1
          WOBS(I,J)=0.
          EMCA(I,J)=0.
        ENDDO
      ENDDO

      XQ=Xsi(2)+Xsi(3)

      call POPsi_i(RS,TE,ZEL,XQ,PHHE,IFPOP)

      COSI1=0.
      SI1TOT=0.
      DO I=1,NSI1
         DO J=1,NSI1
c            WOBSI1(J,I)=WOBS(J,I)
            SI1TOT=SI1TOT+WOBS(J,I)/DEL(IK)
            COSI1=COSI1+EMCA(J,I)
         enddo
      enddo

      do k=1,kmaxp(10,1)
         wlix(10,1,k)=wlin(k)
         cinx(10,1,k)=weh(k)/(del(ik))
         ilabcfx(10,1,k)=1001
      enddo
      line_cool(10,1)=cosi1
            
C     *************************************************************
C     *****
C     CALCULATE POPULATIONS AND LINE EMISSION FROM FE I
C     *****
C     *************************************************************
      ZEL=DEL(IK)-ABN(14)*(XF(2)+2.*XF(3))
      IFPOP=0
c      DO I=1,NFEL
c        SRED(81?+I)=0.
c        SR(I)=0.
c        IF(I.LE.200) WEH(I)=0.
c      ENDDO
      DO I=2,NFEi
        DO J=1,I-1
          WOBS(I,J)=0.
          EMCA(I,J)=0.
        ENDDO
      ENDDO
      XFEI=1.-XF(2)-XF(3)-XF(4)
      TEFEI=TE
c!!
c      IF(TE.LT.500.) TEFEI=500.
      IF(TE.LT.100.) TEFEI=100.
c!!
      IF(IFEMUL.LT.0) THEN
        EPSFEi=1.
      ELSE
        EPSFEi=1.E-10
      ENDIF
      XQ=XF(1)+XF(2)
      ifeim=0
      IF((ABN(14)*XF(1).GT.EPSFEi.AND.TE.GE.TEFEi).or.itime.eq.1) THEN
         ifeim=1
         NFEiS=NFEi
c        IF(TE.LT.2.E3.AND.NFEiS.GT.60) NFEi=46
c        IF(TE.LT.1.E3.AND.NFEiS.GT.20) NFEi=20
c        IF(TE.LT.0.5E3.AND.NFEiS.GT.10) NFEi=9
         IF(NFEi.EQ.NFEiS) THEN
            NQ=NFEiS+1
            IQ=1
         ELSEIF(NFEi.EQ.46) THEN
            NQ=47
            IQ=2
         ELSEIF(NFEi.EQ.20) THEN
            NQ=21
            IQ=3
         ELSEIF(NFEi.EQ.9) THEN
            NQ=10
            IQ=4
         ENDIF
         DO K=1,NQ
            IF(INIFEi(IQ).EQ.1) BOLFEi(K)=BOLFEiS(IQ,K)
         ENDDO
         IF(INIFEi(IQ).EQ.0) BOLFEi(NFEi+1)=BOLFEi(NFEiS+1)
         IF(INIFEi(IQ).EQ.1) TOLDFEi=TOLDFEiS(IQ)
         IF(IQ.EQ.1) INIFE(1)=1
         IF(IQ.EQ.2) INIFE(2)=1
         IF(IQ.EQ.3) INIFE(3)=1
         IF(IQ.EQ.4) INIFE(4)=1
      endif

c Fe I
      ilcf=1401
      z=abn(14)*xf(1)
      XQ=xf(2)+xf(3)
      call popchianti_new(14,1,te,cofe1)
      do k=1,kmaxp(14,1)
         if(ipopch.eq.1) then
            cinx(14,1,k)=weh(k)/del(ik)
            taulinex(14,1,k)=taul(k)
            wlix(14,1,k)=wlin(k)
         endif
         ilabcfx(14,1,k)=1401
      enddo
      line_cool(14,1)=cofe1
      if(ifpop.eq.0.and.ifeim.eq.1) then
         
         FEIFS_C=EXP(-0.6499/T3)/(T3**0.358*(DEN*XEL/1.E7)**0.96)
         FEIFS_C=FEIFS_C*ABN(14)*XF(1)
         
         FEIFS_C = cofei/FEIFS_C
         fecon = cofei/feian
      endif

      IF(ABN(14)*XF(1).LT.EPSFEi.or.ifpop.eq.1) COFEI=CIN18
c      IF(ABN(14)*XF(1).GT.EPSFEi.and.ifpop.eq.0) CIN(18)=FEITOT
C     *************************************************************
C     *****
C     CALCULATE POPULATIONS AND LINE EMISSION FROM FE II
C     *****
C     *************************************************************


c Fe II
      ilcf=1402
      z=abn(14)*xf(2)
      call popchianti_new(14,2,te,cofe2)
      do k=1,kmaxp(14,2)
         if(ipopch.eq.1) then
            cinx(14,2,k)=weh(k)/del(ik)
            taulinex(14,2,k)=taul(k)
            wlix(14,2,k)=wlin(k)
         endif
         ilabcfx(14,2,k)=1402
      enddo
      line_cool(14,2)=cofe2
      nlinetot=ll

      if(ifpop.eq.0) then

         FEiIFS_C=EXP(-0.6499/T3)/(T3**0.358*(DEN*XEL/1.E7)**0.96)
         FEIiFS_C=FEiIFS_C*ABN(14)*XF(2)
         
         FEIiFS_C = cofeii/FEIiFS_C
      endif
C     FOR TE < 500 K SCALE POPFE AT 500 K WITH THE TEMPERATURE DEP.
C           FR!makOM AXELROD
c!!
      IF(TE.LT.TEFE) THEN
         TFE3=0.5
         FEII500=EXP(-0.5077/TFE3)/TFE3**0.162
         FEIIT=EXP(-0.5077/T3)/T3**0.162
         COFEIIs=COFEII*FEIIT/FEII500
      ENDIF
      IF(ABN(14)*XF(2).LT.EPSFE) COFEII=CIN19
C     **************************************************************
C     *****
C     PHOTOIONIZATION HEATING
C     *****
C     **************************************************************
8372  HYOLD=HY

      CALL EMISS(TE,del(ik),den)
      iqq=5
      call opac_source(iqq,0,dr,den,TE,del(ik))
c      CALL SPEC_pl(TE)
      ZOLD=ZQ
      ZQ=DNEW
c      call rateaug(1,1)
c      call rateaug(2,2)
      call rateaug(6,6)
      call rateaug(7,7)
      call rateaug(8,8)
      call rateaug(10,10)
      call rateaug(12,3)
      call rateaug(14,14)
      call rateaug(16,16)
      call rateaug(18,7)
      call rateaug(26,26)

      
C     CALCULATE NON-THERMAL IONIZATION
      CALL NONTH(XEL,FHEAT,FIH,FEH,FIHE,FEHE)



c     TOTAL HEATING
      GATOT=0.

      do iel=3,14
         do iz=1,nionel(iel)
            dga = gea(iel,iz)*abn(iel)*xion(ik,iel,iz)
            gatot = gatot + gea(iel,iz)*abn(iel)*xion(ik,iel,iz)
         enddo            
      enddo            
      do iel=3,14
         do iz=1,nionel(iel)

            ge(iel,iz) = geion(iel,iz)/den
            ze(iel,iz) = zea(iel,iz)/den
            zk(iel,iz) = zka(iel,iz)/den

         enddo       
      enddo            

      MQ=MQ+1
C      IF(ZQ.NE.0.) DZ=ABS(ZQ-ZOLD)/ZQ
      MI=1


C     OTS
      IF(OTSP(4).LE.1) POX(5)=0.

      DO L=1,6
         PC(L)=FHEAT*ABN(3)*xion(ik,3,l)*GE(3,l)
         dheat(3,l)=FHEAT*ABN(3)*xion(ik,3,l)*GE(3,l)
      enddo

C     OTS
      IF(OTSP(5).LE.1) PC(5)=0.

      DO L=1,7
         Pn(L)=FHEAT*ABN(4)*xion(ik,4,l)*GE(4,l)
         dheat(4,l)=FHEAT*ABN(4)*xion(ik,4,l)*GE(4,l)
      enddo

C     OTS
      IF(OTSP(6).LE.1) PN(1)=0.

      PO=0.
      POT=0.
      DO J=1,9
         POJ=FHEAT*ABN(5)*XNQ(2,J)*PHEAT(2,J)/DEN
         POTJ=FHEAT*ABN(5)*XNQ(2,J)*PHEATT(2,J)/DEN
         PO=PO+POJ
         POT=POT+POTJ
      ENDDO
      DO L=1,8
         POX(L)=FHEAT*ABN(5)*xion(ik,5,l)*GE(5,L)
         dheat(5,l)=FHEAT*ABN(5)*xion(ik,5,l)*GE(5,l)
      enddo

      pne=0.
      do l=1,10
         PNEI(l)=FHEAT*ABN(6)*Xne(L)*GE(6,L)
         dheat(6,l)=FHEAT*ABN(6)*Xne(L)*GE(6,L)
         PnE=pnei(l) + PnE
      enddo

      PNAT=0.
      DO L=1,3
         PNAT=PNAT+FHEAT*ABN(7)*xion(ik,7,l)*GE(7,l)
         dheat(7,l)=FHEAT*ABN(7)*xion(ik,7,l)*GE(7,l)
      enddo

      PMG=0.
      DO L=1,3
         PMG=PMG+FHEAT*ABN(8)*xion(ik,8,l)*GE(8,l)
         dheat(8,l)= FHEAT*ABN(8)*xion(ik,8,l)*GE(8,l)
      enddo

      PAL=0.
      DO L=1,4
         PAL=PAL+FHEAT*ABN(9)*Xion(ik,9,L)*GE(9,l)
         dheat(9,l)=FHEAT*ABN(9)*Xion(ik,9,L)*GE(9,l)
      enddo

      PSI=0.
      DO L=1,14
         PSI=PSI+FHEAT*ABN(10)*xion(ik,10,L)*GE(10,l)
         dheat(10,l)=FHEAT*ABN(10)*xion(ik,10,L)*GE(10,l)
      enddo
      psu=0.
      DO L=1,16
         PSU=PSU+FHEAT*ABN(11)*xion(ik,11,l)*GE(11,L)
         dheat(11,l)=FHEAT*ABN(11)*xion(ik,11,l)*GE(11,L)
      enddo

      par = 0.
      do l=1,8
         par=par + FHEAT*ABN(12)*xar(L)*GE(12,l)
         dheat(12,l) = FHEAT*ABN(12)*xar(L)*GE(12,l)
      enddo
      PCAL=0.
      DO L=1,3
         PCAL=PCAL+FHEAT*ABN(13)*xion(ik,13,L)*GE(13,L)
         dheat(13,l)=FHEAT*ABN(13)*xion(ik,13,L)*GE(13,L)
      enddo
      PFE=0.
      DO  L=1,15
         PF(L)=FHEAT*ABN(14)*XF(L)*GE(14,l)
         dheat(14,l)=FHEAT*ABN(14)*XF(L)*GE(14,l)
         PFE=PF(L)+PFE
      enddo

c calculate the highest heating rate
      
      dhmax=0.
      do iel=3,14
         do iz=1,nionel(iel)

            if(dheat(iel,iz).gt.dhmax) then
               dhmax=dheat(iel,iz)
               izmax=iz
               ielmax=iel
            endif
         enddo
      enddo

C     **************************************************************
C     *****
C     TOTAL COOLING AND HEATING RATES
C     *****
C     **************************************************************

       COCAR=C211+C212+coc2+C231+C232+C241 + coc5 + coc6
       coolx(3,1)=c211+c212
       coolx(3,2)=coc2
       coolx(3,3)=c231+c232
       coolx(3,4)=c241
       coolx(3,5)=coc5
       coolx(3,6)=coc6

      CON=con2+con3+con4+C351+C311+ con6 + con7
      coolx(4,1)=c311
      coolx(4,2)=con2
      coolx(4,3)=con3
      coolx(4,4)=con4
      coolx(4,5)=c351
      coolx(4,6)=con6
      coolx(4,7)=con7            

      COO=C161+C152+C151+C131+c141+C142+C143+C111+coo2+C132+coo7+coo8
      coolx(5,1)=c111
      coolx(5,2)=coo2
      coolx(5,3)=c131+c132
      coolx(5,4)=c141+c142+c143
      coolx(5,5)=c151+c152
      coolx(5,6)=c161
      coolx(5,7)=coo7            
      coolx(5,8)=coo8

      
      CONE=C1322+coneiii+C1351+C1341+cone6+cone7+cone8+cone9+cone10
      coolx(6,2)=c1322
      coolx(6,3)=coneiii
      coolx(6,4)=c1341
      coolx(6,5)=c1351
      coolx(6,6)=cone6
      coolx(6,7)=cone7
      coolx(6,8)=cone8
      coolx(6,9)=cone9
      coolx(6,10)=cone10
      CONA=C1011
      coolx(7,1)=c1011
      COMG=C511+C512 + comg2 + comg10 + comg11 + comg12
      coolx(8,1)=c511+c512
      coolx(8,2)=comg2
      coolx(8,10)=comg10
      coolx(8,11)=comg11
      coolx(8,12)=comg12
      COSI=C411+C412+C421+C422+C431+C441 + cosi8 + cosi9 + 
     &     cosi10 + cosi11 + cosi12 + cosi13 + cosi14
      coolx(10,1)=c411+c412
      coolx(10,2)=c421+c422
      coolx(10,3)=c431
      coolx(10,4)=c441
      coolx(10,8)=cosi8
      coolx(10,9)=cosi9
      coolx(10,10)=cosi10
      coolx(10,11)=cosi11
      coolx(10,12)=cosi12
      coolx(10,13)=cosi13
      coolx(10,14)=cosi14

      COS=C1211+C1261+C1241+C1242 + cos2 +cos3 +cos13 +
     &     cos14 + cos15 + cos16
      coolx(11,1)=c1211
      coolx(11,2)=cos2
      coolx(11,3)=cos3
      coolx(11,4)=c1241+c1242
      coolx(11,5)=c1251
      coolx(11,6)=c1261
      coolx(11,13)=cos13
      coolx(11,14)=cos14
      coolx(11,15)=cos15
      coolx(11,16)=cos16
      coar = c1222 + coariii + coariv + coarv + coarvi
      coolx(12,2)=c1222
      coolx(12,3)=coariii
      coolx(12,4)=coariv
      coolx(12,5)=coarv
      coolx(12,6)=coarvi
      cocal = cocal + coca9 + coca10 + coca11 + coca12 + coca13 + coca14
      coolx(13,1)=0.
      coolx(13,2)=cocal
      coolx(13,9)=coca9
      coolx(13,10)=coca10
      coolx(13,11)=coca11
      coolx(13,12)=coca12
      coolx(13,13)=coca13
      coolx(13,14)=coca14
c Fe I NOT included! cin18 includes Axelrod's estimate
      COFE=CIN18+COFEII+cofe3+FEIV+FEV + cofe7 + cofe10 + 
     &     cofe11 + cofe12 + cofe13 + cofe14 + cofe15 + 
     &     cofe16 + cofe17 + cofe18
      coolx(14,1)=cofei
      coolx(14,2)=cofeii           
      coolx(14,3)=cofe3
      coolx(14,4)=cofeiv
      coolx(14,5)=cofev
      coolx(14,6)=0.
      coolx(14,7)=cofe7
      coolx(14,8)=0.
      coolx(14,9)=0.
      coolx(14,10)=cofe10
      coolx(14,11)=cofe11
      coolx(14,12)=cofe12
      coolx(14,13)=cofe13
      coolx(14,14)=cofe14
      coolx(14,15)=cofe15
      coolx(14,16)=cofe16
      coolx(14,17)=cofe17
      coolx(14,18)=cofe18

      coolt=0.
      do iel=3,14
         do ion=1,27
            coolt=coolt+line_cool(iel,ion)
c            write(6,9267)iel,ion,line_cool(iel,ion),coolt
 9267       format('line cooling ',2i5,1pe12.3,10e12.3)
         enddo
      enddo


      COOL=FF+RECCA+RECOX+RECSI+RECS+RECFE + coolt
      POXY=PO
c     obs! O I SEPARATELY
      DO L=2,8
         POXY=POXY+POX(L)
      ENDDO
      PCA=0.0E0
      DO L=1,6
         PCA=PCA+PC(L)
      ENDDO
      PNA=0.
      DO L=1,7
         PNA=PNA+PN(L)
      ENDDO

      HEAT=CO+POXY+PCA+PNA+PFE+PMG+PAL+PSI+pne+psu+
     &     PCAL+PNAT+PAR


      temin=100.
      stop_tmin=0
      if(te.le.temin.and.rcomm.le.0.) then
         write(6,*)'obs!!!  stop cooling at T= ',te,cool
         cool = 0.
         stop_tmin=1
      endif


      RAD=(HEAT-(HY+COOL)*DEL(IK))
9323  format(5e12.4)
C     APPROX. REC. EMISSION FROM C, O, S, SI, FE
      RECAPPC=DEL(IK)*RECCA*1.307e5/(1.5*TE)
      RECAPPO=DEL(IK)*RECOX*1.581e5/(1.5*TE)
      RECAPPSI=DEL(IK)*RECSI*0.946e5/(1.5*TE)
      RECAPPS=DEL(IK)*RECS*1.194e5/(1.5*TE)
      RECAPPFE=DEL(IK)*RECFE*0.9136e5/(1.5*TE)

      if(ik.le.300) ioterm=1
      if(ik.lt.300) ioterm=1
      iout=6

      ioterm=1

      XELEC=DEL(IK)

      COOLQ=COOL*DEL(IK)
      tion=time*den
      WRITE(iout,7362)kradius,TE,DEL(IK),r(kradius),den,HEAT,COOL,cool*del(ik),
     &     tion,RAD
           

8278   IF(IFPOP.EQ.1) WRITE(6,*)' NO CONVERGENCE IN POP-ROUTINE'
 7362 FORMAT(1X,'T,X,R(k),DEN,H,C,C*X,TION,R',i5,1pe14.6,e14.6,e17.9,10E14.6)
 7366 FORMAT(1X,'HE',5E12.4)

 7365 FORMAT(1X,'CO',1pe12.4,8E12.4)
9276  FORMAT('CF ',7E10.3)
      l=0

      RETURN
       END

      subroutine read_coll
      implicit real*8(a-h,o-z)
      common/ca4/teca4(50),coll_ca4(3,50)
      open(11,file='./ATDAT/Ca_IV_coll_strengths_Nahar_v2.dat',status='old')
      do i=1,50
         read(11,*)ii,teca4(i),(coll_ca4(k,i),k=1,3)
      enddo
      close(11)
      return
      end
