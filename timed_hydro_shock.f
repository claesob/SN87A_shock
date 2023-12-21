      subroutine trans
c     index in trans:
c     i = index of shell. here i = 2
c     k = index of shell in rad. transf. calc. k can be up to md = 350
c     ik = k

      IMPLICIT REAL*8(A-H,O-Z)
      character*72 text1,text2,text3,TEXT4,TEXT5,INMOD
      CHARACTER*8 LAB(200)
      parameter (nlp=30000)
      CHARACTER*18 LINELAB(300),linelabex(nlp)
      character*7 labmol
      REAL*8 MTOT,MNI56
      include 'param'
      PARAMETER (NL=340,NLP1=NL+1)
      PARAMETER (NFEL=3000)
      common/inputv/RMCEN,DELM,R1,T00,A1IN,B1IN,dminit,rminner,revmass,
     &     EJDENS,TECOOL,totcolumn,FILLING,MODNR,NMAX,MAX,ICEN,IREV,
     &     NPRINT,IPAR,MQMAX,NIND,IDEP,IRINNER,ISTEPION
      COMMON/TEXT/TEXT1,TEXT2,TEXT3,TEXT4,TEXT5
      COMMON/MOD/INMOD
      COMMON/IOINF/IOSPEC,IOLEVEL,IOEMISS,IOTERM,IOION,IOELITER
      COMMON/THBAL/TMIN,TMAX,IBAL,IKBAL(1000)
      COMMON/CONT/CONR(NL),CONT(NL)
      COMMON/SECX/CSEC(20),DCSDR(20),CISEC(20),DCISDR
      COMMON/DENS/DEN0,R0,RN
      COMMON/TAUFBO/TAFB(30,3)
      COMMON/QSOM/qso,IAGN,ISTAT,ISOBOL
      COMMON/INSPEC/IBLACKB,IPULSSP,IFERNET,INLR,IFERMAT,ICRAB,
     &     IFREEFREE,imewe,iobstar
      common/initmw/initmw,initr2
      COMMON/AGNPARA/DECON,GAMMION,IPRESS,IDEN,ICDEN
      COMMON/PULSAR/RLTOT,ALFA,EMIN,EMAX
      COMMON/FNORM/FNORM
      COMMON/LOWION/ILOWION
      COMMON/SPH/ISPH
      COMMON/INUT/IUY
      COMMON/RQW/TEFF,RQ
      COMMON/TPAR/RIN,DRQ,R1Q,TDAYS
      COMMON/INI/INITH,INIT,INITCA,INITHE,INITFE,initfei
      COMMON/copop/xco(NL),initcop
      common/croinit/initcr
      COMMON/ITER/ITE
      COMMON/MPAR/MTOT,MNI56,VEXP
      COMMON/TERMA/TERMAL
      COMMON/QSPEC/GAMMA,ALQSO
      COMMON/LITER/N
      COMMON/IND/Ik
      COMMON/mdIND/kradius
      COMMON/FRE/NINQ,JMIN,JJ
      COMMON/ITSPH/NITSPH
      COMMON/THER/A1,B1,TIN,E10,E20
      COMMON/SPOT/OTSP(16)
      COMMON/COL/RTE(4),FF,HY,HE,C131,C142,C151,C161,C231,C241,COH,
     &     COHE,C351,C321,C331,C341,C222,C232,C221,C332,C441
      COMMON/EQUIH/FRH(2,401),WEH(401),SR(NFEL+nlp),WOBS(NL,NL)
      COMMON/GAMPAR/GAMLUM,GAELH,GAHE
      COMMON/T/TES,SS
      COMMON/REHEL/REHE21,REHE22,REC31,REN41
      COMMON/CHG/CHGH,CHGHE,CHTE
      COMMON/GSREC/ALGS(100),ALEX(100),ALTOT(100),RECEX(100),RECGS(100) 
      COMMON/HEA/CO,PHEO,PHEI,PO(8),PC(6),PN(7),PMG,PSI,PFE,
     &     pnei(10)
      COMMON/INT/FL(2,NE1:NE2),SI(14,27,NE1:NE2),E1(NE1:NE3),E(NE1:NE3)
      COMMON/DIFFP/FH0(NE1:NE2),FHD(NE1:NE2)
      COMMON/TRES/ EL(14,27),EK(14,27)
      COMMON/SIK/SK(14,27,NE1:NE2)
      COMMON/PHQ/ZE(14,27),GE(14,27),ZK(14,27)
      common/chec/gec(14,27)
      COMMON/OTS/NOTS
      COMMON/EDDFLUX/EDDFLUX(NE1:NE2),SURFJ(NE1:NE2)
      COMMON/DIF/EM(MD,NE1:NE2),TAU(MD,NE1:NE2),TAUTOT(MD,NE1:NE2),
     &     EMC(MD,NE1:NE2)
      COMMON/REC/AL2(16)

      parameter (mh=2,mhe=3,mc=7,mn=8,mo=9,mne=11,mna=12,mmg=13,
     &     mal=14,msi=15,ms=17,mar=19,mca=21,mfe=27,mni=29,
     &     mmat=mh+mhe+mc+mn+mo+mne+mna+mmg+mal+msi+ms+mar+mca+mfe+mni,
     &     ispecy=15)
      common/abc1/alc(mc),aln(mn),alo(mo),alne(mne),alna(mna),
     &     almg(mmg),alal(mal),alsi(msi),alsu(ms),alar(mar),alca(mca),
     &     alfe(mfe),alni(mni)
      COMMON/ABC/AL(16)
      COMMON/HYDROS/HSCALE,DEN1,XQL,AMEAN
      COMMON/RADIE/R(0:MD)
      COMMON/DXA/DR(MD)
      COMMON/ELEC/DEL(MD)
      COMMON/PHY/DEN(MD)
      COMMON/MASSES/RM(MD)
      common/abl/abn(15)
      COMMON/TEM/TE(MD)
      COMMON/NHY/NH
      COMMON/NLEV/e00,NIONI,NHY,NP1H,NMA,NMI
      COMMON/BOLD/TOLD,TOLDH,TOLDCA,TOLDO,TOLDHE,TOLDFE,toldfei,BB(NL),
     &     BOLH(NL),BOL(NL),BOLCA(NL),BOLHE(NL),BOLFE(NL),bolfei(nl)
      COMMON/HYPOP/XNH(NL)
      COMMON/HPOP/XN(6,NL),XN1,XN2,XN3
      COMMON/A5/TAUE,ALFAQ,EN,TSN,XL40,TEXA,RMAX
      COMMON/A7/C3,C33
      COMMON/A10/AV
      COMMON/A11/R11,CV,FLUX
      COMMON/A12/ITER
      COMMON/PHEAT/PHE(NL)
      COMMON/A17/RECNET(NL),RECNDT(NL),RECCO(NL),REC(NL),
     &     PHOT(NL),DRDT(NL),RECT(6,NL),PHET(NL)
      COMMON/HRAT/RHYD,ZEL,XEL,HYR,HEAT,COOL
      common/ctime/tcool,trad
      COMMON/A14/CQ(NL,NL),CI(NL),G(NL),EI(NL),AQ(NL,NL),WL(NL,NL)
     &     ,DCDT(NL,NL),DCIDT(NL)
      COMMON/COLH/DXI,COLH(6,NL),COLHVT(6,NL),COLHVTT(6,NL),TAUB
      COMMON/A19/EMIS(NL,NL),ESC(NL,NL),TTOT(NL,NL),TOP(NL,NL)
      COMMON/A31/TOPI(10,10,10),EMHI(NL,NL)
      COMMON/COLD/COLTOT(14,27),COLTOTT(14,27),SRED(NFEL+nlp)
      COMMON/SPECT/TEL,FD(MD,NE1:NE2),F0(NE1:NE2),IPARA
      COMMON/DTAU/FLUXED(NE1:NE2)
      COMMON/IOPAR/IOSHELL,IORAD1,IORAD2,IOFLUX
      COMMON/CINOUT/INOUT,IPULS
      COMMON/CICS/RSHOCK,ics,ics87a
      common/revs/terev,tecs
      COMMON/FE3/WLFE3(30,30),AFE3(30,30),EMFE3(30,30)
      COMMON/FE4/WLFE4(30,30),AFE4(30,30),EMFE4(30,30)
      COMMON/PROFFE/WLFEIII(300),WLFEIV(300),FEIII(300,MD),FEIV(300,MD)
      COMMON/FELEV/NFEII,nfei
      COMMON/WFE/WOBFE(NL,NL),wobfei(nl,nl),IQPOPH
      COMMON/EMHY/RECEM(6,NL),TWOPH,TWOPHHEI,COHI,PHQ,PHT,POQ,POT
      common/heiires/x2he,he2coll12,he2jla,he2j2g,he2jrec(4),
     &     he2jem(7,4)
      common/test1/qa(10,20,20),rateph(10,20),tsav(10,20,20),esav(10,20,20)
      common/escfe/esfe2(nl,nl),taufe2(nl,nl)
      common/initstr/initstruc
      common /initcoll/initcoll
      common/ikoldc/ikold
      common/hydphot/phhi(nl)
      common/coolants/cl(nlp)
      common/pescpd/pdbe(nlp,2)
      COMMON/LYBFLUOR/OIFLUOR,OIFLUORH,OIFLUORO,OILYB,OIREC
      COMMON/A26/SIGO(20,NE1:NE2),SIGCA(10,NE1:NE2),SIGHE(16,NE1:NE2),
     &     SIGFE(NL,NE1:NE2),sigfei(nl,ne1:ne2),SIGH(50,NE1:NE2),SIGHE2(50,NE1:NE2)
      COMMON/FESAVE/INIFE(4),inifei(4),BOLFES(4,NL),TOLDFES(4),
     &     bolfeis(4,nl),toldfeis(4)
      common/deltar/drscgs
      common/linepoint/nlines,jline(nlp),iion(nlp)
      common/linephoto/phline(nlp)
      common/escapemax/be00
      common/contpoint/otscon(100),ncont,jcont(100),iioncont(100)
      common/oiiilines/em4959,em5007
      common/multisol/imultih,imultihe,imultio,imultica,imultifei,
     &     imultifeii
      parameter(nspec=91)
      common/initwl/iwl,initprl
      common/inidi/initdi
      common/emind/iemiss
      common/rinner/rmin,rmaxq
      COMMON/EDDFL/FH(MD,NE1:NE2),FHC(MD,NE1:NE2)
      dimension totfl(ne1:ne2),contfl(ne1:ne2),tototsl(ne1:ne2)
      DIMENSION COLD(14,27),COLDN(14,27),coltottemp(14,27),FBI(110,MD)
      DIMENSION DLO(50)
      DIMENSION FL0(2,NE1:NE2),WREC(30),WRHE(10),HER(10),OXR(10) 
      common/totemiss2/RECHETOT,TWOPHHEIT,FREEFREE,TWOPHTOT,HEATTOT, 
     &     COOLTOT,GAMTOT,REHE1,REN1,REC1,REHE2,REHE3,REHE4
      common/tot_line_em/tot_intens(md,14,26,401)
      common/taushok/taush(md)
      DIMENSION EMTOT(NL,NL),CNRTOT(NL),
     &     DWFEII(NL,NL),CITOTFE(5000),WLFET(9999),
     &     DWFEI(NL,NL),CITOTFEi(5000),WLFEiT(5000)
      DIMENSION DRECHI(NL),BX(NL)
      DIMENSION DWFEIII(30,30),DWFEIV(30,30)            
      DIMENSION DCIT(nlp),DFBQ(110),RMO(MD),TEINT(MD),
     &     SIRED(nlp+NFEL,MD),deint(md),xelint(md),rint(md)
      DIMENSION BN(MD,5),THA(MD),THB(MD),ABUN(MD,20)
      DIMENSION BQH(NL),BQ(NL),BQCA(NL),BQHE(NL),BQFE(NL),BQFEi(NL)
      DIMENSION DLW(100),OLDFLU(NE1:NE2)
      DIMENSION BBOL(10,MD,NL),axsi(20),gaion(20)
      DIMENSION WLF1(300),IABU(20),ABOLD(20),KLI(300),KKMAX(6),
     &     SPL(100)
      DIMENSION WLSP(100),TREC(100)
      dimension wlsf1(5000),wlsf2(5000),dwsfei(5000),dwsfeii(5000)
      DIMENSION EJMIN(7),EJMAX(7),SION(7),SIONO(7),SIONn(7)
      dimension jres(10),wlres(10),phres(10),fej0(2,ne1:ne2),
     &     fejd(2,ne1:ne2)
      dimension escfe2(10000),tfe2(10000),febin(100),flfebin(100)
      dimension wlphot(nlp)
      common/lineem/rcgs,wlp(300),totl(300)
      common/colines/wlcor(200),cincor(200),covib1h,covib1h2,covib2h2,
     &     covib1e
      common/ionx/xion(md,14,27)
      common/ionxold/te_old(md),xion_old(md,14,27)
      common/timecheck/time,itime
      common/nperm/nlinetot
      common/opacity/TA(NE1:NE2),S(MD,NE1:NE2),copac(md,ne1:ne2)
      common/diffsh/fdshock(ne1:ne2)
      common/preion/ipre
      common/lambit/lambda
      common/heitau/tauheiold,tauhei,dtauhei
      common/fbappt/fb6300,fb5007
      common/coolfu/cfu_old,cfu,radfu_old,radfu,hfu_old,hfu 
      common/initr/initrec
      common/initox/initoi
      common/column/coltr(2,md,14,27)
      common/trq/rtr(2,md),drtr(2,md),dentr(2,md)
      common/velgrad/dr_dv  
      common/numbow/nbow
      COMMON/contspec/fc(MD,NE1:NE2)
      common/nhion/np1hs
      common/position/rcomm
      common/debug/ideb
      common/intcosio/initco_sio

      common/levels/nlevco,nlevsio

      parameter(ncosi=400)
      common/cosioem/wl_co(ncosi,ncosi),wl_sio(ncosi,ncosi),
     &     cin_co(ncosi,ncosi),cin_sio(ncosi,ncosi)
      common/citot_cosio/citot_co(ncosi,ncosi),citot_sio(ncosi,ncosi)
      common/colines/wli_co(1200),dcit_co(1200),citot_cot(1200)
      common/tauelec/tau_tot_es,tau_elec_sh,tau_elec(md)
      common/initchianti/initch
      integer initcexh
      common/initcexh/initcexh
      integer inithe2
      common/inithe/inithe2
      common/init_crh/initcrh,initcrhe1,initcrhe2
      common/init_h_hepop/inithpop,inithe2pop
      integer nz,nions,nshell,i,k
      parameter(nz=30,nions=27,nshell=10)
      integer nsh,kmax
      real*8 fr_aug,eionaug,en_aug,en_augi,eioni
      common/augerfrac/eionaug(nz,nions,nshell),en_aug(nz,nions,nshell),
     &     fr_aug(nz,nions,nshell,10),kmax(nz,nions,nshell),nsh(nz,nions),
     &     init_augfrac
      common/rec_coll/alrec(30,30),collion(30,30),ction(14,27)
      integer stop_tmin
      common/mi_te/stop_tmin
c     nion = number density of ions
c     (don't mix with nions = max # of ionization stages))      
      real*8 nion1,nion2,ni0,ni,nion(md),ne(md),nel
      real*8 kb
      real*8 hold(md),rcgsi(md)
      data kb/1.38e-16/
cxx      dimension emtr(md,ne1:ne2),emtrc(md,ne1:ne2)
      dimension cor(200),dcor(200)
      dimension photok(14,26),photolm(14,26)
      dimension xtot(md),dist_fr_sh(md)
      dimension xe(md),dens(md),rfd(md)
      dimension taus(ne1:ne2),fls(ne1:ne2),tauhi(md)
      dimension tautots(md,ne1:ne2),sources(md,ne1:ne2),ems(md,ne1:ne2),
     &     copacs(md,ne1:ne2),flaan(md,ne1:ne2),flaanc(md,ne1:ne2),
     &     emcs(md,ne1:ne2)
c      dimension flaan(md,ne1:ne2),flaanc(md,ne1:ne2)
      dimension rs(md),drs(md),densave(md),coltrs(md,14,27)
      integer nionel(14)
      data nionel/1,2,6,7,8,10,11,12,13,14,16,18,20,26/
      dimension chi(30,30),nionelz(30),vel(30000),dist(30000),
     &     tempk(30000)

      dimension vbin(1000),dfdnu(1000),xioni_old(14,27)

      DATA EJMIN/3.4d0,13.5987d0,24.5d0,54.4d0,200.d0,1.d3,1.d4/
      DATA EJMAX/13.5987d0,24.5d0,54.4d0,200.d0,1.d3,1.d4,1.d6/
      DATA IABU/1,1,1,1,0,1,1,1,0,0,1,0,0,0,6*0/
      DATA KKMAX/5,9,12,0,5,0/
      DATA PI/3.1415926d0/,ELCH/1.60219d-12/,AMU/1.660531d-24/
c     WAVELENGTHS FOR LINES IN ARRAY SPL
C     HE II 304, TWO-G, 1640.,4686, H I TWO-G, LY-C, BA-C, PA-C, BR-C
C     PF-C, FREE-FREE
C     
      DATA WLSP/304.d0,0.d0,1640.d0,4686.d0,0.d0,911.5d0,
     &     3646.0d0,8203.6d0,14584.d0,
     &     22788.d0,0.d0,89*0.d0/
      DATA SOLMA/1.989d33/
      DATA WLEV/1.239854e4/
      real*8 atw(15)
      data atw/1.,4.,12.,14.,16.,20.,22.,24.,26.,
     &     28.,32.,36.,40.,56.,58./
      real*8 cpu_start, cpu_tsec
      common/cpu_t/cpu_start,cpu_tsec      

c !! output to fort.61 = 1?
      iout61=0
      ideb = 0

      TMINI=500.
c     minimum temperature in preionization region
      temin_pre=2.5e2
      
      KKMAX(6)=NFEI
      KKMAX(4)=NFEII
      initch=1
      iemiss=1
      ifeesc=0
      initstruc=0
      initip = 1
      initlabel=0
      initcoll=0
      initprl=0
      initcr=0
      initoi=1
      inife1=0
      inife2=0
      initphd=0
      initrec=0
      initdi=0
      ikold=0
      iwl=1
      ii=0
      imultio=1
      imultica=1
      imultifei=1
      imultifeii=1
      initcop = 1
      nlinetot=186
      initco_sio = 1
      initcrh=1
      initcrhe1=1
      initcrhe2=1
      inithpop=0
      inithe2pop=0


C     ***********************************************************
C     *****
C     INITIALIZE PARAMETERS
C     *****
C     ***********************************************************

c     calculate ionization potentials
      tq=1.e6

      call ion_pot(nionelz,chi)

      OPEN(32,FILE='./ATDAT/labnum.dat',status='old')
      open(24,file='./ATDAT/line_pointers_v2.dat',status='old')
      do k=1,nlp
         read(24,*,end=99)ikl,wlqq,iline
         if(wlqq.gt.373.9.and.wlqq.lt.374.1.and.ikl.eq.503) then
            nbow = iline
         endif
         wlphot(iline)=wlqq
      enddo
 99   continue
      close(24)
      DO K=1,300
         READ(32,9537,END=9538)KLI(K),LINELAB(K)
      ENDDO
 9537 FORMAT(I5,A18)
 9538 KLINE=K-1


C     ***********************************************************
C     *****
C     INITIALIZE PARAMETERS
C     *****
C     ***********************************************************

      INI=0
C     
C     ELECTRON DENSITY EXCEPT FOR H
C     
      ZEL=ABN(2)+ABN(5)+ABN(3)+ABn(10)+ABn(8)+ABn(14)+ABn(9)
      TSN=1.E-10
      DO J=JMIN,JJ
         DO  KK=1,2
            FL0(KK,J)=0.
            FL(KK,J)=0.
         enddo
         DO iel=1,14
            do iz=1,nionel(iel)
               SI(iel,iz,J)=0.
               SK(iel,iz,J)=0.
            enddo
         enddo
         F0(J)=0.
         DO I1=1,MD
            FLUXED(J)=0.
            EM(I1,J)=0.
            TAU(I1,J)=0.
            TAUTOT(I1,J)=0.
            FD(I1,J)=0.
         enddo
      enddo

      IF(IUY.NE.1) GOTO 3862
      DO IA=1,NL
         BB(IA)=1.E-10
         BOLCA(IA)=1.
         BOLfe(IA)=1.
         do k=1,4
            bolfes(k,ia)=1.
            bolfeis(k,ia)=1.
         enddo
         BOL(IA)=1.+0.1*FLOAT(IA-1)
      enddo
 3862 BB(1)=1.E0
      NPR=NPRINT
      IPARA=IPAR
      IOIND=1
      gahe=1.d-33

c     inner radius 

      rmin = 1.0d15

C     
C     R15=PHOTOSPHERIC RADIUS IN 1E15 CM
C     
      R15=R1
      R11=1.E4*R15
C     
C     CV = VELOCITY OF LIGHT / VELOCITY AT THE POSITION OF THE  BOUNDARY
C     
      CV=3.E10/VEXP
      RMAXCGS=RMAX*1.E15	
      RMCV=RMAX*1.E15*CV
C     
C     OPTICAL DEPTH CONSTANT
C     
      C33=1.3263E-21*CV*RMAX
      VTERM=1.2856E4

      IF(ISTAT.EQ.1) then

c     a21*wl**3 * c/8 * pi * vterm         

         C33=2.2448E-26/VTERM

      endif

C     
C     COUNTER OF ITERATIONS IN RAD. TRANSFER CALC.
C     
      N=1
      NITSPH=N
      IF(INOUT.EQ.0) NITSPH=2
      KQ=1
C     
C     NUMBER OF LEVELS IN H-MODEL ATOM
C     
      NH=3
      NP1=NH+1
      C4=1.44E8
      NUM=28
      NOTS=1
      do iel=1,14
         do iz=1,nionel(iel)
            ZK(iel,iz)=0.
            ZE(iel,iz)=0.0E0
            GE(iel,iz)=0.0E0
         enddo
      enddo

C     
C     INITIAL GUESSES OF IONIC ABUNDANCIES
C     
      XION(1,1,1)=0.05
      XION(1,1,2)=0.95
      XION(1,2,1)=0.01
      XION(1,2,2)=0.05
      XION(1,2,3)=0.94
      XION(1,3,4)=1.
      XION(1,4,5)=1.
      XION(1,5,6)=1.
      XION(1,6,3)=1.
      XION(1,7,3)=1.
      XION(1,8,3)=1.
      XION(1,9,3)=1.
      XION(1,10,4)=1.
      XION(1,11,3)=1.
      XION(1,13,4)=1.
      XION(1,14,7)=1.

c     Ionized initial conditions
      do iel=1,14
         do iz=1,nionel(iel)+1
            xion(1,iel,iz)=0.
         enddo
      enddo
      XION(1,1,1)=0.001
      XION(1,1,2)=0.999
      XION(1,2,1)=0.001
      XION(1,2,2)=0.001
      XION(1,2,3)=0.998
      XION(1,3,7)=1.
      XION(1,4,8)=1.
      XION(1,5,9)=1.
      XION(1,6,10)=1.
      XION(1,7,11)=1.
      XION(1,8,12)=1.
      XION(1,9,4)=1.
      XION(1,10,6)=1.
      XION(1,11,8)=1.
      XION(1,12,9)=1.
      XION(1,13,7)=1.
      XION(1,14,9)=1.


      
C     

C     
C     ASSUME CASE B INITIALLY
C     
      OTSP(1)=0.
      OTSP(2)=0.
      OTSP(3)=1.
      DEL(1)=ABN(1)+ABN(2)+ABN(5)+ABN(3)+ABN(4)+ABn(10)+ABn(14)
      DEL(1)=0.1
      I=1
C     
C     DENSITY AT INNER BOUNDARY
C     
C     **************************************************************
C     *****
C     R(1)=  BOUNDARY RADIUS IN TERMS OF THE PHOTOSPHERIC RADIUS
C     R(I) = R(CGS)/R(PHOT)
C     DR(I) = SHELL THICKNESS IN CM
C     DRA(I) = SHELL THICKNESS RELATIVE TO THE  BOUNDARY RADIUS.
C     DEN(I) = NUMBER DENSITY IN SHELL I.
C     *****
C     **************************************************************
C     
C     ISTRU=1 IF DENSITY FROM WOOSLEY WR-MODEL, 0 IF UNIFORM OX & HE
C     
      ISTRU=1
      if(icsn87a==1) then
         istru=0
      endif
      RSUN=MTOT/SOLMA
      RM(1)=RSUN/1.000001
c     starting from center specify the inner mass
      IF(INOUT.EQ.1) RM(1)=RMINNER
      IF(IREV.EQ.1) THEN
         RM(1)=RSUN/1.000001
c!!!  
         RM(1)=1.e-6
         rm(1)=rsun
         RI15=RSHOCK
      ENDIF
cww   
c!!   RM(1)=1.e-6
      IF(IRINNER.EQ.1) THEN
         RI15=RIN
      ENDIF
      RI15=RSHOCK
      RM1CGS=RI15*1.E15
      i=2
c!!!  
      ri15=0.
c!!   assume distance to shock = 1e15 cm
      ri15=1.

      CALL STRUC(IREV,ISTRU,RM(1),DENR,RI15,RIN15,RM1CGS)
      DENR=DENR/FILLING
      IF(IRINNER.EQ.1) THEN
         RI15=RIN
      ENDIF
      icsm=0
      denr=1.e4
      write(0,*)' Number density in pre-shock gas?'
      read(5,*)denr
      denr0=denr

      if(ics87a==1) then
         densi=denr
      endif

      R(1)=RI15/R15
      RQ=R(1)
      RCGS=1.E15*R(1)*R15
      rshock_cgs=rcgs
      RMO(1)=RM(1)
      RMASS=0.
C     CALL ABUND(INI,RMASS)
      DEN(1)=DENR
      TE(1)=TIN
      IF(ICS.NE.1.and.ics87a.ne.1) THEN
C     TOTAL COLUMN DENSITY
         COTOT=0.
         COTOTN=0.
         RMI=RSUN/1.0000001
         RM1CGS=RIN*1.E15
         CALL STRUC(IREV,ISTRU,RMI,DENR,RI15,RIN15,RM1CGS)
         DENR=DENR/FILLING
         DO I1=1,1000
            ROLD=RI15
            RMI=-I1*RSUN/(1000)+RSUN+1.E-5
            RM1CGS=RI15*1.E15
            CALL STRUC(IREV,ISTRU,RMI,DENR,RI15,RIN15,RM1CGS)
            IF(ICDEN.EQ.1.or.ipress.eq.1) THEN
               DENR=DECON
               DENR=DENR/FILLING
               DEN1=DENR
            else
               DENR=DENR/FILLING
            ENDIF
            COTOTN=1.E15*DENR*(-RI15+ROLD)+COTOT
            COTOT=1.E15*DENR*(-RI15+ROLD)*AMEAN*AMU+COTOT
         enddo
      ENDIF
C     ***********************************************************
C     *****
C     PRINT OUT INPUT PARAMETERS
C     *****
C     ***********************************************************
      INI=0

      call pointers


c     read atomic data for x-ray lines

c     call atdat_mewe_gs

C     *************************************************************
C     EVALUATTION OF CROSSECTIONS FOR EACH ENERGY
C     *************************************************************
      E00=11.871
      NHY=NH
      NMA=NH
      NIONi=1
      CALL CROSSSECT
      CALL ATDATO
      CALL CROXY
      CALL ATDATFEI
      CALL CRFEI
      CALL ATDATFE
      CALL CRFEII
      CALL IONLABEL(LAB)
      NIq=1
      CALL ATDAT
      CALL CRCA
      call read_coll
      iqsx=1
 9753 format(i4,30e12.5)
 9754 format('tres ',2i4,10e12.5)
      iqsxc=1
 1534 continue
c     Ly a, C IV, C III, Mg II 2800 pointers
      wlres(1)=1216.
      wlres(2)=1550.
      wlres(3)=1909.
      wlres(4)=2800.
      do k=1,4
         emk=12400./wlres(k)
         do j=jmin,jj
            if(emk.ge.e(j).and.emk.lt.e(j+1)) then
               jres(k)=j
            endif
         enddo      
      enddo      
C     ***********************************************************
C     *****
C     INFALLING SPECTRUM
C     *****
C     ***********************************************************
C     FL=MEAN INTENSITY PER UNIT OF ENERGY (ERG/(CM**2 STER EV))
      IF(IPAR.EQ.1) GEOM=4.
      IF(IPAR.NE.1) GEOM=2.


      DO J=JMIN,JJ
C     
C     TOTAL MEAN INTENSITY
C     
         FL(1,J)=FMEAN(1,E1(J))
         FL(2,J)=FL(1,J)
         FLUXED(J)=FL(1,J)
         OLDFLU(J)=FLUXED(J)
         F0(J)=FL(1,J)
      ENDDO

      
C     CALCULATE FLUX
      iini=2
      ik = 1
      DO J=JMIN,JJ
         CALL RADTRANPULS(1,J,TAU(iini,J),S)
      ENDDO
      DO IA=1,NL
         DO K=1,6
            COLHVTT(K,IA)=2.
            IF(IA.EQ.1) COLHVTT(K,IA)=1.D33
         ENDDO
      ENDDO
C     
C     START OF ITERATION LOOP FOR THE DIFFUSE FIELD
      ICONTI=0
 555  CONTINUE
      INIT=1
      INITCA=1
      INITFE=1
      INITMG=1

      DXI=DR(2)
      REHE1=0.
      REN1=0.
      REC1=0.
      REHE2=0.
      rehe3=0.
      rehe4=0.
      REC31=0.
      REN41=0.
      DO  J=JMIN,JJ
         TAU(1,J)=0.
         FPRIM=FMEAN(1,E1(J))
         FL(2,J)=FPRIM+FD(1,J)
         FL(1,J)=FL(2,J)
      enddo
      TEL=0.
      COLHYD=0.
      DO iel=1,14
         do iz=1,nionel(iel)
C     ASSUME INFINITE OPTICAL DEPTHS IN THE DIRECTION OUTWARDS
C     FROM THE SOURCE
            COLTOTT(iel,iz)=1.D100
            COLTOT(iel,iz)=0.
            COLD(iel,iz)=0.
            COLDN(iel,iz)=0.
         ENDDO
      ENDDO

      do iel=3,14
         do ion=1,nionel(iel)
            do line=1,401
               do iqq=1,md
                  tot_intens(iqq,iel,ion,line)=0.
               enddo
            enddo
         enddo
      enddo

      do k=1,200
         cor(k)=0.
      enddo
      covibv1=0.
      covibv2=0.
      DO IU=2,NFEI
         DO IL=1,IU-1
cq            WFEI(IU,IL)=0.
         ENDDO
      ENDDO
      DO IU=2,NFEII
         DO IL=1,IU-1
cq            WFEII(IU,IL)=0.
         ENDDO
      ENDDO
      NFEIII=30
      DO IU=2,NFEIII
         DO IL=1,IU-1
cq            WFEIII(IU,IL)=0.
         ENDDO
      ENDDO
      NFEIV=22
      DO IU=2,NFEIV
         DO IL=1,IU-1
cq            WFEIV(IU,IL)=0.
         ENDDO
      ENDDO
      RECHETOT=0.
      DO K=1,5
cq         RECHITOT(K)=0.
      ENDDO
      TWOPHTOT=0.
      TWOPHHEIT=0.
      FREEFREE=0.
      HEATTOT=0.
      COOLTOT=0.
      GAMTOT=0.
      DO I1=1,NL
c     CONTOT(I1)=0.
         CNRTOT(I1)=0.
         PHOT(I1)=0.
         REC(I1)=0.
         DO K=1,6
            COLH(k,I1)=0.
            COLHVT(k,I1)=1.
         ENDDO
         DO J=1,NL
            ESC(I1,J)=0.
            TTOT(I1,J)=1.E20
            TOP(I1,J)=0.
            EMTOT(I1,J)=0.
            EMIS(I1,J)=0.
         ENDDO
      ENDDO
      DO K=1,30
         WREC(K)=0.
      enddo
      DO K=1,40
         WEH(K)=0.
      ENDDO
      TAUB=0.
      citt1=0.
      citt2=0.
      citt3=0.
      citt4=0.
      citt5=0.
      citt6=0.
      citt7=0.
      citt8=0.
      citt9=0.
      citt10=0.

C     ***********************************************************
C     *****
C     CALCULATE PHOTOIONIZATION RATES FOR SLAB 1
C     *****
C     ***********************************************************
      DO J=JMIN,JJ
         TAU(1,J)=0.0E0
      ENDDO
      RMASS=1.E-6
c     RMASS=RM(1)
      COLUMN=COTOT
c     from center set column density = 0
      IF(INOUT.EQ.1) COLUMN=0.
      IF(ICS.EQ.1) COLUMN=0.
      COLUMNN=COTOTN
      IF(INOUT.EQ.1) COLUMNN=0.
      IF(ICS.EQ.1) COLUMNN=0.
      TAUELE=0.
      M=2
      MQW=0
      MOXY=0
      IFE=0
      ISI=0
      IOMG=0
      IOC=0
      IF(IUY.EQ.1) TOLD=TIN
      IF(IUY.EQ.1) TOLDO=TIN
      IF(IUY.EQ.1) TOLDCA=TIN
      IF(IUY.EQ.1) TOLDFE=TIN
      IF(IUY.EQ.1) TOLDFEi=TIN
      IF(IUY.EQ.1) TOLDMG=TIN
      IF(IUY.EQ.1) TOLDQ=TIN
      IF(IUY.EQ.1) TOLDHQ=TIN
      IF(IUY.EQ.1) TOLDHEQ=TIN
      IF(IUY.EQ.1) TOLDOQ=TIN
      IF(IUY.EQ.1) TOLDCAQ=TIN
      IF(IUY.EQ.1) TOLDFEQ=TIN
      IF(IUY.EQ.1) TOLDFEiQ=TIN
      IF(IUY.EQ.1) TOLDMGQ=TIN
      IF(IUY.EQ.1) then
         do k=1,4
            toldfeis(k)=TIN
            toldfes(k)=TIN
         enddo
      endif
      tot_num_rec=0.

C     ***********************************************************
C     *****
C     MAIN LOOP IN DEPTH FOR EACH ITERATION
C     *****
C     ***********************************************************
      lambda = 1

      ipre = 0

      do j=jmin,jj

c     mean int. at shock
         
         fdshock(j) = 0.

      enddo

      write(6,9277)(abn(ia),ia=3,14)
 9277 format(' Abund. ',1pe12.3,14e12.3)

      fb6300=0.
      fb5007=0.

      itemp = 1

      ipp = 0

      ishock = 1

 777  continue

c$$$      if(ipre==1) then
c$$$         xtot(1)=-1.e8
c$$$         xtot(2)=-1.e8
c$$$         istat=1
c$$$c!!!  231030
c$$$         dt=-1.e2
c$$$         write(6,*)' Pre-ionization put xtot slightly neg.',ns,xtot(ns),dt
c$$$      else
c$$$         istat=0
c$$$      endif
      l_list= 25 
      
      IF(ISTAT.EQ.1) then
C     OPTICAL DEPTH CONSTANT
c     a21*wl**3 * c/8 * pi * vterm
         VTERM=1.2856E4
         C33=2.2448E-26/VTERM
      else   
         C33=1.3263E-21*CV*RMAX
      endif

      
c$$$      write(6,*)' New iteration: lambda= ',lambda,' k, ipre= ',ipre,k
c$$$      write(75,*)' New iteration: lambda= ',lambda

      if(icsm==1.and.ipre==1) then
         r(k)=rshock_cgs
         iss=1
      elseif(icsm==1.and.ipre==0) then
         r(k)=rshock_cgs
      endif
      r(k)=rshock_cgs

      iteit = 0

      initisp=1

      isp=1

      ishkp=0

      dt=1.
      
      ns = 1

      nc = 0
      
      do kk=1,nlp
cq         CITOT(kk)=0.
cq         te_av(kk)=0.
      enddo

      do iel=1,14
         do ion=1,nionel(iel)
            do line=1,401
               do iqq=1,md
                  tot_intens(iqq,iel,ion,line)=0.
               enddo
            enddo
         enddo
      enddo


      if(lambda.ge.2) then
         DO J=JMIN,JJ
            DO I1=1,MD
               EM(I1,J)=0.
               TAU(I1,J)=0.
               TAUTOT(I1,J)=0.
            enddo
         enddo
      endif

c     set column density = 0 for first step

c     for post-shock ok for first iteration

c     for pre-shock ok for first iteration since integrating from 
c     shock and back

c     for second iteration ok since coltot etc inverted after iteration s
c     so that column densit = 0 at shock and increasing in each 
c     direction 

      DO iel=1,14
         do iz=1,nionel(iel)
            if(lambda.le.1) then
C     ASSUME INFINITE OPTICAL DEPTHS IN THE DIRECTION OUTWARDS
C     FROM THE SOURCE
               COLTOTT(iel,iz)=1.D100
            endif
            COLTOT(iel,iz)=0.
            COLTOTtemp(iel,iz)=0.
            COLD(iel,iz)=0.
            COLDN(iel,iz)=0.
         ENDDO
      ENDDO

      if(initip.eq.1) then

         initip = 0

c     ion density in front of shock

         dnum0=denr0

         write(0,*)' V_s (km/s) '

         read(5,*)vs

         write(0,*)' ionization time (n0 x t in cgs)'
         read(5,*)t_ioniz

c     vbin in cm/s     

         ivmax = 801

         deltav = vs*1.e5/(4.*real(ivmax-1))

         uold = vs*1.e5/4.

         do iv = 1,ivmax
            vbin(iv) =  deltav*real(iv-1)
            dfdnu(iv) = 0.
         enddo

c     magnetic field in gauss
         b0=1.e-7*sqrt(dnum0)
c!!!  .
c         b0=0.

      endif

      tlim_max = 1.e6
      tlim_min = 1.e4

      xh=abn(1)
      xhe = abn(2)
      xo = abn(5)

      xe0 = 1.*xh + 2.*xhe + 6.*xo
      if(lambda.eq.2) then
         
         xe0 = xeold

      elseif(lambda.gt.2) then
         
         xe0 = xeold

      endif

      amean=0.
      do iel=3,14
         amean=amean+abn(iel)*atw(iel)
      enddo

      nel=0.
      do iel=3,14
         do ion=1,nionel(iel)
            nel=nel + abn(iel)*xion(1,iel,ion)*real(ion-1)
         enddo                    
 9182    format('iel,ion,abn,xion,nel ',2i4,1pe12.3,10e12.3)
      enddo

      rmus=amean/(1.+nel)

      te0=2.269e9*rmus*(vs/1.e4)**2

c     vs in km/s

      u0=vs*1.e5
      
      d00=amu*dnum0*amean
      
      te0=2.269e9*rmus*(vs/1.e4)**2

c     first denm = pre-shock density
      denm = d00

      ni0=denm/(amu*amean)

c     assume temperature in front of the shock is constant = t00

      te00 = 1.0e4

      write(6,9246)vs,te0,dnum0,b0,t_ioniz
 9246 format('V_s = ',1pe12.3,' T_s = ',e12.3,' n_0 = ',e12.3,
     &     ' B_0 = ',e12.3,' t_ion = ',e12.3)

c!!!  skip next 4 lines

      te00 = 1.e4

      if(ipre.eq.1.and.lambda==1) then
         te00 = 1.e4
      endif



c!!!  modif 030807

c     te00 = 0.7e3


c     enthalpy in front of the shock
c     

c!!!  index?
      ik = 1

c     if(ipre==0) then
      eion = 0.

      do iel=3,14
         do iz=1,nionel(iel)
            eion = eion + abn(iel)*elch*chi(iel,iz)*
     &           xion(ik,iel,iz+1)
         enddo
      enddo

      eion=0
      
      h0 = b0**2/(4.*pi*ni0) + amean*amu*u0**2/2. + 
     &     2.5*(1.+xe0)*kb*te00 + eion
      
      call ion_energy(ik,chi,nionel,eion)
      
      h0prim = h0-eion
      
c     mass density flux
      denfl = u0*d00
      
c     ram pressure
      pram0 = d00 * u0**2
      eb0 = b0**2/(8.*pi)
      ptherm0 = dnum0*(1.+xe0)*kb*te00
      p0 = pram0 + ptherm0 + eb0


      write(6,9189)d00,u0,b0,dnum0,xe0,te00,pram0,ptherm0,eb0,p0
 9189 format(' pram ',1pe12.4,10e12.4)
      xc = 4.1


      call compress(p0,eb0,ni0,u0,denfl,h0prim,xc)

      time=0.

      dt=1.e3
      

c     initial guesses for compression


      if(ipre.eq.0) then

c     first iteration only post shock

         xc=4.1

         dt = 1.e1


         do j=jmin,jj
            
            fl(2,j)=0.
            
         enddo


      else

c     second and higher pre-shock

         if(xtot(ns).le.0.) then
            xc=1.01            
         else
            xc = 4.01
         endif

         if(lambda.eq.1) then
            dt = -1.e3
            if(ipre==1) then
               dt = -1.e2
            endif
         else
            dt = 1.e2
         endif

         do j=jmin,jj
            
            fl(2,j)=fdshock(j)
            
         enddo

         write(6,*) ' preionization '

      endif

      ni = xc*ni0

c     index?
      write(6,*)'icsm,ipre ',icsm,ipre
      if(icsm==0) then

         ns = 1
      
         if(lambda.eq.1.and.ipre==0) then

            do ij=1,md

               xtot(ij) = 1.e-10

            enddo

         else

            do ij=1,md

               xtot(ij) = rstart

            enddo

         endif

c$$$      elseif(icsm==1.and.ipre==1) then
c$$$
c$$$         xtot(2)=r(k-1)
c$$$
c$$$         xtot(2) = rshock_cgs
c$$$
c$$$         write(6,*)' xtot(2),r(k-1),r(k) ',xtot(2),r(k-1),r(k)
         
      endif

      tdiff = 0.

      tshell = 2.e6

      iint=0

      if(lambda.eq.1.and.ipre==0) then
c index?         
         k = 1
         r(1) = 0.
c!!! next statement questionable 211020         
      elseif(lambda.ge.2) then
         k = 0
c     f!!!!! modif 020807
         k = 1
      endif
c??!!

      totcool = 0.
      totcoolt = 0.
      totcoolet = 0.
      totheat = 0.
      totemt=0.
      totemit=0.

      totha=0.
      totha2=0.
      ispec=1

      itimediff=0

      do itime=1,100000


c     follow the shell in time withh constant index ns
c     when optical depth of shell is larger than a ceratin limit add a new shell and save
c     the old with index k
         
c     loop:

c     1. r(k) = xtot(i)

c     2. dr(k) = r(k) - r(k-1)

c     3. emiss -> em(k)

c     4. opac -> ta(k)*den(i)*dr(k) = dtau(k), tautot(k)

c     5. coolf -> rad(Te(i)) -> h

c     6. tempcalc -> xc -> den(i), T(i), u(i)

c     7. xtot(i) = u(i)*dt + xtot(i)

c     r(k) = xtot(ns)

c     goto 1 

c     Note: 

c     cooling fcn calc. with old temp. Should really have been between the 
c     time steps (in the mid point)

c     radiation transport calc. should use emissivities etc at midpoints for calc. 
c     of mean intensity at the grid points at r(k).


c     First iteration: r(1) = 0 at shock
c     r(2) = 

c     r(1)=0            r(2)              r(3)                              r(imax)

c     dr(2)            dr(3)                              dr(imax)                note! dr(1)=0

c     after inversion

c     r(nshock)=0   r(nshock-1)       r(nshock-2)             r(2)            r(1)
         
c     dr(nshock)    dr(nshock-1)                             dr(2)

c     preshock

c     r(nshock+2)     r(nshock+1)       r(nshock)=0           
         
c     dr(nshock+2)     dr(nshock+1)                    


c     Second and higher

c     r(1) = rstart (at upstream boundary)


         if(itime.eq.10) then
            iwl=1
         endif

         if(itime.eq.1) then
            newshell = 1
            dt2 = dt
            ns=1
         endif

         old_time = time

         te_oldi=te(i)

         time=old_time + dt

         if(itime.le.20) then

            dr_dv = time

         endif

         tday = time/8.64e4

c!!   not always true, but is probably correct for shocks

         tdays = tday

c     this is used only for the soboloev length

         drscgs = xtot(ns)

c!!!  use sobolev everywhere!!!!

         drscgs = 1.e33

c     rcgs0 = 1.e18
         
c     r(1) = rcgs0/(1.e15*r15)

         if(lambda.eq.1.and.ipre==0) then

            r(1)=0.

         else

            if(icsm==0.and.ipre.ne.1) then
               r(1) = rstart
            endif
         endif

c         write(6,*)' step 1 start r',k,lambda,r(k),dt

         if(lambda.ge.2) then

            rtr(2,1) = r(1)

         endif

c     rshock=rcgs0/(1.e15*r15)

         eflux_old = eflux 
         etot_old = etot 

c     if(xtot(ns).ge.rfshock.and.lambda.ge.2) then
c     xc=4.1


         if(itime.ge.2) then
            do is=2,ns
               te_old(is)=te(is)
               do iel=3,14
                  do j=1,nionel(iel)+1
                     xion_old(is,iel,j)=xion(is,iel,j)
                  enddo
               enddo
            enddo
         endif

c     check if new shell should be added

c     tdiff = tdiff + dt
c     if(tdiff.gt.tshell) then
c     tdiff = 0.
c     newshell = 1
c     endif
c!!!! ????  shoild (itime==1.and.ipre==0) be added?
         if(itime.eq.1) then
c     
c     add new shell at shock (first post-shock shell)
c     

            newshell = 0
c     note that ns is always 2 for the shell being calculated!
            ns = ns+1
c     rcgsi(ns)=rcgs0
c     r(ns) = rcgsi(i)/(1.e15*r15)
c     dr(ns) = rcgsi(ns-1)-rcgsi(ns)
            xe(ns)=xe0
            te(ns)=1.e5

c$$$            if(ipre==0) then
               hold(ns) = h0

               write(6,*)' entalpi hold a ',h0

c               write(6,*)' step 1 r(k) init ',k,ns,r(k)

               h = hold(ns)

               write(6,*)' entalpi a ',h

               call ion_energy(ns,chi,nionel,eion)

               hprim = h - eion

               write(6,*)' iss,h,eion,hprim ',iss,h,eion,hprim

               call compress(p0,eb0,ni0,u0,denfl,hprim,xc)

               write(6,*)'xxx itime,i,ik,ns,k,te(ik),den(ik),del(ik) b ',itime,i,ik,k,ns,te(ik),den(ik),del(ik)

               dencgs = xc*d00

               u = denfl/dencgs
               
               p = p0-dencgs*u**2-eb0*xc**2

               nion(ns)=dencgs/(amu*amean)
               ne(ns)=nion(ns)*xe(ns)
               den(ns)=nion(ns)

               dens(ns)=den(ns)
               den1 = den(ns)
               te(ns)=p/(1.38e-16*(ne(ns)+nion(ns)))

c$$$            elseif(ipre==1.and.icsm==1) then
c$$$c does not apply to 87a!               
c$$$c     den(ns)=dmdt/(4.*pi*r(k)**2*uwind*amean*amu)
c$$$c     dencgs=dmdt/(4.*pi*r(k)**2*uwind)
c$$$c     if(r(k) > r_precursor ) then
c$$$c     den(ns)=den(ns)/prec_boost
c$$$c     dencgs=dencgs/prec_boost
c$$$c     endif
c$$$               x=(r(k)-r_precursor)/r_precursor
c$$$               dblog=prec_boost*(atan(x/ramp)/pi+1/2.)
c$$$               denslog=densi_rs_log-2.*log10(r(k)/r_precursor)-dblog
c$$$               den(ns)=10.**denslog
c$$$               dencgs=den(ns)*amean*amu
c$$$               dens(ns)=den(ns)
c$$$               nion(ns)=dencgs/(amu*amean)
c$$$               ne(ns)=nion(ns)*xe(ns)
c$$$               write(6,9261)k,dmdt,r(k),amean,
c$$$     &              den(ns),ne(ns)
c$$$ 9261          format(' new wind shell a ',i6,1pe12.3,10e12.3)
c$$$
c$$$            elseif(ipre==1.and.icsm==0) then
c$$$c does not apply to 87a!               
c$$$c     den(ns)=dmdt/(4.*pi*r(k)**2*uwind*amean*amu)
c$$$c     dencgs=dmdt/(4.*pi*r(k)**2*uwind)
c$$$c     if(r(k) > r_precursor ) then
c$$$c     den(ns)=den(ns)/prec_boost
c$$$c     dencgs=dencgs/prec_boost
c$$$c     endif
c$$$
c$$$               den(ns)=dnum0
c$$$               dencgs=den(ns)*amean*amu
c$$$               dens(ns)=den(ns)
c$$$               nion(ns)=dencgs/(amu*amean)
c$$$               ne(ns)=nion(ns)*xe(ns)
c$$$               write(6,9563)k,dmdt,r(k),amean,
c$$$     &              den(ns),ne(ns)
c$$$ 9563          format(' new preion shell a ',i6,1pe12.3,10e12.3)
c$$$            endif

            des2 = nion(ns)
            des1 = des2

            tes2 = te(ns)
            tes1 = tes2

            xes2 = xe(ns)
            xes1 = xes2

            write(6,*)amean,nion(ns),ne(ns),te(ns)
            if(te(ns).le.0.) then
               te(ns) = te00
            endif

            tempold = te(ns)

            write(6,*)' te den ',ns,te(ns),den(ns),dens(ns)

            do iel=3,14

               sumx=0.

c     if(lambda.eq.1) then
c     if(lambda.ge.1) then
               if(ipre.le.1.and.lambda.eq.1) then

                  do j=1,nionel(iel)+1
                     xion(ns,iel,j)=0.
                     if(j.eq.3) then
                        xion(ns,iel,j)=1.0
                     elseif(j.eq.1) then
                        xion(ns,iel,j)=1.e-5
                     else 
                        xion(ns,iel,j)=1.e-10
                     endif
c     sumx=sumx+xion(ns,iel,j)
                  enddo
                  xion(ns,1,2)=1.0
                  xion(ns,1,1)=1.e-3
                  xion(ns,2,1)=1.e-3
                  xion(ns,2,2)=1.e-2
                  xion(ns,2,3)=1.0

               elseif(ipre.eq.-1) then

c     in pre-shock gas start with H II, He II
c     C IV, N V, O VI, Ne V, 
c     Na II, Mg II, Si II, S II, Ca II, Fe II
                  do j=1,nionel(iel)+1
                     xion(ns,iel,j)=1.e-3
                  enddo

                  xion(ns,1,2)=1.
                  xion(ns,2,3)=1.0
                  xion(ns,3,5)=1.
                  xion(ns,4,5)=1.
                  xion(ns,5,6)=1.
                  xion(ns,6,6)=1.
                  xion(ns,7,5)=1.
                  xion(ns,8,5)=1.
                  xion(ns,9,5)=1.
                  xion(ns,10,5)=1.
                  xion(ns,11,5)=1.
                  xion(ns,12,5)=1.
                  xion(ns,13,5)=1.
                  xion(ns,14,6)=1.

                  do j=1,nionel(iel)+1
                     sumx = sumx+xion(ns,iel,j)
                  enddo
                  write(6,*)' sumx pre ion ',sumx,xion(ns,iel,3)
               elseif(ipre.eq.2) then

c     in pre-shock gas start with neutral medium with exception to C II,
c     Na II, Mg II, Si II, S II, Ca II, Fe II
                  do j=1,nionel(iel)+1
                     xion(ns,iel,j)=0.
                     if(j.eq.1) then
                        xion(ns,iel,j)=1.0
                     elseif(j.eq.2) then
                        xion(ns,iel,j)=1.e-1
                     else 
                        xion(ns,iel,j)=1.e-10
                     endif
                     sumx=sumx+xion(ns,iel,j)
                  enddo
                  xion(ns,1,2)=1.e-4
                  xion(ns,2,1)=1.e0
                  xion(ns,2,2)=1.e-3
                  xion(ns,2,3)=1.0e-5
                  xion(ns,3,1)=1.e-4
                  xion(ns,3,2)=1.
                  xion(ns,7,1)=1.e-4
                  xion(ns,7,2)=1.
                  xion(ns,8,1)=1.e-4
                  xion(ns,8,2)=1.
                  xion(ns,9,1)=1.e-4
                  xion(ns,9,2)=1.
                  xion(ns,10,1)=1.e-4
                  xion(ns,10,2)=1.
                  xion(ns,11,1)=1.e-4
                  xion(ns,11,2)=1.
                  xion(ns,13,1)=1.e-4
                  xion(ns,13,2)=1.
                  xion(ns,14,1)=1.e-4
                  xion(ns,14,2)=1.
               elseif(ipre.eq.1.and.lambda.ge.2) then


c     in pre-shock gas start with neutral medium with exception to C II,
c     Na II, Mg II, Si II, S II, Ca II, Fe II
                  do j=1,nionel(iel)+1
                     xion(ns,iel,j)=0.
                     if(j.eq.1) then
                        xion(ns,iel,j)=1.0
                     elseif(j.eq.2) then
                        xion(ns,iel,j)=1.e-1
                     else 
                        xion(ns,iel,j)=1.e-10
                     endif
                     sumx=sumx+xion(ns,iel,j)
                  enddo
                  xion(ns,1,2)=1.e-8
                  xion(ns,2,1)=1.e0
                  xion(ns,2,2)=1.e-8
                  xion(ns,2,3)=1.0e-6
                  xion(ns,3,1)=1.e-8
                  xion(ns,3,2)=1.
                  xion(ns,7,1)=1.e-8
                  xion(ns,7,2)=1.
                  xion(ns,8,1)=1.e-8
                  xion(ns,8,2)=1.
                  xion(ns,9,1)=1.e-8
                  xion(ns,9,2)=1.
                  xion(ns,10,1)=1.e-8
                  xion(ns,10,2)=1.
                  xion(ns,11,1)=1.e-8
                  xion(ns,11,2)=1.
                  xion(ns,13,1)=1.e-8
                  xion(ns,13,2)=1.
                  xion(ns,14,1)=1.e-8
                  xion(ns,14,17)=1.


                  write(6,9165)iel,(xion(ns,iel,j),j=1,nionel(iel)+1)
                  write(6,9165)iel,(xion_old(ns,iel,j),
     &                 j=1,nionel(iel)+1)
 9165             format(' init xion ',i5,1pe12.3,30e12.3)
               endif
               
            enddo

            do iel=3,14
               sumx=0.
               do j=1,nionel(iel)+1
                  sumx = sumx+xion(ns,iel,j)
               enddo

               do j=1,nionel(iel)+1
                  xion(ns,iel,j) = xion(ns,iel,j)/sumx
               enddo

               do j=1,nionel(iel)+1
                  xion_old(ns,iel,j)=xion(ns,iel,j)
               enddo                  

            enddo

            call enum_cf(ns)
            
            call elecdens(ns,xe(ns),zelec)

            ne(ns)=nion(ns)*xe(ns)

            if(lambda.eq.1) then
               do j=jmin,jj
                  fd(ns,j)=0.
               enddo
            endif

            ik = k

            kk = k

            if(itime.eq.1.and.lambda.ge.2) then
               ik = 2
               kk = 1
            endif

            iqq=1
            call opac_source(iqq,0,dr(kk),den(ns),te(ns),xe(ns))
         endif
c     end of first shell         
         ipri=0
         iint=iint+1
         if(iint.eq.100) then
            ipri=1
            iint=0
         endif
         if(xtot(ns).lt.rfshock.and.lambda.ge.2222) then
c     if(xtot(ns).lt.rfshock.and.lambda.lt.2) then
c     New PRE-SHOCK
c     for pre-shock shells use the old shell structure

            if(xtot(ns).gt.r(k)) then
               
               write(6,*)' new pre-shock shell ',k,ik,xtot(ns),r(k),
     &              r(k+1),dr(k),den(k)
 9374          format(' new pre-shock shell ',2i5,1pe13.5,10e12.3)
               
               k = k+1
               ispec = 1
               write(6,*)'k, r(k) c ',k, r(k)
            endif

         else
c     new POST_SHOCK shell?
c     
c     for post-shock shells invoke a new shell structure (k > 1)
c!!   use this also for pre shock after first iteration

            temax = 0.05
            temax = 0.1
            tediff = abs(te(ns) - tempold)/tempold
            ispec = 0

c     change in He optical depth

            if(itime.gt.3) then
               taudiff=log10(tauhei/tauheiold)
            else
               taudiff=0.
c     calculate spectrum for itime < 3               
               ispec=1
            endif

            taudmax=0.05
            taudmax=0.025
            taudmax=0.015
            
c     if(tdiff.gt.tdmax) then
            itimediff=itimediff+1

c     check if new spectrum should be calculated for large enough
c     change in the  optical depth

c     for lambda = 1 use all shells calculated. for lambda >= 2 only those 
c     with delta tau > taudmax

            if(lambda.ge.3.and.xtot(ns).gt.0.) then

c     new shell if thickness larger than 0.02 x total 

               deltark = xtot(ns)-r(k)

               ratio_r = abs(deltark/rtotsh)

               write(6,9163)xtot(ns),dxtot,rtr(2,imax-1),ratio_r
 9163          format('ratio_r',1pe15.6,10e12.3)
               if(ratio_r.gt.0.02) then
                  write(6,*)' ratio_r large ',ratio_r
               endif
            elseif(lambda.ge.2.and.xtot(ns).lt.0.) then

c     new shell if thickness larger than 0.02 x total 

               deltark = xtot(ns)-r(k)

               ratio_r = abs(deltark/r(1))

               write(6,9263)xtot(ns),dxtot,r(1),ratio_r
 9263          format('xtot(ns),dxtot,r(1),ratio_r ',1pe12.3,10e12.3)
               if(ratio_r.gt.0.02) then
                  write(6,*)' ratio_r large ',ratio_r
               endif

            else
               
               ratio_r = 0.

            endif
c     check if optical depth change > taudmax or other conditions for a new shell!

            if(taudiff.gt.taudmax.or.itimediff.gt.20.or.tediff.gt.temax.
     &           or.(ipre.eq.1.and.lambda.eq.1).or.ratio_r.gt.0.02.
     &           or.iss.eq.1) then
c$$$               write(6,*)'taudiff/taudmax ',taudiff/taudmax            
c$$$               write(6,*)'itimediff ',itimediff
c$$$               write(6,*)'tediff,temax ',tediff,temax
c$$$               write(6,*)'ratio_r ',ratio_r
               
               ITIMEDIFF=0
               TEMPOLD = TE(NS)
C               WRITE(6,*)' TE(NS) Q1 ',NS,TE(NS)

               IF(IPRE.EQ.1.AND.ITIME.EQ.1.AND.LAMBDA.EQ.1) THEN
C     SHOCK!!
                  K = NSHOCK
                  WRITE(6,*)' K = NSHOCK ',K

               ENDIF
C     SAVE OLD RADIUS 

C     INCREASE NUMBER OF SHELLS  BY ONE SO THAT K+1 IS THE SHELL CALCULATED
               K = K+1
c               WRITE(6,*)'NEW SHELL AT LINE ',2252,K,ITIME
c               WRITE(6,*)'NEW K ',K,ITIME
c               WRITE(47,*)'NEW K ',K,ITIME
C     IK=K
c               WRITE(6,*)'XTOT,R,TE AFT SHOCK AND BEF PRE-SHOCK ',NSHOCK
               IF(K >= L_LIST) THEN
                  WRITE(6,*)'STRONGEST LINES FOR K= ',K, IK
                  CALL SORT(IK,100,0)
                  L_LIST=L_LIST+5
               ENDIF
               DO II=1,NSHOCK
                  WRITE(6,9981)II,XTOT(II),R(II),TE(II)
 9981             FORMAT(I5,1PE20.12,E20.12,E12.3)
               ENDDO
               
               IF(K.GE.MD) THEN
                  
                  WRITE(6,*)' STOP MD ',K,MD

C     STOP

               ENDIF
C     NOTE THAT K HAS INCREASED BY ONE. FIRST SHELL FOR PREIONIZED IS K=3!

C     SAVE RADII, DENSITIES AND IONIZATION FRACTIONS AT K

C     STEP 1 IN LOOP

               R(K) = XTOT(NS)
               KRADIUS=K

               RCOMM = XTOT(NS)

C     DR IS NOT DEFINED FOR THE FIRST SHELL FOR THE SECOND AND HIGHER ITERATION

C     STEP 2 IN LOOP
               IF(LAMBDA.GE.2.AND.K.EQ.1) THEN
                  DR(K) = U*DT
               ELSE
                  DR(K) = R(K)- R(K-1)

               ENDIF

c               WRITE(6,*)' STEP 2 DR(K) ',K,DR(K)
C     SETS DENSITY TO THE SAME AS AT NS=2 (NS DOES NOT CHANGE)
               DEN(K) = DENS(NS)

C     FOR PRE-SHOCK  USE A DENSITY DETERMINED BY WIND  IF ICSM=1

               RTR(2,K) = R(K)
               DRTR(2,K) = DR(K)
               DENTR(2,K) = DEN(K)

               IF(K.EQ.2) THEN
                  DRTR(2,1) = RTR(2,2)-RTR(2,1)
                  DENTR(2,1) = DEN(K)
               ENDIF

               DO IEL=3,14
                  DO IZ=1,NIONEL(IEL)+1
                     COLTR(2,K,IEL,IZ) = COLTOTTEMP(IEL,IZ)
                  ENDDO
               ENDDO

               IF(LAMBDA.GE.2.AND.K.LE.3) THEN
                  XEOLD = XE(NS)
               ENDIF


C     SAVE XION FOR SHELL K, SINCE NS=2 AND A NEW XION(NS WILL BE
C     CALCULATED. MEANS THAT XION(2, IS ALWAYS UPDATED
C     AND NOT FOR K=2
               DO IEL=3,14
                  DO J=1,NIONEL(IEL)+1
                     XION(K,IEL,J)=XION(NS,IEL,J)
                  ENDDO
               ENDDO

               ISPEC = 1
C               WRITE(6,*)' NEW SHELL ',K,TIME,TE(NS),DEN(K),R(K),DR(K)
               IF(IDEB.EQ.1) WRITE(0,*)' NEW SHELL ',K,TIME,TE(NS),
     &              DEN(K),R(K),DR(K)

C     TEST!!! 211102
               DEN(I)=DEN(K)
               DENS(I)=DEN(K)


            ENDIF

         ENDIF


C     ONLY FOLLOW ONE SHELL  => NS = 2

C     DO I=2,NS

C     SKIP LOOP IN I. ONLY USED FOR RAD TRANSFER AT EACH TIME STEP. USE IK INSTEAD = K         
         I=2

         DENEL=XE(I)*DENS(I)
         
C     MOD 21
         IK=I
         
         DEL(IK)=XE(I)
         DEN(IK) = DENS(I)
         
         IF(ISPEC.EQ.1) THEN

C     WHEN SPECTRUM IS CALCULATED CALCULATE EMISSIVITY AND OPACITIES


            IK = K
C     STOP HERE! IS NEXT OK??
C!!!! CHANGE ABOVE TO OK?
C     IK=I
            IF(ITIME.EQ.1.AND.(IPRE.EQ.0.OR.LAMBDA.GE.2)) THEN
               IK = 2
            ENDIF

            IF(ITIME.EQ.1.AND.IPRE.EQ.1) THEN
               IK = 2
            ENDIF

            DEL(IK)=XE(I)
            DEN(IK) = DENS(I)

            IF(IK.EQ.2.AND.DEN(IK).EQ.0.) THEN
               WRITE(6,*)'I,DENS(I) ',IK,DENS(IK)
               DEN(IK) = DEN(1)
               DENS(I) = DEN(IK)
            ENDIF

C     STEP 3 IN LOOP
            
C     CALCULATE THE EMISSION FOR A SHELL WITH INDEX IK   
            CALL EMISS(TE(I),XE(I),DENS(I))

            CALL EMISS_CONT(TE(I),XE(I),DENS(I),COOLI)

            JX = 3

C     STEP 4 IN LOOP

            IQQ=2
            CALL OPAC_SOURCE(IQQ,1,DR(K),DEN(I),TE(I),XE(I))

C     SAVE OPACITIES

            DO J=JMIN,JJ
               COPAC(K,J) = TA(J)
            ENDDO


C     ONLY SOLVE RADIATION TRANSPORT IN THE POST SHOCK DIRECTION FIRST TIME

            IF(K.GE.3.AND.LAMBDA.EQ.1) THEN
               IMODE = 1
               EDDINGFLUX_OLD = EDDINGFLUX 
               ICONT = 0
               CALL DIFFUSE_PL(IMODE,ICONT,K,ISHOCK,TE(I),EDDINGFLUX) 

               ICONT = 1
               CALL DIFFUSE_PL(IMODE,ICONT,K,ISHOCK,TE(I),
     &              EDDINGFLUXC) 

               EBAND = 0.
               ENTOT = 0.
               EMIN=0.5E3
               EMAX = 2.E3

               DO J=JMIN,JJ

                  IF(E1(J).GT.EMIN.AND.E1(J).LE.EMAX) THEN

                     EBAND = EBAND + FD(K,J)*(E(J+1)-E(J))

                  ENDIF
                  
                  ENTOT = ENTOT + FD(K,J)*(E(J+1)-E(J))

               ENDDO

               DEDD=EDDINGFLUX_OLD-EDDINGFLUX    
               DXTOT=XTOT(IK)-XTOT_OLD
               DEDD_DX=DEDD/DXTOT
               XTOT_OLD=XTOT(I)

            ENDIF

            IK=I
            
            IOUT=6

            IF(ITIME.LE.3) THEN
               TAUHEI=TA(102)*ABS(DELTAX)*DENS(I)
               TAUHEIOLD = TAUHEI                  
            ELSE
               TAUHEI=TAUTOT(K,102)
            ENDIF
            TAUHEIOLD = TAUHEI

         ENDIF

C     FOR SECOND ITERATIONJ INTERPOLATE THE MEAN INTENISITY AT XTOT(I)
C     
         IF(LAMBDA.GE.2) THEN

            WRITE(6,*)' TO FLUX_INTERPOL NSPEC,XTOT(I)',ISP,
     &           NDIFF,XTOT(I)
C!!!! 

            ISP = 1

            IF(ISP.EQ.1) THEN

C     ISHOCK = NPRE +3

C     OPTICAL DEPTH FROM SHOCK IN CURRENT ITERATION

               TAUHS = TAUTOT(K,3)

               WRITE(6,9276)K,IK,ISHOCK,TAUHS,XTOT(I)
 9276          FORMAT('ISHOCK,TAU ',3I5,1PE12.3,10E12.3)

               CALL FLUX_INTERPOL(IMAX,ISHOCK,XTOT(I),RFSHOCK,TAUHS,
     &              TAUHI,RFD,ISP)

            ENDIF

            IF(ISP.EQ.0.AND.LAMBDA.EQ.1) THEN

               IF(INITISP.EQ.1) THEN
                  INITISP=0

C     SAVE MEAN INTENSITY AND OPTICAL DEPTH FOR SIMPLE CALC.

                  WRITE(6,*)' SAVED SPECTRUM ',K
                  IF(IDEB.EQ.1) WRITE(0,*)' SAVED SPECTRUM ',K

                  DO J=JMIN,JJ
                     
                     FLS(J)=FL(2,J)
                     TAUS(J) = TAUTOT(K,J)
                     WRITE(6,9288)J,TAUS(J),FLS(J)
 9288                FORMAT(' SAVED SP ',I5,1PE12.3,10E12.3)

                  ENDDO

               ENDIF

               WRITE(6,*)' RADSIMP I,K ', IK,K

               CALL RADSIMP(K,TAUS,FLS)

            ENDIF


         ENDIF


         CF = COOLF(I,K,ISHOCK,ISS,DT,TE(I),DENS(I),XE(I),
     &        PHOTOK,PHOTOLM)

         DO IV = 1,IVMAX
            IF(U.GT.VBIN(IV).AND.U.LE.VBIN(IV+1)) THEN
               DELTAU = UOLD - U
 9157          FORMAT('DFDNU ',I5,F10.2,1PE12.3,10E12.3)
            ENDIF
         ENDDO

         UOLD = U

         IF(IDEB.EQ.1) WRITE(0,*)' TO TREC ',IK,DENEL


c     do this for the post shock always, and for the pre shock gas after 
c     the first shell

            irk = 0


c     this applies now

c     change in enthalpy dh/dt = n * lambda

               cfu_old = cfu

               hfu_old = hfu 

               teold = te(i)

               write(6,*)' xe(i),xeold ',i,xe(i),xeold,teold
               xeold = xe(i)

               deold = dens(i)

 111           time=old_time + dt

               xe(i) = xeold

               te(i) = teold

               te_oldi=teold

               dens(i) = deold

               do iel=3,14
                  do iz=1,nionel(iel)+1
                     xion(i,iel,iz)=xion_old(i,iel,iz)
                  enddo
               enddo
               xel=xe(i)
               zeq = xel -abn(1)*xion(i,1,2) -abn(2)*(xion(i,2,2)+
     &              2.d0*xion(i,2,3))

c     step 5 in loop

c     start of iteration?                  

c     if temp too large do a shoryter time step withh same energy loss cool1

               iter=0
c!!!  test this!               
c     te_oldi=te(i)

 7777          iter=iter+1
               if(itime>= 10) then
c                  write(6,*)' itime > 10  use old xions ',itime
                  do iel=3,14
                     do ion=1,iel+1
                        xion(i,iel,ion)=xioni_old(iel,ion)
                     enddo
                  enddo
                  te(i)=te_oldi
               endif

               cool1 = coolf(i,k,ishock,iss,dt,te(i),dens(i),xe(i),
     &              photok,photolm)     

               h = hold(i) + dt* cool1 * nion(i)

               write(6,9167)iter,dt,hold(i),nion(i),cool1,h,te(i)
 9167          format(' New iteration for T ',i5,1pe12.3,10e12.3)

               if(iss.eq.1) then

                  h = hold(i)

               endif

               ishort = 1

c!!!! skipt step shortening 140627

               ishort = 0
c     OK with ns? not i??
               if(xtot(ns).gt.0..and.ishkp.eq.1) then

                  ishkp=0

                  ishort = 0

                  xc = 4.1

c     ishock = k

               endif

c     step 6 in loop

               nc = nc + 1

               vel(nc) = u

               dist(nc) = xtot(i)

               tempk(nc) = te(i)

               iss = 0

               if(nc.ge.51) then

                  xm = 0.
                  tem = 0.
                  nstep = 50
                  rnstep = real(nstep)
                  do ii=1,nstep
                     xm = xm + dist(nc-ii+1)
                     tem = tem + tempk(nc-ii+1)
                  enddo
                  xm = xm/rnstep
                  tem = tem/rnstep                 
                  syy=0.
                  sxx=0.
                  sxy=0.
                  ss = 0.
                  
                  do ii=1,nstep

                     sxx=dist(nc-ii+1)**2 + sxx

                     syy=tempk(nc-ii+1)**2 + syy

                     sxy=tempk(nc-ii+1)*dist(nc-ii+1) + sxy

                  enddo

                  sxx = -rnstep*xm**2 + sxx
                  syy = -rnstep*tem**2 + syy
                  sxy= -rnstep*tem*xm + sxy                 
                  sigx2 = sxx/rnstep
                  sigy2 = syy/rnstep
                  covxy = sxy/rnstep
                  b = sxy/sxx
                  a= tem - b * xm
                  
                  chisq2=0.
                  do ii=1,nstep
                     teinter=b*dist(nc-ii+1)+a
                     chisq2 = (tempk(nc-ii+1)-tei)**2/teinter**2 + 
     &                    chisq2
c     write(6,9175)nc-ii+1,dist(nc-ii+1),tempk(nc-ii+1),
c     &                       teinter,chisq2
 9175                format(' chi2 te ',i5,1pe20.12,10e12.3)                     
                  enddo

                  sig2=chi2/(nstep-1.)
                  tflow = tempk(nc)/(b*vel(nc))
                  tcool = 3. *kb *tempk(nc)/(cool*dens(i))
c                  write(6,9178)nc,time,tempk(nc),dens(i),vel(nc),cool
c     &                 ,b,tflow, tcool,tcool/tflow
 9178             format(' tflow, tcool ',i5,1pe15.5,10e12.3)
                  

                  write(6,9282)te(i),heat,heat-xel*cool,(heat-xel*cool)/heat
 9282             format('te,H,rad,rad/heat ',1pe12.3,10e12.3)
                  if(abs((heat-xel*cool)/heat) < 1.e-2) then         
                     write(6,*)' steady state ',te(i)  
                     iss = 1
c!!!  no steady state!!!                        
                     iss = 0
                     write(6,*)
     &                    ' transit to steady state because tcool< tflow/100'
                  endif
               endif
c!!!  no steady state!!!                        
               iss=0

               if(iss.eq.1) then
                  write(6,*)'steady state assumed!'

                  dtold=dt

                  dt=1.d33

                  tlim_min = 0.9 * te(i)

                  tlim_max = 1.05 * te(i)

                  write(6,9289)(xion(i,14,kk),kk=1,26)
 9289             format(' Fe ion ',1pe12.3,10e12.3)

                  write(6,*)' den1 a ',den1

                  call biseq(i,k,ishock,iss,dt,dens(i),xe(i),
     &                 photok,photolm,tlim_min,tlim_max,
     &                 d00,p0,eb0,ni0,u0,denfl,amean,h,
     &                 chi,xc,nion(i),ne(i),u,p,te(i))
                  write(6,9272)k,xtot(i),dens(i),xe(i),nion(i),te(i)
 9272             format(' after biseq ',i5,1pe20.15,10e13.5)

                  h = hold(i) 
                  write(6,*)' entalpi d ',h
                  dt = dtold
                  dens(i)=nion(i)
                  den1=nion(i)
                  if(itime.eq.1.and.ipre.eq.1) then                     
c     deltax = -1.e15/nion(i)
c     dt = deltax/u
                  endif
                  write(6,*)' nion(i),ne(i),amean,dencgs ',itime,
     &                 nion(i),ne(i),amean,dencgs 
                  deltax = u * dt

c     new x
                  rcgsi(i)=rcgsi(i)+deltax

                  xtot(i) = xtot(i) + deltax

                  write(6,*)' xtot,deltax,u,dt a ',k,xtot(i),deltax,u,dt
c                  write(47,*)' xtot,deltax,u,dt a ',k,xtot(i),deltax,u,dt

               else

                  call tempcalc(i,iss,d00,p0,eb0,ni0,u0,denfl,amean,
     &                 xe(i),chi,h,xc,nion(i),ne(i),u,p,te(i))

                  write(6,*)'Temp change ',itime,iter,te(i),te_oldi,abs(te(i)-te_oldi)/te_oldi,dt
                  if(abs(te(i)-te_oldi)/te_oldi > 0.01 .and. itime>=10)  then
                  dt=dt/2.
                  write(6,*)'Temp change too large ',iter,te(i),te_oldi,abs(te(i)-te_oldi)/te_oldi,dt
c!!   change from 5 -> 10 230626
                  if(iter > 20) then
                     goto 3777
                  else
                     goto 7777
                  endif
               endif
 3777          continue
               do iel=3,14
                  do ion=1,iel+1
                     xioni_old(iel,ion)=xion(i,iel,ion)
                  enddo
               enddo
            endif
c     te_oldi=te(i)
c            write(6,*)' step 6  ',k,i,te(i)               

            tediff = abs(teold-te(i))/teold

c!!!  mod 050807   if(tediff.gt.0.02.and.ishort.eq.1.and.k.ge.1) then
            if(tediff.gt.0.05.and.ishort.eq.1.and.k.ge.1.
     &           and.iss.eq.0) then
               
               dt = dt/2.

c               write(6,*)' shorter time step!',dt,teold,te(i)

               goto 111

            endif

            cfu = cool

            hfu =heat

         dens(i)=nion(i)
         den1=nion(i)
         

         if(lambda.ge.2) then

         endif

         if(ishkp.eq.1) then

            write(67,*)' flux from preionized region',itime
            write(68,*)' flux from preionized region',itime
            iu=2
            call emlines_new(iu,k,kline,1,kli,linelab,lab,filling,
     &           dens(i),xe(i),te(i),xtot(i),u,rabscm,dcit)


            do iel=3,14
               do ion=1,nionel(iel)
                  do line=1,401
                     do iqq=1,md
                        tot_intens(iqq,iel,ion,line)=0.
                     enddo
                  enddo
               enddo
            enddo
         else
            iu=2
            call emlines_new(iu,k,kline,1,kli,linelab,lab,filling,
     &           dens(i),xe(i),te(i),xtot(i),u,rabscm,dcit)
c     call sort(ik,100,0)
         endif

c     step 7 in loop

         deltax = u * dt

c     new x
         rcgsi(i)=rcgsi(i)+deltax
         
         xtot(i) = xtot(i) + deltax

 545  continue
      call coldens(k,xtot(i),te(i),del(i),filling,deltax,cold,coldn,coltottemp)


c     note that coltottemp is not used for the interpolation, only coltot
c     first iteration it can be used directly, but for higher iterations
c     it need to be calculated from shock first, since coltottemp is 
c     calculated from the preshock boundary

      dtauhei=ta(102)*abs(deltax)*nion(i)
      tauhei=tauhei+dtauhei

c     energy lost by cooling (minus because it is cooling)

      decool = cool2*nion(i)*ne(i)*abs(deltax)

      t4=te(i)/1.e4

      totcool = totcool - decool
      shockenergy = 0.5*d00*u0**3

c     energy lost by emission

      totem=0.
      toteml=0.
      totemc=0.
      totemi=0.
      totemin=0.
 9823 format(' delta em ',i5,1pe12.3,10e12.3)

      totemt=totemt + totem*dt

      totemit=totemit + totemi*dt

      totemitn=totemitn + totemin*dt

      ekin= 0.5 * dencgs * u**2 

      eterm = 1.5*(nion(i)+ne(i))*kb*te(i) 

      eion = 0.

      do iel=3,14
         do iz=1,nionel(iel)
            eion = eion + abn(iel)*elch*chi(iel,iz)*
     &           xion(i,iel,iz+1)
         enddo
      enddo

      eion=0

      eion = eion*den(i)

      etot = ekin + eterm + eion

      eflux = u*etot

      delta_etot = etot -etot_old

      dheat = heat * nion(i) * deltax

      totheat = totheat + heat *nion(i) * dt

      dcool = cool2 *nion(i)*ne(i) * deltax

      dcool2 = cool2 *nion(i)*ne(i) * dt

      totcoolt = totcoolt + cool *nion(i)* dt

      totcoolet = totcoolet + cool *ne(i)* dt

      delta_eflux = eflux - eflux_old

      dc=cool2 *nion(i)*ne(i)
      dh = heat *nion(i)**2

c     calculate velocity gradient
      
      if(istat.eq.1) then

         ivgrad = 0

      else

         ivgrad = 1

      endif

      if(ivgrad.eq.1) then
         nc = nc + 1

         vel(nc) = u

         dist(nc) = xtot(i)


         if(nc.ge.11) then

            dr_dv1 = (dist(nc)-dist(nc-10))/(vel(nc)-vel(nc-10))

            nsig = 0
            sig2 = vel(nc-5)**2

            nstep = 10
            dr_dv_old = dr_dv

 1444       nsig = nsig+1
            
            sx=0.
            sy=0.
            sxx=0.
            sxy=0.
            ss = 0.
            
            do ii=1,nstep
               ss = ss + 1./sig2
               sx=dist(nc-ii+1)/sig2 + sx
               sxx=dist(nc-ii+1)**2/sig2 + sxx
               sy=vel(nc-ii+1)/sig2 + sy
               sxy=vel(nc-ii+1)*dist(nc-ii+1)/sig2 + sxy
            enddo
            
            delta=ss*sxx-sx**2
            
            if (delta.eq.0.) then
               nstep = 2*nstep
               goto 1444
            endif
            a=(sxx*sy-sx*sxy)/delta
            b=(ss*sxy-sx*sy)/delta
            
            dr_dv = 1./b

            chi2=0.
            do ii=1,nstep
               vint=b*dist(nc-ii+1)+a
               chi2=(vel(nc-ii+1)-vint)**2/vint**2+chi2
               
            enddo
            
            if (chi2.gt.1..and.nstep.ge.40) then
               dr_dv = dr_dv_old
               goto 1445
            endif

            if (chi2.gt.1.) then
               nstep = 2*nstep
               goto 1444
            endif
            
            if(nsig.le.2) then
               sig2=chi2/(nstep-1.)
               goto 1444
            endif
            
 1445       continue
            
            dr_dv = abs(dr_dv)
            

         else
            
            dr_dv = time
            
         endif

c Obs! Sobolev evrerywhere!!!         
         dr_dv = time
         
         if(dr_dv.gt.time) then
            
            dr_dv = time
            
         endif

      endif

c$$$      write(6,*)'Ionization balance at i= ',i,k
c$$$      do iel=3,14
c$$$         write(6,9138)iel,(xion(ik,iel,iz),iz=1,nionel(iel)+1)
c$$$ 9138    format(' xion ',i5,1pe12.3,28e12.3)
c$$$      enddo
      dist_fr_sh(k)=xtot(i)
c$$$      write(6,91)itime,k,i,dist_fr_sh(k),time,dr(k),xc,u,p,nion(i),ne(i),
c$$$     &     xe(i),tautot(k,3),cool1,te(i),
c$$$     &     (xion(i,1,iz),iz=1,2),(xion(i,2,iz),iz=1,3),
c$$$     &     (xion(i,5,iz),iz=1,8),
c$$$     &     (xion(i,14,iz),iz=1,18)
      if(icsm==1.and.ipre==0) then
         rabs=rshock_cgs-xtot(i)
      elseif(icsm==1.and.ipre==1) then
         rabs=xtot(i)
      endif
c     Fe I=28, XXVI=53, Si I=54, XIV=67, S I=68, XVI=83, Ar I=84, Ar XX=101, Ca I=102            

      write(63,9291)k,xtot(i),xc,u,p,nion(i),ne(i),xe(i),
     &     tautot(k,3),cool1,te(i),
     &     (xion(i,5,ion),ion=1,8),(xion(i,10,ion),ion=1,14),
     &     (xion(i,11,ion),ion=1,16),(xion(i,12,ion),ion=1,18),
     &     (xion(i,13,ion),ion=1,10),(xion(i,14,ion),ion=1,14)
 9291 format(i4,1pe15.7,3e11.3,131e11.3)
c$$$      write(56,9562)itime,k,ipre,rabs,xtot(i),time,xc,u,p,nion(i),ne(i),
c$$$     &     xe(i),tautot(k,3),cool1,te(i),
c$$$     &     (xion(i,1,iz),iz=1,2),(xion(i,2,iz),iz=1,3),
c$$$     &     (xion(i,5,iz),iz=1,8),(xion(i,14,iz),iz=1,26),
c$$$     &     (xion(i,10,iz),iz=1,14),
c$$$     &     (xion(i,11,iz),iz=1,16),(xion(i,12,iz),iz=1,18),
c$$$     &     (xion(i,13,iz),iz=1,20)
      iz=2
c$$$      write(6,*)'T = ,H =15-16,, He=17-19, O=20-27, Fe=28-63, Si=64-77, S=78-93, Ar=94-111, Ca=112-131'
c$$$      write(6,9481)xion(i,1,iz),xion(i,2,iz),
c$$$     &     xion(i,5,iz),xion(i,14,iz),
c$$$     &     xion(i,10,iz),xion(i,11,iz),xion(i,12,iz),xion(i,13,iz)
c$$$ 9481 format(' First ion',1pe11.3,20e11.3)
c$$$ 9562 format(3i7,1pe22.14,2e18.9,1pe11.3,150e11.3)
c$$$c     if(ideb.eq.1) write(0,91)itime,k,xtot(i),time,dr(k),xc,u,p,nion(i),ne(i),
c$$$c     &           xe(i),tautot(k,3),cool1,te(i),(xion(i,5,iz),iz=1,8),
c$$$c     &           (xion(i,1,iz),iz=1,2),(xion(i,2,iz),iz=1,3)
c$$$
c$$$c     if(ideb.eq.1) write(0,9991)itime,k,xtot(i),time,xc,nion(i),te(i),fl(2,3),
c$$$c     &           fl(2,105),fl(2,142),fl(2,211),fl(2,256),fl(2,315)
c$$$      write(6,9122)ik,decool,totem
c$$$ 91   format('x,xc ',3i7,1pe18.10,1e15.7,450e12.3)
c$$$ 9991 format('x,sp ',2i7,1pe16.7,1e15.7,20e12.3)
c$$$ 9122 format('decool,emtot',i7,1pe15.7,1e15.7,10e12.3)
c$$$      write(6,*)'xion g ',(xion(i,1,iz),iz=1,2)
      hold(i)= h

      if(ns.ge.4) then
         
         xtot(1) = xtot(2) + (xtot(2)-xtot(3))
         do ix=2,ns
            dr(ix) = xtot(ix) - xtot(ix)
         enddo

      endif

      if(xtot(i).lt.-1.e16) then

         write(6,*)' x < -1.e16 cm',xtot(i)

c     stop

      endif

c!!!  stop  integration of shocked region if te < 1000 K and goto 2321 for renumbering (after writing spectrum etc)
      temin_shock=2.5e2
c      temin_shock=2.e4
      if(te(ns).lt.temin_shock.and.itime.ge.25.and.ipre.eq.0) then
         write(6,*)' shock calculation done! te_min_shock',temin_shock
         goto 2321
      endif
      if(lambda.eq.1) then

         tsec_num0 = time * dnum0

          if(tsec_num0.gt.t_ioniz) then


 9635       format(i5,f12.3,1pe12.3,20e12.3)

            
c     if t x n > t_ion goto pre-shock

            write(6,*)' te x n > t_ion ',tsec_num0,t_ioniz

            goto 2321

         endif

         if(te(ns).lt.100.) then 

            do iv = 1,ivmax
               vmid = (vbin(iv)+vbin(iv+1))/2.
c     q               deltawl = wli(63)*(1.+vmid/3.e10)-wli(63)
               write(6,9057)iv,deltawl,vbin(iv),vbin(iv+1),dfdnu(iv)
 9057          format('lineprof',i5,1pe12.3,10e12.3)
            enddo

         endif

      else
         if(te(ns).lt.temin_pre) goto 2321
      endif

      if(itime.ge.2) then
         
         veli = u

         if(lambda.ge.3) then
            
            dxtotsh = rtotsh

         else

            dxtotsh = 1.e33

         endif

         call timestep(iss,k,ns,lambda,xtot(ns),veli,den(ns),
     &        dxtotsh,dt,dtnew)

         if(k > 5 .and. te(i) >= 3.0e3) then
c     use 0.02 x cooling time as time step
            dtnew=0.02*tcool
            dtnew_old=dtnew
c     dtnew=0.04*tcool
            write(6,*)' time step determined by cooling ',k,itime,dtnew
         elseif(k>5.and.te(i) < 3.e3) then
c 501
            dtnew=dtnew_old*(3.e3/te(i))**4
c 502
            dtnew=min(dtnew,3.e4)
c 503
            dtnew=min(dtnew,1.e5)            
c            dtnew=0.04*tcool
            write(6,*)' time step NOT determined by cooling ',k,itime,te(i),dtnew_old,dtnew
         else
            dtnew=2.e3
         endif


         dt = dtnew

      endif

      if(stop_tmin==1) then
         write(6,*)' T_min reached ! ',i,k,te(i)
      endif

c     main time step loop stops just here!      
      enddo


 2321 imode = 2

c     save last point

      k = k+1

      if(k.ge.md) stop

      r(k) = xtot(ns)

      rtr(2,k) = r(k)
      drtr(2,k) = r(k)- r(k-1)
      dentr(2,k) = dens(ns)

c     extra step for last point 
      do iel=1,14
         do iz=1,nionel(iel)+1
            coltr(2,k,iel,iz) = coltottemp(iel,iz)
         enddo
      enddo

      ik = k
c???
      ik = i

      do iel=1,14
         do j=1,nionel(iel)+1
            xion(k,iel,j)=xion(ns,iel,j)
         enddo
      enddo

      call elecdens(k,xe(k),zelec)

      call emiss(te(ns),xe(ns),dentr(2,k))

      call emiss_cont(te(ns),xe(ns),dentr(2,k),cooli)

      iqq=3
      call opac_source(iqq,1,drtr(2,k),dentr(2,k),te(i),xe(i))

c     save opacities

      do j=jmin,jj
         copac(k,j) = ta(j)
      enddo

      if(lambda.eq.1.and.ipp.eq.0) then

         ipp = 1

         nshock = k


         write(6,*)' before renum ',k

         do ij=1,nshock
            do j=jmin,jj
               tautots(ij,j) = tautot(ij,j)
               sources(ij,j) = s(ij,j)
               ems(ij,j) = em(ij,j)
               emcs(ij,j) = emc(ij,j)
               if(ij.eq.2) then                  
               endif
               copacs(ij,j) = copac(ij,j)
            enddo
            rs(ij) = rtr(2,ij)
            drs(ij) = drtr(2,ij)
            densave(ij) = dentr(2,ij)

            do iel=1,14
               do iz=1,nionel(iel)+1
                  coltrs(ij,iel,iz) = coltr(2,ij,iel,iz)
               enddo
            enddo

         enddo

c     now renumber with i=1 at post shock boundary and i=nshock at shock
         call renum(k,nshock)
         
         write(6,*)' after renum ',nshock
         write(6,*)' em f ',k,emc(k,3),em(k,3)

         jx=-28

         write(6,*)' after inversion in trans'
         do ij=1,nshock
c     if(i.lt.5.or.i.gt.nshock-5) then
            write(6,9232)ij,r(ij),xtot(ij),rtr(2,ij),drtr(2,ij),
     &           dentr(2,ij),tautot(ij,jx),em(ij,jx),copac(ij,jx),
     &           s(ij,jx),coltr(2,ij,1,1)
 9232       format('i,r,xtot,rtr,drtr,dentr,tautot,em,cop,s,N ',i5,1pe15.7,10e12.3)
c            endif
         enddo


         

         

         write(6,*)' To emlines_new at ',k,ik         
         write(6,*)'Strongest lines for post-shock ',k, ik
c         write(101,*)'Strongest lines for post-shock ',k, ik
         iu=2
         call emlines_new(iu,k,kline,1,kli,linelab,lab,filling,
     &        dens(i),xe(i),te(i),xtot(i),u,rabscm,dcit)
         
         call sort(k,500,1)


      endif


 9100 FORMAT(1X,'*******************************************************
     &******************')
 9899 FORMAT(A)
 112  FORMAT(12(1X,1PE8.2))
 9776 FORMAT(' INITIAL VALUES OF A1 AND B1 : ',1PE12.5,E12.5)
 7778 FORMAT(1X,'NUMBER OF IONIZING PHOTONS: S(13.6-54 EV)=',1PE13.6
     &     ,' S( >54 EV)=',E13.6)
 7779 FORMAT(1X,'IONIZATION PARAMETER: U1(13.6-54 EV)=',1PE13.6
     &     ,' U2( >54 EV)=',E13.6)
 2758 FORMAT(1X,'N(O)=',5E10.3)
 2759 FORMAT(1X,'N(CA)=',5E10.3)
 2757 FORMAT(1X,'N(HE)=',7E10.3,/,7E10.3)
 2756 FORMAT(1X,'N(H)=',7E10.3)
 9477 format(1x,'Balm. ion. Ly a:',f7.2,' CIV:',f7.2,' CIII:',f7.2,
     &     ' MgII:',f7.2,' %')
 9478 format(1x,'P(H)=     ',1pe10.2,4e10.2)
 9479 format(1x,'X(H)*P(H)=',1pe10.2,10e10.2)
 9264 format(2i6,f10.1,f10.3,1pe11.3,10e11.3)
 2754 FORMAT(1X,'N(FE I)=',7E10.3,/,7E10.3)
 2755 FORMAT(1X,'N(FE II)=',7E10.3,/,7E10.3)
 9572 FORMAT(' H, HE ',1PE10.2,9E10.2)
 9573 FORMAT(' C     ',1PE10.2,9E10.2)
 9574 FORMAT(' N     ',1PE10.2,9E10.2)
 9575 FORMAT(' O     ',1PE10.2,9E10.2)
 9576 FORMAT(' SI    ',1PE10.2,9E10.2)
 9577 FORMAT(' SI    ',1PE10.2,9E10.2)
 9578 FORMAT(' MG    ',1PE10.2,9E10.2)
 9579 FORMAT(' FE    ',1PE10.2,9E10.2)
 9580 FORMAT(' AL    ',1PE10.2,9E10.2)
 9589 FORMAT(' AR    ',1PE10.2,9E10.2)
 9581 FORMAT(' CA    ',1PE10.2,9E10.2)
 9582 FORMAT(' NA    ',1PE10.2,9E10.2)
 9583 FORMAT(' S     ',1PE10.2,9E10.2)
 9584 FORMAT(' NE    ',1PE10.2,9E10.2)
 9532 FORMAT(' H2,OH,CO,CO+,H2O ',1PE10.2,9E10.2)
 525  FORMAT(' ')
 9721 FORMAT(1X,'ABUND ',7E10.3)
 9722 FORMAT(1X,1PE12.3,17E12.3)
 900  FORMAT(1X,'I=',I3,1X,'V=',F9.3,1X,'DX=',1PE12.3,1x,
     &     'N=',E9.3,1X,'X=',E9.3,1X,'T=',0PF9.1)
 9910 FORMAT(1X,'I=',I3,1X,'DR=',1PE11.5,' CM',1X,
     &     'DX=',e12.3,'DENS=',E10.4,1X,
     &     'EL=',E10.4,1X,'TEMP=',E10.4)
 9911 FORMAT(1X,'I=',I3,1X,'DR=',1PE11.5,' CM',
     &     1x,'DENS=',E10.4,1X,
     &     'EL=',E10.4,1X,'TEMP=',E10.4)
 9458 FORMAT(' MINIMUM AND MAXIMUM VELOCITY OF THE OXYGEN ZONE',
     &     1PE12.4,E12.4)
 9283 FORMAT(1X,'M(R)=',1PE10.3,1X,'N=',E10.3,1X,'TAU=',E10.3
     &     ,' DEN=',E10.3,' <A>=',0PF6.2)
 903  FORMAT(1X,'ENERGY (TOT) ',1PE12.4,' (PRIM)',E12.4,
     &     ' 13.6-54 ',E12.5,' > 54 EV '
     &     ,E12.5,' DE ',E12.5,' DEDV ',E12.4,' DEDV(LINES) ',E12.4,
     &     ' R=',E13.6)
 904  FORMAT(1X,'COOLING=',8E12.5)
 912  FORMAT(1X,'FLUX AT INNER BOUNDARY BELOW 13.6 EV=',1PE12.5,
     &     ' ABOVE 13.6 EV',E12.5)
 914  FORMAT(1X,'FLUX AT SURFACE BELOW 13.6 EV=',1PE12.5,
     &     ' ABOVE 13.6 EV',E12.5)
 913  FORMAT(1X,'DIFFUSE EMISSION BELOW 13.6 EV=',1PE12.5,
     &     ' ABOVE 13.6 EV',E12.5)
 3966 FORMAT(1X,'T(13.6)=',1PE12.4,1X,'T(24.5)=',E12.4,1X,
     &     'T(54.0)=',E12.4,1X,'E(TAU=1)',E12.4)
 1008 FORMAT(1X,'HE II (1640) =',1pE10.3,'  HE II (4686)  =',E10.3)
 1098 FORMAT(1X,'He II (304)  =',1pe10.3,'  HE II (Two-g) =',E10.3)
c     &'  C III (1909) REC.=',E10.3,'  N IV (1486) REC.=',E10.3)
 7468 FORMAT(' REC=',8E11.4)
 7469 FORMAT(' RES=',8E11.4)
 7467 FORMAT(1X,7E12.4)
 7278 format(1x,7f10.1)
 901  FORMAT(1X,'RADIUS= ',1PE14.7)
 9372 FORMAT(' RECCO',8E11.4)
 9373 FORMAT(' REC',8E11.4)
 9420 FORMAT(1X,5E12.5)
 699  FORMAT(1X,'FINAL RESULTS')
 1001 FORMAT(1X,'C II  (1334) =',1PE9.2,' (2326) =',E9.2
     &     ,'  C III (977) =',E9.2,/,' C III (1909) =',E9.2,
     &     '  C IV (1550)  =',E9.2)
 1002 FORMAT(1X,'N III (1750) =',1PE9.2,'  N III (990)  =',E9.2,/,
     &     ' N IV (1486)  =',E9.2,'  N V (1240)   =',E9.2)
 1003 FORMAT(1X,'O III (1664) =',1PE9.2,'  O IV (1400)  =',E9.2,/,
     &     ' O V (1216)   =',E9.2,'  O VI (1034)  =',E9.2)
 1004 FORMAT(1X,' SI II (2335) =',1PE9.2,
     &     ' SI III (1892) =',E9.2,' SI IV (1400) =',E9.2)
 1005 FORMAT(1X,'LA=',1PE9.2,'  LB=',E9.2,'  LG=',
     &     E9.2,'  LD=',E9.2)
 1006 FORMAT(1X,'HA=',1PE9.2,'  HB=',E9.2,'  HD=',E9.2,'  PA=',
     &     E9.2,'  PB=',E9.2,'  BA=',E9.2)
 1009 FORMAT(1X,'N I (3468) =',1PE9.2,' (5201) =',E9.2,' (10406) ='
     &     ,E9.2,'  N II (3063) =',E9.2,' (5755) =',E9.2,' (6548-83) ='
     &     ,E9.2)
 1010 FORMAT(1X,'O I   (2964) =',1PE9.2,' (5581) =',E9.2,' (6302-64) ='
     &     ,E9.2)
 1034 FORMAT(' O II  (2470) =',1PE9.2,' (3726-29) =',E9.2,' (7320-30) ='
     &     ,E9.2)
 9171 FORMAT(1X,'O I (1302) =',1PE9.2,' (1356) =',E9.2,' (7774) ='
     &     ,E9.2,' (8446) =',E9.2)
 1011 FORMAT(1X,'O III (2321) =',1PE9.2,' (4363) =',E9.2,
     &     ' (4959-5007) =',E9.2)
 9116 FORMAT(' MG I  (3834) =',1PE9.2,' (4571) =',E9.2,
     &     ' (5176) =',E9.2,' (8806) = ',E9.2)
 9101 FORMAT(' MG I  (4571) =',1PE9.2,' (2853) =',E9.2,
     &     ' MG II (2800) =',E9.2)
 9108 FORMAT(' CA II (7300) =',1PE9.2,' (3950) =',E9.2,
     &     ' (8500) =',E9.2)
 9115 FORMAT(' NA I (5889) CE =',1PE9.2,' (5889) REC =',E9.2)
 9118 FORMAT(' SI I (10995) =',1PE9.2,' (16360) =',E9.2,
     &     ' S II (2335) =',E9.2)
 9119 FORMAT(' S II (6716-31) =',E9.2,' (4069-76) =',E9.2,
     &     ' (10286-370) =',E9.2)
 9151 FORMAT(' FE I   =',1PE9.2,' FE II =',E9.2,' FE III =',E9.2,
     &     ' FE IV =',E9.2)
 1013 FORMAT(1X,'C I   (2967) =',1PE9.2,' (4619) =',E9.2,' (8729) ='
     &     ,E9.2,'  (9812) =',E9.2)
 2011 FORMAT(' C I   (370.4 M) =',1PE9.2,' (609.1 M) =',E9.2,
     &     ' C II  (157.7 M) =',E9.2)
 1012 FORMAT(' O I   (63.17 M) =',1PE9.2,' (145.6 M) =',E9.2)
 2012 FORMAT(' NE II (12.81 M) =',1PE9.2)
 2013 FORMAT(' SI I  (68.49 M) =',1PE9.2,' (129.7 M) =',E9.2,' SI II ',
     &     '(34.81 M) =',E9.2)
 2014 FORMAT(' S  I  (25.25 M) =',1PE9.2,' (56.31 M) =',E9.2,' S  II ',
     &     '(314.5 M) =',E9.2)
 2015 FORMAT(' A  II (6.985 M) =',1PE9.2)
 1014 FORMAT(' EPOCH =',F8.1,' DAYS  MAX RADIUS = ',1PE11.4)
 1015 FORMAT(' TAU(GAMMA) =',1PE11.4,' GAMMA LUM =',E11.4,
     &     ' TAU(EL) =',E11.4)
 1016 FORMAT(' M(TOT) = ',F7.2,' M(NI) = ',F8.4,
     &     ' V(EXP) = ',1PE11.4)
 9745 FORMAT(' COL(N=1, N=2), TAU(BALMER)',1PE12.3,4E12.3)
 475  FORMAT(I4,E11.4,F11.2,E10.3,10F7.3)
 477  FORMAT(4E11.4,E10.3,14F7.3)
 476  FORMAT(I4,2E11.4,17F7.3)
 9430 FORMAT(1X,E11.5,14F7.3)
      RETURN
      END

      subroutine coldens(k,xtot,tei,deli,filling,dx,cold,coldn,coltot)
      IMPLICIT REAL*8(A-H,O-Z)
      include 'PARAM'   
      PARAMETER (NFEL=3000)
      parameter (nlp=30000)
      COMMON/IND/I
      COMMON/DXA/DR(MD)
      COMMON/PHY/DEN(MD)
      common/abl/abn(15)
      common/ionx/xion(md,14,27)
      parameter(nspec=91)
      DIMENSION COLD(14,27),COLDN(14,27),coltot(14,27)
      real*8 atw(15)
      data atw/1.,4.,12.,14.,16.,20.,22.,24.,26.,
     &28.,32.,36.,40.,56.,58./
      dimension nionel(14)
      data nionel/1,2,6,7,8,10,11,12,13,14,16,18,20,26/
      common/atmolw/atmol(nspec)
      common/preion/ipre
      common/tauelec/tau_tot_es,tau_elec_sh,tau_elec(md)

C     ***********************************************************
C     *****
C     CALCULATE COLUMN DENSITY * SQRT(M/TE) TO BE USED FOR LINE DEPTH
C     *****
C     ***********************************************************

      DO iel=3,14
         do iz=1,nionel(iel)

            COLDN(iel,iz)=ABN(iel)*xion(i,iel,iz)*DEN(I)*FILLING*abs(dx)
     &           + COLDN(iel,iz)

            COLD(iel,iz)=ABN(iel)*xion(i,iel,iz)*SQRT(atw(iel)/TEI)

            COLTOT(iel,iz)=COLD(iel,iz)*DEN(I)*abs(dx)*FILLING+
     &           COLTOT(iel,iz)
         enddo
      enddo

      return
      end


      subroutine ion_pot(nion,chi)
      implicit real*8(a-h,o-z)
      dimension chiry(30,30),chi(30,30),nion(30)


      open(11,file='./ATDAT/ionization_pot_z.dat',status='old')

      do iz=1,30
         read(11,*)izq,(chiry(iz,kk),kk=1,iz)
         if(iz.le.2) then
            iel=iz
            incl = 1
         elseif(iz.ge.6.and.iz.le.8) then
            iel=iz-3
            incl = 1
         elseif(iz.ge.10.and.iz.le.14) then
            iel=iz-4
            incl = 1
         elseif(iz.eq.16) then
            iel=iz-5
            incl = 1
         elseif(iz.eq.18) then
            iel=iz-6
            incl = 1
         elseif(iz.eq.20) then
            iel=iz-7
            incl = 1
         elseif(iz.eq.26) then
            iel=iz-12
            incl = 1
         elseif(iz.eq.27) then
            iel=iz-12
            incl = 1
         else
            incl = 0
         endif
         if(incl.eq.1) then
            nion(iel)=iz
            do kk=1,iz
               chi(iel,kk)=chiry(iz,kk)*13.5987
            enddo
            write(6,9)iz,iel,(chi(iel,kk),kk=1,iz)
 9          format(2i5,30f10.3)
         endif
      enddo

      return
      end


      subroutine compress(p0,eb0,ni0,u0,denfl,hprim,xc)
c 
c newton-raphson for the shock compression factor x
c
      implicit real*8(a-h,o-z)
      real*8 ni0

      eps=1.e-10
      x1 = xc
      xcold = xc
      do i=1,20
         call fent(p0,eb0,ni0,u0,denfl,hprim,x1,funx,dfdx)
         x2 = x1 -funx/dfdx
c         write(6,98)i,x1,x2,funx,dfdx
 98      format('fent',i5,1pe12.3,10e12.3)
         if(abs((x2-x1)/x1).lt.eps) then
            goto 1
         endif
         x1 = x2
      enddo
      
 1    xc = x2

      diff = (xc-xcold)

      if(xcold.gt.1.e2) then

         eps2 = 0.5e-4
      else
         eps2 = 1.e-2
      endif

      if(abs(diff/xcold).gt.eps2) then

c         xc = xcold + diff/4.

c         write(6,97)i,x1,x2,funx,dfdx,xcold,xc
 97      format('damp xc',i5,1pe12.3,10e12.3)

      endif

      return
      end

      subroutine fent(p0,eb0,ni0,u0,denfl,hprim,x,funx,dfdx)
c
c equation for shock compression (see cox 1972)
c
c hprim = h - eion
c
      implicit real*8(a-h,o-z)
      real*8 ni0

c      write(6,*)'rk p0,eb0,ni0,denfl,hp,x ',p0,eb0,ni0,denfl,hprim,x

      funx = eb0 *x**2/5. + 2.*hprim*ni0*x/5. + 
     &     4.*denfl*u0/(5.*x)-p0

      dfdx = 2.*eb0 *x/5. + 2.*hprim*ni0/5. - 
     &     4.*denfl*u0/(5.*x**2)

c      write(6,9)eb0,denfl,h,ni0,u0,x,funx,dfdx
c      write(6,9)eb0,denfl,hprim,ni0,u0,p0,x,funx,dfdx
 9    format(' rk comp ',1pe12.3,10e12.3)
      return
      end





      subroutine emlines_new(iu,ik,kline,ioutp,kli,linelab,lab,filling,
     &              densi,xeli,tei,xdist,velo,rabscm,dcit)

      IMPLICIT REAL*8(A-H,O-Z)
      character*72 text1,text2,text3,TEXT4,TEXT5,INMOD
      CHARACTER*8 LAB(200)
      parameter (nlp=30000)
      CHARACTER*18 LINELAB(300),linelabex(nlp)
      character*7 labmol
      REAL*8 MTOT,MNI56
      include 'param'
      PARAMETER (NL=340,NLP1=NL+1)
      PARAMETER (NFEL=3000)
      common/ionx/xion(md,14,27)
      COMMON/TEXT/TEXT1,TEXT2,TEXT3,TEXT4,TEXT5
      COMMON/MOD/INMOD
      COMMON/IOINF/IOSPEC,IOLEVEL,IOEMISS,IOTERM,IOION,IOELITER
      COMMON/THBAL/TMIN,TMAX,IBAL,IKBAL(1000)
      common/sorc/sor(NE1:NE2),KM(NE1:NE2)
      COMMON/CONT/CONR(NL),CONT(NL)
      COMMON/SECX/CSEC(20),DCSDR(20),CISEC(20),DCISDR
      COMMON/DENS/DEN0,R0,RN
      COMMON/TAUFBO/TAFB(30,3)
      COMMON/QSOM/qso,IAGN,ISTAT,ISOBOL
      COMMON/INSPEC/IBLACKB,IPULSSP,IFERNET,INLR,IFERMAT,ICRAB,
     &     IFREEFREE,imewe,iobstar
      common/initmw/initmw,initr2
      COMMON/AGNPARA/DECON,GAMMION,IPRESS,IDEN,ICDEN
      COMMON/PULSAR/RLTOT,ALFA,EMIN,EMAX
      COMMON/FNORM/FNORM
      COMMON/LOWION/ILOWION
      COMMON/SPH/ISPH
      COMMON/INUT/IUY
      COMMON/RQW/TEFF,RQ
      COMMON/TPAR/RIN,DRQ,R1Q,TDAYS
      COMMON/INI/INITH,INIT,INITCA,INITHE,INITFE,initfei
      common/croinit/initcr
      COMMON/ITER/ITE
      COMMON/MPAR/MTOT,MNI56,VEXP
      COMMON/TERMA/TERMAL
      COMMON/QSPEC/GAMMA,ALQSO
      COMMON/LITER/N
c      COMMON/IND/Ik
      COMMON/FRE/NINQ,JMIN,JJ
      COMMON/ITSPH/NITSPH
      COMMON/THER/A1,B1,TIN,E10,E20
      COMMON/SPOT/OTSP(16)
      COMMON/COL/RTE(4),FF,HY,HE,C131,C142,C151,C161,C231,C241,COH,
     &COHE,C351,C321,C331,C341,C222,C232,C221,C332,C441
      COMMON/EQUIV/RADF,FR(2,nlp),W(nlp),CIN(nlp),WLI(nlp+NFEL),FB(300),tauline(nlp+nfel)
      COMMON/EQUIH/FRH(2,401),WEH(401),SR(NFEL+nlp),WOBS(NL,NL)
      COMMON/GAMPAR/GAMLUM,GAELH,GAHE
      COMMON/T/TES,SS
      COMMON/REHEL/REHE21,REHE22,REC31,REN41
      COMMON/CHG/CHGH,CHGHE,CHTE
      COMMON/GSREC/ALGS(100),ALEX(100),ALTOT(100),RECEX(100),RECGS(100) 
      COMMON/HEA/CO,PHEO,PHEI,PO(8),PC(6),PN(7),PMG,PSI,PFE,
     &     pnei(10)
      COMMON/INT/FL(2,NE1:NE2),SI(14,27,NE1:NE2),E1(NE1:NE3),E(NE1:NE3)
      COMMON/DIFFP/FH0(NE1:NE2),FHD(NE1:NE2)
      COMMON/OTS/NOTS
      COMMON/EDDFLUX/EDDFLUX(NE1:NE2),SURFJ(NE1:NE2)
      COMMON/DIF/EM(MD,NE1:NE2),TAU(MD,NE1:NE2),TAUTOT(MD,NE1:NE2),
     &           EMC(MD,NE1:NE2)
      COMMON/REC/AL2(16)

      parameter (mh=2,mhe=3,mc=7,mn=8,mo=9,mne=11,mna=12,mmg=13,
     &     mal=14,msi=15,ms=17,mar=19,mca=21,mfe=27,mni=29,
     &     mmat=mh+mhe+mc+mn+mo+mne+mna+mmg+mal+msi+ms+mar+mca+mfe+mni,
     &     ispecy=15)

      COMMON/RADIE/R(0:MD)
      COMMON/DXA/DR(MD)
      COMMON/ELEC/DEL(MD)
      COMMON/PHY/DEN(MD)
      COMMON/MASSES/RM(MD)
      COMMON/TEM/TE(MD)
      COMMON/NHY/NH
      COMMON/NLEV/e00,NIONI,NHY,NP1H,NMA,NMI
      COMMON/EMHY/RECEM(6,NL),TWOPH,TWOPHHEI,COHI,PHQ,PHT,POQ,POT
      common/coolants/cl(nlp)
      common/pescpd/pdbe(nlp,2)
      COMMON/LYBFLUOR/OIFLUOR,OIFLUORH,OIFLUORO,OILYB,OIREC
      COMMON/A26/SIGO(20,NE1:NE2),SIGCA(10,NE1:NE2),SIGHE(16,NE1:NE2),
     &SIGFE(NL,NE1:NE2),sigfei(nl,ne1:ne2),SIGH(50,NE1:NE2),SIGHE2(50,NE1:NE2)
      COMMON/FESAVE/INIFE(4),inifei(4),BOLFES(4,NL),TOLDFES(4),
     &     bolfeis(4,nl),toldfeis(4)
      common/deltar/drscgs
      common/linepoint/nlines,jline(nlp),iion(nlp)
      common/linephoto/phline(nlp)
      common/escapemax/be00
      common/contpoint/otscon(100),ncont,jcont(100),iioncont(100)
      common/oiiilines/em4959,em5007
      common/multisol/imultih,imultihe,imultio,imultica,imultifei,
     &     imultifeii
      parameter(nspec=91)
      common/label/labmol(100)
      common/molpop/xmol(nspec)
      common/abl/abn(15)
      common/timecheck/time,itime
      common/line_em/cinx(14,26,401),taulinex(14,26,401),wlix(14,26,401)
     &     ,ilabcfx(14,26,401)
      common/tot_line_em/tot_intens(md,14,26,401)
      real*8 tot_line_lum
      common/linelum/tot_line_lum(14,27,401)
      COMMON/mdIND/kradius
      dimension nionel(14)
      data nionel/1,2,6,7,8,10,11,12,13,14,16,18,20,26/
      character*2 elid(14)
      data elid/'H ','He','C ','N ','O ','Ne','Na','Mg','Al','Si','S ','Ar','Ca','Fe'/
      dimension totfl(ne1:ne2),contfl(ne1:ne2),tototsl(ne1:ne2)
      DIMENSION DRA(MD),COLD(14,27),COLDN(14,27),FBI(110,MD)
      DIMENSION DLO(50)
      DIMENSION FL0(2,NE1:NE2),WREC(30),WRHE(10),HER(10),OXR(10)
      DIMENSION DRECHI(NL),BX(NL)
      DIMENSION DWFEIII(30,30),DWFEIV(30,30)
      DIMENSION DCIT(nlp),DFBQ(110),RMO(MD),TEINT(MD),
     &          SIRED(nlp+NFEL,MD),deint(md),xelint(md),rint(md)
      DIMENSION BN(MD,5),THA(MD),THB(MD),ABUN(MD,20)
      DIMENSION BQH(NL),BQ(NL),BQCA(NL),BQHE(NL),BQFE(NL),BQFEi(NL)
      DIMENSION DLW(100),OLDFLU(NE1:NE2)
      DIMENSION BBOL(10,MD,NL),axsi(20),gaion(20)
      DIMENSION WLF1(300),IABU(20),ABOLD(20),KLI(300),KKMAX(6),
     &          SPL(100)
      DIMENSION WLSP(100),TREC(100)
      dimension wlsf1(5000),wlsf2(5000),dwsfei(5000),dwsfeii(5000)
      DIMENSION EJMIN(7),EJMAX(7),SION(7),SIONO(7),SIONn(7)
      dimension jres(10),wlres(10),phres(10),fej0(2,ne1:ne2),
     &     fejd(2,ne1:ne2)
      dimension escfe2(10000),tfe2(10000),febin(100),flfebin(100)
      dimension wlphot(nlp),wlpr(1000),citpr(1000),dcitpr(1000)
      dimension dint(14,26,401)
      dimension ilabpr(1000)
      common/lineem/rcgs,wlp(300),totl(300)
      common/colines/wlcor(200),cincor(200),covib1h,covib1h2,covib2h2,
     &     covib1e
      common/abc1/alc(mc),aln(mn),alo(mo),alne(mne),alna(mna),
     &     almg(mmg),alal(mal),alsi(msi),alsu(ms),alar(mar),alca(mca),
     &     alfe(mfe),alni(mni)
      common/initwl/iwl,initprl
      common/fbapp/fbapp6300,fbapp5007
      common/fbappt/fb6300,fb5007
      common/debug/ideb
      common/emind/iemiss
      common/rinner/rmin,rmax

      common/levels/nlevco,nlevsio
      parameter(ncosi=400)
      common/cosioem/wl_co(ncosi,ncosi),wl_sio(ncosi,ncosi),
     &     cin_co(ncosi,ncosi),cin_sio(ncosi,ncosi)

      common/citot_cosio/citot_co(ncosi,ncosi),citot_sio(ncosi,ncosi)
      parameter(nco=279,nsio=400)
      common/comoldata/eco(nco),aco(nco,nco),vco(nco),jco(nco)
     &     ,nmco            
      
      common/kmaxpop/kmaxp(14,27)
      common/colines/wli_co(1200),dcit_co(1200),citot_cot(1200)
c      common/emtest/emlines,emtota
      dimension wlspec(50000),flspec(50000)

      save ik_old,itout,kpr,ilabpr,wlpr,xdist_old


      DATA EJMIN/3.4d0,13.5987d0,24.5d0,54.4d0,200.d0,1.d3,1.d4/
      DATA EJMAX/13.5987d0,24.5d0,54.4d0,200.d0,1.d3,1.d4,1.d6/
      DATA IABU/1,1,1,1,0,1,1,1,0,0,1,0,0,0,6*0/
      DATA KKMAX/5,9,12,0,5,0/
      DATA PI/3.1415926d0/,ELCH/1.60219d-12/,AMU/1.660531d-24/
c     WAVELENGTHS FOR LINES IN ARRAY SPL
C     HE II 304, TWO-G, 1640.,4686, H I TWO-G, LY-C, BA-C, PA-C, BR-C
C     PF-C, FREE-FREE
C
      DATA WLSP/304.d0,0.d0,1640.d0,4686.d0,0.d0,911.5d0,
     &     3646.0d0,8203.6d0,14584.d0,
     &     22788.d0,0.d0,89*0.d0/
      iout61=0
      initlabel=0  
      iout=6

         
      dx=xdist-xdist_old

      if(itime.eq.1) then

         rabscm = rmin

      else

         rabscm = rabscm + dx

      endif
      if(iout61==1) then
         write(61,9487)iu,ik,xdist,xdist_old,dx,dr(ik),rabscm,densi,xeli,tei
 9487    format('ik,xdist,xdist_old,dx ',2i5,1pe15.7,10e12.3)
      endif

      xdist_old=xdist

c     New calc. for each ion

      dcitot=0.
      citott=0.      
      if(iout61==1) then
         write(61,*)'Intensities '
      endif
      
      totn=4.*pi*rabscm**2*dx*densi**2*xeli
c      write(62,*)kradius,rabscm,dx,totn,densi,xeli

      tot_line_int=0.
      do iel=3,14
         do ion=1,nionel(iel)
            ich=0
            do line=1,401
               dint(iel,ion,line)=4.*pi*rabscm**2*densi**2*xeli*
     &              cinx(iel,ion,line)
               if(ik.ge.2) then
                  tot_intens(ik,iel,ion,line)=tot_intens(ik-1,iel,ion,line)+
     &                 dint(iel,ion,line)*FILLING*abs(dx)
                  if(tot_intens(ik,iel,ion,line) > 0.) then
c                     write(6,934)iel,ion,line,wlix(iel,ion,line),dint(iel,ion,line), tot_intens(ik,iel,ion,line)
 934                 format('tot_intens ',3i4,f12.2,1pe12.3,10e12.3)
                  endif
               else
                  tot_intens(ik,iel,ion,line)=dint(iel,ion,line)*FILLING*abs(dx)
               endif
               tot_line_int=tot_line_int + dint(iel,ion,line)*FILLING*abs(dx)
               if(tot_intens(ik,iel,ion,line) > 0.) then
                  ich=1
                  nel=int(ilabcfx(iel,ion,line)/100.)
c                  write(6,*)'elid(nel),nel,iel,ion',elid(nel),nel,iel,ion
c                  write(61,*)'elid(nel),nel,iel,ion',elid(nel),nel,iel,ion

c for output to 62. 75. and to be read like in PWN program
                  if(wlix(iel,ion,line) > 900. .and. line < 20) then
                     tot_line_lum(iel,ion,line)=totn*cinx(iel,ion,line) +
     &                    tot_line_lum(iel,ion,line)
c                     write(62,*)iel,ion,line,wlix(iel,ion,line),totn*cinx(iel,ion,line),
c     &                    tot_intens(ik,iel,ion,line)
                  endif


                  if(iout61==1) then
                     write(61,2928)elid(nel),ion,line,
     &                    wlix(iel,ion,line),cinx(iel,ion,line),
     &                    taulinex(iel,ion,line),dint(iel,ion,line),
     &                    tot_intens(ik,iel,ion,line),tot_intens(ik,iel,ion,line)
 2928                format(a2,i3,i5,f15.3,1pe12.3,10e12.3)
                  endif
               endif
            enddo
            if(iout61==1) then
               if(ich==0) then
                  write(61,*)'ich=0',iel,ion
               elseif(ich==2) then
                  write(61,*)'ich=2',iel,ion
               endif
            endif
         enddo
      enddo
      tot_em=0.
      tot_em_vol=0.
      do j=jmin,jj
         tot_em=tot_em + 4.*pi*em(ik,j) * 4.*pi*rabscm**2*densi**2*abs(dx)*
     &        (E(J)-E(J-1))
         tot_em_vol=tot_em_vol + 4.*pi*em(ik,j) * (E(J)-E(J-1))
      enddo

      tot_cool=4.*pi*rabscm**2*densi**2*xeli*cool*abs(dx)
      tot_heat=4.*pi*rabscm**2*densi**2*heat*abs(dx)
      emlinesi=4.*pi*rabscm**2*densi**2*emlines*abs(dx)      
      emtoti=4.*pi*rabscm**2*densi**2*emtota*abs(dx)    

      if(iout61==1) then
         write(61,*)' '
         write(61,*)'Strongest lines '
      endif

c      call sort(ik,100,0)

c     hbeta=tot_intens(ik,1,1,5)
c!!!
      hbeta=0.
      do iel=3,14
         do ion=1,nionel(iel)

            do line=1,kmaxp(iel,ion)
               if(tot_intens(ik,iel,ion,line) > 2.e-2*hbeta) then
                  nel=int(ilabcfx(iel,ion,line)/100.)
c                  write(6,*)'qq',iel,ion,line,ilabcfx(iel,ion,line),nel
                  if(iout61==1) then
                     write(61,2928)elid(nel),ion,line,
     &                    wlix(iel,ion,line),cinx(iel,ion,line),
     &                    taulinex(iel,ion,line),dint(iel,ion,line),
     &                    tot_intens(ik,iel,ion,line),tot_intens(ik,iel,ion,line)
                  endif
                  if(tot_intens(ik,iel,ion,line) > 5.e-2*hbeta) then
c     write(6,2928)elid(nel),ion,line,
c     &                    wlix(iel,ion,line),cinx(iel,ion,line),
c     &                    taulinex(iel,ion,line),dint(iel,ion,line),
c     &                    tot_intens(ik,iel,ion,line),tot_intens(ik,iel,ion,line)
                  endif
               endif
            enddo
         enddo
      enddo
      ik_old=ik

      return

      return
      end


      subroutine flux_interpol(imax,ishock,r,rfshock,tauh,
     &     tauhi,rfd,isp)
c tauh  = optical depth of hydrogen from shock in current iteration
c tauhi = optical depth of hydrogen from last iteration (total)
      implicit real*8(a-h,o-z)
      include 'param'
      PARAMETER (NFEL=3000)
      parameter (nlp=30000)
      common/trq/rtr(2,md),drtr(2,md),dentr(2,md)
      COMMON/FRE/NINQ,JMIN,JJ
      COMMON/SPECT/TEL,FD(MD,NE1:NE2),F0(NE1:NE2),IPARA
      common/opacity/TA(NE1:NE2),S(MD,NE1:NE2),copac(md,ne1:ne2)
      common/column/coltr(2,md,14,27)
      COMMON/INT/FL(2,NE1:NE2),SI(14,27,NE1:NE2),E1(NE1:NE3),E(NE1:NE3)
      COMMON/COLD/COLTOT(14,27),COLTOTT(14,27),SRED(NFEL+nlp)
      common/taushok/taush(md)
      common/timecheck/time,itime
      common/debug/ideb
      dimension rfd(md),tauhi(md)
      dimension nionel(14)
      data nionel/1,2,6,7,8,10,11,12,13,14,16,18,20,26/
      is = 0
      isp = 0
      idebq=1
      if(ideb.eq.1) write(6,*)' r, rfshock ',r, rfshock,rtr(1,1),
     &     rtr(1,imax)
      write(6,*)'rtr ',rtr(1,1),rtr(1,2),rtr(1,imax-1),
     &     rtr(1,imax),rtr(1,imax+1)
      if(ideb.eq.1) write(0,*)' r, rfshock ',r, rfshock,rtr(1,1),rtr(1,imax)
      if(ideb.eq.1) write(0,*)'rtr ',rtr(1,1),rtr(1,2),rtr(1,imax-1),
     &     rtr(1,imax),rtr(1,imax+1)

c r increases with i

      if(rtr(1,2).lt.rtr(1,imax-1)) then

         if(r.lt.rtr(1,1)) then

            dx=abs(r-rtr(1,1))

            do j=jmin,jj

               fl(2,j)=exp(-dx*dentr(1,1)*ta(j))*fd(1,j)

            enddo

            ii=1
c            write(6,938)ii,r,rtr(1,1),dentr(1,1),dx*dentr(1,1)*ta(3),
c     &           fd(1,3),fl(2,3)

            goto 333

         elseif(r.gt.rtr(1,imax)) then

            dx=abs(r-rtr(1,imax))

            if(idebq.eq.1) write(6,*)'dx ',dx,dentr(1,imax-1),
     &           dentr(1,imax)

            do j=jmin,jj

               flold = fl(2,j)

               fl(2,j)=exp(-dx*dentr(1,imax)*ta(j))*fd(imax,j)

c            write(6,938)j,ta(j),
c     &           dx*dentr(1,imax)*ta(j),fd(imax-1,j),fd(imax,j),
c     &              flold,fl(2,j)
            enddo

            if(idebq.eq.1)write(6,*)' imax ',imax
            if(idebq.eq.1)write(6,938)imax,r,rtr(1,imax),dentr(1,imax),
     &           dx*dentr(1,imax)*ta(3),fd(imax,3),fl(2,3)

            goto 333

         endif

 938     format(' Extrapolation ',i5,1pe15.7,e15.7,10e12.4)

c r decreases with i

      elseif(rtr(1,2).gt.rtr(1,imax-1)) then

         if(r.gt.rtr(1,1)) then
            
            dx=abs(r-rtr(1,1))

            if(idebq.eq.1)write(6,*)'dx ',dx,dentr(1,1)

            do j=jmin,jj

               flold = fl(2,j)

               fl(2,j)=exp(-dx*dentr(1,1)*ta(j))*fd(2,j)
               

               if(j.eq.3.or.j.eq.-200) then 
                  if(idebq.eq.1)write(6,938)j,ta(j),
     &                 dx*dentr(1,1)*ta(j),fd(2,j),fd(3,j),
     &                 flold,fl(2,j)
               endif
            enddo
            
            ii=1
            if(idebq.eq.1)write(6,938)ii,r,rtr(1,1),dentr(1,1),
     &           dx*dentr(1,1)*ta(3),fd(1,3),fl(2,3)

            goto 333

         elseif(r.lt.rtr(1,imax)) then

            dx=abs(r-rtr(1,imax))

            do j=jmin,jj

               fl(2,j)=exp(-dx*dentr(1,imax)*ta(j))*fd(imax,j)

            enddo

c            write(6,938)imax,r,rtr(1,imax),dentr(1,imax),
c     &           dx*dentr(1,imax)*ta(3),fd(imax,3),fl(2,3)

            goto 333

         endif
         
      endif


      if(r.lt.rfshock) then
         
         icase = 1

         do i=1,imax

c note that fd(i,j) is calculated at midpoint between r(i-1) and r(i)

            if(rtr(1,imax).gt.rtr(1,1)) then

               if(r.ge.rtr(1,i).and.r.lt.rtr(1,i+1)) then

                  is = i
                  isp = 1

               endif

            elseif(rtr(1,imax).lt.rtr(1,1)) then

               if(r.ge.rtr(1,i+1).and.r.lt.rtr(1,i)) then

                  is = i
                  isp = 1

               endif

            endif

         enddo

         if(idebq.eq.1)write(6,*)' imax, isp, is, r, rtr(1,is), 
     &        rtr(1,is+1) ',imax,isp,is,r,rtr(1,is),rtr(1,is+1),
     &        rtr(1,imax)

      else

         icase = 2

         if(idebq.eq.1)write(6,928)imax,ishock,tauh,taush(1),
     &        taush(ishock),taush(imax)
 928     format('tauh i inter',2i5,1pe12.4,10e12.4)

         do i=ishock-1,imax

            tauhsh = taush(i)

            tauhshp1 = taush(i+1)

            if(tauh.ge.tauhsh.and.tauh.lt.tauhshp1) then
               is = i
               isp=1
               goto 11
            endif

         enddo
 11      if(idebq.eq.1) write(6,98)is,isp1,tauh,tauhsh,tauhshp1
 98      format(' tauh ',2i5,1pe12.4,10e12.4)

      endif



      if(isp.eq.0) then

         if(idebq.eq.1)write(6,*)' isp = 0 ',r,rtr(1,1),rtr(1,imax)

         if(r.gt.rtr(1,imax)) then

            do j=jmin,jj
               
               fl(2,j) = fd(imax-1,j) 
               
            enddo

         endif

      else

         if(icase.eq.1) then

            do j=jmin,jj

c note that fd(i,j) is midpoint between r(i-1) and r(i)

               fl(2,j) = fd(is+1,j) + (r-rtr(1,is))*
     &              (fd(is+2,j)-fd(is+1,j))/(rtr(1,is+1)-rtr(1,is))



c               if(is.eq.174) then
c                  write(6,926)is,j,r,rtr(1,is),rtr(1,is+1),
c     &                 (r-rtr(1,is)),fd(is+1,j),fd(is+2,j),fl(2,j)
c               endif

               if(itime.le.2) then
                  if(idebq.eq.1) write(6,926)is,j,r,rtr(1,is),
     &                 rtr(1,is+1),
     &                 (r-rtr(1,is)),fd(is+1,j),fd(is+2,j),fl(2,j)
 926              format(' fl i int ',2i5,1pe15.7,10e15.7)
              endif
            enddo

            do iel = 3,14
               do iz = 1,nionel(iel)+1
                  coltot(iel,iz) = coltr(1,is,iel,iz) + 
     &              (r-rtr(1,is))*(coltr(1,is+1,iel,iz)-
     &                 coltr(1,is,iel,iz))/(rtr(1,is+1)-rtr(1,is))
               enddo
            enddo

         elseif(icase.eq.2) then

            do j=jmin,jj

               fl(2,j) = fd(is,j) + (tauh-tauhsh)*
     &              (fd(is+1,j)-fd(is,j))/(tauhshp1-tauhsh)
                  
            enddo

            do iel = 3,14
               do iz = 1,nionel(iel)+1
                  coltot(iel,iz) = coltr(1,is,iel,iz) + 
     &              (tauh-tauhsh)*(coltr(1,is+1,iel,iz)-
     &                 coltr(1,is,iel,iz))/(tauhshp1-tauhsh)
               enddo
            enddo
         endif
      endif
 333  continue
      return
      end

      subroutine radsimp(i,taus,fls)
c
c simple extinction assuming plane geometry from a point with optical 
c    depth taus and mean intensity fls
c
      implicit real*8(a-h,o-z)
      include 'param'
c      PARAMETER (NE1=-100,NE2=250,NE3=NE2+1)
      COMMON/INT/FL(2,NE1:NE2),SI(14,27,NE1:NE2),E1(NE1:NE3),E(NE1:NE3)
c      COMMON/IND/I
      COMMON/FRE/NINQ,JMIN,JJ
      COMMON/DIF/EM(MD,NE1:NE2),TAU(MD,NE1:NE2),TAUTOT(MD,NE1:NE2),
     &     EMC(MD,NE1:NE2)
      common/debug/ideb
      dimension taus(ne1:ne2),fls(ne1:ne2)


      do j=jmin,jj

c optical depth from end point of last iteration to the current point

         taudiff = tautot(i,j) - taus(j)

         fl(2,j) = fls(j)*e2(taudiff)

      enddo

      return
      end




      real*8 function coolf(i,k,ishock,iss,dt,te,dens,xe,
     &     photok,photolm)
      implicit real*8(a-h,o-z)
      include 'PARAM'
      integer stop_tmin
      common/ctime/tcool,trad
      common/mi_te/stop_tmin
      common/lambit/lambda
      common/preion/ipre
      common/ionx/xion(md,14,27)
      common/abl/abn(15)
      COMMON/IND/Ik
      COMMON/PHY/DEN(MD)
      COMMON/DXA/DR(MD)
c      PARAMETER (NE1=-100,NE2=250,NE3=NE2+1)
      COMMON/FRE/NINQ,JMIN,JJ
      COMMON/DIF/EM(MD,NE1:NE2),TAU(MD,NE1:NE2),TAUTOT(MD,NE1:NE2),
     &     EMC(MD,NE1:NE2)
      common/opacity/TA(NE1:NE2),S(MD,NE1:NE2),copac(md,ne1:ne2)
      common/debug/ideb
      COMMON/ELEC/DEL(MD)
      COMMON/HRAT/RHYD,ZEL,XEL,HYR,HEAT,COOL
      common/chem_rec_ph/alchh,alchhe,alchheii,alchc,alchn,alcho,alchne,
     &     alchna,alchmg,alchal,alchsi,alchs,alchar,alchca,alchfe,
     &     phionh,phionhe,phionheii,phionc,phionn,phiono,phionne,
     &     phionna,phionmg,phional,phionsi,phions,phionar,phionca,
     &     phionfe
      common/chem_rec_phii/alchcii,alchnii,alchoii,alchneii,
     &     alchnaii,alchmgii,alchalii,alchsiii,alchsii,alcharii,
     &     alchcaii,alchfeii,
     &     phioncii,phionnii,phionoii,phionneii,
     &     phionnaii,phionmgii,phionalii,phionsiii,phionsii,phionarii,
     &     phioncaii,phionfeii
      parameter(nspec=91,neq=nspec-3)
      common/number/iel,igam,icr
      COMMON/PHQ/ZEA(14,27),GEA(14,27),ZKA(14,27)
      COMMON/INT/FL(2,NE1:NE2),SI(14,27,NE1:NE2),E1(NE1:NE3),E(NE1:NE3)
      common/cophoto/cophr
      common/rec_coll/alrec(30,30),collion(30,30),ction(14,27)
      dimension photok(14,26),photolm(14,26)
      dimension nionel(14)
      data nionel/1,2,6,7,8,10,11,12,13,14,16,18,20,26/

      ideb = 1

      den(ik) = dens

      if((ipre.eq.1.and.lambda.eq.1).or.iss.eq.1) then

         call emiss(te,xe,dens)

         call emiss_cont(te,xe,dens,cooli)

c add continuum to line emissivity (already in cont_emiss.f)

c         do j=jmin,jj
c            em(i,j)=em(i,j) + emc(i,j)
c         enddo

         iqq=4

         call opac_source(iqq,1,dr(k),den(2),te,xe)
c     last stop

c save opacities

         do j=jmin,jj
            copac(k,j) = ta(j)
         enddo

c only solve radiation transport in the post shock direction first time

         icont = 0

         call diffuse_pl(1,icont,k,ishock,te,eddingflux1)
         
         icont = 1

         call diffuse_pl(1,icont,k,ishock,te,eddingflux2)
         jx=100

      endif

      iter = 0

 11   continue

      xeold = xe

      tq=1.e6
      call coion_fit(1,1,tq,cov)
      call rec_colli(i,te,dens,xe)

      call rate(te,xe,photok,photolm)

      call ctionrate(i,te,ction)

      denel=xe*dens

      denh = abn(1)*dens

      dent = dens

      ze = xe -abn(1)*xion(i,1,2) - 
     &     abn(2)*(xion(i,2,2) + 2.*xion(i,2,3))

      if(ze.lt.0.) then

         write(6,9281)xe,(xion(i,iq,1),xion(i,iq,2),xion(i,iq,3),
     &        iq=1,14),ze
         write(6,9281)abn(1),abn(2)
 9281    format(' stop! ',1pe13.5,10e13.5)
c         stop

      endif

      call elecdens(i,xe,zelec)

      call ionization_h_he(i,dt,dent,denel,
     &     photok,xe,zelec)
      
      ze = xe -abn(1)*xion(i,1,2) - 
     &     abn(2)*(xion(i,2,2) + 2.*xion(i,2,3))
      zel=ze
      xel=xe

      call elecdens(i,xe,zelec)
      
      ze = xe -abn(1)*xion(i,1,2) - 
     &     abn(2)*(xion(i,2,2) + 2.*xion(i,2,3))
      zel=ze
      xel=xe
      
      if(ze.lt.-1.e-15) then

         write(6,9281)xe,xion(i,1,2),xion(i,2,2),xion(i,2,3),ze

c!!!  

         ze = 0.

c     stop
         
      endif

      denel=xe*dens


c     call ionization(i,dt,denel,denh,dent,alrec,collion,photok,
c     &        photolm,ction)

      call ionization_aug(i,dt,denel,denh,dent)

      call elecdens(i,xe,zelec)

      ze = xe -abn(1)*xion(i,1,2) - 
     &     abn(2)*(xion(i,2,2) + 2.*xion(i,2,3))
      zel=ze
      xel=xe

      if(ze.lt.0.) then

         write(6,9281)xe,xion(i,1,2),xion(i,2,2),xion(i,2,3),ze
         stop
         
      endif

      deltaxe = abs(xe-xeold)/xe

      do iel=1,14
         do j=1,nionel(iel)+1
c     xion(k,iel,j)=xion(2,iel,j)
         enddo
      enddo
      
      if(ipre.eq.1.and.iter.le.20) then

         iter = iter + 1

         if(deltaxe.lt.1.e-3) goto 12

         goto 11

      endif

 12   call enum_cf(i)
           
      del(ik)=xe
      
      coolf = rad(te,xe,dens,dr(k),ifp)

      trad=3.*1.38e-16*te/(dens*xe*coolf)
      tcool=3.*1.38e-16*te/(dens*xe*cool)

      write(6,9781)i,k,te,xe,dens,coolf,trad,tcool
 9781 format('Cooling time ',2i5,1pe12.3,10e12.3)

      ideb = 0

      return

      end

      subroutine ion_energy(i,chi,nionel,eion)
      implicit real*8(a-h,o-z)
      include 'PARAM'
      common/ionx/xion(md,14,27)
      common/abl/abn(15)
      common/debug/ideb
      dimension chi(30,30)
      dimension nionel(14)
c      data nionel/1,2,6,7,8,10,11,12,13,14,16,18,20,26/
      DATA ELCH/1.60219d-12/
      eion = 0.

      do iel=1,14
         do iz=1,nionel(iel)
            eion = eion + abn(iel)*elch*chi(iel,iz)*
     &           xion(i,iel,iz+1)
         enddo
      enddo

c!!! eion=0

      eion=0

      return

      end

      subroutine tempcalc(i,iss,d00,p0,eb0,ni0,u0,denfl,amean,
     &     xe,chi,h,xc,nion,ne,u,p,te)
c calculate temperature and densities from new enthalphy h
      implicit real*8(a-h,o-z)
      real*8 nion,ne,ni0
      COMMON/THER/A1,B1,TIN,E10,E20
      dimension chi(30,30),nionelz(30)
      dimension nionel(14)
      data nionel/1,2,6,7,8,10,11,12,13,14,16,18,20,26/

      data amu/1.660531d-24/

      call ion_energy(i,chi,nionel,eion)
      
      if(iss.eq.1) then
         
         h = 0.

         a1 = 0.95 * te 
         
         b1 = 1.05 * te

c     try this! 170924

         a1 = 0.90 * te 
         
         b1 = 1.2 * te

c e10 not used!

         e10 =  0.

         e20 =  1.e-3

         tmin = 0.

         call bis(te,xe,tmin,ifail,ifpop,iconv,istop)

         write(0,92)te,xe,tmin,ifail,ifpop,iconv,istop
         write(6,92)te,xe,tmin,ifail,ifpop,iconv,istop
 92      format(' from bis ',1pe12.3,10e12.3)

      endif

      hprim = h - eion

      call compress(p0,eb0,ni0,u0,denfl,hprim,xc)

      dencgs = xc*d00
      
      u = denfl/dencgs
         
      p = p0-dencgs*u**2-eb0*xc**2
      
      nion=dencgs/(amu*amean)
      ne=nion*xe

      if(iss.eq.0) then

         te=p/(1.38e-16*(ne+nion))

      endif

      return
      end

      subroutine biseq(i,k,ishock,iss,dt,dens,xe,
     &     photok,photolm,tlim_min,tlim_max,
     &     d00,p0,eb0,ni0,u0,denfl,amean,h,
     &     chi,xc,nion,ne,u,p,te)

      implicit real*8(a-h,o-z)
      real*8 nion,ne,ni0

      COMMON/HRAT/RHYD,ZEL,XEL,HYR,HEAT,COOL

      dimension photok(14,26),photolm(14,26)
      common/rec_coll/alrec(30,30),collion(30,30),ction(14,27)

      dimension chi(30,30)
      dimension nionel(14)
      data nionel/1,2,6,7,8,10,11,12,13,14,16,18,20,26/

      data amu/1.660531d-24/

      write(6,93)i,k,ishock,dens,xe,
     &     tlim_min,tlim_max,
     &     d00,p0,eb0,ni0,u0,denfl,amean,h,
     &     xc,nion,ne,u,p,te

 93   format(' in biseq ',3i5,1pe12.3,20e12.3)

      call ion_energy(i,chi,nionel,eion)

c first coolf at upper and lower te-limit

      te1 = tlim_min

      mit=10
      do m=1,mit
         cool1 = coolf(i,k,ishock,iss,dt,te1,dens,xe,
     &        photok,photolm)
      enddo

      te2 = tlim_max

      do m=1,mit
         cool2 = coolf(i,k,ishock,iss,dt,te2,dens,xe,
     &        photok,photolm)
      enddo
      
      write(6,94)te1,te2,cool1,cool2
 94   format('irst bracket te1,te2,cool1,cool2 ',1pe12.3,10e12.3)

      if(cool1*cool2.gt.0.) then

         write(6,*)' same sign for pre-temp',
     &        te1,te2,cool1,cool2

c find zero point by starting at high temp and gradually decreasing it 
c     until sign shifts                  

         te2 = 1.2* tlim_max

         do m=1,mit
            cool2 = coolf(i,k,ishock,iss,dt,te2,dens,xe,
     &           photok,photolm)
         enddo
                  
         te1 = te2

         do kk=1,100

            te1 = te1/1.1

            do m=1,mit
               cool1 = coolf(i,k,ishock,iss,dt,te1,dens,xe,
     &              photok,photolm)
               write(6,*)'scan of cool ',te1,cool1,cool2
            enddo

            if(cool1*cool2.lt.0.) then
               write(6,*)' shift of sign for net cooling',te1,te2,cool1,cool2
               goto 322

            endif

            te2 = te1

         enddo

         write(6,*)' could not find zero point of rad!!!'
         
         stop

c     322     te2 = 1.01*te2
 322     continue
         do m=1,mit
            cool2 = coolf(i,k,ishock,iss,dt,te2,dens,xe,
     &           photok,photolm)
         enddo
         
         te1 = 0.99*te1
         do m=1,mit
            cool1 = coolf(i,k,ishock,iss,dt,te1,dens,xe,
     &           photok,photolm)
         enddo
         tem = (te1 + te2)/2.

         write(6,*)' aft scan ',te1,te2,tem,cool1,cool2

      endif

      tem = (te1 + te2)/2.

 444  continue

c now coolf at midpoint
      do m=1,mit
         coolm = coolf(i,k,ishock,iss,dt,tem,dens,xe,
     &        photok,photolm)               
      enddo
      
      if(cool1*coolm.lt.0.) then

         te2 = tem

         tec = tem

         tem = (te1 + tem)/2.

      elseif(cool2*coolm.lt.0.) then

         te1 = tem

         tec = tem

         tem = (te2 + tem)/2.

      else

         write(6,*)' same sign for pre-temp',
     &        te1,te2,tem,cool1,cool2,coolm

      endif

 9183 format('te1,te2,cool1,cool2 ',i5,1pe12.3,10e12.3)

      tediff = abs(te1-te2)/te2

      epste = 1.e-3

      epscool=1.e-3

      coolbal = abs(coolm/cool)

      write(6,9488)te1,te2,tem,cool,coolm
 9488 format('te1,te2,tem,cool,coolm ',1pe12.3,10e12.3)

      if(tediff.gt.epste.and.coolbal.gt.epscool) then

         goto 444   

      else

         if(tediff.lt.epste) then

            te = tem
         
         elseif(coolbal.lt.epscool) then

            te = tec         

         endif

c cool = heat -> h = const


         hprim = h - eion

         write(6,9378)h,eion,hprim
 9378    format('h,eion,hprim ',1pe12.3,10e12.3)

         call compress(p0,eb0,ni0,u0,denfl,hprim,xc)
            
         dencgs = xc*d00
      
         u = denfl/dencgs
         
         p = p0-dencgs*u**2-eb0*xc**2
      
         nion=dencgs/(amu*amean)

         ne=nion*xe

      endif

      write(6,9)i,d00,p0,eb0,ni0,u0,denfl,amean,
     &     xe,h,nion,ne,xc,te
 9    format(' out of biseq ',i5,1pe13.5,20e13.5)

      return

      end

      
      subroutine sort(ik,nmax,istr)
      implicit real*8 (a-h,o-z)
      include 'param'
c      common/kmax/kmax(14,27)
      common/kmaxpop/kmaxp(14,27)
c     number of ionization stages on ordered scale 1-14      
      integer nionstage,done(14,27,500),hit    
c      real*8 tot_line_lum
c      common/linelum/tot_line_lum(14,27,500)
      common/tot_line_em/tot_intens(md,14,26,401)
      common/line_em/cinx(14,26,401),taulinex(14,26,401),wlix(14,26,401)
     &     ,ilabcfx(14,26,401)
      dimension nionstage(14)
      data nionstage/1,2,6,7,8,10,11,12,13,14,16,18,20,26/
      data pi/3.14159/
      character*2 elid(14)
      data elid/'H ','He','C ','N ','O ','Ne','Na','Mg','Al','Si','S ','Ar','Ca','Fe'/
      character*5 ionid(26)
      data ionid/'I   ','II   ','III  ','IV   ','V    ','VI   ',
     &     'VII  ','VIII ','IX   ','X    ','XI   ','XII  ','XIII ',
     &     'XIV  ','XV   ','XVI  ','XVII ','XVIII','XIX  ','XX   ',
     &     'XXI  ','XXII ','XXIII','XXIV ','XXV  ','XXVI '/

      write(6,*)' Strongest lines at ',ik
      do iel=1,14
         do ion=1,nionstage(iel)
            do k=1,500
               done(iel,ion,k)=0
            enddo
         enddo
      enddo
      CIMAX=1.D100
      kx=1
      DO KL=1,nmax
         CMAX=0.
         hit=0
         do iel=1,14
            do ion=1,nionstage(iel)
c               write(6,*)'kl,iel,ion ',kl,iel,ion,nionstage(iel),kmaxp(iel,ion)
               do k=1,kmaxp(iel,ion)
c     if(iel.eq.11)
cc                   write(6,92)iel,ion,k,done(iel,ion,k),wlix(iel,ion,k),tot_intens(ik,iel,ion,k),cmax
 92               format('Str ',4i5,f10.1,1pe12.3,10e12.3)
                  if(tot_intens(ik,iel,ion,k) > cmax .and. done(iel,ion,k).eq.0) then
                    cmax=tot_intens(ik,iel,ion,k)
                    ielx=iel
                    ionx=ion
                    kx=k
c                    write(6,9)kl,elid(ielx),ionid(ionx),wlix(ielx,ionx,kx),cmax
                    hit=1
                  endif
                enddo
            enddo
         enddo
         if(kx > 500) then
            write(0,*)' kx > 500 ',ielx,ionx,kx
            kx=1
         endif
c         write(6,*)'ielx,ionx,kx ',ielx,ionx,kx
c     if(ielx.ne.1.and.ionx.ne.1.and.kx.ne.1) then
         if(hit==1) then
            done(ielx,ionx,kx)=1
c            write(6,*)'done ',ielx,ionx,kx
            if(wlix(ielx,ionx,kx) > -900.) then
               write(6,97)elid(ielx),ionid(ionx),wlix(ielx,ionx,kx),cmax
               if(istr==1) then
                  write(101,97)elid(ielx),ionid(ionx),wlix(ielx,ionx,kx),cmax
 9                format('strong ',i5, 2a6,f12.2,1pe12.3)
               endif
 97            format(2a6,f12.2,1pe12.3)
            endif
         endif            
      enddo
c      stop
      return
      end



      
      subroutine renum(k,nshock)
c      include 'crenum.inc'
c      implicit none
      implicit real*8 (a-h,o-z)
      integer*4 k,nshock,nionel
      include 'parameters.h'
      COMMON/DIF/EM(MD,NE1:NE2),TAU(MD,NE1:NE2),TAUTOT(MD,NE1:NE2),
     &     EMC(MD,NE1:NE2)
      common/opacity/TA(NE1:NE2),S(MD,NE1:NE2),copac(md,ne1:ne2)
      common/trq/rtr(2,md),drtr(2,md),dentr(2,md)
      common/column/coltr(2,md,14,27)
      COMMON/RADIE/R(0:MD)
      COMMON/FRE/NINQ,JMIN,JJ
      dimension nionel(14)
      data nionel/1,2,6,7,8,10,11,12,13,14,16,18,20,26/
      dimension rs(md),drs(md),densave(md),coltrs(md,14,27)
      dimension emtr(md,ne1:ne2),emtrc(md,ne1:ne2)
      dimension tautots(md,ne1:ne2),sources(md,ne1:ne2),ems(md,ne1:ne2),
     &     copacs(md,ne1:ne2),flaan(md,ne1:ne2),flaanc(md,ne1:ne2),
     &     emcs(md,ne1:ne2)
      write(6,*)' before renum ',k
      write(6,*)'jmin,jj ',jmin,jj

c!!   Seg problem with this section
      do ij=1,nshock
         do j=jmin,jj
            tautots(ij,j) = tautot(ij,j)
            sources(ij,j) = s(ij,j)
            ems(ij,j) = em(ij,j)
            emcs(ij,j) = emc(ij,j)
            if(ij.eq.2) then                  
cxx   ems(ij,j) = emtr(ij,j)
cxx   emcs(ij,j) = emtrc(ij,j)
            endif
            copacs(ij,j) = copac(ij,j)
         enddo
         rs(ij) = rtr(2,ij)
         drs(ij) = drtr(2,ij)
         densave(ij) = dentr(2,ij)

         do iel=1,14
            do iz=1,nionel(iel)+1
               coltrs(ij,iel,iz) = coltr(2,ij,iel,iz)
            enddo
         enddo

         if(ij.lt.-5.or.ij.gt.nshock-5) then
cxx   write(6,9332)ij,rtr(2,ij),drtr(2,ij),dentr(2,ij),
cxx   &              tautot(ij,3),em(ij,3),emtr(ij,3),copac(ij,3),s(ij,3)
cxx   &              ,coltr(2,ij,1,1)
 9332       format('ij,r,tau,em,s,N ',i5,1pe15.7,10e12.3)
         endif
      enddo

c     now renumber with i=1 at post shock boundary and i=nshock at shock


      write(0,*)'nshock',nshock
      do ij=1,nshock
         do j=jmin,jj
            tautot(ij,j) = tautots(nshock-ij+1,j)
            s(ij,j) = sources(nshock-ij+1,j)
            em(ij,j) = ems(nshock-ij+1,j)
            emc(ij,j) = emcs(nshock-ij+1,j)
            emtr(ij,j) = em(ij,j)
            emtrc(ij,j) = emc(ij,j)
            copac(ij,j) = copacs(nshock-ij+1,j)
         enddo
         rtr(2,ij) = rs(nshock-ij+1)
         drtr(2,ij) = drs(nshock-ij+1)
c!!!  
         drtr(2,ij) = drs(nshock-ij+2)
         dentr(2,ij) = densave(nshock-ij+1)

         do iel=1,14
            do iz=1,nionel(iel)+1
               coltr(2,ij,iel,iz) = coltrs(nshock-ij+1,iel,iz)
            enddo
         enddo

         jx=-28

         if(ij.lt.5.or.ij.gt.nshock-5) then
            write(6,9332)ij,rtr(2,ij),drtr(2,ij),dentr(2,ij),
     &              tautot(ij,jx),em(ij,jx),copac(ij,jx),s(ij,jx),
     &              coltr(2,ij,1,1)
         endif
      enddo

      drtr(2,1) = abs(rtr(2,2)-rtr(2,1))
      dentr(2,1) = dentr(2,2)

      write(6,*)' rtr(1) ',k,rtr(2,1),drtr(2,1),dentr(2,1)

c     save i = 1 & 2 for radiative transfer calc. and interpolation 
c     since these are overwritten 

c     calculate optical depth from post shock and upstreams

      do j= jmin,jj
         do i=2,nshock
            tautot(i,j)=tautot(1,j) - tautot(i,j)
         enddo
         tautot(1,j) = 0.
      enddo
      
      write(6,*)' after optical depth inversion'
      do ij=1,nshock
c     if(i.lt.5.or.i.gt.nshock-5) then
         write(6,9232)ij,r(ij),rtr(2,ij),drtr(2,ij),
     &           dentr(2,ij),tautot(ij,jx),em(ij,jx),copac(ij,jx),
     &           s(ij,jx),coltr(2,ij,1,1)
 9232    format('i,r,rtr,drtr,dentr,tautot,em,copac,s,coltr ',i5,1pe15.7,10e12.3)
c     endif
      enddo

      imax = k

      if(ideb.eq.1) write(0,*)' k,imax',k,imax

      rstart = rtr(1,imax-1)
      nshock = k-1

c     for forward shock (FS)= and csm
      if(icsm==1.and.ipre==0) then
         write(6,*)'after shift in radial direction for a FS ',
     &        rshock_cgs/1.e15,' 1e15 cm'
         do ij=1,nshock
            rtr(2,ij)=rshock_cgs-rtr(2,ij)
            r(ij)=rtr(2,ij)
            drtr(2,ij)=drtr(2,ij)
            write(6,9332)ij,rtr(2,ij),r(ij),drtr(2,ij),dentr(2,ij),
     &              tautot(ij,jx),em(ij,jx),copac(ij,jx),s(ij,jx),
     &              coltr(2,ij,1,1)
         enddo
      endif

      do j=jmin,jj
c     em(i,j) = emtr(i,j)
c     emc(i,j) = emtrc(i,j)
      enddo

      return
      end
