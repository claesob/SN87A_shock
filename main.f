
C***********************************************************************
C**********
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 MTOT,MNI,MSUN
      character*72 text1,text2,text3,TEXT4,
     &     INMOD,infil,MODNRSTR
      character*9 
     &     FILE5,FILE6
      character*6 FILE1,FILE2,FILE4
      character*7 FILE3
      CHARACTER*25 LABEL
      character*11 file57,file67,file68,file69,file75,file77,
     &     file78,file96,file56,file59,file81,file62
      character*20 file101
      character*21 file63
c     include 'param'
      include 'parameters.h'
c      PARAMETER (NE1=-100,NE2=250,NE3=NE2+1)
      COMMON/MOD/INMOD
      common/isupp/isupp
      COMMON/TEXT/TEXT1,TEXT2,TEXT3,TEXT4,TEXT6
      COMMON/ITER/ITE
      COMMON/A5/TAUE,ALFA2,EN,TSN,XL40,TEXA,RMAX
      COMMON/PULSAR/RLTOT,ALFA,EMIN,EMAX
      common/obstar/rltotob,teffob
      COMMON/DENS/DEN0,R0,RN
      COMMON/MPAR/MTOT,MNI,VEXP
      COMMON/RQW/TEFF,RQ
      COMMON/INUT/IUY
      COMMON/ELDEN/IDENS
      COMMON/NBA/NBACK
      COMMON/TPAR/RIN,DRQ,R1Q,TDAY
      COMMON/FRE/NINQ,JMIN,JJ
      COMMON/INT/FL(2,NE1:NE2),SI(14,27,NE1:NE2),E1(NE1:NE3),E(NE1:NE3)
      COMMON/AUG/AUG
      COMMON/QSOM/QSO,IAGN,ISTAT,ISOBOL
      COMMON/INSPEC/IBLACKB,IPULSSP,IFERNET,INLR,IFERMAT,ICRAB,
     &     IFREEFREE,imewe,iobstar
      common/initmw/initmw,initr2
      COMMON/AGNPARA/DECON,GAMMION,IPRESS,IDEN,ICDEN
      COMMON/MULTI/IHMUL,IHEMUL,IOMUL,ICAMUL,IFEMUL
      COMMON/UVBL/IUVBLANK
      COMMON/SPH/ISPH
      COMMON/GAMPAR/GAMLUM,GAELH,GAHE
      COMMON/HYDROS/HSCALE,DEN1,XEL,AMEAN
      COMMON/CHG/CHGH,CHGHE,CHTE
      COMMON/THER/A1,B1,TIN,E10,E20
      COMMON/ITSPH/NITSPH
      COMMON/A12/ITER
      COMMON/A10/AV
      COMMON/IOPAR/IOSHELL,IORAD1,IORAD2,IOFLUX
      COMMON/CICS/RSHOCK,ics,ics87a
      common/revs/terev,tecs
      COMMON/CINOUT/INOUT,IPULS
      COMMON/LOWION/ILOWION
      COMMON/FELEV/NFEII,nfei
      COMMON/IOINF/IOSPEC,IOLEVEL,IOEMISS,IOTERM,IOION,IOELITER
      COMMON/CHEMCOMP/ABIN(12),ISOLAR,ISNCOMP,INPUTCOMP
      COMMON/THBAL/TMIN,TMAX,IBAL,IKBAL(1000)
      common/fill/fillingq
      integer nz,nion,nshell,i,k
      parameter(nz=30,nion=27,nshell=10)
      integer ns,kmax
      real*8 fr_aug,eion,en_aug,en_augi,eioni
      common/augerfrac/eion(nz,nion,nshell),en_aug(nz,nion,nshell),
     &     fr_aug(nz,nion,nshell,10),kmax(nz,nion,nshell),ns(nz,nion),
     &     init_augfrac
      common/inputv/RMCEN,DELM,R1,T00,A1IN,B1IN,dminit,rminner,revmass,
     &     EJDENS,TECOOL,totcolumn,FILLING,MODNR,NMAX,MAX,ICEN,IREV,
     &     NPRINT,IPAR,MQMAX,N,IDEP,IRINNER,ISTEPION

      DATA PI/3.1415926E0/,ELCH/1.60219E-12/,AMU/1.660531E-24/
      DATA SOLMA/1.989E33/
c     dielectronic supression on = 1 
      isupp=1
      INOUT=2
      IPULS=1
C
C     INITIALIZE ALL PARAMETERS
C
C
C	OPEN FILES
C
      open(27,file='./ATDAT/feIV_22.dat',status='old')
      open(28,file='./ATDAT/feIII_30.dat',status='old')
C     INPUT PARAMETERS IN THIS FILE
      OPEN(99,FILE='shock_input.txt',access='sequential',form='formatted')
      READ(99,9899)TEXT1
      READ(99,9899)TEXT2
      READ(99,9899)TEXT4
      READ(99,9899)TEXT3
      READ(99,9899)TEXT5      
      READ(99,*)MODNR
      MN1=IFIX(real(modnr)/100)
      MN2=MODNR-100*MN1
      MN2=IFIX(real(MN2)/10)
      MN3=MODNR-MN1*100-MN2*10
      MN3=IFIX(real(MN3))
      MODNRSTR(1:1)=CHAR(MN1+48)
      MODNRSTR(2:2)=CHAR(MN2+48)
      MODNRSTR(3:3)=CHAR(MN3+48)

c      file81(1:8)='fort_81_'
c      file81(9:11)=modnrstr(1:3)
c      open(81,file=file81)
c      file56(1:8)='fort_56_'
c      file56(9:11)=modnrstr(1:3)
c      open(56,file=file56)
      file101(1:13)='strong_lines_'
      file101(14:16)=modnrstr(1:3)
      file101(17:20)='.txt'
      open(101,file=file101)
c      file59(1:8)='fort_59_'
c      file59(1:8)='spectrum'      
c      file59(9:11)=modnrstr(1:3)
c      open(59,file=file59)
      FILE63(1:5)='struc'
      FILE63(6:8)=MODNRSTR(1:3)
      OPEN(63,FILE=FILE63)
   
c      file75(1:8)='fort_75_'
c      file75(9:11)=modnrstr(1:3)
c      open(75,file=file75)
c      file77(1:8)='fort_77_'
c      file77(9:11)=modnrstr(1:3)
c      open(77,file=file77)
c      open(77,file='fort_77_390')
c      file62(1:8)='line_em_'
c      file62(9:11)=modnrstr(1:3)
c      open(62,file=file62)
c      file67(1:8)='fort_67_'
c      file67(9:11)=modnrstr(1:3)
c      open(67,file=file67)
c      open(67,file='fort_67_390')
c      file68(1:8)='fort_68_'
c      file68(9:11)=modnrstr(1:3)
c      open(68,file=file68)
c      open(68,file='fort_68_390')
c      file69(1:8)='fort_69_'
c      file69(9:11)=modnrstr(1:3)
c      open(69,file=file69)
c      open(69,file='fort_69_390')
c      file57(1:8)='fort_57_'
c      file57(9:11)=modnrstr(1:3)
c      open(57,file=file57)
c      open(57,file='fort_57_390')
c      file78(1:8)='fort_78_'
c      file78(9:11)=modnrstr(1:3)
c      open(78,file=file78)
c      file96(1:8)='fort_96_'
c      file96(9:11)=modnrstr(1:3)
c      open(96,file=file96)


c      FILE5(1:6)='slltot'
c      FILE5(7:9)=MODNRSTR(1:3)
c      OPEN(33,FILE=FILE5)
c      FILE6(1:6)='sltab.'
c      FILE6(7:9)=MODNRSTR(1:3)
c      OPEN(37,FILE=FILE6)
9899  FORMAT(A)
C     SAVE FLUXES ETC. IN THIS
c      OPEN(14,FILE='slutin2',FORM='UNFORMATTED')
C     FILE1= GENERAL OUTPUT FILE = UNIT 6
      FILE1(1:3)='slt'
      FILE1(4:6)=MODNRSTR(1:3)
C     FILE2= LINE PROFILES = UNIT 16
c      FILE2(1:3)='sll'
c      FILE2(4:6)=MODNRSTR(1:3)
C     FILE3= EMISSION RATES IN THIS FILE = UNIT 17
c      FILE3(1:4)='slem'
c      FILE3(5:7)=MODNRSTR(1:3)
C     FILE4= FINAL MODEL = UNIT 15
c      FILE4(1:3)='slf'
c      FILE4(4:6)=MODNRSTR(1:3)
      OPEN(6,FILE=FILE1)
c      OPEN(15,FILE=FILE4)
c      OPEN(16,FILE=FILE2)
c      OPEN(17,FILE=FILE3)
c      OPEN(22,FILE='slasktx2')
c      OPEN(21,FILE='slrecem2')
c      OPEN(18,FILE='slabund2')
c      OPEN(31,FILE='sleddflux')
c$$$      WRITE(6,9899)TEXT1
c$$$      WRITE(6,9899)TEXT2
c$$$      WRITE(6,9899)TEXT4      
c$$$      WRITE(6,9899)TEXT3
c$$$      WRITE(6,9899)TEXT5      
c$$$      WRITE(6,*)MODNR
      JMIN=-249
c      jmin=-499
c      jmin=-699
      JJ=247
      NINQ=2
      CALL ENINT
      init_augfrac=1

      call auger_fr
      call atdatar5
c     Auger data      
      call readferec
c data for recomb. of Si, S and Ar      
      call recomb_adas
c call recomb data from Badnell et al for Na, Mg, Al, Si, P sequencies      
      call badnell_et_al
C     DEFAULT PARAMETERS 
      MAX=459
      ICEN=35
      MNI=1.E-20
      A1IN=3000.
      B1IN=10000.
c  cf:s old
      NFEII=116
c  johns Fe II atom
      NFEII=121
c  cf new
      nfeii=191
      nfei=121
      TDAY=1.E20
      TOTCOLUMN=1.E33
      RIN=1.E14/1.E15
      FILLING=1.
      CGGH=0.1
      CHTE=0.2
      ISTEPION=1
      VEXP=1.00E9
      IOSHELL=0
      IORAD1=0
      IORAD2=0
      IOFLUX=0
      IOSPEC=0
      iobstar=0
      IOLEVEL=0
      IOTERM=0
      IOION=0
      IOELITER=0
      IDEP=0
      IHMUL=1
      IHEMUL=1
      IOMUL=1
      ICAMUL=1
      IFEMUL=1
C     IF IDENS < 0 USE A CONSTANT ELECTRON FRACTION IN POPULATION CALC.
      IDENS=1
C     IF NBACK = 0 PUT BACKGROUND FLUX IN MULTI LEVEL CALC. TO ZERO.
      NBACK=1
C     IF AUG>0 INCLUDE AUGER IONIZATIONS
      AUG=999.
C     QSO=INDEX SPECIFYING IF THE CALCULATION APPLIES TO A QUASAR
C         OR NOT. QSO>0 : NO BACKGROUND RADIATION FIELD, QSO<0 :
       QSO=-999.
      DO K=1,100
        READ(99,9)LABEL
9       FORMAT(A21)
        write(*,*)'label ',label

        IF(LABEL(1:7).EQ.'INMODEL') THEN
          READ(99,9899)INMOD
        ELSEIF(LABEL(1:5).EQ.'INOUT') THEN
C         IF INOUT = 1 START FROM CENTER, 0 from outer boundary
          read(99,*)INOUT
        ELSEIF(LABEL(1:7).EQ.'TEMPLIM') THEN
          READ(99,*)A1IN,B1IN
        ELSEIF(LABEL(1:7).EQ.'BALANCE') THEN
C         SCAN THE HEATING AND COOLING CURVES FROM TMIN TO TMAX FOR
C         ALL I = IBAL(K)
          READ(99,*)TMIN,TMAX
          IBAL=1
          DO KQ=1,1000
            READ(99,*)IKBAL(KQ)
            IF(IKBAL(KQ).EQ.0) GOTO 569
          ENDDO
569       CONTINUE
        ELSEIF(LABEL(1:4).EQ.'MTOT') THEN
C         MTOT = TOTAL MASS IN G
          READ(99,*)MTOT
        ELSEIF(LABEL(1:6).EQ.'COLTOT') THEN
C         TOTCOLUMN = TOTAL COLUMN DENSITY IN CM-2
          READ(99,*)TOTCOLUMN
        ELSEIF(LABEL(1:5).EQ.'RIN15') THEN
C         RIN= RADIUS OF INNER BOUNDARY IN 1E15 CM
          READ(99,*)RIN
          IRINNER=1
        ELSEIF(LABEL(1:7).EQ.'FILLING') THEN
C         FILLING FACTOR
          READ(99,*)FILLING
        ELSEIF(LABEL(1:4).EQ.'TDAY') THEN
C         TDAY=TIME SINCE EXPLOSION IN DAYS
          READ(99,*)TDAY
        ELSEIF(LABEL(1:9).EQ.'DIFF ITER') THEN
C         NMAX=NUMBER OF ITERATIONS FOR RADIATION FIELD. EG. NMAX=1 MEANS
C          NO ITERATION , NMAX=2 , ONE ITERATION ETC.
          READ(99,*)NMAX
        ELSEIF(LABEL(1:6).EQ.'SHELLS') THEN
C         MAX=TOTAL NUMBER OF RADIAL SHELLS AND NUMBER OF SHELL IN CENTRE
          READ(99,*)MAX
        ELSEIF(LABEL(1:7).EQ.'CSHELLS') THEN
C         ICEN=NUMBER OF SHELLS IN CENTRE
          READ(99,*)ICEN
        ELSEIF(LABEL(1:5).EQ.'RMCEN') THEN
C         MAX=TOTAL NUMBER OF RADIAL SHELLS AND NUMBER OF SHELL IN CENTRE
          READ(99,*)RMCEN
        ELSEIF(LABEL(1:6).EQ.'DMINIT') THEN
C         DMINIT = FRACTION OF STROEMGREN THICKNESS IN FIRST MASS-STEP      
          READ(99,*)dminit
        ELSEIF(LABEL(1:7).EQ.'RMINNER') THEN
          READ(99,*)rminner
        ELSEIF(LABEL(1:9).EQ.'UNIF STEP') THEN
          ISTEPION=0
        ELSEIF(LABEL(1:9).EQ.'CH LOGOPT') THEN
C         CHGH=MAXIMUM LOGARITMIC CHANGE IN OPTICAL DEPTH AT J=2
          READ(99,*)CHGH
        ELSEIF(LABEL(1:7).EQ.'CH LOGT') THEN
C         CHTE=RATIO BETWEEN ALLOWED CHANGE IN THE LOG(TEMPERATURE) TO THE
          READ(99,*)CHTE
        ELSEIF(LABEL(1:11).EQ.'ENERGY BINS') THEN
C         JMIN=INDEX OF MINIMUM ENERGY BIN
C         JJ= D:O OF MAXIMUM BIN
          READ(99,*)ENMIN,ENMAX
        ELSEIF(LABEL(1:5).EQ.'FE II') THEN
C         NUMBER OF FE II LEVELS
          READ(99,*)NFEII
        ELSEIF(LABEL(1:5).EQ.'FE  I') THEN
C         NUMBER OF FE I LEVELS
          READ(99,*)NFEI
        ELSEIF(LABEL(1:19).EQ.'ITERATIVE MULTI H I') THEN
C         SIMPLE MULTILEVEL CALCN. FOR H I
          IHMUL=0
        ELSEIF(LABEL(1:20).EQ.'ITERATIVE MULTI HE I') THEN
C         SIMPLE MULTILEVEL CALCN. FOR HE I
          IHEMUL=0
        ELSEIF(LABEL(1:19).EQ.'ITERATIVE MULTI O I') THEN
C         SIMPLE MULTILEVEL CALCN. FOR O I
          IOMUL=0
        ELSEIF(LABEL(1:21).EQ.'ITERATIVE MULTI CA II') THEN
C         SIMPLE MULTILEVEL CALCN. FOR CA II
          ICAMUL=0
        ELSEIF(LABEL(1:21).EQ.'ITERATIVE MULTI FE II') THEN
C         SIMPLE MULTILEVEL CALCN. FOR FE II
          IFEMUL=0
        ELSEIF(LABEL(1:18).EQ.'NO MULTI FOR FE II') THEN
C         SIMPLE MULTILEVEL CALCN. FOR FE II
          IFEMUL=-1
        ELSEIF(LABEL(1:5).EQ.'DECON') THEN
C         AGN DENSITY AT INNER BOUNDARY
          READ(99,*)DECON
          ICDEN=1
        ELSEIF(LABEL(1:6).EQ.'IONPAR') THEN
C         AGN IONIZATION PARAMETER
          READ(99,*)GAMMION
        ELSEIF(LABEL(1:9).EQ.'BLACKBODY') THEN
C         BLACK-BODY SPECTRUM
          READ(99,*)TEFFBB
          IBLACKB=1
        ELSEIF(LABEL(1:19).EQ.'PHOTOSPHERIC RADIUS') THEN
C         PHOTOSPHERIC RADIUS IN CM
          READ(99,*)RPH
          IPH=1
        ELSEIF(LABEL(1:11).EQ.'PULSAR SPEC') THEN
C         PULSAR SPECTRUM
          READ(99,*)IPULSSP
        ELSEIF(LABEL(1:9).EQ.'FREE FREE') THEN
C         Free-Free SPECTRUM
          IFREEFREE=1
          read(99,*)imewe
          initmw=0
          initr2=0
        ELSEIF(LABEL(1:4).EQ.'CRAB') THEN
C         CRAB SPECTRUM
          ICRAB = 1
        ELSEIF(LABEL(1:14).EQ.'FERLAND NETZER') THEN
C         FERLAND & NETZER AGN SPECTRUM
          READ(99,*)IFERNET
        ELSEIF(LABEL(1:15).EQ.'FERLAND MATHEWS') THEN
C         FERLAND & MATHEWS AGN SPECTRUM
          READ(99,*)IFERMAT
        ELSEIF(LABEL(1:3).EQ.'NLR') THEN
C         NLR AGN SPECTRUM
          READ(99,*)INLR
        ELSEIF(LABEL(1:3).EQ.'NIC') THEN
C         56 NI MASS
          READ(99,*)MNI
        ELSEIF(LABEL(1:11).EQ.'COMPOSITION') THEN
C         SPECIFY CHEMICAL COMPOSITION
           READ(99,9)LABEL
           write(6,*)'label ',label
          IF(LABEL(1:5).EQ.'SOLAR') THEN
C           SOLAR COMPOSITION
            ISOLAR=1
          ELSEIF(LABEL(1:6).EQ.'SNCOMP') THEN
C           SUPERNOVA COMPOSITION
            ISNCOMP=1
          ELSEIF(LABEL(1:5).EQ.'INPUT') THEN
C           SPECIFY CHEMICAL COMPOSITION
            INPUTCOMP=1
            READ(99,*)(ABIN(L),L=1,12)      
          ELSEIF(LABEL(1:5).EQ.'INP19') THEN
C           SPECIFY CHEMICAL COMPOSITION
            INPUTCOMP=1
c     ring 87A
            open(19,file='./EXTRAS/shockspec_ring.dat',status='old')
            READ(19,*)(ABIN(L),L=1,12)      
            write(*,*)' ab ',abin(1),abin(2),abin(3),abin(5)
          ENDIF
        ELSEIF(LABEL(1:3).EQ.'VEL') THEN
C         EXPANSION VELOCITY AT THE OUTER BOUNDARY IN CM/S
          READ(99,*)VEXP
        ELSEIF(LABEL(1:7).EQ.'NO BACK') THEN
C         NO BACKGROUND CONT. IN LINE EXC.
          NBACK=0
        ELSEIF(LABEL(1:6).EQ.'OUTPUT') THEN
C         READ PARAMETERS CONTROLLING THE OUTPUT
          READ(99,*)IOSHELL,IORAD1,IORAD2,IOFLUX
        ELSEIF(LABEL(1:10).EQ.'PRINT SPEC') THEN
C         PRINT SPECTRUM FOR EACH ZONE
          READ(99,*)IOSPEC
        ELSEIF(LABEL(1:9).EQ.'DEPARTURE') THEN
C         PRINT DEPARTURE COEFFICIENTS RATHER THEN NUMBER FRACTIONS
          IDEP=1
        ELSEIF(LABEL(1:9).EQ.'PRINT LEV') THEN
C         PRINT LEVEL INFORMATION FOR H I AND He I
          READ(99,*)IOLEVEL
        ELSEIF(LABEL(1:11).EQ.'PRINT EMISS') THEN
C         PRINT EMISSIVITIES FOR EACH ZONE
          READ(99,*)IOEMISS
        ELSEIF(LABEL(1:12).EQ.'PRINT TERMAL') THEN
C         PRINT COOLING AND HEATING FOR EACH ITERATION
          READ(99,*)IOTERM
        ELSEIF(LABEL(1:16).EQ.'PRINT IONIZATION') THEN
C         PRINT IONIZATION RATES FOR EACH ZONE
          READ(99,*)IOION
        ELSEIF(LABEL(1:19).EQ.'PRINT ELECTRON ITER') THEN
C         PRINT IONIZATION RATES FOR EACH ZONE
          READ(99,*)IOELITER
        ENDIF
C       CONSTANT PRESSURE IPRESS = 1
        IF(LABEL(1:3).EQ.'PRE') IPRESS=1
C       CONSTANT DENSITY IDEN = 1
        IF(LABEL(1:3).EQ.'DEN') IDEN=1
C       AGN SPECTRUM 
        IF(LABEL(1:3).EQ.'AGN') IAGN=1
C       PULSAR 
        IF(LABEL(1:4).EQ.'PUL ') THEN
C         SPECTRUM PARAMETERS FOR PULSAR AND FREE FREE SPECTRA
          READ(99,*)RLTOT,ALFA,EMIN,EMAX
        ENDIF
        if(label(1:6).eq.'OBSTAR') then
           iobstar=1
           read(99,*)rltotob,teffob
        endif
C
C       IPULS=0 ; only low ionized ions   
C       IPULS=1 ; higher ionized ions included
C
        IF(LABEL(1:3).EQ.'PUL') IPULS=1
C       CIRCUMSTELLAR:
C         SET ICS = 1
C         INOUT = 0
C         SPECIFY RSHOCK, MTOT
C       FOR REVERSE SHOCK 
C         ICS = 1
C         INOUT=1
C         IREV=1
        IF(LABEL(1:3).EQ.'CS ') ICS=1
        IF(LABEL(1:5).EQ.'CS87A') ICS87A=1
        irev=0
        IF(LABEL(1:3).EQ.'REV') THEN
          IREV=1
          READ(99,*)TECOOL
        ENDIF
C       ONLY LOW IONIZATION IONS INCLUDED: ILOWION = 1
        IF(LABEL(1:3).EQ.'LOW') ILOWION=1
C       STATIC MEDIUM
        IF(LABEL(1:3).EQ.'STA') ISTAT=1
C       SOBOLEV APP. 
        IF(LABEL(1:3).EQ.'SOB') ISOBOL=1
        write(6,*)'isobol= ',isobol
C       ISPH=INDEX SPECIFYING WHETER A SPHERICAL (ISPH=1)
C         OR PLANE (ISPH=0) GEOMETRY SHOULD BE USED
        IF(LABEL(1:3).EQ.'SPH') ISPH=1
C       IPAR=INDEX SPECIFYING WHETHER A PARALLELL (IPAR=1) OR ISOTROPIC
C       BEAM OF IONIZING RAD. SHOULD BE USED.
        IF(LABEL(1:3).EQ.'PAR') IPAR=1
C       IF IUVBLANK = 1 NO UV-CONTINUUM BELOW 3000 A
        IF(LABEL(1:3).EQ.'UVB') IUVBLANK=1
        IF(LABEL(1:4).EQ.'STOP') GOTO 33
      ENDDO
33    CONTINUE

      IPAR=1
      ISPH=1
      IF(ISTAT.EQ.0.AND.ISOBOL.EQ.0) WRITE(0,*)' Sobolev or static?'
      IF(ISTAT.EQ.0.AND.ISOBOL.EQ.0) STOP
      IF(INOUT.EQ.2) WRITE(0,*)' You must specify in or out!'
      IF(INOUT.EQ.2) STOP
      WRITE(6,93)INMOD
93    FORMAT(' INMOD = ',A72)
      CALL PRINTMODELI(' INOUT    = ',INOUT)
      CALL PRINTMODEL (' A1IN     = ',A1IN)
      CALL PRINTMODEL (' B1IN     = ',B1IN)
      CALL PRINTMODEL (' MTOT     = ',MTOT)
      CALL PRINTMODEL (' TOTCOLUMN= ',TOTCOLUMN)
      CALL PRINTMODEL (' RIN      = ',RIN)
      CALL PRINTMODEL (' TDAY     = ',TDAY)
      CALL PRINTMODELI(' NMAX     = ',NMAX)
      CALL PRINTMODELI(' MAX      = ',MAX)
      CALL PRINTMODELI(' ICEN     = ',ICEN)
      CALL PRINTMODEL (' RMCEN    = ',RMCEN)
      CALL PRINTMODEL (' dminit   = ',dminit)
      CALL PRINTMODEL (' rminner  = ',rminner)
      CALL PRINTMODEL (' CHGH     = ',CHGH)
      CALL PRINTMODEL (' CHTE     = ',CHTE)
      CALL PRINTMODEL (' ENMIN    = ',ENMIN)
      CALL PRINTMODEL (' ENMAX    = ',ENMAX)
      CALL PRINTMODELI(' NFEII    = ',NFEII)
      CALL PRINTMODEL (' DECON   = ',DECON)
      CALL PRINTMODEL (' GAMMION  = ',GAMMION)
      CALL PRINTMODEL (' RLTOT    = ',RLTOT)
      CALL PRINTMODEL (' ALFA     = ',ALFA)
      CALL PRINTMODEL (' EMIN     = ',EMIN)
      CALL PRINTMODEL (' EMAX     = ',EMAX)
      CALL PRINTMODEL (' MNI      = ',MNI) 
      CALL PRINTMODEL (' VEXP     = ',VEXP)
      CALL PRINTMODELI(' ISPH     = ',ISPH)
      CALL PRINTMODELI(' IPAR     = ',IPAR)
      CALL PRINTMODELI(' IPULS    = ',IPULS)
      CALL PRINTMODEL (' FILLING    = ',FILLING)
      IF(IPRESS.EQ.1)  WRITE(6,*)' Constant pressure'
      IF(IDEN.EQ.1)  WRITE(6,*)' Constant density'
      IF(IAGN.EQ.1) WRITE(6,*)' AGN-spectrum'
      IF(ICS.EQ.1) WRITE(6,*)' Circumstellar interaction'
      IF(IREV.EQ.1) WRITE(6,*)' Reverse shock'
      IF(ILOWION.EQ.1) WRITE(6,*)' Only low ionization ions'
      IF(ISTAT.EQ.1) WRITE(6,*)' Static geometry'
      IF(ISOBOL.EQ.1) WRITE(6,*)' Sobolev line transfer'
      IF(IUVBLANK.EQ.1) WRITE(6,*)' UV-blanking!'
C     ************************************************************
C     *****
C     ENERGY INTERVALS
C     *****
C     ************************************************************
      JMIN0=JMIN
      JJ0=JJ
      DO J=JMIN0,JJ0
        IF(E(J).LT.ENMIN) JMIN=J+1
        IF(E(J).LT.ENMAX) JJ=J+1
      ENDDO
      WRITE(6,*)'JMIN, JJ,     ENMIN,     ENMAX'
      WRITE(6,92)JMIN,JJ,ENMIN,ENMAX
92    FORMAT(2I5,1PE12.3,E12.3)
C     NINQ= NUMBER OF INTERVALS -1 BETWEEN 13.6 AND 15.3 EV
      NINQ=2
      JJ=jj+NINQ
      IF(INOUT.EQ.0.AND.ICS.EQ.0) CHGH=0.25
C     IF ITE=1 SOLV THERMAL EQ. IF ITE=+1 SKIP THIS.
      ITE=+1
      DO 100 IQ=1,1
      IF(ICS.EQ.1.AND.IPULSSP.NE.1.and.ifreefree.eq.0) THEN
        write(*,*)' to shockspec'
        read(19,*)tday
        write(*,*)'tday ',tday
      ENDIF
      DO 100 IU=1,1
      IUY=1
c!!
C     ITER=NUMBER OF ITERATIONS OF STATISTICAL EQUATIONS
      ITER=200
      IGAMM=1
C     AV=POWERLAW OF VELOCITY LAW, V(R)=R**AV
      AV=1.
C     R1= PHOTOSPHERIC RADIUS IN 1E15 CM
C     TEFF = EFFECTIVE TEMPERATURE OF  BLACKBODY CONTINUUM
      CALL  SNPAR(TDAY,R1,TEFF,RLUM)
      IF(IBLACKB.EQ.1) THEN
        TEFF=TEFFBB
      ENDIF
      IF(IPH.EQ.1) THEN
        R1=RPH
      ENDIF
C     CIRCUMSTELLAR INTERACTION
      IF(ICS.EQ.1.AND.IPULSSP.NE.1.and.ifreefree.eq.0) THEN
C       SHOCK RADIUS IN 1E15 CM
        read(19,*)rshock
        r1=rshock
        IF(IREV.EQ.1) THEN
          R1=RSHOCK
        ENDIF
      ENDIF
      R1=1.E-15*R1
      TIME=8.64E4*TDAY
C     RADIUS OF THE OUTER BOUNDARY IN CM
      RMAX=VEXP*TIME
      IF(RIN.LT.R1) RIN=R1+1.E-5
C
      ejdens=0.
      IF(ICS.EQ.1.AND.IPULSSP.NE.1.and.ifreefree.eq.0) THEN
        read(19,*)mtot
c       temperature of reverse and c-s shocks
        read(19,*)terev,tecs
        read(19,*)revmass
        read(19,*)EJDENS
        write(6,*)' Trev = ',terev
        write(6,*)' Rev. mass: ',revmass
        write(6,*)' Ej. dens. ',ejdens
      ENDIF   
      IF(IAGN.EQ.1) MTOT=1.D66
      IF(IAGN.EQ.1) REVMASS=1.D66
      MSUN=MTOT/SOLMA
      RMASS=MTOT/SOLMA
C     DELM = INITIAL MASS INTERVAL
      DELM=(MSUN-RMCEN)/FLOAT(MAX-1-ICEN)
      RMEN=MTOT*2./3.
      RMCO=RMEN/3.
      R0=RMAX
C     GAMMA RAY LUMINOSITY DUE TO THE DECAY OF 56-COBOLT
C     56   CO DECAYS
      GAMLUM=1.36D43*MNI*EXP(-TDAY/111.26)
C     57 CO DECAYS
C     ASSUME SOLAR RATIO OF 56CO/57CO = 0.0243
      GAM57=3.47D38*(MNI/0.075)*EXP(-TDAY/391.)
      GAMLUM=GAMLUm+GAM57
      RMAX=1.E-15*RMAX
      WRITE(6,9023)R1,RMAX
 9023 FORMAT(' R1= ',E11.4,' RMAX= ',E11.4)
C     NPRINT=INTERVALS BETWEEN PRINTOUTS
      NPRINT=0
C     A1=LOWER LIMIT OF START TEMPERATURE
C     B1=UPPER        D:O
      A1=A1IN
      B1=B1IN
C     TIN=INITIAL GUESS OF TEMPERATURE
      TIN=0.55E4
C     E10=TOLERANCE IN TEMPERATURE
      E10=2.E1
C     E20 = RELATIVE ERROR IN FUNCTION RAD=HEAT-XEL*COOL
      E20=0.01
      MQMAX=MAX
      R1Q=R1
      fillingq=filling

      call trans

 100  CONTINUE
      END

      SUBROUTINE PRINTMODEL(LABEL,VAR)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*(*) LABEL
      WRITE(6,9)LABEL,VAR
9     FORMAT(A12,1E12.4)
      RETURN
      END

      SUBROUTINE PRINTMODELI(LABEL,IVAR)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*12 LABEL
      WRITE(6,9)LABEL,IVAR
9     FORMAT(A12,I5)
      RETURN
      END

      DOUBLE PRECISION FUNCTION DE(R)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/DENS/D0,R0,RN
      COMMON/TPAR/RIN,DR,R1,TDAY
      COMMON/HYDROS/HSCALE,DEN1,XEL,AMEAN
      RCGS=1.E15*R1*R
C     POLYTROPIC DENSITY FOR WHITE DWARF
      IF(RCGS.LT.R0) DE=D0/((1.+10.13*(RCGS/R0)**3)*AMEAN)
C     DE=D0/AMEAN
C     IF(RCGS.GE.R0) DE=0.
C     SHELL
      RINCGS=1.E15*RIN
      IF(RCGS.LT.RINCGS) DE=0.
C     IF(RCGS.GT.RINCGS) DE=D0/AMEAN
C     DENSITY FROM STRUC
      DE=DEN1
      RETURN
      END

