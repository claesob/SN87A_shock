

      SUBROUTINE PROFILE(IMAX,R1CGS,VELI,RMCV,FBIA,SIREDA,citot,fbq)
C     ***********************************************************
C     *****
C     CALCULATE THE LINE PROFILE FOR THE FORBIDDEN LINES
C
C     INTEGRATE FROM INNER RADIUS R(JF), CORRESPONDING TO  FREQUENCY
C     FREK(JF), OUT TO R(IMAX)= OUTER RADIUS. 
C     FTOT IS THE LUMINOSITY PER ANGSTROM (ERG/S A)
C     *****
C     ***********************************************************
C     ENUMERATION OF LINES
C     W(I)
C      1 = C II  1334,  2 = C II 2326,  3 = C III  977,  4 = C III 1909
C      5 = C IV  1550,  6 = O VI 1030,  7 = O V   1220,  8 = O IV  1400
C      9 = O III 1664, 10 = N II 2143, 11 = N III 1750, 12 = N III  990
C     13 = N IV  1486, 14 = N V 1240,  15 = MG II 2800, 16 = SI IV 1400
C     17 = SI II 2335, 18 = FE I    ,  19 = FE II     , 20 = MG I  4572
C     21 = MG I  2852, 22 = NA I 5889, 23 = CA II 7300, 24 = CA II 3950
C     25 = CA II 8500, 26 = C II 158M, 27 = C I   2966, 28 = MG  I 5176
C     29 = MG  I 8806, 30 = MG I 34, 31 = O I 6300-64 32 = O I   2972
C     33 = O I   5572, 34 = O I  1356, 35 = O I   1302, 36 = O I   7774
C     37 = O I   8446, 38 = O I  1027, 39 = O I    990, 40 = MG I  4571(
C     41 = NA I 5889(R)42 = FE III     43 = S II 314.5M 44 = SI II 34.81M
C     45 = A II 6.985M 46 = NE II 12.81M 47 = HE I
C
C     61 =             62 = HE I  63 = H I   ...72 = H I 73 SI III 1892
C
C     74 = S IV 749.6,
C     75 = NE VIII 773.7,  76 = S V 786.5,     77 = S VI 937.1, 
C     78 = S IV 1069.6,    79 = S III 1197.5,  80 = S IV 1393.4            
C     81 = FE IV
C
C     FEII : 82 - (NFEL+81)
C
CK NEW LINES : 74 - 80
C
C
C     FB(I)
C      1 = O III 4959-5007,  2 = O III     2321,  3 = O III       4363
C      4 = O II       3726,  5 = O II      2470,  6 = O II        7320
C      7 = N II    6548-83,  8 = N II      3063,  9 = N II        5755
C     10 = O I     6300-64, 11 = O I       2964, 12 = O I         5581
C     13 = N I        3468, 14 = N I       5201, 15 = N I        10406
C     17 = O I        63 M, 19 = O I      145 M,
C     20 = C I        9818, 21 = C I       4619, 22 = C I         8729
C     23 = S II  6718-6733, 24 = S II      4071, 25 = S II       10320
C     26 = S I      25.25M, 27 = S I     56.31M, 28 = SI I      129. M
C     29 = SI I 68.5 M    , 30 = SI I    1.636M, 31 = SI I     1.0995M
C     32 = C I     370.4 M, 33 = C I     609.1M,
C     34 = S III      3722, 35 = S III     6312, 36 = S III  9069-9532
C     37 = NE III     1815, 38 = NE III    3342, 39 = NE III 3869-3968
C     40 = NE IV      1602, 41 = NE IV  2423-25, 42 = NE IV    4714-24      
C     43 = NE V       1575, 44 = NE V      2975, 45 = NE V   3346-3426
C
C    "EMPTY" : 46-50  ,  O : 51-54  ,  FEIII : 55-354 , FEIV : 355-654
C
CK NEW FB : 34 - 45
C
C     DATA WLHEI/584.,10830.,3889.,3188.,20581.,7065.,5876.,4713.,
C    &4471.,42944.,12547.,186233.,21120.,17002.,19543.,18648./
C
      implicit real*8(a-h,o-z)
      REAL*8 MTOT,MNI
c      PARAMETER (MD=350,MDP1=MD+1)
      include 'PARAM'
      PARAMETER (NL=340,NLP1=NL+1)
      PARAMETER (NFEL=3000)
      COMMON/RADIE/RA(0:MD)
      COMMON/PHY/DENA(MD)
      COMMON/ELEC/DELA(MD)
      COMMON/MPAR/MTOT,MNI,VEXP
      parameter (nlp=30000)
      COMMON/EQUIV/RADF,FR(2,nlp),W(nlp),CIN(nlp),WLI(nlp+NFEL),FB(300),tauline(nlp+nfel)
      COMMON/PROFFE/WLFEIII(300),WLFEIV(300),FEIIIA(300,MD),
     &                 FEIVA(300,MD)
      COMMON/CINOUT/INOUT,IPULS
change ck : sired(100,md) -> sired(250,md),sireda(100,md) -> sireda(250,md)
      DIMENSION FBIA(110,MD),FTOTAL(300+NFEL,MDP1),FREK(MDP1),
     &                                 SIREDA(300+NFEL,MD)
      DIMENSION WLF1(300),WLF(1000),citot(nlp),fbq(100)
      DIMENSION R(MD),DEN(MD),DEL(MD),SIRED(300+NFEL,MD),FBI(900,MD)
      DATA PI/3.1415926E0/
      OPEN(26,FILE='sltspec')
      rewind 26
      WRITE(26,*)IMAX,R1CGS,VELI,RMCV,INOUT
      WRITE(26,*)FBIA,SIREDA,citot,fbq
      WRITE(26,*)RA,DENA,DELA
      WRITE(26,*)RADF,FR,W,CIN,WLI,FB
      WRITE(26,*) WLFEIII,WLFEIV,FEIIIA,FEIVA
      close (26)
      goto 1234
      IF(INOUT.EQ.1) GOTO 8456
C     
C     CHANGE DIRECTION OF VARIABLES FROM CENTER TO SURFACE
C
      nforb=110
      ntforb=600+nforb
      DO 180 I=1,IMAX
      R(I)=RA(IMAX+1-I)
      DEN(I)=DENA(IMAX+1-I)
      DEL(I)=DELA(IMAX+1-I)
      DO 181 K=1,ntforb
      IF(K.LE.nforb) FBI(K,I)=FBIA(K,IMAX+1-I)
      IF(K.GT.nforb.AND.K.LT.300+nforb) 
     &     FBI(K,I)=FEIIIA(K-nforb,IMAX+1-I)
      IF(K.GT.300+nforb.AND.K.LT.600+nforb)
     &     FBI(K,I)=FEIVA(K-300-nforb,IMAX+1-I)
 181  CONTINUE
c$$$      DO K=1,nlp+NFEL
c$$$         SIRED(K,I)=SIREDA(K,IMAX+1-I)
c$$$      enddo
180   continue
8456  CONTINUE
C
C     FIRST FORBIDDEN LINES
C
      CALL WLFORB(WLF1)
      nforb=110
      DO K=1,600+nforb
      IF(K.LE.nforb)WLF(K)=WLF1(K)
      IF (K.GT.nforb.AND.K.LE.300+nforb)WLF(K)=WLFEIII(K-nforb)
      IF (K.GT.300+nforb.AND.K.LE.600+nforb)WLF(K)=WLFEIV(K-300-nforb)
      ENDDO
      ntforb=600+nforb
      WRITE(16,*)VELI
      WRITE(16,*)(IMAX-1)
      FREK(IMAX)=1.
      DO 3287 K=1,ntforb
      IF(WLF(K).LE.0.) GOTO 3287
      DO JF=1,IMAX
      FTOTAL(K,JF)=0.
      FREK(JF)=R(JF)/R(IMAX)
      DO 1100 I=JF,IMAX
C      IF(JF.GT.1)DFTOT=(FBI(K,I)*DEL(I)+FBI(K,I+1)*DEL(I+1))*
C     &                              DEN(I+1)**2*(R(I+1)**2-R(I)**2)
C      IF(JF.EQ.1)DFTOT=FBI(K,2)*DEL(2)*DEN(2)**2*(R(2)**2-R(1)**2)
      FBIDEN2=FBI(K,I)*DEN(I)**2*DEL(I)
      IF(I.GT.JF.AND.I.LT.IMAX)DFTOT=FBIDEN2*R(I)*(R(I+1)-R(I-1))/2.
      IF(I.EQ.JF) DFTOT=FBIDEN2*R(I)*(R(I+1)-R(I))/2.
      IF(I.EQ.IMAX) DFTOT=FBIDEN2*R(I)*(R(I)-R(I-1))
 1100 FTOTAL(K,JF)=FTOTAL(K,JF)+2.*PI*R1CGS**2*RMCV*DFTOT/WLF(K)
      ENDDO
 3287 CONTINUE
C
      FMAX=0.D0
      DO K=1,ntforb
      IF (FTOTAL(K,2).GT.FMAX)FMAX=FTOTAL(K,2)
      ENDDO
C
C
C     PRINT RESULT TO FILE 16
C
      DO 3286 K=1,ntforb
      IF (FTOTAL(K,2).LE.1.D-5*FMAX) GOTO 3286
      WRITE(16,*)K,WLF(K)
      VC=VELI/3.E10
      FLUM=0.
      DO 7372 IL=1,IMAX-1
      IF(IL.GT.1) WL1=WLF(K)*(1.-VC*FREK(IL-1))
      IF(IL.EQ.1) WL1=WLF(K)
      WL2=WLF(K)*(1.-VC*FREK(IL+1))
      DWL=(WL1-WL2)/2.
      FLUM=FLUM+FTOTAL(K,IL)*DWL
7372  WRITE(16,9274)FREK(IL),FTOTAL(K,IL)
	FLUM=2.*FLUM
	if(flum.le.0.) flum=1.
      IF(FBQ(K).GT.0..and.k.le.100) RATFF=abs(log(fbq(k)/flum))
c      if(RATFF.gt.0.5.and.k.le.100) WRITE(6,9483)K,WLF(K),FLUM,
c     &                                                      fbq(k)
9483  FORMAT(' LARGE DIFF. IN PROFILE!',I5,4E13.5)
9274  FORMAT(1PE13.5,E13.5)
 3286 CONTINUE
      write(16,*)'nuedeslut'
C
C     TWO LEVEL LINES FROM RLOSS
C
      DO 287 K=1,81+NFEL
      DO 287 JF=1,IMAX
      FTOTAL(K,JF)=0.
      FREK(JF)=R(JF)/R(IMAX)
      DO 100 I=JF,IMAX
      IF(I.GT.JF.AND.I.LT.IMAX)DFTOT=SIRED(K,I)*R(I)*(R(I+1)-R(I-1))/2.
      IF(I.EQ.JF) DFTOT=SIRED(K,I)*R(I)*(R(I+1)-R(I))/2.
      IF(I.EQ.IMAX) DFTOT=SIRED(K,I)*R(I)*(R(I)-R(I-1))
c      IF(JF.GT.1)DFTOT=(SIRED(K,I)+SIRED(K,I+1))*(R(I+1)**2-R(I)**2)
c     IF(JF.EQ.1)DFTOT=SIRED(K,2)*(R(2)**2-R(1)**2)
  100 FTOTAL(K,JF)=FTOTAL(K,JF)+4.*PI*2.*PI*R1CGS**2*DFTOT
  287 CONTINUE
C
      FMAX=0.D0
      DO K=1,81+NFEL
      IF (FTOTAL(K,2).GT.FMAX)FMAX=FTOTAL(K,2)
      ENDDO
C
C
C     PRINT RESULT TO FILE 16
C
      DO 286 K=1,nlp+NFEL
      IF (FTOTAL(K,2).LE.1.D-5*FMAX) GOTO 286
      WRITE(16,*)K,WLI(K)
      FLUM=0.
      DO 372 IL=1,IMAX-1
      IF(IL.GT.1) WL1=WLI(K)*(1.-VC*FREK(IL-1))
      IF(IL.EQ.1) WL1=WLI(K)
      WL2=WLI(K)*(1.-VC*FREK(IL+1))
      DWL=(WL1-WL2)/2.
      FLUM=FLUM+FTOTAL(K,IL)*DWL
 372  WRITE(16,9274)FREK(IL),FTOTAL(K,IL)
	FLUM=2.*FLUM
	if(flum.le.0.) flum=1.
      IF(CITOT(K).GT.0.) RATCF=abs(log(CITOT(k)/flum))
c      WRITE(6,9489)K,WLI(K),FLUM,citot(k),ratcf
c      if(RATCF.gt.0.5)WRITE(6,9483)K,WLI(K),FLUM
c      wRITE(6,9483)K,WLI(K),FLUM
c     &                                               ,CITOT(k)
  286 CONTINUE
      write(16,*)'nuedeossoslut'
1234  continue
      RETURN
      END

      SUBROUTINE IONLABEL(LAB1)
      CHARACTER*8 LAB(200),LAB1(200)
      DATA LAB/'H I     ','HE I    ','HE II   ','O VI    ','O VII   '
     &        ,'O VIII  ','O V     ','C III   ','C IV    ','C V     '
     &        ,'C VI    ','C I     ','C II    ','O I     ','O II    '
     &        ,'O III   ','O IV    ','N I     ','N II    ','N III   '
     &        ,'N IV    ','N V     ','N VI    ','N VII   ','SI I    '
     &        ,'SI II   ','SI III  ','SI IV   ','SI V    ','SI VI   '
     &        ,'SI VII  ','SI VIII ','SI IX   ','SI X    ','SI XI   '
     &        ,'SI XII  ','SI XIII ','SI XIV  ','MG I    ','MG II   '
     &        ,'MG III  ','FE I    ','FE II   ','FE III  ','FE IV   '
     &        ,'AL I    ','AL II   ','AL III  ','AL IV   ','CA I    '
     &        ,'CA II   ','CA III  ','NA I    ','NA II   ','NA III  '
     &        ,'S I     ','S II    ','S III   ','S IV    ','S V     '
     &        ,'S VI    ','S VII   ','S VII   ','S IX    ','S X     '
     &        ,'S XI    ','S XII   ','S XIII  ','S XIV   ','S XV    '
     &        ,'S XVI   ','NE I    ','NE II   ','NE III  ','NE IV   '
     &        ,'NE V    ','NE VI   ','NE VII  ','NE VIII ','NE IX   '
     &        ,'NE X    ','AR I    ','AR II   ','FE V    ','FE VI   '
     &        ,'FE VII  ','FE VIII ','FE IX   ','FE X    ','FE XI   '
     &        ,'FE XII  ','FE XIII ','FE XIV  ','FE XV   ','O III 1D'
     &      ,5*'        ','BALMER  '
     &        ,'PASCHEN ','BRACKETT','PFUND   ','O I 4   ','O I 5   '
     &        ,'O I 6   ','O I 7   ','O I 8   ','O I 9   ','FF      '
     &        ,89*' '/
      DO K=1,111
            LAB1(K)=LAB(K)
c            write(0,9)k,lab(k)
 9          format('lab ',i5,a8)
      ENDDO
      RETURN
      END

      SUBROUTINE WLFORB(WLF)
      IMPLICIT REAL*8(A-H,O-Z)
      common/forbwl/wlft(300)
      common/wlforbk/wlfm(300)
      DIMENSION WLF(300),WLF1(71)
C     WAVELENGTHS OF ALL FORBIDDEN LINES
      DATA WLF1/5007.,2321.,4363.,3726.,2470.,7320.,6548.,3063.,
     &         5755.,6300.,2964.,5581.,5201.,3468.,10406.,0.,631800.,
     &         440560.,1455600.,9818.,4619.,8729.,6718.,4071.,10320.,
     &         252460.,563060.,1297016.,684932.,16360.,10995.,
     &         3704000.,6091000.,9069.,3722.,6312.,3869.,1815.,3342.,
     &         2423.,1602.,4714.,1575.,2975.,3346.,155500.,108640.,
     &         360200.,883300.,518000.,11287.,7002.,9264.,6157.,
     &         89910.,2.184e5,7135.8,7751.1,3006.1,3109.1,5191.8,
     &         4711.3,4740.2,2853.7,2868.2,7.741e5,564721.,7237.3,
     &         7170.6,7333.4,7262.8/
      nforb=110
      do k=1,nforb
         wlf(k)=wlfm(k)
      enddo
      DO K=1,71
            WLF(K)=WLF1(K)
      ENDDO
c!!!
      goto 1111
c  All lines which should not be included in summation
      do k=1,71
         if(k.ge.1.and.k.le.6) then
            wlf(k)=-1.
         elseif(k.ge.13.and.k.le.16) then
            wlf(k)=-1.
         elseif(k.ge.23.and.k.le.25) then
            wlf(k)=-1.
         elseif(k.ge.40.and.k.le.45) then
            wlf(k)=-1.
         elseif(k.ge.49.and.k.le.50) then
            wlf(k)=-1.
         endif
      enddo    
 1111 continue
      RETURN
      END

      subroutine printbalance(ion,te,den,xel,xion,ab,b)
      implicit real*8(a-h,o-z)
      PARAMETER (NL=340,NLP1=NL+1)
      common/test1/qa(6,20,20),rateph(6,20),tsav(6,20,20),esav(6,20,20)
      COMMON/NLEV/e00,NION,NHY,NP1H,NMA,NMI
      COMMON/HEIION/ZHEICO,ZHEILA,ZOTSHI,ZBALM,ZHBALM
      COMMON/A14/C(NL,NL),CI(NL),G(NL),E(NL),A(NL,NL),WL(NL,NL)
     &,DCDT(NL,NL),DCIDT(NL)
      common/qradcoll/qcoll(nl,nl),qrad(nl,nl)
      dimension qcollt(nl,nl),qradt(nl,nl)
      dimension poplte(nl),b(nl),qap(6,20,20)
      if(ion.eq.1) nlq=5
      if(ion.eq.5) nlq=5
      if(ion.eq.3) nlq=16
      if(ion.eq.4) nlq=191
      if(ion.eq.6) nlq=121
      if(nlq.gt.20) nlq=20
      if(ion.eq.1) write(6,*)' Ca II balance'
      if(ion.eq.5) write(6,*)' Hydrogen balance'
      if(ion.eq.3) write(6,*)' Helium balance'
      if(ion.eq.6) write(6,*)' Fe I balance'
      if(ion.eq.4) write(6,*)' Fe II balance'
      if(ion.eq.3) then
c$$$        call atdathe
      elseif(ion.eq.1) then
        call atdat
      elseif(ion.eq.5) then
        call atdath
      elseif(ion.eq.4) then
        call atdatfe
      elseif(ion.eq.6) then
        call atdatfei
      endif
      tev=te/1.1609e4
      denf=den*xel*xion
      do i=1,20
        poplte(i)=2.07e-16*g(i)*denf*exp((e00-e(i))/tev)/te**1.5
        do k=1,12
          qap(ion,i,k)=qa(ion,i,k)*poplte(i)
        enddo
        if(ion.eq.3) then
           do j=1,nlq
              qcollt(i,j)=qcoll(i,j)*poplte(i)
              qradt(i,j)=qrad(i,j)*poplte(i)
           enddo
        endif
      enddo
      write(6,*)' LTE populations'
      write(6,9345)(poplte(i),i=1,nlq+1)
      write(6,*)' NLTE populations'
      write(6,9345)(poplte(i)*b(i),i=1,nlq+1)
      write(6,*)' departure coeff.'
      write(6,9345)(b(i),i=1,nlq+1)
      write(6,*)'optical depths'
      do j=2,8
        write(6,9345)(tsav(ion,i,j),i=1,j-1)
      enddo
      write(6,*)'escape prob.'
      do j=2,8
        write(6,9345)(esav(ion,i,j),i=1,j-1)
      enddo
      write(6,*)'net collisional + radiative rates'
      write(6,9345)(qap(ion,i,1),i=1,nlq)
      write(6,*)'net radiative rates to and from lower levels'
      write(6,9345)(qap(ion,i,9),i=1,nlq)
      write(6,*)'net radiative rates to and from upper levels'
      write(6,9345)(qap(ion,i,8),i=1,nlq)
      write(6,*)'net coll. rates to and from lower levels'
      write(6,9345)(qap(ion,i,11),i=1,nlq)
      write(6,*)'net coll. rates to and from upper levels'
      write(6,9345)(qap(ion,i,10),i=1,nlq)
      if(ion.eq.3) then
         write(6,*)' ind. coll. rates'
         do i=1,nlq
            write(6,9345)(qcollt(i,j),j=1,nlq)
         enddo
         write(6,*)' ind. rad. rates'
         do i=1,nlq
            write(6,9345)(qradt(i,j),j=1,nlq)
         enddo
      endif
      write(6,*)'photoionization rates'
      write(6,9345)(qap(ion,i,4),i=1,nlq)
      write(6,*)'absolute photoionization rates'
      write(6,9345)(rateph(ion,i),i=1,nlq)
      if(ion.eq.5) then
        write(6,*)' He I 584 and 504 ionization fraction, OTS, N=2'
        write(6,9345)ZHEILA,ZHEICO,ZOTSHI,ZHBALM
        write(6,*)'non-th excitation'
        write(6,9345)(qap(ion,i,12),i=1,nlq)
      endif
      write(6,*)'recombination rates'
      rect=0.
      do i=1,nlq
         rect=rect+qap(ion,i,2)
      enddo
      if(rect.gt.0.) then
         rec1=qap(ion,1,2)/rect
      else
         rect1=0.
      endif
      write(6,9345)(qap(ion,i,2),i=1,nlq),rect,rec1
      write(6,*)'coll. ionization rates'
      write(6,9345)(qap(ion,i,5),i=1,nlq)
      write(6,*)'gamma ionization rates'
      write(6,9345)qap(ion,1,6)
      write(6,*)'total net rate (should be small)'
      write(6,9345)(qap(ion,i,7),i=1,nlq)
9345  format(1pe12.3,5e12.3)
      return
      end

      real*4 function etime(tarr)
      real*4 tarr(2)
      return
      end

c23456789012345678901234567890123456789012345678901234567890123456789012
