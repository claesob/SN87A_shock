      DOUBLE PRECISION FUNCTION FX(X)
      IMPLICIT REAL*8(A-H,O-Z)
      include 'param'
      PARAMETER (NL=340,NLP1=NL+1)
      COMMON/IND/I
      COMMON/DIF/EM(MD,NE1:NE2),TAUQ(MD,NE1:NE2),TAUTOT(MD,NE1:NE2),
     &           EMC(MD,NE1:NE2)
      COMMON/NLEV/e00,NION,NI,NP1H,NHMAX,NHMIN
      COMMON/A5/TAU,ALFA,EN,TSN,XL40,TXEV,R15
      COMMON/NION/ION
      COMMON/A14/C(NL,NL),CI(NL),G(NL),E(NL),A(NL,NL),
     &WL(NL,NL),DCDT(NL,NL),DCIDT(NL)
      COMMON/A16/SIG(NL),GA(NL),BP(NL)
      COMMON/A21/N
      COMMON/A22/R,TE,MQ
      COMMON/QSOM/qso,IAGN,ISTAT,ISOBOL
      COMMON/SPECT/TEL,FD(MD,NE1:NE2),F0(NE1:NE2),IPARA
      DATA PI/3.1415926E0/,ELCH/1.60219E-12/,AMU/1.660531E-24/
      FX=0.
      IF(TSN.LE.0.) WRITE(6,*)'TS',TSN
      TBEV=TSN/1.1609D4
      TEEV=TE/1.1609D4
      DEN=DE(R)
C     CHI_N = IONIZATION POTENTIAL OF LEVEL N
      CHI_N=E00-E(N)
      ELOG0=LOG10(CHI_N)
      E1=CHI_N*10.**X+1.D-10
C
C     FJ = MEAN INTENSITY OF RADIATION ( IN CGS UNITS )
C
      
C!    CORRECT FOR ABSORBTION IN THE BALMER CONTINUUM OF THE INCOMING RAD. (NOT DONE!)
c      IF(E1.GT.3.4.AND.E1.LT.13.6) THEN
c      FJ=EXP(-TAUTOT(I-1,-13)*(3.4/E1)**3)*FMEAN(2,E1)
c      ELSE 
      FJ=FMEAN(2,E1)
c     ENDIF
      
      IF(E1.LT.CHI_N) GOTO 100
      EX3=(E1-CHI_N)/TEEV
C
C     CROSS-SECTION AND GAUNT FACTOR FROM SEATON(1959)
C
      IF(ION.NE.2) GOTO 3887
C
C     CROSSECTIONS OF O I FROM DAWIES AND LEWIS
C
      IF(N.EQ.4) ERY=(E1-4.481)/13.61
      IF(N.EQ.5) ERY=(E1-4.110)/13.61
      IF(N.EQ.6) ERY=(E1-2.886)/13.61
      IF(N.EQ.7) ERY=(E1-2.642)/13.61
      IF(N.EQ.8) ERY=(E1-1.545)/13.61
      IF(N.EQ.9) ERY=(E1-1.095)/13.61
      IF(ERY.LE.0.) GOTO 100
      IF(N.EQ.4.AND.ERY.LT..4) SL=(-1.034*ERY+1.365-.0491/ERY)-2.
      IF(N.EQ.4.AND.ERY.GT..4) SL=-1.20-2.71*(ERY/.4)
      IF(N.EQ.5.AND.ERY.LT..5) SL=(-.4038*ERY+1.783-.0449/ERY)-2.
      IF(N.EQ.5.AND.ERY.GT..5) SL=-.483*ERY-.286
      S=10.**SL
      IF(N.EQ.6) S=10.**(-7.204*ERY+1.657-1.)+10.**(-.3748*ERY-.5877)
      IF(N.EQ.7) S=10.**(-5.509*ERY+.2457)+10.**(-.54*ERY-.1715)
      IF(N.EQ.8) S=20.5*(1.545/E1)**3.56
      IF(N.EQ.9) S=20.5*(1.545/E1)**3.56
      SIGM=S
 3887 CONTINUE
      GAUNT=1.
      if(ion.ne.5) GAUNT=1.
      if(ion.eq.5.and.n.le.6) then
         GAUNT=gbf(n,e1)
      elseif(ion.eq.5) then
         gaunt = 1.
      endif
      IF(ION.EQ.1) then
         SIGM=CRPHCA(N,E1)
      ELSEIF(ION.eq.2) then
         SIGM=SIGOX(N,E1)
      ELSE
         SIGM=SIG(N)*GAUNT*(CHI_N/E1)**GA(N)
      ENDIF
      
      if(mq==1.or.mq==3.or.mq==4.or.mq==6.or.mq==7) then
C
C     RECOMBINATION RATE (mq=1)
C
         IF(EX3.GT.700.) GOTO 100
         FX=4.36685D9*SIGM*(2.0845D-4*E1**3.+FJ)*EXP(-EX3)
         

C     NET RECOMBINATION COOLING (MQ=3)C
         IF(MQ.EQ.3) then
            FX=ELCH*(E1-CHI_N)*FX
         endif

C     DR/DT (MQ=4)            
         IF(MQ.EQ.4) then
            FX=FX*EX3/TE
         endif
         
C     TOTAL REC COOLING RATE
         IF(MQ.EQ.7) then
            FX=ELCH*E1*FX
         endif

C     DRECCOL/DT (MQ=6)
         IF(MQ.EQ.6) then
            FX=ELCH*(E1-CHI_N)*FX*EX3/TE
         endif

      elseif(MQ.eq.2) then
C     
C     HEATING RATE (MQ=2)
C 
C     FOR NET HEATING (= electron heating rate) use
         FX=ELCH*4.36685E9*SIGM*(E1-CHI_N)*FJ
         GOTO 100

C     IONIZATION RATE (MQ=5)
C     4*PI*LN(10.)*1.E-18/H = 4.36685E9
C
      elseIF(MQ.EQ.5) then
         FX=4.36685E9*SIGM*FJ
         GOTO 100

      elseIF(MQ.EQ.8) then
C     FOR TOTAL HEATING RATE (including ionization energy)
         FX=ELCH*E1*4.36685D9*SIGM*FJ
      endif
      
100   CONTINUE

      RETURN
      END

      SUBROUTINE HORATE(ionel,R,TE)
c split of 2s and 2p for H      
      IMPLICIT REAL*8(A-H,O-Z)
c      PARAMETER (MD=350,MDP1=MD+1)
      include 'param'
      PARAMETER (NL=340,NLP1=NL+1)
c      PARAMETER (NE1=-100,NE2=250,NE3=NE2+1)
C
C     CALCULATE RECOMBINATION AND PHTOOIONIZATION RATES 
C
      COMMON/IND/I
      COMMON/A14/C(NL,NL),CI(NL),G(NL),EN(NL),A(NL,NL),
     &WL(NL,NL),DCDT(NL,NL),DCIDT(NL)
      COMMON/A17/RECNET(NL),RECNDT(NL),RECCO(NL),RTE(NL),
     &PH(NL),DRDT(NL),RECT(6,NL),PHET(NL)
      COMMON/NLEV/e00,NION,N,NP1H,NHMAX,NHMIN
      COMMON/PHEAT/PHE(NL)
      COMMON/PHOTOHEAT/PHEAT(10,NL),PHEATT(10,NL)
      COMMON/A26/SIGO(20,NE1:NE2),SIGCA(10,NE1:NE2),SIGHE(16,NE1:NE2),
     &SIGFE(NL,NE1:NE2),sigfei(nl,ne1:ne2),SIGH(50,NE1:NE2),SIGHE2(50,NE1:NE2)
      COMMON/FRE/NINT,JMIN,JJ
      COMMON/INT/FL(2,NE1:NE2),SI(14,27,NE1:NE2),E1(NE1:NE3),E(NE1:NE3)
      COMMON/DTAU/FLU(NE1:NE2)
      COMMON/DIF/EM(MD,NE1:NE2),TAU(MD,NE1:NE2),TAUTOT(MD,NE1:NE2),
     &           EMC(MD,NE1:NE2)
      COMMON/HEIION/ZHEICO,ZHEILA,ZOTSHI,ZBALM,ZHBALM
      COMMON/RECAL/RECO(NL)
      common/hydphot/phhi(nl)
      common/ctionfe/phfe(2)
      common/fill/filling
      parameter(nlp=30000)
      common/linepoint/nlines,jline(nlp),iion(nlp)
      common/lineopac/totopl(nlp)
      common/linephoto/phline(nlp) 
      common/timecheck/time,itime
      common/preion/ipre
      DATA PI/3.1415926E0/,ELCH/1.60219E-12/,AMU/1.660531E-24/
      TS=TE
c      write(6,*)' horate ',itime,ion,n

      do j=1,n
         ph(j)=0.
      enddo

      if(ionel==5) then
         nj=16
      else
         nj=n
      endif
 
      DO 500 J=1,Nj
C
C     FIRST CALCULATE THE RECOMBINATION RATE FOR ALL LEVELS
C
         CALL RECOMB(IONel,J,TS,AL,RI)
         RECO(J)=AL
C
C     THEN THE IONIZATION AND HEATING RATE
C     
         E0=E00-EN(J)
         ZA=0.
         HEAT=0.
         HEATT=0.
         phtot=0.
         DO J1=JMIN,JJ
            IF(E1(J1).GT.E0) then
               if(ionel==5) then
               elseif(ionel==4) then
                  SIGM=SIGFE(J,J1)
               elseif(ionel==6) then
                  SIGM=SIGFEI(J,J1)
               elseif(ionel==2) then
                  SIGM=SIGO(J,J1)
               endif

               if(sigm.gt.0.) then
                  ZC=4.*PI*FL(2,J1)*SIGM*(E(J1+1)-E(J1))/(ELCH*E1(J1))
                  ZA=ZA+ZC
                  HEAT=ELCH*ZC*(E1(J1)-E0)+HEAT
                  dHEAT=ELCH*ZC*(E1(J1)-E0)
                  HEATT=ELCH*ZC*E1(J1)+HEATT
               endif
            endif            
         enddo
         PHEAT(IONel,J)=HEAT

         PHEATT(IONel,J)=HEATT
         PH(J)=ZA
         
c     include photoionization due to continous absorption of line
c     photons from rlossd
c     
         if(ionel.eq.5.and.j.eq.1) then
c     ground state
            k=1
         elseif(ionel.eq.5.and.j.ge.2) then
c     n=2 and higher continua k=101 = Balmer etc
            k=99+j   
         elseif(ionel.eq.2.and.j.eq.1) then
            k=14
         elseif(ionel.eq.3.and.j.eq.1) then
            k=2
         elseif(ionel.eq.4.and.j.eq.1) then
            k=43
         elseif(ionel.eq.6.and.j.eq.1) then
            k=42
         else
            goto 333
         endif
         do l=1,nlines
            dphl=0.
            dphhl=0.
            if(jline(l).gt.jmin) then
               IF(IONel.EQ.4)  then
                  SIGM=SIGFE(J,jline(l))
               elseIF(IONel.EQ.6)  then
                  SIGM=SIGFEI(J,jline(l))
               elseIF(IONel.EQ.2)  then
                  SIGM=SIGO(J,jline(l))
               endif
               dphl=phline(l)*sigm
               dphhl=dphl*elch*(e1(jline(l))-e0)
               if(ionel.eq.5) then
c                  write(6,*)' dphl ',l,jline(l),sigm,phline(l),dphl
               endif
            endif
            ph(j)=ph(j)+dphl
           
           pheat(ionel,j)=pheat(ionel,j)+dphhl
         enddo
 333     continue
c     Include charge transfer ionization of Fe I and Fe II for first 5 
c     gs levels
         if(ionel.eq.6.and.j.le.5) ph(j)=phfe(1)
         if(ionel.eq.4.and.j.le.5) ph(j)=phfe(2)
         RTE(J)=RI
 500  CONTINUE

      RETURN
      END

      

      double precision function wd(r)
      implicit real*8(a-h,o-z)
      COMMON/QSOM/qso,IAGN,ISTAT,ISOBOL
      COMMON/SPH/ISPH
      COMMON/CICS/RSHOCK,ics,ics87a
      IF(ISPH.EQ.1) then
         wd=0.5*(1.-sqrt(1.-1./r**2))
      else
         wd=0.5
      endif
      wd=1.d-30
      return
      end


      SUBROUTINE ATDATn2(te)
      IMPLICIT REAL*8(A-H,O-Z)
      character*1 dum
      PARAMETER (NL=340,NLP1=NL+1)
      COMMON/A14/C(NL,NL),CI(NL),G(NL),E(NL),A(NL,NL),
     &     WL(NL,NL),DCDT(NL,NL),DCIDT(NL)
      COMMON/A27/OM(NL,NL)
      COMMON/NBA/NBACK
      COMMON/A16/SIG(NL),GAP(NL),BP(NL)
      COMMON/NLEV/e00,NION,NH,NP1H,NHMAX,NHMIN
      COMMON/HEIRECION/AL0(16),BE(16),B0(16),GA(16)
      common/collfit/omt(12),omfit(5,20,20,15),tefit(5,15)
      TEV=Te/1.1609E4
      
      ion=1
      n=12
      rewind (44)

      do i=1,n
         g(i)=0.
         e(i)=0.
         do j=1,n
            a(j,i)=0.
            c(i,j)=0.
         enddo
      enddo
      open(44,file='ATDAT/nii_lev_coll_rad.dat',status='old')
      read(44,987)  dum
      read(44,987)  dum
987   format(a)
      do i=1,n
         read(44,*)iq,jlev,wn
         g(i)=2.*jlev+1.
         e(i)=wn/8065.46d0
      enddo

      do i=1,1000
         read(44,*,err=11,end=11)il,iu,wlq,wll,wlu,gl,gu,a21
         a(iu,il)=a21
      enddo

 11   continue
      read(44,987)  dum
      nte=12
      read(44,*)(tefit(ion,k),k=1,nte)
      do i=1,1000
         read(44,*,err=12,end=12)ii,il,iu,(omt(k),k=1,nte)
         do k=1,nte
            omfit(ion,iu,il,k)=omt(k)
         enddo
      enddo
 12   continue

      telog=log10(te)
      do i=1,n
         ip1=i+1
         DO J=IP1,N
            wl(j,i)=12398.54/(e(j)-e(i))
            wl(i,j)=wl(j,i)
            do k=1,nte-1 
               if(telog.gt.tefit(ion,nte)) then
                  omint = omfit(ion,j,i,nte)
               elseif(telog.lt.tefit(ion,1)) then
                  omint = omfit(ion,j,i,1)
               elseif(telog.gt.tefit(ion,k).and.
     &                 telog.le.tefit(ion,k+1)) 
     &                 then
                  omint = omfit(ion,j,i,k) + (telog-tefit(ion,k))*
     &                 (omfit(ion,j,i,k+1)-omfit(ion,j,i,k))/
     &                 (tefit(ion,k+1)-tefit(ion,k))
               endif
            enddo 
            EIJ=ABS(E(J)-E(I))
            ET=EIJ/TEV
            IF(EIJ/TEV.GT.700.) ET=700.
            C(I,J)=8.63E-6*omint*EXP(-ET)/(G(I)*SQRT(Te))
            C(j,i)=G(i)*EXP(ET)*C(i,j)/G(j)
         enddo
      enddo
c      if(te.gt.980.and.te.lt.1005) then
c         do i=1,n
c            do j=1,n
c               write(6,9)i,j,wl(i,j),a(i,j),c(i,j),te
c 9             format(' niidata ',2i5,f15.2,1pe12.3,10e12.3)
c            enddo
c         enddo
c      endif
      return
      end


      SUBROUTINE ATDATar5(te)
      IMPLICIT REAL*8(A-H,O-Z)
      character*1 dum
      PARAMETER (NL=340,NLP1=NL+1)
      COMMON/A14/C(NL,NL),CI(NL),G(NL),E(NL),A(NL,NL),
     &     WL(NL,NL),DCDT(NL,NL),DCIDT(NL)
      COMMON/A27/OM(NL,NL)
      COMMON/NBA/NBACK
      COMMON/A16/SIG(NL),GAP(NL),BP(NL)
      COMMON/NLEV/e00,NION,NH,NP1H,NHMAX,NHMIN
      COMMON/HEIRECION/AL0(16),BE(16),B0(16),GA(16)
      common/collfit/omt(12),omfit(5,20,20,15),tefit(5,15)
      TEV=Te/1.1609E4
      
      ion=2
      n=5
      rewind (45)

      do i=1,n
         g(i)=0.
         e(i)=0.
         do j=1,n
            a(j,i)=0.
            c(i,j)=0.
         enddo
      enddo
      open(45,file='./ATDAT/ar5_lev_coll_rad.dat',status='old')
      read(45,987)  dum
      read(45,987)  dum
987   format(a)
      do i=1,n
         read(45,*)iq,jlev,wn
         g(i)=2.*jlev+1.
         e(i)=wn/8065.46d0
      enddo

      do i=1,1000
         read(45,*,err=11,end=11)il,iu,wlq,el,eu,gl,gu,a21
         a(iu,il)=a21
      enddo

 11   continue
      read(45,987)  dum
      nte=11
      read(45,*)(tefit(ion,k),k=1,nte)
      do i=1,1000
         read(45,*,err=12,end=12)iu,il,(omt(k),k=1,nte)
         do k=1,nte
            if(iu.eq.4.or.iu.eq.5.and.il.ne.4) then
c  fine structure splitting for 4th and 5th levels
               omfit(ion,iu,1,k)=1.*omt(k)/9.
               omfit(ion,iu,2,k)=3.*omt(k)/9.
               omfit(ion,iu,3,k)=5.*omt(k)/9.
            else
               omfit(ion,iu,il,k)=omt(k)
            endif
         enddo
      enddo
 12   continue
c      close(45)
      telog=log10(te)
      do i=1,n
         ip1=i+1
         DO J=IP1,N
            wl(j,i)=12398.54/(e(j)-e(i))
            wl(i,j)=wl(j,i)
            do k=1,nte-1 
               if(telog.gt.tefit(ion,nte)) then
                  omint = omfit(ion,j,i,nte)
               elseif(telog.lt.tefit(ion,1)) then
                  omint = omfit(ion,j,i,1)
               elseif(telog.gt.tefit(ion,k).and.
     &                 telog.le.tefit(ion,k+1)) 
     &                 then
                  omint = omfit(ion,j,i,k) + (telog-tefit(ion,k))*
     &                 (omfit(ion,j,i,k+1)-omfit(ion,j,i,k))/
     &                 (tefit(ion,k+1)-tefit(ion,k))
               endif
            enddo 
            EIJ=ABS(E(J)-E(I))
            ET=EIJ/TEV
            IF(EIJ/TEV.GT.700.) ET=700.
            C(I,J)=8.63E-6*omint*EXP(-ET)/(G(I)*SQRT(Te))
            C(j,i)=G(i)*EXP(ET)*C(i,j)/G(j)
c            write(6,9)i,j,wl(j,i),a(j,i),omint,c(i,j),c(j,i)
 9          format(' Ar 5 ',2i5,f15.2,1pe12.3,10e12.3)
         enddo
      enddo
c      do i=1,n
c         do j=1,n
c         enddo
c      enddo
      return
      end

      SUBROUTINE ATDATfe7(te)
      IMPLICIT REAL*8(A-H,O-Z)
      character*1 dum
      PARAMETER (NL=340,NLP1=NL+1)
      COMMON/A14/C(NL,NL),CI(NL),G(NL),E(NL),A(NL,NL),
     &     WL(NL,NL),DCDT(NL,NL),DCIDT(NL)
      COMMON/A27/OM(NL,NL)
      COMMON/NBA/NBACK
      COMMON/A16/SIG(NL),GAP(NL),BP(NL)
      COMMON/NLEV/e00,NION,NH,NP1H,NHMAX,NHMIN
      COMMON/HEIRECION/AL0(16),BE(16),B0(16),GA(16)
      common/collfit/omt(12),omfit(5,20,20,15),tefit(5,15)
      TEV=Te/1.1609E4
      
      ion=3
      n=9
      rewind (46)

      do i=1,n
         g(i)=0.
         e(i)=0.
         do j=1,n
            a(j,i)=0.
            c(i,j)=0.
         enddo
      enddo
      open(46,file='./ATDAT/fe7_lev_coll_rad.dat',status='old')
      do i=1,n
         read(46,*)iq,jlev,wn
         g(i)=2.*jlev+1.
         e(i)=wn/8065.46d0
      enddo
      read(46,987)  dum
      read(46,987)  dum
      read(46,987)  dum
987   format(a)
      do i=1,1000
         read(46,*,err=11,end=11)il,iu,wlq,wll,a21
         a(iu,il)=a21
      enddo

 11   continue
      read(46,987)  dum
      read(46,987)  dum
      read(46,987)  dum
      read(46,987)  dum
      read(46,987)  dum
      read(46,987)  dum
      nte=8
      read(46,*)(tefit(ion,k),k=1,nte)
      do i=1,1000
         read(46,*,err=12,end=12)il,iu,aqq,(omt(k),k=1,nte)
         do k=1,nte
            omfit(ion,iu,il,k)=omt(k)
         enddo
      enddo
 12   continue
      telog=log10(te)
      do i=1,n
         ip1=i+1
         DO J=IP1,N
            wl(j,i)=12398.54/(e(j)-e(i))
            wl(i,j)=wl(j,i)
            do k=1,nte-1 
               if(telog.gt.tefit(ion,nte)) then
                  omint = omfit(ion,j,i,nte)
               elseif(telog.lt.tefit(ion,1)) then
                  omint = omfit(ion,j,i,1)
               elseif(telog.gt.tefit(ion,k).and.
     &                 telog.le.tefit(ion,k+1)) 
     &                 then
                  omint = omfit(ion,j,i,k) + (telog-tefit(ion,k))*
     &                 (omfit(ion,j,i,k+1)-omfit(ion,j,i,k))/
     &                 (tefit(ion,k+1)-tefit(ion,k))
               endif
            enddo 
            EIJ=ABS(E(J)-E(I))
            ET=EIJ/TEV
            IF(EIJ/TEV.GT.700.) ET=700.
            C(I,J)=8.63E-6*omint*EXP(-ET)/(G(I)*SQRT(Te))
            C(j,i)=G(i)*EXP(ET)*C(i,j)/G(j)
c            write(6,9)i,j,wl(j,i),a(j,i),omint,c(i,j),c(j,i)
 9          format(' fe7 ',2i5,f15.2,1pe12.3,10e12.3)
         enddo
      enddo
      return
      end





      SUBROUTINE ATDATo_new
      IMPLICIT REAL*8(A-H,O-Z)
      character*8 dum
      PARAMETER (NL=340,NLP1=NL+1)
      COMMON/A14/C(NL,NL),CI(NL),G(NL),E(NL),A(NL,NL),
     &WL(NL,NL),DCDT(NL,NL),DCIDT(NL)
      COMMON/A27/OM(NL,NL)
      COMMON/NBA/NBACK
      COMMON/A16/SIG(NL),GAP(NL),BP(NL)
      COMMON/NLEV/e00,NION,NH,NP1H,NHMAX,NHMIN
      common/initox/initoi
      common/oicoll/omoi(13,13),alfaoi(13,13)
      dimension eoi(13),goi(14),wloi(13,13),aoi(13,13)
      save eoi,goi,wloi,aoi
      N=NHMAX
      N=13
      E00=13.618
      if(initoi.eq.1) then         
         initoi=0
         
         do i=1,n
            do j=1,n
               aoi(i,j)=0.
            enddo
         enddo

         rewind 25
         OPEN(25,FILE='./ATDAT/OI_kb.dat',status='old')
         read(25,987)dum
         read(25,987)dum
         read(25,987)dum
 987     format(a)
         do i=1,n
            read(25,*)nn,wn,goi(i)
c     read(25,*)nn,wn,g(i),al0(i),be(i),b0(i),ga(i)
            eoi(i)=wn/8065.46d0
         enddo


C     O II

         goi(N+1)=4.

         read(25,987)dum
         read(25,987)dum

         do i=1,n-1
            do j=i+1,n

               read(25,*)i1,i2,aoi(j,i),OM5,OM10
               if(om5.gt.0.) then
                  alfaoi(i,j) = log(om10/om5)/log(2.d0)
                  omoi(i,j) = om10
               else
                  omoi(i,j) = 1.d-20
                  alfaoi(i,j) = 0.
               endif
c     omoi(i,j) = om10*(te/1.e4)**alfa

            enddo
         enddo

      endif

      DO I=1,N
         e(i) = eoi(i)
         g(i) = goi(i)
      enddo


      DO 5395 I=1,N
         DO 5394 J=1,N
            WL(I,J)=0.0
            IF(I.EQ.J) GOTO 5394
            IF(E(I).EQ.E(J)) GOTO 5394
            WL(I,J)=ABS(12398.54/(E(I)-E(J)))
 5394    CONTINUE
 5395 CONTINUE

      DO  I=1,N
         DO J=1,N
            OM(I,J)=0.
            A(I,J)=aoi(i,j)            
         enddo
      enddo

C
      SIG(1)=7.91
      SIG(2)=15.
      DO I=3,N
         SIG(I)=7.91*REAL(I)
      enddo
      SIG(8)=20.5
      SIG(9)=20.5
      GAP(1)=2.99
      GAP(2)=2.5
      DO I=3,N
         GAP(I)=3.
      enddo
      GAP(8)=3.56
      GAP(9)=3.56
      DO I=1,N
         C(I,I)=0.
      enddo
      RETURN
      END



      SUBROUTINE ATDATO
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NL=340,NLP1=NL+1)
      COMMON/A14/C(NL,NL),CI(NL),G(NL),E(NL),A(NL,NL),
     &WL(NL,NL),DCDT(NL,NL),DCIDT(NL)
      COMMON/NBA/NBACK
      COMMON/A16/SIG(NL),GA(NL),BP(NL)
      COMMON/NLEV/e00,NION,NH,NP1H,NHMAX,NHMIN
      NH=9
      N=NHMAX
      N=NH
      E00=13.618
      E(1)=0.
      E(2)=1.961
      E(3)=4.172
      E(4)=9.137
      E(5)=9.508
      E(6)=10.732
      E(7)=10.976
      E(8)=12.073
      E(9)=12.523
      G(1)=9.
      G(2)=5.
      G(3)=1.
      G(4)=5.
      G(5)=3.
      G(6)=15.
      G(7)=9.
      G(8)=15.
      G(9)=15.
      G(N+1)=4.
      DO 5395 I=1,N
      DO 5394 J=1,N
      WL(I,J)=0.0
      IF(I.EQ.J) GOTO 5394
      IF(E(I).EQ.E(J)) GOTO 5394
      WL(I,J)=ABS(12398.54/(E(I)-E(J)))
5394  CONTINUE
5395  CONTINUE
      DO I=1,N
         DO J=1,N
            C(I,J)=0.
            A(I,J)=0.
         enddo
      enddo
c     kestner & bhatia
      A(2,1)=7.43E-3
      A(3,1)=7.63E-2
      A(4,1)=2.246e4
      A(5,1)=6.15E8
      A(8,1)=.748E8
c     Christenson and Cunningham, J. Geophys. Res. 83, 4393 -78
      A(9,1)=2.28E8
      A(3,2)=1.22E0
c 1641 
      a(5,2)=1.84e3
      a(5,3)=4.6
      A(5,4)=0.
      A(6,4)=0.356E8
      A(7,4)=4.00e2
      A(6,5)=6.47e1
      a(7,5)=3.35e7
      a(8,5)=1.39e2
      a(8,6)=1.17e2
      A(8,7)=.309E8
c     Christenson and Cunningham, J. Geophys. Res. 83, 4393 -78
      A(9,7)=9.59E4
7362  CONTINUE
C
      SIG(1)=7.91
      SIG(2)=15.
      DO I=3,N
         SIG(I)=7.91*REAL(I)
      enddo
      SIG(8)=20.5
      SIG(9)=20.5
      GA(1)=2.99
      GA(2)=2.5
      DO I=3,N
         GA(I)=3.
      enddo
      GA(8)=3.56
      GA(9)=3.56
      DO I=1,N
         C(I,I)=0.
      enddo
      RETURN
      END


      SUBROUTINE CROXY
      IMPLICIT REAL*8(A-H,O-Z)  
      PARAMETER (NL=340,NLP1=NL+1)
      include 'param'
c      PARAMETER (NE1=-100,NE2=250,NE3=NE2+1)
      COMMON/A26/SIGO(20,NE1:NE2),SIGCA(10,NE1:NE2),SIGHE(16,NE1:NE2),
     &SIGFE(NL,NE1:NE2),sigfei(nl,ne1:ne2),SIGH(50,NE1:NE2),SIGHE2(50,NE1:NE2)
      COMMON/FRE/NINT,JMIN,JJ
      COMMON/INT/FL(2,NE1:NE2),SI(14,27,NE1:NE2),E1(NE1:NE3),E(NE1:NE3)
      DO J1=JMIN,JJ
         EN=E1(J1)
C     O I (MCALPINE CORRECTED TO AGREE WITH REILMAN AND MANSON AT
C     HIGH ENERGIES).
         SIGO(1,J1)=SI(5,1,J1)
         DO  J=4,9
            S=0.
            S=SIGOX(J,EN)
            SIGO(J,J1)=S*1.E-18
         enddo
      enddo
      END

      DOUBLE PRECISION FUNCTION SIGOX(J,EN)
C
C     CROSSECTIONS OF O I FROM DAWIES AND LEWIS IN 1.E-18 CM
C
      IMPLICIT REAL*8(A-H,O-Z)
      IF(J.EQ.1) ERY=(EN-13.618)/13.618
      IF(J.EQ.4) ERY=(EN-4.481)/13.618
      IF(J.EQ.5) ERY=(EN-4.110)/13.618
      IF(J.EQ.6) ERY=(EN-2.886)/13.618
      IF(J.EQ.7) ERY=(EN-2.642)/13.618
      IF(J.EQ.8) ERY=(EN-1.545)/13.618
      IF(J.EQ.9) ERY=(EN-1.095)/13.618
      S=0.
      IF(ERY.LE.0.) GOTO 200
      IF(J.EQ.4.AND.ERY.LT..3488) SL=(-1.034*ERY+1.365-.0491/ERY)-2.
      IF(J.EQ.4.AND.ERY.GT..3488) SL=-0.5481-0.6746*(ERY/.4)
      IF(J.EQ.5.AND.ERY.LT..5) SL=(-.4038*ERY+1.783-.0449/ERY)-2.
      IF(J.EQ.5.AND.ERY.GT..5) SL=-.483*ERY-.2672
      S=10.**SL
      IF(J.EQ.6) S=10.**(-7.204*ERY+1.657-1.)+10.**(-.3748*ERY-.5877)
      IF(J.EQ.7.and.ery.gt.0.25) S=10.**(-.54*ERY-.1715)
      IF(J.EQ.7.and.ery.lt.0.25) S=10.**(-.54*0.25-.1715)
      IF(J.EQ.7.and.ery.lt.0.12) S=10.**(-6.*ERY+.28)
      IF(J.EQ.8) S=20.5*(1.545/EN)**3.56
      IF(J.EQ.9) S=20.5*(1.545/EN)**3.56
      IF(J.EQ.1) THEN
            S=9.05*(4.378*(13.618/EN)**1.5-3.378*(13.618/EN)**2.5)
            S=2.94*(2.66*(13.618/EN)**1.-1.66*(13.618/EN)**2.)
            IF(EN.GT.100.) 
     &             S=1.778*(3.65*(100./EN)**3-2.65*(100./EN)**4)
      ENDIF
200   SIGOX=S
      RETURN
      END


      SUBROUTINE CRFEI
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NL=340,NLP1=NL+1)
      include 'param'
c      PARAMETER (NE1=-100,NE2=250,NE3=NE2+1)
      COMMON/A14/C(NL,NL),CI(NL),G(NL),EN(NL),AS(NL,NL),
     &WL(NL,NL),DCDT(NL,NL),DCIDT(NL)
      COMMON/NLEV/e00,NION,N,NP1H,NMAX,NMIN
      COMMON/A26/SIGO(20,NE1:NE2),SIGCA(10,NE1:NE2),SIGHE(16,NE1:NE2),
     &SIGFE(NL,NE1:NE2),sigfei(nl,ne1:ne2),SIGH(50,NE1:NE2),SIGHE2(50,NE1:NE2)
      COMMON/FRE/NINT,JMIN,JJ
      COMMON/INT/FL(2,NE1:NE2),SI(14,27,NE1:NE2),E1(NE1:NE3),E(NE1:NE3)
      COMMON/SIK/SK(14,27,NE1:NE2)
      DO J1=JMIN,JJ
c  same crossection to all sublevels of the ground state
         do j=1,36
            ethresh=e00-en(j)
            if(e(j1).gt.ethresh) then
               SIGFEI(j,J1)=SI(14,1,J1)
            endif
         enddo
         DO J=37,121
            S=0.
            SIGFEI(J,J1)=S*1.E-18
         enddo
      enddo
      RETURN
      END

      SUBROUTINE CRFEII
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NL=340,NLP1=NL+1)
      include 'param'
c      PARAMETER (NE1=-100,NE2=250,NE3=NE2+1)
      COMMON/A14/C(NL,NL),CI(NL),G(NL),EN(NL),AS(NL,NL),
     &WL(NL,NL),DCDT(NL,NL),DCIDT(NL)
      COMMON/NLEV/e00,NION,N,NP1H,NMAX,NMIN
      COMMON/A26/SIGO(20,NE1:NE2),SIGCA(10,NE1:NE2),SIGHE(16,NE1:NE2),
     &SIGFE(NL,NE1:NE2),sigfei(nl,ne1:ne2),SIGH(50,NE1:NE2),SIGHE2(50,NE1:NE2)
      COMMON/FRE/NINT,JMIN,JJ
      COMMON/INT/FL(2,NE1:NE2),SI(14,27,NE1:NE2),E1(NE1:NE3),E(NE1:NE3)
      COMMON/SIK/SK(14,27,NE1:NE2)
      DO J1=JMIN,JJ
c  same crossection to all sublevels of the ground state
         do j=1,36
            sigfe(j,j1)=0.
            ethresh=e00-en(j)
            if(e(j1).gt.ethresh) then
               SIGFE(j,J1)=SI(14,2,J1)
            endif
         enddo
         DO J=37,191
            S=0.
            SIGFE(J,J1)=S*1.E-18
         enddo
      enddo
      RETURN
      END

      SUBROUTINE HPCA(R,TE)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NL=340,NLP1=NL+1)
      include 'param'
c      PARAMETER (NE1=-100,NE2=250,NE3=NE2+1)
C
C     CALCULATE RECOMBINATION AND PHOTOIONIZATION RATES FOR CALCIUM
C
      COMMON/IND/I
      COMMON/A14/C(NL,NL),CI(NL),G(NL),EN(NL),A(NL,NL),
     &WL(NL,NL),DCDT(NL,NL),DCIDT(NL)
      COMMON/A17/RECNET(NL),RECNDT(NL),RECCO(NL),RTE(NL),
     &PH(NL),DRDT(NL),RECT(6,NL),PHET(NL)
      COMMON/NLEV/e00,NION,N,NP1H,NHMAX,NHMIN
      COMMON/PHEAT/PHE(NL)
      COMMON/A26/SIGO(20,NE1:NE2),SIGCA(10,NE1:NE2),SIGHE(16,NE1:NE2),
     &SIGFE(NL,NE1:NE2),sigfei(nl,ne1:ne2),SIGH(50,NE1:NE2),SIGHE2(50,NE1:NE2)
      COMMON/FRE/NINT,JMIN,JJ
      COMMON/INT/FL(2,NE1:NE2),SI(14,27,NE1:NE2),E1(NE1:NE3),E(NE1:NE3)
      COMMON/DTAU/FLU(NE1:NE2)
      DATA PI/3.1415926E0/,ELCH/1.60219E-12/,AMU/1.660531E-24/
      TS=TE
      DO J=1,N
C     
C     FIRST CALCULATE THE RECOMBINATION RATE FOR ALL LEVELS
C     
C     CALL CARECQ(J,TS,RI)
C     
C     THEN THE IONIZATION AND HEATING RATE
C     
         E0=E00-EN(J)
         ZA=0.
         HEAT=0.
         DO J1=JMIN,JJ
            IF(E1(J1).gT.E0) then
               SIGM=SIGCA(J,J1)
               if(flu(j1).lt.0.) WRITE(6,9292)i,J,J1,FLU(J1),E1(J1),SIGM
 9292          FORMAT('CA le0 ',3I5,11E11.4)
               ZC=4.*PI*FLU(J1)*SIGM*(E(J1+1)-E(J1))/(ELCH*E1(J1))
               ZA=ZA+ZC
               HEAT=ELCH*ZC*(E1(J1)-E0)+HEAT
            endif
         enddo
         PH(J)=ZA
C     RTE(J)=RI
      enddo
      RETURN
      END

      SUBROUTINE CRCA
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NL=340,NLP1=NL+1)
      include 'param'
c      PARAMETER (NE1=-100,NE2=250,NE3=NE2+1)
      COMMON/A26/SIGO(20,NE1:NE2),SIGCA(10,NE1:NE2),SIGHE(16,NE1:NE2),
     &SIGFE(NL,NE1:NE2),sigfei(nl,ne1:ne2),SIGH(50,NE1:NE2),SIGHE2(50,NE1:NE2)
      COMMON/A14/C(NL,NL),CI(NL),G(NL),EN(NL),A(NL,NL),
     &WL(NL,NL),DCDT(NL,NL),DCIDT(NL)
      COMMON/NLEV/e00,NION,N,NP1H,NHMAX,NHMIN
      COMMON/FRE/NINT,JMIN,JJ
      COMMON/INT/FL(2,NE1:NE2),SI(14,27,NE1:NE2),E1(NE1:NE3),E(NE1:NE3)
      DO J1=JMIN,JJ
         DO J=1,3
            SIGCA(J,J1)=1.E-18*CRPHCA(J,E1(J1))
         enddo
         SIGCA(1,J1)=SI(13,2,J1)
      enddo
      RETURN
      END

      DOUBLE PRECISION FUNCTION CRPHCA(J,E1)
C
C     CROSSECTIONS OF CA II FROM SHINE
C
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NL=340,NLP1=NL+1)
      COMMON/NLEV/e00,NION,N,NP1H,NHMAX,NHMIN
      COMMON/A14/C(NL,NL),CI(NL),G(NL),EN(NL),A(NL,NL),
     &WL(NL,NL),DCDT(NL,NL),DCIDT(NL)
      E0=E00-EN(J)
      S=0.
      IF(E1.LE.E0) GOTO 200
      IF(J.EQ.1) THEN
C     APPROX. FIT TO REILMAN AND MANSON
        IF(E1.LT.35.) S=0.20*(E0/E1)**0.0001
        IF(E1.GE.35..AND.E1.LT.45.) S=0.10*(35./E1)**1.
        IF(E1.GE.45..AND.E1.LT.160.) S=1.5
        IF(E1.GE.160..AND.E1.LT.390.) S=1.15*(160./E1)**1.66
        IF(E1.GE.390..AND.E1.LT.4000.) S=2.5*(390./E1)**2.18
        IF(E1.GT.4000.) S=6.85E-2*(4000./E1)**2.43
      ELSE
      IF(J.EQ.2.or.j.eq.3) THEN
         S=6.15*(E0/E1)**1.
C!!! skip Lyman alpha ionization of Ca II
c         IF(E1.LT.10.5) S=0.
      ENDIF
      IF(J.EQ.4.or.j.eq.5) S=2.38*(E0/E1)**3.65
      ENDIF
 200  CRPHCA=S
      RETURN
      END

      SUBROUTINE ATDAT
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NL=340,NLP1=NL+1)
      COMMON/A14/C(NL,NL),CI(NL),G(NL),E(NL),A(NL,NL),
     &WL(NL,NL),DCDT(NL,NL),DCIDT(NL)
      COMMON/NBA/NBACK
      COMMON/A16/SIG(NL),GA(NL),BP(NL)
      COMMON/NLEV/e00,NION,NH,NP1H,NHMAX,NHMIN
      N=NHMAX
      E00=11.871
      E(1)=0.
      E(2)=7323.9
      E(3)=7291.5
      e(4)=3968.5
      e(5)=3933.7
      do i=2,5
         e(i)=12398.54/e(i)
      enddo
      G(1)=2.
      G(2)=4.
      G(3)=6.
      g(4)=2.
      g(5)=4.
c     CA III 1S0
      G(N+1)=1.
      DO  I=1,N
         DO  J=1,N
            WL(I,J)=0.0
            IF(I.ne.J) then
               IF(E(I).ne.E(J)) then
                  WL(I,J)=ABS(12398.54/(E(I)-E(J)))
               endif
            endif
         enddo
      enddo
      DO I=1,N
         DO  J=1,N
            C(I,J)=0.
            A(I,J)=0.
         enddo
      enddo
      A(2,1)=1.3
      a(3,1)=1.3
      A(4,1)=1.46e8
      a(5,1)=1.50e8
      A(3,2)=0.
      a(4,2)=1.06e7
      a(5,2)=1.11e6
      a(4,3)=0.
      a(5,3)=9.9e6
      a(5,4)=0.
7362  CONTINUE
C
      SIG(1)=0.20
      SIG(2)=6.15
      sig(3)=sig(2)
      SIG(4)=2.38
      sig(5)=sig(4)
      GA(1)=0.
      GA(2)=1.0
      ga(3)=ga(2)
      GA(4)=3.65
      ga(5)=ga(4)
      DO 165 I=1,N
165   C(I,I)=0.
      RETURN
      END

      SUBROUTINE ATDATFE
      IMPLICIT REAL*8(A-H,O-Z)
      character*1 dum
      character*2 ch2
      character*4 ch4
      character*5 ch5
      SAVE
      PARAMETER (NL=340,NLP1=NL+1)
      COMMON/A14/C(NL,NL),CI(NL),G(NL),E(NL),A(NL,NL),
     &WL(NL,NL),DCDT(NL,NL),DCIDT(NL)
      COMMON/A27/OM(NL,NL)
      COMMON/NBA/NBACK
      COMMON/A16/SIG(NL),GA(NL),BP(NL)
      COMMON/NLEV/e00,NION,NH,NP1H,NHMAX,NHMIN
      COMMON/FELEV/NFEII,nfei
	common/initat/initfeii,initfei
      DIMENSION GS(NL),ES(NL),AS(NL,NL),WLS(NL,NL),OMS(NL,NL)
      ntot=191
      nfeii=191
      N=NFEII
      E00=16.19791
      IF(INITFEii.EQ.0) THEN
      INITFEii=1
      OPEN(23,FILE='./ATDAT/feIInew2.dat',status='old')
      rewind 23
      do i=1,10
         read(23,987)dum
      enddo
987   format(a)
      do i=1,ntot
         read(23,935)nn,wn,gs(i),ch2,ch4,ig
 935     format(i3,f13.3,f6.1,1x,a2,1x,a4,i4)
         es(i)=wn/8065.46d0
      enddo
      DO 5395 I=1,N
      DO 5394 J=1,N
      WLS(I,J)=0.0
      IF(I.EQ.J) GOTO 5394
      IF(ES(I).EQ.ES(J)) GOTO 5394
      WLS(I,J)=ABS(12398.54/(ES(I)-ES(J)))
5394  CONTINUE
5395  CONTINUE
      DO I=1,N
         DO  J=1,N
            OMS(I,J)=0.
            AS(I,J)=0.
         enddo
      enddo
c      read(23,987)dum
      read(23,987)dum
      read(23,987)dum
      do i=2,nfeii
         do j=1,i-1
            read(23,937)i1,i2,ch5,g2,ch5,g1,wlq,as(i,j),OMs(J,I)
c     read(23,937)i1,i2,ch5
c     937  format(2i4,2x,a5,i2,2x,a5,i2,2x,f12.2,1pe12.4,e12.4)
 937        format(2i4,2x,a5,f2.0,2x,a5,f2.0,2x,f12.2,1pe12.4,e12.4)
c     1234123412123451212123451212123456789012123456789012123456789012
c     2   1  a6De  8  a6De 10     259895.52  2.1300E-03  5.7280E+00
         enddo
      enddo
      close (23)
      do iu=2,191
         do il=1,iu-1
c            write(62,955)iu,il,wls(iu,il),as(iu,il),oms(il,iu)
 955        format(2i4,f10.1,1pe12.4,e12.4)
            enddo
         enddo
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
      DO I=3,N
         SIG(I)=0.
      enddo
      GA(1)=2.99
      GA(2)=2.5
      DO I=3,N
         GA(I)=3.
      enddo
      GA(8)=3.56
      GA(9)=3.56
      DO I=1,N
         C(I,I)=0.
      enddo
      RETURN
      END

      SUBROUTINE ATDATFEI
      IMPLICIT REAL*8(A-H,O-Z)
      character*80 dum
      SAVE
      PARAMETER (NL=340,NLP1=NL+1)
      COMMON/A14/C(NL,NL),CI(NL),G(NL),E(NL),A(NL,NL),
     &WL(NL,NL),DCDT(NL,NL),DCIDT(NL)
      COMMON/A27/OM(NL,NL)
      COMMON/NBA/NBACK
      COMMON/A16/SIG(NL),GA(NL),BP(NL)
      COMMON/NLEV/e00,NION,NH,NP1H,NHMAX,NHMIN
      COMMON/FELEV/NFEII,nfei
	common/initat/initfeii,initfei
      DIMENSION GS(NL),ES(NL),AS(NL,NL),WLS(NL,NL),OMS(NL,NL)
c??
      NTOT=121
      N=NFEI
      E00=7.870
      IF(INITFEi.EQ.0) THEN
      INITFEi=1
      OPEN(23,FILE='./ATDAT/feI.dat',status='old')
      rewind 23
      read(23,987)dum
      read(23,987)dum
987   format(a)
      do 1 i=1,ntot
      read(23,*)nn,wn,gs(i)
1     es(i)=wn/8065.46d0
      DO 5395 I=1,N
      DO 5394 J=1,N
      WLS(I,J)=0.0
      IF(I.EQ.J) GOTO 5394
      IF(ES(I).EQ.ES(J)) GOTO 5394
      WLS(I,J)=ABS(12398.54/(ES(I)-ES(J)))
5394  CONTINUE
5395  CONTINUE
      DO I=1,N
         DO J=1,N
            OMS(I,J)=0.
            AS(I,J)=0.
         enddo
      enddo
      read(23,987)dum
      read(23,987)dum
      do i=2,nfei
         do j=1,i-1
            read(23,*)i1,i2,wlq,as(i,j),OMs(J,I)
         enddo
      enddo
	close (23)
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
      DO I=3,N
         SIG(I)=0.
      enddo
      GA(1)=2.99
      GA(2)=2.5
      DO I=3,N
         GA(I)=3.
      enddo
      GA(8)=3.56
      GA(9)=3.56
      DO I=1,N
         C(I,I)=0.
      enddo
      RETURN
      END



      SUBROUTINE CION(ION,N,N0,T,CI)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NL=340,NLP1=NL+1)
C
C     COLLISIONAL IONIZATION RATES
C
      COMMON/A14/CQ(NL,NL),CIQ(NL),G(NL),E(NL),AQ(NL,NL),WL(NL,NL)
     &,DCDT(NL,NL),DCIDT(NL)
       COMMON/NLEV/e00,NION,NHY,NP1H,NHMAX,NHMIN
      CI=0.
      IF(ION.NE.1) GOTO 3
C     COLLISIONAL IONIZATION RATES FOR CA II FROM SHINE
      TEV=T/1.1609E4
      EI=E00-E(N)
      EITEV=DMIN1(700.D0,EI/TEV)
      CI=1.14E-10*SQRT(T)*EXP(-EITEV)/EI**2
      GOTO 500
C
C     COLLISIONAL IONIZATION RATES FOR O I FROM SUMMERS
C
 3    IF(ION.eq.2) then
         T4=T/1.E4
         ci=0.
         IF(N.eq.1) then
            IF(T.gT.2.E3) then
               CI=1.09E-10*SQRT(T)*EXP(-1.58E5/T)/(1.+0.1*T/1.58E5)
            else
               ci=0.
            endif
         else
            IF(N.EQ.2) OM=1.68*T4**.956
            IF(N.EQ.3) OM=.496*T4**.952
            IF(N.EQ.4) OM=2.54*T4**.920
            IF(N.EQ.5) OM=1.80*T4**.915
            IF(N.EQ.6) OM=17.4*T4**.892
            IF(N.EQ.7) OM=12.3*T4**.885
            IF(N.EQ.8) OM=54.6*T4**.840
            IF(N.EQ.9) OM=101.*T4**.806
            EITEV=DMIN1(700.D0,(-E(N)+E00)*11609./T)
            CI=8.63E-6*OM*EXP(-EITEV)/(SQRT(T)*G(N))
         endif
      endif
      IF(ION.eq.4) then
         CALL CIONFE(N,E(N),T,CI)
      endif
 500  continue
      RETURN
      END

      SUBROUTINE CEX(ION,T,C)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NL=340,NLP1=NL+1)
C
C     COLLISIONAL EXCITATION RATES 
C
      COMMON/NLEV/e00,NION,N,NP1H,NMAX,NMIN
      COMMON/A14/CQ(NL,NL),CI(NL),G(NL),E(NL),AQ(NL,NL),WL(NL,NL)
     &,DCDT(NL,NL),DCIDT(NL)
      COMMON/A27/OM(NL,NL)
      COMMON/FELEV/NFEII,nfei
      common/oicoll/omoi(13,13),alfaoi(13,13)
      DIMENSION C(NL,NL)
      T4=T/1.E4
      T3=T/1.E3
      TEV=T/1.1609E4
      IF(ION.NE.1) GOTO 3
C     COLLISIONAL EXCITATION RATES FOR CA II FROM SHINE as given
c     by Li & McCray 1993
      c(1,2)=5.67*t3**(-0.01)
      c(1,3)=8.51*t3**(-0.01)
      c(2,3)=20.60*t3**0.01
      c(1,4)=3.99*t3**0.12
      c(2,4)=1.68*t3**0.23
      c(3,4)=3.93*t3**0.01
      c(1,5)=7.99*t3**0.12
      c(2,5)=6.27*t3**0.05
      c(3,5)=6.49*t3**0.29
      c(4,5)=10.08*t3**(-0.01)
3     IF(ION.NE.2) GOTO 4
c coll. excit. for O I
c 1 - 3 pradhan compil.
      IF(T4.GT.0.5) GOTO 10
      C(1,2)=.0151*T3**1.31
      C(1,3)=.00184*T3**1.32
      C(2,3)=.031*T3**.534
      GOTO 50
 10   IF(T4.GT.1.) GOTO 20
      C(1,2)=.275*T4**1.10
      C(1,3)=.034*T4**1.08
      C(2,3)=.105*T4**.52
      GOTO 50
 20   C(1,2)=.266*T4**.91
      C(1,3)=.0324*T4**.91
      C(2,3)=.105*T4**.50
c bathia & kastner
 50   C(1,4)=4.16e-2*T4**1.88
c bathia & kastner
      C(1,5)=.345*T4**.61
      C(1,6)=.046*T4**.82
      C(1,7)=.041*T4**.50
      C(1,8)=.0017*T4**2.16
      C(1,9)=.0088*T4**2.21
      C(4,6)=18.47*TEV**.64
      C(4,7)=2.5E-2*TEV**.64
      C(5,7)=9.63*TEV**.64
      C(7,8)=20.7*TEV**.64
      C(7,9)=5.3E-2*TEV**.64



C     MEWE FORMULA
      C(1,4)=1.28
      C(1,5)=1.59
      C(1,6)=2.69
      C(1,7)=2.76
      C(1,8)=.51
      C(1,9)=1.34
      C(4,6)=148.*TEV**.64
      C(5,7)=94.*TEV**.64
      C(7,8)=447.*TEV**.64
      C(7,9)=93.*TEV**.64

      do i=1,12
         do j=i,13
            c(i,j) = omoi(i,j)*t4**alfaoi(i,j)
         enddo
      enddo

4     IF(ION.NE.3) GOTO 5
C     HE I 
      DO I=1,12
         DO J=I+1,12
            C(I,J)=OM(I,J)
         enddo
      enddo
5     IF(ION.NE.4) GOTO 7
C     FE II
      DO I=1,NFEII
         DO  J=I+1,NFEII
           C(I,J)=OM(I,J)
         enddo
      enddo
7     IF(ION.NE.6) GOTO 32
C     FE I
      DO  I=1,NFEI
	DO  J=I+1,NFEI	
           c(I,J)=OM(I,J)
        enddo
      enddo
 32   continue
9     CONTINUE
      DO I=1,N
         IP1=I+1
         DO J=IP1,N
            EIJ=ABS(E(J)-E(I))
            TEV=T/1.1609E4
            ET=EIJ/TEV
            IF(EIJ/TEV.GT.700.) ET=700.
            C(I,J)=8.63E-6*C(I,J)*EXP(-ET)/(G(I)*SQRT(T))
         enddo
      enddo
15    CONTINUE
      RETURN
      END


      SUBROUTINE CIONFE(N,EN,T,CI)
      IMPLICIT REAL*8(A-H,O-Z)
C
C     COLLISIONAL IONIZATION RATES FOR FE II
C
c!! ?
      a1=10.
c!! ?
      eion=16.1
      CON=5.465E-11
      Q=0.
      T4=T/1.E4
      YN=1.1609*EION/T4
      CI=CON*SQRT(T)*Q*A1
      RETURN
      END

      SUBROUTINE CIONFEI(N,EN,T,CI)
      IMPLICIT REAL*8(A-H,O-Z)
C
C     COLLISIONAL IONIZATION RATES FOR FE I
C
c!! ?
      a1=10.
c!! ?
      eion=7.870
      CON=5.465E-11
      Q=0.
      T4=T/1.E4
      YN=1.1609*EION/T4
      CI=CON*SQRT(T)*Q*A1
      RETURN
      END



      SUBROUTINE OXREC(OXR,TE)
C     ************************************************************
C     ******
C     O I RECOMBINATION EMISSION CALCULATED FROM JULIENNE ET AL
C     AND SCALED TO 1.E4 K BY A POWER T4**-0.7
C     WAVELENGTHS IN ORDER 1304, 8446, 11287, 7002, 1356, 7774,
C                          9264, 6157
C     ******
C     *************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION OXR(10)
      T4=TE/1.E4
      OXR(1)=9.44E-25/T4**.7
      OXR(2)=.143
      OXR(3)=.070
      OXR(4)=.0133
      OXR(5)=1.68
      OXR(6)=0.288
      OXR(7)=0.135
      OXR(8)=.033
      DO I=2,8
         OXR(I)=OXR(I)*OXR(1)
      enddo
      RETURN
      END

      SUBROUTINE RECOMB(ION,I,TE,A,RI)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NL=340,NLP1=NL+1)

      parameter (mh=2,mhe=3,mc=7,mn=8,mo=9,mne=11,mna=12,mmg=13,
     &     mal=14,msi=15,ms=17,mar=19,mca=21,mfe=27,mni=29,
     &     mmat=mh+mhe+mc+mn+mo+mne+mna+mmg+mal+msi+ms+mar+mca+mfe+mni,
     &     ispecy=15)
      common/abc1/alc(mc),aln(mn),alo(mo),alne(mne),alna(mna),
     &     almg(mmg),alal(mal),alsi(msi),alsu(ms),alar(mar),alca(mca),
     &     alfe(mfe),alni(mni)

      COMMON/ABC/AL(16)

      COMMON/NLEV/e00,NION,N,NP1H,NHMAX,NHMIN
      COMMON/A14/C(NL,NL),CI(NL),G(NL),E(NL),AS(NL,NL),
     &WL(NL,NL),DCDT(NL,NL),DCIDT(NL)
      COMMON/HEIRECION/AL0(16),BE(16),B0(16),GA(16)
      DATA CSAHA/2.0708E-16/
      A=0.
      TEV=TE/1.1609E4
      TEV0=TEV*10.
      T4=TE/1.E4
      IF(ION.NE.2) GOTO 1
C     O I
c      IF(I.EQ.1) A=1.36E-13/TEV**.7
c      IF(I.EQ.4) A=1.6E-14/TEV0**.7
c      IF(I.EQ.5) A=.14E-13/TEV0**.7
c      IF(I.EQ.6) A=4.77E-13/TEV0**.7
c      IF(I.EQ.7) A=0.47E-13/TEV0**.7
c      IF(I.EQ.8) A=1.68E-13/TEV0**.7
c      IF(I.EQ.9) A=0.52E-13/TEV0**.7

c rec. rates from Julienne, P.S., Davies, J. & Oran, E. 1974, 
c J. Geophys. Res. 79, 2540  
c 3P divided according to statistical weights. Total 4.4e-13
      IF(I.EQ.1) A=2.44
      IF(I.EQ.2) A=1.466
      IF(I.EQ.3) A=0.4888
      IF(I.EQ.4) A=0.
      IF(I.EQ.5) A=0.
      IF(I.EQ.6) A=0.10
      IF(I.EQ.7) A=0.04
      IF(I.EQ.8) A=2.1
      IF(I.EQ.9) A=0.11
      IF(I.EQ.10) A=0.
      IF(I.EQ.11) A=0.
      IF(I.EQ.12) A=2.7
      IF(I.EQ.13) A=1.2

c correction to give same total rate as Chung, Lin & Lee 
c + same temp. dep.      
c     temp in 0.1 eV
      t01ev=te/1164.
      a=(1.47/1.065)*a*1.e-13/t01ev**0.8433


      GOTO 2
1     IF(ION.NE.3) GOTO 2
C     HE I
      T5000=TE/5.E3
      IF(I.EQ.1) A=2.23E-13/T5000**.488
      IF(I.EQ.2) A=1.98E-14/T5000**.42
      IF(I.EQ.3) A=3.6E-14/T5000**.784
      IF(I.EQ.4) A=1.135E-13/T5000**.610
      IF(I.EQ.5) A=7.2E-14/T5000**.784
      IF(I.EQ.6) A=8.45E-15/T5000**.45
      IF(I.EQ.7) A=3.92E-14/T5000**.63
      IF(I.EQ.8) A=4.894E-14/T5000**.85
      IF(I.EQ.9) A=0.343E-14/T5000**.478
      IF(I.EQ.10) A=1.77E-14/T5000**.620
      IF(I.EQ.11) A=3.14E-14/T5000**.851
      IF(I.EQ.12) A=4.87e-14/T5000**1.173
C     FROM ALMOG & NETZER
      A=AL0(I)*T4**BE(I)
2     IF(ION.NE.4) GOTO 32
C     FE II 
c  all recombinations to the a6D ground state
      IF(I.EQ.1) A=10.*alfe(3)/30.
      IF(I.EQ.2) A=8.*alfe(3)/30.  
      IF(I.EQ.3) A=6.*alfe(3)/30.  
      IF(I.EQ.4) A=4.*alfe(3)/30.  
      IF(I.EQ.5) A=2.*alfe(3)/30.  
      IF(I.ge.6) A=0.
32     IF(ION.NE.6) GOTO 3
C     FE I 
c all recombinations to the a5D ground state
      IF(I.EQ.1) A=9.*alfe(2)/25.
      IF(I.EQ.2) A=7.*alfe(2)/25.  
      IF(I.EQ.3) A=5.*alfe(2)/25.  
      IF(I.EQ.4) A=3.*alfe(2)/25.  
      IF(I.EQ.5) A=1.*alfe(2)/25.  
      IF(I.ge.6) A=0.
c!!!! stryk!
c      IF(I.EQ.1) A=alfe(2)
c      if(i.ge.2) a=0.
3     RI=A*EXP((E(I)-E00)/TEV)*TE**1.5/(CSAHA*G(I))
      RETURN
      END

      subroutine readferec
      implicit real*8 (a-h,o-z)
      common/cferec_f1/tef1(81),rectf1(81)
      common/cferecl_f1/rcmpf1(13,81),rcsumf1(10,81)
      common/cferec_f2/tef2(81),rectf2(81),rcmpf2(42,81),rcsumf2(17,81)
      dimension tv(81)
c
c Read total Fe I recombination coefficients
      open(22,file='./ATDAT/rec_tot_fe1.dat',status='old')
      do i=1,4
         read(22,*) 
      enddo
      do i=1,81
         read(22,*) tef1(i),rectf1(i)
      enddo
      close(22)
c Read total Fe II recombination coefficients
      open(22,file='./ATDAT/rec_tot_fe2.dat',status='old')
      do i=1,4
         read(22,*) 
      enddo
      do i=1,81
         read(22,*) tef2(i),rectf2(i)
      enddo
      close(22)
c Read Fe I recombination coefficients to specific multiplets
      open(22,file='./ATDAT/rec_lev_fe1.dat',status='old')
      read(22,*)
      read(22,701) (tv(i),i=1,81)
c
      do imp=1,13
         read(22,*)
         read(22,701) (rcmpf1(imp,it),it=1,81)
      enddo
c Read Fe I sum of recombination coefficients to higher multiplets
      read(22,*)
      do imps=1,10
         read(22,702) (rcsumf1(imps,it),it=1,81)
      enddo
 701  format(19x,85e10.2)
 702  format(11x,85e10.2)
      close(22)
c Read Fe II recombination coefficients to specific multiplets
      open(22,file='./ATDAT/rec_lev_fe2.dat',status='old')
      read(22,*)
      read(22,701) (tv(i),i=1,81)
      do imp=1,42
         read(22,*)
         read(22,701) (rcmpf2(imp,it),it=1,81)
      enddo
c Read Fe II sum of recombination coefficients to higher multiplets
      read(22,*)
      do imps=1,17
         read(22,702) (rcsumf2(imps,it),it=1,81)
      enddo
      close(22)
c
      return
      end

      subroutine ferecomb(ife,te)
      implicit real*8 (a-h,o-z)
c
      integer mie,mion
      integer nlf1,nlf2,nlf3,nlf4
      parameter(nlf1=122,nlf2=192,nlf3=113,nlf4=46)
      parameter(mie=15,mion=5)
c
      data ieh/1/,iehe/2/,iec/3/,ien/4/,ieo/5/,iene/6/,iena/7/,
     &     iemg/8/,iesi/9/,ies/10/,iear/11/,ieca/12/,iefe/13/,
     &     ieco/14/,ieni/15/,iemax/15/
c
      COMMON/NLEV/e00,NION,N,NP1H,NHMAX,NHMIN
      common/cnlev/nlev(mie,mion)
      common/cferec_f1/tef1(81),rectf1(81)
      common/cferecl_f1/rcmpf1(13,81),rcsumf1(10,81)
      common/cferec_f2/tef2(81),rectf2(81),rcmpf2(42,81),rcsumf2(17,81)
      common/cfedat/gf1(nlf1),gf2(nlf2),gf3(nlf3),gf4(nlf4)
      common/ferecombi/recof1(nlf1),recof2(nlf2)
      write(6,*)'in ferecomb ',ife,iefe,n
      if(ife==1) then
c     
c     Fe I
         do il=1,nlf1
            recof1(il)=0.d0
         enddo
         if (te.lt.tef1(1)) then
            tconst=0.d0
            it1=1
            it2=1
         elseif (te.ge.tef1(81)) then
            tconst=0.d0
            it1=81
            it2=81
         else
            do i=1,80
               if (te.ge.tef1(i).and.te.lt.tef1(i+1)) then
                  it1=i
                  it2=i+1
                  goto 9
               endif
            enddo
 9          continue
            tconst=(te-tef1(it1))/(tef1(it2)-tef1(it1))
         endif
c     
         rctfe1=rectf1(it1)+tconst*(rectf1(it2)-rectf1(it1))
c     
c     a 5De ; imp=1, il=1(g=9),2(7),3(5),4(3),5(1); gsum=25
c     Add 5De-sum; imps=2
         rcmp=rcmpf1(1,it1)+tconst*(rcmpf1(1,it2)-rcmpf1(1,it1))
         rcsum=rcsumf1(2,it1)+tconst*(rcsumf1(2,it2)-rcsumf1(2,it1))
         rcmp=rcmp+rcsum
         recof1(1)=rcmp*0.36
         recof1(2)=rcmp*0.28
         recof1(3)=rcmp*0.20
         recof1(4)=rcmp*0.12
         recof1(5)=rcmp*0.04
c     a 5Fe ; imp=2, il=6(11),7(9),8(7),9(5),10(3); gsum=35
c     Add 5Fe-sum; imps=3
         rcmp=rcmpf1(2,it1)+tconst*(rcmpf1(2,it2)-rcmpf1(2,it1))
         rcsum=rcsumf1(3,it1)+tconst*(rcsumf1(3,it2)-rcsumf1(3,it1))
         rcmp=rcmp+rcsum
         recof1(6)=rcmp*0.3143
         recof1(7)=rcmp*0.2571
         recof1(8)=rcmp*0.2
         recof1(9)=rcmp*0.1429
         recof1(10)=rcmp*0.0857
c     a 5Pe ; imp=3, il=14(7),15(5),16(3); gsum=15
c     Add 5Pe-sum; imps=1
         rcmp=rcmpf1(3,it1)+tconst*(rcmpf1(3,it2)-rcmpf1(3,it1))
         rcsum=rcsumf1(1,it1)+tconst*(rcsumf1(1,it2)-rcsumf1(1,it1))
         rcmp=rcmp+rcsum
         recof1(14)=rcmp*0.4667
         recof1(15)=rcmp*0.3333
         recof1(16)=rcmp*0.2
c     z 7Do ; imp=4, il=18(11),21(9),23(7),25(5),26(3); gsum=35
c     Add 7Do-sum; imps=9
         rcmp=rcmpf1(4,it1)+tconst*(rcmpf1(4,it2)-rcmpf1(4,it1))
         rcsum=rcsumf1(9,it1)+tconst*(rcsumf1(9,it2)-rcsumf1(9,it1))
         rcmp=rcmp+rcsum
         recof1(18)=rcmp*0.3143
         recof1(21)=rcmp*0.2571
         recof1(23)=rcmp*0.2
         recof1(25)=rcmp*0.1429
         recof1(26)=rcmp*0.0857
c     z 7Fo ; imp=5, il=34(13),36(11),38(9),40(7),41(5),42(3),43(1); gsum=49
c     Add 7Fo-sum; imps=10
         rcmp=rcmpf1(5,it1)+tconst*(rcmpf1(5,it2)-rcmpf1(5,it1))
         rcsum=rcsumf1(10,it1)+tconst*(rcsumf1(10,it2)-rcsumf1(10,it1))
         rcmp=rcmp+rcsum
         recof1(34)=rcmp*0.2653
         recof1(36)=rcmp*0.2245
         recof1(38)=rcmp*0.1837
         recof1(40)=rcmp*0.1429
         recof1(41)=rcmp*0.1020
         recof1(42)=rcmp*0.0612
         recof1(43)=rcmp*0.0204
c     z 7Po ; imp=6, il=44(9),47(7),50(5); gsum=21
c     Add 7Po-sum; imps=8
         rcmp=rcmpf1(6,it1)+tconst*(rcmpf1(6,it2)-rcmpf1(6,it1))
         rcsum=rcsumf1(8,it1)+tconst*(rcsumf1(8,it2)-rcsumf1(8,it1))
         rcmp=rcmp+rcsum
         recof1(44)=rcmp*0.4286
         recof1(47)=rcmp*0.3333
         recof1(50)=rcmp*0.2381
c     z 5Do ; imp=7, il=54(9),56(7),58(5),61(3),62(1); gsum=25
         rcmp=rcmpf1(7,it1)+tconst*(rcmpf1(7,it2)-rcmpf1(7,it1))
         recof1(54)=rcmp*0.36
         recof1(56)=rcmp*0.28
         recof1(58)=rcmp*0.20
         recof1(61)=rcmp*0.12
         recof1(62)=rcmp*0.04
c     z 5Fo ; imp=8, il=65(11),66(9),67(7),69(5),70(3); gsum=35
         rcmp=rcmpf1(8,it1)+tconst*(rcmpf1(8,it2)-rcmpf1(8,it1))
         recof1(65)=rcmp*0.3143
         recof1(66)=rcmp*0.2571
         recof1(67)=rcmp*0.2
         recof1(69)=rcmp*0.1429
         recof1(70)=rcmp*0.0857
c     z 5Po ; imp=9, il=73(7),78(5),79(3); gsum=15
         rcmp=rcmpf1(9,it1)+tconst*(rcmpf1(9,it2)-rcmpf1(9,it1))
         recof1(73)=rcmp*0.4667
         recof1(78)=rcmp*0.3333
         recof1(79)=rcmp*0.2
c     y 5Do ; imp=10, il=89(9),91(7),94(5),96(3),98(1); gsum=25
c     Add 5Do-sum; imps=5
         rcmp=rcmpf1(10,it1)+tconst*(rcmpf1(10,it2)-rcmpf1(10,it1))
         rcsum=rcsumf1(5,it1)+tconst*(rcsumf1(5,it2)-rcsumf1(5,it1))
         rcmp=rcmp+rcsum
         recof1(89)=rcmp*0.36
         recof1(91)=rcmp*0.28
         recof1(94)=rcmp*0.20
         recof1(96)=rcmp*0.12
         recof1(98)=rcmp*0.04
c     y 5Fo ; imp=11, il=92(11),97(9),99(7),101(5),104(3); gsum=35
c     Add 5Fo-sum; imps=6
         rcmp=rcmpf1(11,it1)+tconst*(rcmpf1(11,it2)-rcmpf1(11,it1))
         rcsum=rcsumf1(6,it1)+tconst*(rcsumf1(6,it2)-rcsumf1(6,it1))
         rcmp=rcmp+rcsum
         recof1(92)=rcmp*0.3143
         recof1(97)=rcmp*0.2571
         recof1(99)=rcmp*0.2
         recof1(101)=rcmp*0.1429
         recof1(104)=rcmp*0.0857
c     z 5Go ; imp=12, il=105(11),106(13),107(9),109(7),111(5); gsum=45
c     Add 5Go-sum; imps=7
         rcmp=rcmpf1(12,it1)+tconst*(rcmpf1(12,it2)-rcmpf1(12,it1))
         rcsum=rcsumf1(7,it1)+tconst*(rcsumf1(7,it2)-rcsumf1(7,it1))
         rcmp=rcmp+rcsum
         recof1(105)=rcmp*0.2444
         recof1(106)=rcmp*0.2889
         recof1(107)=rcmp*0.2
         recof1(109)=rcmp*0.1556
         recof1(111)=rcmp*0.1111
c     y 5Po ; imp=13, il=114(7),115(5),117(3); gsum=15
c     Add 5Po-sum; imps=4
         rcmp=rcmpf1(13,it1)+tconst*(rcmpf1(13,it2)-rcmpf1(13,it1))
         rcsum=rcsumf1(4,it1)+tconst*(rcsumf1(4,it2)-rcsumf1(4,it1))
         rcmp=rcmp+rcsum
         recof1(114)=rcmp*0.4667
         recof1(115)=rcmp*0.3333
         recof1(117)=rcmp*0.2
c     
c     Add the remaining recombinations to the 8 highest levels 
c     (114 - 121; y5P, y3F, y3D) weighted
c     with their statistical weights so that the total recombination 
c     coefficient is correct ; gsum=42 (7,5,7,3,5,7,5,3)
c     
         recsum=0.d0
         do il=1,n
            recsum=recsum+recof1(il)
c            write(6,*)'il,recof1(il),recsum,rctfe1',il,recof1(il),recsum
c     &           ,rctfe1
         enddo
         recadd=rctfe1-recsum
         if (recadd.lt.0.) then
            write(50,*) 'Error in calc. rec. coeff. for Fe I !!'
         endif
c     add the remaining recomb. contribution to
c     all the levels, not only the highest
         gsum=0.d0
         do il=1,n
            gsum=gsum+gf1(il)
         enddo
         do il=1,n
            recof1(il)=recof1(il)+recadd*gf1(il)/gsum
c            write(6,928)il,recadd,gf1(il),gsum,recof1(il)
 928        format('recadd ',i5,1pe12.3,10e12.3)
         enddo

      elseif(ife==2) then
c     
c     Fe II
c     
         do il=1,nlf2
            recof2(il)=0.d0
         enddo
         if (te.lt.tef2(1)) then
            tconst=0.d0
            it1=1
            it2=1
         elseif (te.ge.tef2(81)) then
            tconst=0.d0
            it1=81
            it2=81
         else
            do i=1,80
               if (te.ge.tef2(i).and.te.lt.tef2(i+1)) then
                  it1=i
                  it2=i+1
                  goto 10
               endif
            enddo
 10         continue
            tconst=(te-tef2(it1))/(tef2(it2)-tef2(it1))
         endif
c     
         rctfe2=rectf2(it1)+tconst*(rectf2(it2)-rectf2(it1))
c     
c     a 6De ; imp=1, il=1(10),2(8),3(6),4(4),5(2); gsum=30
c     Add 6De-sum; imps=3
         rcmp=rcmpf2(1,it1)+tconst*(rcmpf2(1,it2)-rcmpf2(1,it1))
         rcsum=rcsumf2(3,it1)+tconst*(rcsumf2(3,it2)-rcsumf2(3,it1))
         rcmp=rcmp+rcsum
         recof2(1)=rcmp*0.3333
         recof2(2)=rcmp*0.2667
         recof2(3)=rcmp*0.2
         recof2(4)=rcmp*0.1333
         recof2(5)=rcmp*0.0667
c     a 4Fe ; imp=2, il=6(10),7(8),8(6),9(4); gsum=28
         rcmp=rcmpf2(2,it1)+tconst*(rcmpf2(2,it2)-rcmpf2(2,it1))
         recof2(6)=rcmp*0.3571
         recof2(7)=rcmp*0.2857
         recof2(8)=rcmp*0.2143
         recof2(9)=rcmp*0.1429
c     a 4De ; imp=3, il=10(8),11(6),12(4),13(2); gsum=20
         rcmp=rcmpf2(3,it1)+tconst*(rcmpf2(3,it2)-rcmpf2(3,it1))
         recof2(10)=rcmp*0.4
         recof2(11)=rcmp*0.3
         recof2(12)=rcmp*0.2
         recof2(13)=rcmp*0.1
c     a 4Pe ; imp=4, il=14(6),15(4),16(2); gsum=12
         rcmp=rcmpf2(4,it1)+tconst*(rcmpf2(4,it2)-rcmpf2(4,it1))
         recof2(14)=rcmp*0.5
         recof2(15)=rcmp*0.3333
         recof2(16)=rcmp*0.1667
c     b 4Pe ; imp=5, il=24(6),30(4),31(2); gsum=12
         rcmp=rcmpf2(5,it1)+tconst*(rcmpf2(5,it2)-rcmpf2(5,it1))
         recof2(24)=rcmp*0.5
         recof2(30)=rcmp*0.3333
         recof2(31)=rcmp*0.1667
c     a 4He ; imp=6, il=25(14),27(12),28(10),29(8); gsum=44
c     Add 4He-sum; imps=15
         rcmp=rcmpf2(6,it1)+tconst*(rcmpf2(6,it2)-rcmpf2(6,it1))
         rcsum=rcsumf2(15,it1)+tconst*(rcsumf2(15,it2)-rcsumf2(15,it1))
         rcmp=rcmp+rcsum
         recof2(25)=rcmp*0.3182
         recof2(27)=rcmp*0.2727
         recof2(28)=rcmp*0.2273
         recof2(29)=rcmp*0.1818
c     b 4Fe ; imp=7, il=32(10),33(8),34(6),35(4); gsum=28
         rcmp=rcmpf2(7,it1)+tconst*(rcmpf2(7,it2)-rcmpf2(7,it1))
         recof2(32)=rcmp*0.3571
         recof2(33)=rcmp*0.2857
         recof2(34)=rcmp*0.2143
         recof2(35)=rcmp*0.1429
c     a 6Se ; imp=8, il=36(6) 
c     Add 6Se-sum; imps=1
         rcmp=rcmpf2(8,it1)+tconst*(rcmpf2(8,it2)-rcmpf2(8,it1))
         rcsum=rcsumf2(1,it1)+tconst*(rcsumf2(1,it2)-rcsumf2(1,it1))
         rcmp=rcmp+rcsum
         recof2(36)=rcmp
c     a 4Ge ; imp=9, il=37(12),39(10),40(8),41(6); gsum=36
         rcmp=rcmpf2(9,it1)+tconst*(rcmpf2(9,it2)-rcmpf2(9,it1))
         recof2(37)=rcmp*0.3333
         recof2(39)=rcmp*0.2778
         recof2(40)=rcmp*0.2222
         recof2(41)=rcmp*0.1667
c     b 4De ; imp=10, il=49(4),50(2),51(6),52(8); gsum=20
         rcmp=rcmpf2(10,it1)+tconst*(rcmpf2(10,it2)-rcmpf2(10,it1))
         recof2(49)=rcmp*0.2
         recof2(50)=rcmp*0.1
         recof2(51)=rcmp*0.3
         recof2(52)=rcmp*0.4
c     z 6Do ; imp=11, il=64(10),65(8),66(6),67(4),68(2); gsum=30
c     Add 6Do-sum; imps=4
         rcmp=rcmpf2(11,it1)+tconst*(rcmpf2(11,it2)-rcmpf2(11,it1))
         rcsum=rcsumf2(4,it1)+tconst*(rcsumf2(4,it2)-rcsumf2(4,it1))
         rcmp=rcmp+rcsum
         recof2(64)=rcmp*0.3333
         recof2(65)=rcmp*0.2667
         recof2(66)=rcmp*0.2
         recof2(67)=rcmp*0.1333
         recof2(68)=rcmp*0.0667
c     z 6Fo ; imp=12, il=69(12),70(10),71(8),72(6),73(4),74(2); gsum=42
c     Add 6Fo-sum; imps=5
         rcmp=rcmpf2(12,it1)+tconst*(rcmpf2(12,it2)-rcmpf2(12,it1))
         rcsum=rcsumf2(5,it1)+tconst*(rcsumf2(5,it2)-rcsumf2(5,it1))
         rcmp=rcmp+rcsum
         recof2(69)=rcmp*0.2857
         recof2(70)=rcmp*0.2381
         recof2(71)=rcmp*0.1905
         recof2(72)=rcmp*0.1429
         recof2(73)=rcmp*0.0952
         recof2(74)=rcmp*0.0476
c     z 6Po ; imp=13, il=75(8),76(6),77(4); gsum=18
c     Add 6Po-sum; imps=2
         rcmp=rcmpf2(13,it1)+tconst*(rcmpf2(13,it2)-rcmpf2(13,it1))
         rcsum=rcsumf2(2,it1)+tconst*(rcsumf2(2,it2)-rcsumf2(2,it1))
         rcmp=rcmp+rcsum
         recof2(75)=rcmp*0.4445
         recof2(76)=rcmp*0.3333
         recof2(77)=rcmp*0.2222
c     z 4Fo ; imp=14, il=78(10),80(8),85(6),87(4); gsum=28
         rcmp=rcmpf2(14,it1)+tconst*(rcmpf2(14,it2)-rcmpf2(14,it1))
         recof2(78)=rcmp*0.3571
         recof2(80)=rcmp*0.2857
         recof2(85)=rcmp*0.2143
         recof2(87)=rcmp*0.1429
c     z 4Do ; imp=15, il=79(8),81(6),84(4),86(2); gsum=20
         rcmp=rcmpf2(15,it1)+tconst*(rcmpf2(15,it2)-rcmpf2(15,it1))
         recof2(79)=rcmp*0.4
         recof2(81)=rcmp*0.3
         recof2(84)=rcmp*0.2
         recof2(86)=rcmp*0.1
c     z 4Po ; imp=16, il=88(6),89(4),90(2); gsum=12
         rcmp=rcmpf2(16,it1)+tconst*(rcmpf2(16,it2)-rcmpf2(16,it1))
         recof2(88)=rcmp*0.5
         recof2(89)=rcmp*0.3333
         recof2(90)=rcmp*0.1667
c     c 4Pe ; imp=17, il=93(2),94(4),99(6); gsum=12
         rcmp=rcmpf2(17,it1)+tconst*(rcmpf2(17,it2)-rcmpf2(17,it1))
         recof2(93)=rcmp*0.1667
         recof2(94)=rcmp*0.3333
         recof2(99)=rcmp*0.5
c     c 4Fe ; imp=18, il=95(4),96(6),97(10),98(8); gsum=28
c     Add 4Fe-sum; imps=11
         rcmp=rcmpf2(18,it1)+tconst*(rcmpf2(18,it2)-rcmpf2(18,it1))
         rcsum=rcsumf2(11,it1)+tconst*(rcsumf2(11,it2)-rcsumf2(11,it1))
         rcmp=rcmp+rcsum
         recof2(95)=rcmp*0.1429
         recof2(96)=rcmp*0.2143
         recof2(97)=rcmp*0.3571
         recof2(98)=rcmp*0.2857
c     b 4Ge ; imp=19, il=101(12),102(10),103(6),104(8); gsum=36
c     Add 4Ge-sum; imps=13
         rcmp=rcmpf2(19,it1)+tconst*(rcmpf2(19,it2)-rcmpf2(19,it1))
         rcsum=rcsumf2(13,it1)+tconst*(rcsumf2(13,it2)-rcsumf2(13,it1))
         rcmp=rcmp+rcsum
         recof2(101)=rcmp*0.3333
         recof2(102)=rcmp*0.2778
         recof2(103)=rcmp*0.1667
         recof2(104)=rcmp*0.2222
c     d 4Pe ; imp=20, il=108(6),109(4),110(2); gsum=12
c     Add 4Pe-sum; imps=7
         rcmp=rcmpf2(20,it1)+tconst*(rcmpf2(20,it2)-rcmpf2(20,it1))
         rcsum=rcsumf2(7,it1)+tconst*(rcsumf2(7,it2)-rcsumf2(7,it1))
         rcmp=rcmp+rcsum
         recof2(108)=rcmp*0.5
         recof2(109)=rcmp*0.3333
         recof2(110)=rcmp*0.1667
c     z 4So ; imp=21, il=113(4) 
         rcmp=rcmpf2(21,it1)+tconst*(rcmpf2(21,it2)-rcmpf2(21,it1))
         recof2(113)=rcmp
c     c 4De ; imp=22, il=114(8),115(2),117(4),118(6); gsum=20
c     Add 4De-sum; imps=9
         rcmp=rcmpf2(22,it1)+tconst*(rcmpf2(22,it2)-rcmpf2(22,it1))
         rcsum=rcsumf2(9,it1)+tconst*(rcsumf2(9,it2)-rcsumf2(9,it1))
         rcmp=rcmp+rcsum
         recof2(114)=rcmp*0.4
         recof2(115)=rcmp*0.1
         recof2(117)=rcmp*0.2
         recof2(118)=rcmp*0.3
c     y 4Po ; imp=23, il=116(6),125(2),128(4); gsum=12
         rcmp=rcmpf2(23,it1)+tconst*(rcmpf2(23,it2)-rcmpf2(23,it1))
         recof2(116)=rcmp*0.5
         recof2(125)=rcmp*0.1667
         recof2(128)=rcmp*0.3333
c     z 4Go ; imp=24, il=119(12),120(10),123(8),126(6); gsum=36
         rcmp=rcmpf2(24,it1)+tconst*(rcmpf2(24,it2)-rcmpf2(24,it1))
         recof2(119)=rcmp*0.3333
         recof2(120)=rcmp*0.2778
         recof2(123)=rcmp*0.2222
         recof2(126)=rcmp*0.1667
c     z 4Ho ; imp=25, il=121(14),122(12),124(10),127(8); gsum=44
         rcmp=rcmpf2(25,it1)+tconst*(rcmpf2(25,it2)-rcmpf2(25,it1))
         recof2(121)=rcmp*0.3182
         recof2(122)=rcmp*0.2727
         recof2(124)=rcmp*0.2273
         recof2(127)=rcmp*0.1818
c     z 4Io ; imp=26, il=129(16),130(10),131(14),132(12); gsum=52
c     Add 4Io-sum; imps=17
         rcmp=rcmpf2(26,it1)+tconst*(rcmpf2(26,it2)-rcmpf2(26,it1))
         rcsum=rcsumf2(17,it1)+tconst*(rcsumf2(17,it2)-rcsumf2(17,it1))
         rcmp=rcmp+rcsum
         recof2(129)=rcmp*0.3077
         recof2(130)=rcmp*0.1923
         recof2(131)=rcmp*0.2692
         recof2(132)=rcmp*0.2308
c     y 4Do ; imp=27, il=133(8),138(6),139(2),140(4); gsum=20
         rcmp=rcmpf2(27,it1)+tconst*(rcmpf2(27,it2)-rcmpf2(27,it1))
         recof2(133)=rcmp*0.4
         recof2(138)=rcmp*0.3
         recof2(139)=rcmp*0.1
         recof2(140)=rcmp*0.2
c     y 4Fo ; imp=28, il=134(8),135(6),136(10),137(4); gsum=28
         rcmp=rcmpf2(28,it1)+tconst*(rcmpf2(28,it2)-rcmpf2(28,it1))
         recof2(134)=rcmp*0.2857
         recof2(135)=rcmp*0.2143
         recof2(136)=rcmp*0.3571
         recof2(137)=rcmp*0.1429
c     x 4Do ; imp=29, il=141(8),142(6),143(4),144(2); gsum=20
         rcmp=rcmpf2(29,it1)+tconst*(rcmpf2(29,it2)-rcmpf2(29,it1))
         recof2(141)=rcmp*0.4
         recof2(142)=rcmp*0.3
         recof2(143)=rcmp*0.2
         recof2(144)=rcmp*0.1
c     y 4Go ; imp=30, il=145(12),146(10),147(8),148(6); gsum=36
         rcmp=rcmpf2(30,it1)+tconst*(rcmpf2(30,it2)-rcmpf2(30,it1))
         recof2(145)=rcmp*0.3333
         recof2(146)=rcmp*0.2778
         recof2(147)=rcmp*0.2222
         recof2(148)=rcmp*0.1667
c     x 4Go ; imp=31, il=149(12),150(10),151(8),153(6); gsum=36
         rcmp=rcmpf2(31,it1)+tconst*(rcmpf2(31,it2)-rcmpf2(31,it1))
         recof2(149)=rcmp*0.3333
         recof2(150)=rcmp*0.2778
         recof2(151)=rcmp*0.2222
         recof2(153)=rcmp*0.1667
c     x 4Fo ; imp=32, il=152(10),154(8),157(6),159(4); gsum=28
         rcmp=rcmpf2(32,it1)+tconst*(rcmpf2(32,it2)-rcmpf2(32,it1))
         recof2(152)=rcmp*0.3571
         recof2(154)=rcmp*0.2857
         recof2(157)=rcmp*0.2143
         recof2(159)=rcmp*0.1429
c     y 4Ho ; imp=33, il=155(14),156(12),158(10),160(8); gsum=44
c     Add 4Ho-sum; imps=16
         rcmp=rcmpf2(33,it1)+tconst*(rcmpf2(33,it2)-rcmpf2(33,it1))
         rcsum=rcsumf2(16,it1)+tconst*(rcsumf2(16,it2)-rcsumf2(16,it1))
         rcmp=rcmp+rcsum
         recof2(155)=rcmp*0.3182
         recof2(156)=rcmp*0.2727
         recof2(158)=rcmp*0.2273
         recof2(160)=rcmp*0.1818
c     w 4Po ; imp=34, il=161(6),162(4),164(2); gsum=12
         rcmp=rcmpf2(34,it1)+tconst*(rcmpf2(34,it2)-rcmpf2(34,it1))
         recof2(161)=rcmp*0.5
         recof2(162)=rcmp*0.3333
         recof2(164)=rcmp*0.1667
c     w 4Fo ; imp=35, il=163(4),165(6),166(8),169(10); gsum=28
         rcmp=rcmpf2(35,it1)+tconst*(rcmpf2(35,it2)-rcmpf2(35,it1))
         recof2(163)=rcmp*0.1429
         recof2(165)=rcmp*0.2143
         recof2(166)=rcmp*0.2857
         recof2(169)=rcmp*0.3571
c     w 4Do ; imp=36, il=167(2),168(6),170(8),171(4); gsum=20
         rcmp=rcmpf2(36,it1)+tconst*(rcmpf2(36,it2)-rcmpf2(36,it1))
         recof2(167)=rcmp*0.1
         recof2(168)=rcmp*0.3
         recof2(170)=rcmp*0.4
         recof2(171)=rcmp*0.2
c     v 4Do ; imp=37, il=172(2),173(4),174(6),175(8); gsum=20
         rcmp=rcmpf2(37,it1)+tconst*(rcmpf2(37,it2)-rcmpf2(37,it1))
         recof2(172)=rcmp*0.1
         recof2(173)=rcmp*0.2
         recof2(174)=rcmp*0.3
         recof2(175)=rcmp*0.4
c     w 4Go ; imp=38, il=176(6),177(8),178(10),179(12); gsum=36
c     Add 4Go-sum; imps=14
         rcmp=rcmpf2(38,it1)+tconst*(rcmpf2(38,it2)-rcmpf2(38,it1))
         rcsum=rcsumf2(14,it1)+tconst*(rcsumf2(14,it2)-rcsumf2(14,it1))
         rcmp=rcmp+rcsum
         recof2(176)=rcmp*0.1667
         recof2(177)=rcmp*0.2222
         recof2(178)=rcmp*0.2778
         recof2(179)=rcmp*0.3333
c     y 4So ; imp=39, il=180(4); gsum=4
c     Add 4So-sum; imps=6
         rcmp=rcmpf2(39,it1)+tconst*(rcmpf2(39,it2)-rcmpf2(39,it1))
         rcsum=rcsumf2(6,it1)+tconst*(rcsumf2(6,it2)-rcsumf2(6,it1))
         rcmp=rcmp+rcsum
         recof2(180)=rcmp
c     u 4Po ; imp=40, il=181(2),182(4),183(6); gsum=12
c     Add 4Po-sum; imps=8
         rcmp=rcmpf2(40,it1)+tconst*(rcmpf2(40,it2)-rcmpf2(40,it1))
         rcsum=rcsumf2(8,it1)+tconst*(rcsumf2(8,it2)-rcsumf2(8,it1))
         rcmp=rcmp+rcsum
         recof2(181)=rcmp*0.1667
         recof2(182)=rcmp*0.3333
         recof2(183)=rcmp*0.5
c     t 4Do ; imp=41, il=184(2),185(4),186(6),187(8); gsum=20
c     Add 4Do-sum; imps=10
         rcmp=rcmpf2(41,it1)+tconst*(rcmpf2(41,it2)-rcmpf2(41,it1))
         rcsum=rcsumf2(10,it1)+tconst*(rcsumf2(10,it2)-rcsumf2(10,it1))
         rcmp=rcmp+rcsum
         recof2(184)=rcmp*0.1
         recof2(185)=rcmp*0.2
         recof2(186)=rcmp*0.3
         recof2(187)=rcmp*0.4
c     t 4Fo ; imp=42, il=188(4),189(6),190(10),191(8); gsum=28
c     Add 4Fo-sum; imps=12
         rcmp=rcmpf2(42,it1)+tconst*(rcmpf2(42,it2)-rcmpf2(42,it1))
         rcsum=rcsumf2(12,it1)+tconst*(rcsumf2(12,it2)-rcsumf2(12,it1))
         rcmp=rcmp+rcsum
         recof2(188)=rcmp*0.1429
         recof2(189)=rcmp*0.2143
         recof2(190)=rcmp*0.3571
         recof2(191)=rcmp*0.2857
c     
c     Add the remaining recombinations to the 4 highest levels 
c     (188 - 191; t4F) weighted
c     with their statistical weights so that the total recombination 
c     coefficient is correct ; gsum=28 (4,6,10,8)
c     
         recsum=0.d0
         do il=1,n
            recsum=recsum+recof2(il)
         enddo
         recadd=rctfe2-recsum
         if (recadd.lt.0.) then
            write(50,*) 'Error in calc. rec. coeff. for Fe II !!'
         endif
c     add the remaining recomb. contribution to
c     all the levels, not only the highest
         gsum=0.d0
         do il=1,n
            gsum=gsum+gf2(il)
         enddo
         do il=1,n
            recof2(il)=recof2(il)+recadd*gf2(il)/gsum
         enddo
c
      endif 
      return
      end
      

      SUBROUTINE SIMP(A,B,N,S)
      IMPLICIT REAL*8(A-H,O-Z)
C     CALC. ESCAPE PROB. INTEGRAL
      H=(B-A)/REAL(N-1)
      NA=N-2
      S1=0.
      DO I=2,NA,2
         X1=H*REAL(I-1)
         X2=H*REAL(I)
         S1=S1+16.*FU(X1)+14.*FU(X2)
      enddo
      S1=S1+7.*FU(A)+7.*FU(B)
      S=H*S1/15.+H**2*(FP(A)-FP(B))/15.
      RETURN
      END


      DOUBLE PRECISION FUNCTION FU(X)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/A2/SI,TA,TAUL
      F1=ABS(1.+SI*X*X)
      EX=TA/F1
      IF(EX.LT.30.) GOTO 200
      EX=30.
200   CONTINUE
      FU=F1*(1.-EXP(-EX))/TA
      IF(EX.GT.1.E-5) GOTO 300
      FU=F1*EX/TA
300   CONTINUE
      RETURN
      END

      DOUBLE PRECISION FUNCTION FP(X)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/A2/SI,TA,TAUL
      Q=1.+SI*X*X
      QAB=ABS(Q)
      EX=TA/QAB
      IF(EX.LT.30.) GOTO 200
      EX=30.
200   CONTINUE
      FP=2.*SI*X*(1.-EXP(-EX)*(1.+TA/QAB))/TA
      IF(EX.GT.1.E-5) GOTO 300
      FP=2.*SI*X*(-TA/QAB+EX+EX*TA/QAB)/TA
300   CONTINUE
      IF(Q.GT.0.) GOTO 100
      FP=-FP
  100 CONTINUE
      RETURN
      END

      DOUBLE PRECISION FUNCTION ESCCO(TE,TAU)
      IMPLICIT REAL*8(A-H,O-Z)
      TEV=TE/1.1609E4
      A=13.597/TEV
      IF(TAU.GT.A/3) ESCCO=0.5*(A/(3.*TAU))**.25*EXP(-(A/3.)**.75*
     &TAU**.25)
      IF(TAU.LE.A/3.) ESCCO=0.5*EXP(-TAU)
      RETURN
      END

