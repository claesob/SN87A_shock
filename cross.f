      SUBROUTINE CROSSSECT
      IMPLICIT REAL*8(A-H,O-Z)
      integer skl,skk
      include 'param'
      COMMON/INT/FL(2,NE1:NE2),SI(14,27,NE1:NE2),E1(NE1:NE3),E(NE1:NE3)
      COMMON/FRE/NINT,JMIN,JJ
      COMMON/SIK/SK(14,27,NE1:NE2)
      call enint
      call num_shell
      CALL CRVERN(6,6,3)
      CALL CRVERN(7,7,3)
      CALL CRVERN(8,8,3)
      CALL CRVERN(10,10,3)
      CALL CRVERN(11,3,4)
      CALL CRVERN(12,3,4)
      CALL CRVERN(13,4,5)
      CALL CRVERN(14,14,5)
      CALL CRVERN(16,16,5)
      CALL CRVERN(18,7,5)
      CALL CRVERN(20,3,7)
      CALL CRVERN(26,26,7)
      RETURN
      end

      SUBROUTINE CRVERN(iel,imax,ns)
c     iel = Z (=6 for C)
c gsv2 = cross section for element ielcf (= 3 for C)      
c ns  numer of shells 5 for Si, S, Ar, 7 for Fe       
      IMPLICIT REAL*8(A-H,O-Z)
      integer iel,imax,ns
c      PARAMETER (NE1=-200,NE2=130,NE3=NE2+1)
      include 'parameters.h'
      COMMON/FRE/NINT,JMIN,JJ
c      COMMON/CSI_vern/GSV(14,5,NE1:NE2)
      COMMON/CSI_vern2/GSV2(14,30,7,NE1:NE2)
      COMMON/INT/FL(2,NE1:NE2),SI(14,27,NE1:NE2),E1(NE1:NE3),E(NE1:NE3)
      COMMON/SIK/Sk(14,27,NE1:NE2)
      common/ninn/ninn(30)
      common/ntot/ntot(30)
      common/ph1/ph1(6,30,30,7)
      common/ethresh/ethresh(30,30,7)
      common/nshells/nsh(30,30)
      integer elcf(30),elcfi
      data elcf/1,2,0,0,0,3,4,5,0,6,7,8,9,10,0,11,0,12,0,13,
     &     0,0,0,0,0,14,0,15,0,0/
      elcfi=elcf(iel)
      DO  ion=1,imax
c     # of electrons in ion 
         nel=iel+1-ion
         nout=ntot(nel)
         if(iel.eq.nel.and.iel.gt.18)nout=7
         if(iel.eq.(nel+1).and.(iel.eq.20.or.iel.eq.21.or.iel.eq.22.
     &        or.iel.eq.25.or.iel.eq.26))nout=7
         nsh(iel,ion)=nout
         DO  is=1,nout
            ethresh(iel,ion,is)=ph1(1,iel,nel,is)
            do j=jmin,jj
               GSV2(elcfi,ion,is,j)=0.
               call phfit2(iel,nel,is,e1(j),gsv2(elcfi,ion,is,j))
            enddo
         enddo
      enddo
      DO ion=1,imax
         do j=jmin,jj
            SK(elcfi,ion,J)=GSv2(elcfi,ion,1,J)
            SI(elcfi,ion,J)=0.
            DO MI=2,ns
               SI(elcfi,ion,J)=SI(elcfi,ion,J)+GSv2(elcfi,ion,MI,J)
            enddo
         enddo
      enddo
      RETURN
      END


      SUBROUTINE ENINT
C     ***********************************************************
C     *****
C     ENERGY INTERVALS IN EV
C        E(-1)    E(0)    E(1)    E(
C     *****
C     ***********************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      include 'param'
c      PARAMETER (NE1=-100,NE2=250,NE3=NE2+1)
      COMMON/FRE/NINT,JMIN,JJ
      COMMON/INT/FL(2,NE1:NE2),SI(14,27,NE1:NE2),E1(NE1:NE3),E(NE1:NE3)
      parameter(nlp=30000)
      common/linepoint/nlines,jline(nlp),iion(nlp)
      common/contpoint/otscon(100),ncont,jcont(100),iioncont(100)
      common/specpoint/nspec,jspec(100),iionspec(100)
      PARAMETER(NFEL=3000)      
      COMMON/EQUIV/RADF,FR(2,nlp),W(nlp),CIN(nlp),WLI(nlp+NFEL),FB(300),tauline(nlp+nfel)
      common/numbow/nbow
      dimension wlsp(100),ispec(100),icin(100)
      DATA WLEV/1.239854e4/,WLBAL/3704.9034d0/
c     INTERVALS FROM PASHEN LIMIT TO BALMER
      NPASH=30
      DEP=LOG10(13.598/3.4)/(NPASH+1)
      NBALM=60
      WLX=2500.
      ELX=WLEV/WLX
      DEB=LOG10(13.598/ELX)/(NBALM+1)
      DWL1=100.
      DWL1=25.
      WL=WLX

      DO J=-NBALM,jmin,-1
         IF(WL.GT.10000.) DWL1=0.01*WL
         IF(WL.GT.20000.) DWL1=0.010*WL
         wlP1=WL
         wl=WL+DWL1
         E(J)=WLEV/WL
         EP1=WLEV/WLP1
         E1(J)=SQRT(E(J)*EP1)
      enddo
300   CONTINUE
c     INTERVALS FROM BALMER LIMIT TO LYMAN
      DO J=1-NBALM,1
         E(J)=ELX*10.**(REAL(J+NBALM-1)*DEB)
         E1(J)=ELX*10.**((REAL(J)+NBALM-.5)*DEB)
      enddo
c number of intervals from Lyman to He I edge
      N1=20
      n1=100
      DE1=0.2572/N1
      DO J=2,2+N1-1
         E(J)=13.598*10.**((J-2)*DE1)
         E1(J)=13.598*10.**((J-1.5)*DE1)
      enddo
      J2=2+N1
      n2=25
      N2=40
      DE2=0.34502/N2
      DO J=J2,J2+N2-1
         E(J)=24.587*10.**((J-J2)*DE2)
         E1(J)=24.587*10.**((J-J2+0.5)*DE2)
      enddo
      J3=J2+N2
      N3=25
      DE3=0.4076629/N3
      DO 230 J=J3,J3+N3-1
      E(J)=54.416*10.**((J-J3)*DE3)
  230 E1(J)=54.416*10.**((J-J3+0.5)*DE3)
      J4=J3+N3

      DO 240 J=J4,J4+5
      E(J)=139.12*10.**((J-J4)*0.0526138)
  240 E1(J)=139.12*10.**((J-J4+0.5)*0.0526138)
      J5=J4+5-27
 
      jl = j3

      ns = 50

      call enbin(54.416d0,139.12d0,ns,jl)

      ns = 40

      call enbin(139.12d0,280.d0,ns,jl)

      ns = 5

      call enbin(280.d0,296.d0,ns,jl)

      call enbin(296.d0,317.d0,ns,jl)

      call enbin(317.d0,347.d0,ns,jl)

      call enbin(347.d0,369.d0,ns,jl)

      call enbin(369.d0,392.d0,ns,jl)

      call enbin(392.d0,412.d0,ns,jl)

      call enbin(412.d0,432.d0,ns,jl)

      call enbin(432.d0,459.d0,ns,jl)

      call enbin(459.d0,490.d0,ns,jl)

      call enbin(490.d0,511.d0,ns,jl)

      call enbin(511.d0,533.d0,ns,jl)

      call enbin(533.d0,550.d0,ns,jl)

      call enbin(550.d0,570.d0,ns,jl)

      call enbin(570.d0,595.d0,ns,jl)

      call enbin(595.d0,627.d0,ns,jl)

      call enbin(627.d0,667.d0,ns,jl)

      call enbin(667.d0,702.d0,ns,jl)

      call enbin(702.d0,739.d0,ns,jl)

      call enbin(739.d0,802.d0,ns,jl)

      call enbin(802.d0,870.d0,ns,jl)

      call enbin(870.d0,1.6387d4,100,jl)

      call enbin(1.6387d4,1.d5,20,jl)


      jj = jl-1

      do j=jmin,jj
         wl=wlev/e(j)
         wl1=wlev/e1(j)
      enddo
c      stop
      RETURN
      END


      subroutine enbin(emin,emax,ns,jl)
      IMPLICIT REAL*8(A-H,O-Z)
      include 'param'
      COMMON/FRE/NINT,JMIN,JJ
      COMMON/INT/FL(2,NE1:NE2),SI(14,27,NE1:NE2),E1(NE1:NE3),E(NE1:NE3)

      de1=log10(emax/emin)/(ns-1)
      
      do j=2,ns
         jl=jl+1
         e(jl) = emin*10.**((j-1.)*de1)
         e1(jl-1)=sqrt(e(jl-1)*e(jl))
      enddo

      return
      end


      SUBROUTINE pointers
C     ***********************************************************
C     *****
c calculate pointers for continuum photoionizations
c x-ray lines done in subroutine atdat_mewe_gs
C     *****
C     ***********************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      include 'param'
c      PARAMETER (NE1=-100,NE2=250,NE3=NE2+1)
      COMMON/FRE/NINT,JMIN,JJ
      COMMON/INT/FL(2,NE1:NE2),SI(14,27,NE1:NE2),E1(NE1:NE3),E(NE1:NE3)
      parameter(nlp=30000)
      common/linepoint/nlines,jline(nlp),iion(nlp)
      common/contpoint/otscon(100),ncont,jcont(100),iioncont(100)
      common/specpoint/nspec,jspec(100),iionspec(100)
      PARAMETER(NFEL=3000)      
      COMMON/EQUIV/RADF,FR(2,nlp),W(nlp),CIN(nlp),WLI(nlp+NFEL),FB(300),tauline(nlp+nfel)
      common/numbow/nbow
      dimension wlsp(100),ispec(100),icin(100)
      DATA WLEV/1.239854e4/,WLBAL/3704.9034d0/

c     put all line pointers in their bins
c
c     jline(iline) = energy bin into which line iline from line_pointers.dat
c     is put (between e(jline) and e(jline+1))
c     iion(iline) = ion number of line iline
c
c      goto 2345
      open(24,file='./ATDAT/line_pointers_v2.dat',status='old')
      do i=1,15000
         jline(i)=-999
      enddo
      do i=1,15000
         read(24,*,end=99)ik,wl,iline
         enwl=12398.59/wl
         if(enwl.gt.e(jmin)) then
            do j=jmin,jj
               if(enwl.gt.e(j).and.enwl.lt.e(j+1)) then
                  jline(iline)=j
               endif
            enddo
         else
            jline(iline)=-999
         endif
         iion(iline)=ik
         wli(iline)=wl
      enddo
 92   format('jli ',i6,f10.1,i6,3f10.3)
 99   continue
      nlines=iline
      close (24)

 2345 continue
      open(24,file='./ATDAT/cont_pointers.dat',status='old')
      do i=1,100
         jcont(i)=-999
      enddo
      do i=1,100
         read(24,*,end=98)ik,wl,iline
         enwl=12398.59/wl
         if(enwl.gt.e(jmin)) then
            do j=jmin,jj
               if(enwl.gt.e(j).and.enwl.lt.e(j+1)) then
                  jcont(iline)=j
                  iioncont(iline)=ik
               endif
            enddo
         else
            jcont(iline)=-999
            iioncont(iline)=ik
         endif
         if(jcont(iline).ge.jmin) then
         else
            qqq=-999.
         endif
      enddo
 98   continue
      ncont=iline
      close (24)
c      open(24,file='spec_pointers.dat')

      wl1=911.5
      enwl1=12398.59/wl1
      wl2=500.
      enwl2=12398.59/wl2

      goto 2222
      do i=1,100
         jspec(i)=-999
      enddo
c if icin = 1 include line in cin-array
      do i=1,100
         icin(i)=1
      enddo
c He II Lya. Already included as cin(88)
      wlsp(1)=304.
      icin(1)=0
c He II two-photon ???. Added in "emiss"
      wlsp(2)=500.
      icin(2)=0
c He II Balmer 911.5. Added in "emiss"
      wlsp(3)=911.5
      icin(3)=0
      do i=1,3
         ispec(i)=3
      enddo
c O III Bowen 374 A
      wlsp(4)=374.
c O III Bowen 644 A
      wlsp(5)=644.
c O III Bowen 703 A
      wlsp(6)=703.
c He I 584. already as cin(47)
      wlsp(7)=584.
      icin(7)=0
c O III Bowen 3210
      wlsp(8)=3210.
c O III Bowen 3287
      wlsp(9)=3287.
      do i=4,9
         ispec(i)=16
      enddo
      ispec(7)=2
c N III Bowen 4642
      wlsp(10)=4642.
c N III Bowen 4100
      wlsp(11)=4100.
c N III Bowen 691
      wlsp(12)=691.
c N III Bowen 990
      wlsp(13)=990.
c N III Bowen 452
      wlsp(14)=452.
      do i=10,14
         ispec(i)=20
      enddo
      nspec=14
      do i=1,nspec
c         read(24,*,end=97)ik,wl,iline
         wl=wlsp(i)
         enwl=12398.59/wl
         if(enwl.gt.e(jmin)) then
            do j=jmin,jj
               if(enwl.gt.e(j).and.enwl.lt.e(j+1)) then
                  jspec(i)=j
               endif
            enddo
         else
            jspec(i)=-999
         endif
         if(jspec(i).ge.jmin) then
         else
            qqq=-999.
         endif
      enddo
 97   continue
c add the special lines to the lines in linepointers.dat
      k=nlines
c index of first Bowen line (=374A)      
      nbow=nlines + 1
      do i=1,nspec
         if(icin(i).eq.1) then
            k=k+1
            jline(k)=jspec(i)
            iion(k)=ispec(i)
            wli(k)=wlsp(i)
         endif
      enddo
      nlines=k
 2222 continue
      do iline=1,nlines
         if(jline(iline).ge.jmin) then
         else
            qqq=-999.
         endif
      enddo
c      nspec=iline
c      close (24)
      RETURN
      END

      SUBROUTINE CROSS
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*8 LAB(200)
      include 'param'
c      PARAMETER (NE1=-100,NE2=250,NE3=NE2+1)
      COMMON/INT/FL(2,NE1:NE2),SI(14,27,NE1:NE2),E1(NE1:NE3),E(NE1:NE3)
      COMMON/FRE/NINT,JMIN,JJ
      COMMON/TRES/ EL(14,27),EK(14,27)
      COMMON/SIK/SK(14,27,NE1:NE2)
      common/contpoint/otscon(100),ncont,jcont(100),iioncont(100)
      DIMENSION JLMIN(14,27),JKMIN(14,27)
      DATA PI/3.1415926E0/,ELCH/1.60219E-12/,AMU/1.660531E-24/
      dimension nionel(14)
      data nionel/1,2,6,7,8,10,11,12,13,14,16,18,20,26/
      data rydb/13.5987/
      EK(1,1)=13.598
      EL(1,1)=13.598
      EK(2,1)=24.587
      EL(2,1)=24.587
      EK(2,2)=54.416
      EL(2,2)=54.416

      EL(3,1)=11.26
      EK(3,1)=280.
      EL(3,2)=24.38
      EK(3,2)=296.
      EK(3,3)=317.0
      EL(3,3)=47.9
      EK(3,4)=347.0
      EL(3,4)=64.5
      EK(3,5)=392.
      EL(3,5)=392.
      EK(3,6)=490.
      EL(3,6)=490.

      EL(4,1)=14.53
      EK(4,1)=395.
      EL(4,2)=29.60
      EK(4,2)=412.
      EL(4,3)=47.45
      EK(4,3)=432.
      EL(4,4)=77.47
      EK(4,4)=459.
      EL(4,5)=97.89
      EK(4,5)=496.
      EL(4,6)=552.06
      EK(4,6)=552.06
      EL(4,7)=667.03
      EK(4,7)=667.03

      EL(5,1)=13.62
      EK(5,1)=533.
      EL(5,2)=35.12
      EK(5,2)=550.
      EL(5,3)=54.93
      EK(5,3)=570.
      EL(5,4)=77.41
      EK(5,4)=595.
      EK(5,6)=672.
      EL(5,6)=138.12
      EK(5,7)=739.
      EL(5,7)=739.
      EK(5,8)=870.
      EL(5,8)=870.
      EK(5,5)=627.
      EL(5,5)=113.9

C
CK>  NE I - NE X
C
c from verner routine
      EL(6,1)=21.564
      EK(6,1)=870.1
      EL(6,2)=40.96
      EK(6,2)=883.1
      EL(6,3)=63.74
      EK(6,3)=913.1
      EL(6,4)=97.12
      EK(6,4)=948.0
      EL(6,5)=126.2
      EK(6,5)=987.3
      EL(6,6)=157.9
      EK(6,6)=1031.
      EL(6,7)=207.3
      EL(6,7)=1078.
      EK(6,8)=1125.
      EL(6,8)=239.1
      EK(6,9)=1196.
      EL(6,9)=1196.
      EK(6,10)=1362.
      EL(6,10)=1362.


      EL(7,1)=5.139
      EK(7,1)=1079.
      EL(7,2)=47.286
      EK(7,2)=1097.
      EL(7,3)=71.64
      EK(7,3)=1118.
      EL(7,4)=98.92
      EK(7,4)=1143.
      EL(7,5)=138.4
      EK(7,5)=1185.
      EL(7,6)=172.2
      EK(7,6)=1230.
      EL(7,7)=208.5
      EK(7,7)=1281.
      EL(7,8)=264.2
      EK(7,8)=1335.
      EL(7,9)=299.9
      EK(7,9)=1386.
      EL(7,10)=1465.
      EK(7,10)=1465.
      EL(7,11)=1649.
      EK(7,11)=1649.


      EL(8,1)=7.646
      EK(8,1)=1311.
      EL(8,2)=15.04
      EK(8,2)=1320.
      EL(8,3)=80.14
      EK(8,3)=1336.
      EL(8,4)=109.3
      EK(8,4)=1356.
      EL(8,5)=141.3
      EK(8,5)=1400.
      EL(8,6)=186.5
      EK(8,6)=1449.
      EL(8,7)=224.9
      EK(8,7)=1501.
      EL(8,8)=266.
      EK(8,8)=1558.
      EL(8,9)=328.2
      EK(8,9)=1618.
      EL(8,10)=367.5
      EK(8,10)=1675.
      EL(8,11)=1762.
      EK(8,11)=1762.
      EL(8,12)=1963.
      EK(8,12)=1963.

      EL(9,1)=5.986
      EK(9,1)=1567.
      EL(9,2)=18.83
      EK(9,2)=1571.
      EL(9,3)=28.45
      EK(9,3)=1583.
      EL(9,4)=120.
      EK(9,4)=1604.
      EL(9,5)=153.8
      EK(9,5)=1634.
      EL(9,6)=190.5
      EK(9,6)=1688.
      EL(9,7)=241.4
      EK(9,7)=1739.
      EL(9,8)=284.6
      EK(9,8)=1799.
      EL(9,9)=330.1
      EK(9,9)=1862.
      EL(9,10)=399.4
      EK(9,10)=1929.
      EL(9,11)=442.1
      EK(9,11)=1992.
      EL(9,12)=2086.
      EK(9,12)=2086.
      EL(9,13)=2304.
      EK(9,13)=2304.

c si not verner compatible

      EL(10,1)=8.15
      EK(10,1)=1870.
      EL(10,2)=16.35
      EK(10,2)=1880.
      EL(10,3)=33.49
      EK(10,3)=1910.
      EL(10,4)=45.14
      EK(10,4)=1930.
      EL(10,5)=166.8
      EK(10,5)=1950.
      EL(10,6)=205.1
      EK(10,6)=1980.
      EL(10,7)=246.5
      EK(10,7)=2010.
      EL(10,8)=303.2
      EK(10,8)=2050.
      EL(10,9)=351.
      EK(10,9)=2090
      EL(10,10)=401.
      EK(10,10)=2140.
      EL(10,11)=476.
      EK(10,11)=2210.
      EL(10,12)=523.
      EK(10,12)=2247.
      EL(10,13)=2438.
      EK(10,13)=2438.
      EL(10,14)=2673.
      EK(10,14)=2673.


      EL(11,1)=10.360
      EK(11,2)=2477.
      EL(11,2)=23.33
      EK(11,2)=2478.
      EL(11,3)=34.83
      EK(11,3)=2486.
      EL(11,4)=47.31
      EK(11,4)=2502.
      EL(11,5)=72.68
      EK(11,5)=2522.
      EL(11,6)=88.05
      EK(11,6)=2544.
      EL(11,7)=280.9
      EK(11,7)=2569.
      EL(11,8)=328.2
      EK(11,8)=2641.
      EL(11,9)=379.1
      EK(11,9)=2705.
      EL(11,10)=447.1
      EK(11,10)=2782.
      EL(11,11)=504.8
      EK(11,11)=2859.
      EL(11,12)=564.7
      EK(11,12)=2941.
      EL(11,13)=651.7
      EK(11,13)=3029.
      EL(11,14)=707.2
      EK(11,14)=3107.
      EL(11,15)=3224.
      EK(11,15)=3224.
      EL(11,16)=3494.
      EK(11,16)=3494.

CK<
C
C    AR I - AR II
C
      EL(12,1)=15.759
      EK(12,1)=3203.
      EL(12,2)=27.629
      EK(12,2)=3208.
      EL(12,3)=40.74
      EK(12,3)=3216.
      EL(12,4)=59.81
      EK(12,4)=3228.
      EL(12,5)=75.02
      EK(12,5)=3253.
      EL(12,6)=91.01
      EK(12,6)=3277.
      EL(12,7)=124.3
      EK(12,7)=3303.
      EL(12,8)=143.5
      EK(12,8)=3331.
      EL(12,9)=422.5
      EK(12,9)=3361.
      EL(12,10)=478.7
      EK(12,10)=3446.
      EL(12,11)=539.
      EK(12,11)=3523.
      EL(12,12)=618.3
      EK(12,12)=3613.
      EL(12,13)=686.1
      EK(12,13)=3702.
      EL(12,14)=755.8
      EK(12,14)=3798.
      EL(12,15)=854.8
      EK(12,15)=3898.
      EL(12,16)=918.
      EK(12,16)=3988.
      EL(12,17)=4121.
      EK(12,17)=4121.
      EL(12,18)=4426.
      EK(12,18)=4426.

      EL(13,1)=6.113
      EK(13,2)=4043.
      EL(13,2)=11.871
      EK(13,2)=4047.
      EL(13,3)=50.91
      EK(13,3)=4053.
      EL(13,4)=67.27
      EK(13,4)=4063.
      EL(13,5)=84.51
      EK(13,5)=4078.
      EL(13,6)=108.9
      EK(13,6)=4105.
      EL(13,7)=127.3
      EK(13,7)=4133.
      EL(13,8)=147.2
      EK(13,8)=4163.
      EL(13,9)=188.3
      EK(13,9)=4198.
      EL(13,10)=211.3
      EK(13,10)=4229.
      EL(13,11)=591.9
      EK(13,11)=4265.
      EL(13,12)=657.2
      EK(13,12)=4362.
      EL(13,13)=726.7
      EK(13,13)=4453.
      EL(13,14)=817.7
      EK(13,14)=4555.
      EL(13,15)=894.6
      EK(13,15)=4650.
      EL(13,16)=974.5
      EK(13,16)=4767.
      EL(13,17)=1087.
      EK(13,17)=4880.
      EL(13,18)=1157.
      EK(13,18)=4982.
      EL(13,19)=5129.
      EK(13,19)=5129.
      EL(13,20)=5470.
      EK(13,20)=5470.

C
C  FE I - FE XV
C 
      EL(14,1)=7.902
      EK(14,1)=7124.
      EL(14,2)=16.19
      EK(14,2)=7140.
      EL(14,3)=30.65
      EK(14,3)=7155.
      EL(14,4)=54.8
      EK(14,4)=7169.
      EL(14,5)=75.01
      EK(14,5)=7084.4
      EL(14,6)=99.06
      EK(14,6)=7199.
      EL(14,7)=125.
      EK(14,7)=7217.
      EL(14,8)=151.1
      EK(14,8)=7237.
      EL(14,9)=233.6
      EK(14,9)=7275.
      EL(14,10)=262.1
      EK(14,10)=7316.0
      EL(14,11)=290.2
      EK(14,11)=7359.
      EL(14,12)=330.8
      EK(14,12)=7403.
      EL(14,13)=361.
      EK(14,13)=7450.
      EL(14,14)=392.2
      EK(14,14)=7499.
      EL(14,15)=457.
      EK(14,15)=7553.
      el(14,16) = 35.96*rydb
      el(14,17) = 92.75*rydb
      el(14,18) = 99.81*rydb
      el(14,19) = 107.0*rydb
      el(14,20) = 116.3*rydb
      el(14,21) = 124.1*rydb
      el(14,22) = 132.2*rydb
      el(14,23) = 143.3*rydb
      el(14,24) = 150.4*rydb
      el(14,25) = 648.9*rydb
      el(14,26) = 681.9*rydb


      ek(14,25) = 648.9*rydb
      ek(14,26) = 681.9*rydb

      EL(14,16)=489.3
      EK(14,16)=7599.
      EL(14,17)=1262.
      EK(14,17)=7651.
      EL(14,18)=1358.
      EK(14,18)=7769.
      EL(14,19)=1456.
      EK(14,19)=7918.
      EL(14,20)=1582.
      EK(14,20)=8041
      EL(14,21)=1689.
      EK(14,21)=8184
      EL(14,22)=1799.
      EK(14,22)=8350.
      EL(14,23)=1950.
      EK(14,23)=8484.
      EL(14,24)=2046.
      EK(14,24)=8638.
      EL(14,25)=8829.
      EK(14,25)=8829
      EL(14,26)=9278.
      EK(14,26)=9278



      ELoiii1d=52.42


66    CONTINUE

  300 CONTINUE

C     SET ALL LIMITS OF THE CROSS SECTIONS TO THE NEAREST ENERGY BIN

      do iel=1,14
         do ion=1,nionel(iel)
            JLMIN(iel,ion)=-999
            JkMIN(iel,ion)=-999
            DO J=JMIN,JJ
               IF(EL(iel,ion).GT.E(J).AND.EL(iel,ion).le.E(J+1)) THEN
                  EMEAN=(E(J)+E(J+1))/2.
                  IF(EL(iel,ion).GT.EMEAN) THEN 
                     JLMIN(iel,ion)=J+1
                  ELSE
                     JLMIN(iel,ion)=J
                  ENDIF
               ENDIF
               IF(EK(iel,ion).GT.E(J).AND.EK(iel,ion).Le.E(J+1)) THEN
                  EMEAN=(E(J)+E(J+1))/2.
                  IF(EK(iel,ion).GT.EMEAN) THEN 
                     JKMIN(iel,ion)=J+1
                  ELSE
                     JKMIN(iel,ion)=J
                  ENDIF
               ENDIF
            ENDDO

            DO J=JMIN,JJ
               IF(J.LT.JLMIN(iel,ion)) THEN
                  SI(iel,ion,J)=0.
               ENDIF
               IF(J.LT.JKMIN(iel,ion)) THEN
                  SK(iel,ion,J)=0.
               ENDIF
            ENDDO
            ILQ=0
            IKQ=0

         ENDDO
      enddo


91    FORMAT(a8,2I5,5F10.2)

      call ionlabel(lab)

      RETURN
      END






