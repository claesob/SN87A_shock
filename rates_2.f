      SUBROUTINE EMISS(TE,xel,dens)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 ar(100),BD(100),cr(100),dr(100),er(100)
      real*8 lyem
      include 'param'
      PARAMETER(NFEL=3000)
      PARAMETER (NL=340,NLP1=NL+1)
C     *****************************************************************
C     *****
C     THIS ROUTINE CALC. THE EMISSION RATES FOR GIVEN TEMPERATURE
C     AND ABUNDANCES OF THE RELEVANT IONS.FREE-FREE AND RECOMBI-
C     NATION EMISSION ARE CONSIDERED.
C     ALL RATES ARE EMISSION INTEGRATED OVER 4 PI, EXCEPT EM(I,J)
C     WHICH IS EMISSION PER STERADIAN
c     Only one shell with index i is done
C     *****
C     *****************************************************************
       COMMON/EMHY/RECEM(6,NL),TWOPHH,TWOPHHEI,COHI,PHI,PHIT,PO,POT
      common/bowem/ emb600, emb762, emb3210, emb3287, em374,
     &              emb4642, emb691, emb990, emb4100, emb452
      common/heiires/x2he,he2coll12,he2jla,he2j2g,he2jrec(4),
     &                                                he2jem(7,4)
      COMMON/GSREC/ALGS(100),ALEX(100),ALTOT(100),RECEX(100),RECGS(100)
      COMMON/ELEC/DEL(MD)
      common/ionx/xion(md,14,27)
      COMMON/IND/I
      COMMON/mdIND/kradius
      COMMON/TRES/EL(14,27),EK(14,27)
      parameter (nlp=30000)
      COMMON/EQUIV/RADF,FR(2,nlp),W(nlp),CIN(nlp),WLI(nlp+NFEL),FB(300),tauline(nlp+nfel)
      common/lineioncf/ilabcf(nlp+nfel)
      COMMON/FRE/NINT,JMIN,JJ
      COMMON/INT/FL(2,NE1:NE2),SI(14,27,NE1:NE2),E1(NE1:NE3),E(NE1:NE3)
      COMMON/DIF/EM(MD,NE1:NE2),TAU(MD,NE1:NE2),TAUTOT(MD,NE1:NE2),
     &           EMC(MD,NE1:NE2)
      COMMON/PHY/DEN(MD)
      COMMON/REC/AL2(16)
      COMMON/BOLD/TOLD,TOLDH,TOLDCA,TOLDO,TOLDHE,TOLDFE,toldfei,BB(NL),
     &       BOLH(NL),BOL(NL),BOLCA(NL),BOLHE(NL),BOLFE(NL),bolfei(nl)
      COMMON/A17/RECNET(NL),RECNDT(NL),RECCO(NL),RTEQ(NL),
     &PHN(NL),DRDT(NL),RECT(6,NL),PHET(NL)
      common/abl/abn(15)
      COMMON/SPOT/OTSP(16)
      COMMON/HPOP/XN(6,NL),XN1,XN2,XN3
      COMMON/LINE/XLINE
      COMMON/SIK/SK(14,27,NE1:NE2)
      COMMON/REC1/RE1(7),RHE2B
      COMMON/WFE/WOBFE(NL,NL),wobfei(nl,nl),IQPOPH
      COMMON/FELEV/NFEII,nfei
      COMMON/A14/CQ(NL,NL),CI(NL),G(NL),EI(NL),AQ(NL,NL),WL(NL,NL)
     &     ,DCDT(NL,NL),DCIDT(NL)
      common/twophotdist/twophfh(NE1:NE3),twophfhe(NE1:NE3)
      common/reslpointer/jmgii
      common/nperm/nline
      common/specpoint/nspec,jspec(100),iionspec(100)
      common/contpoint/otscon(100),ncont,jcont(100),iioncont(100)
      common/otsemiss/emissots(100)
      common/xrayem/fxlin(ne1:ne3)
      common/timecheck/time,itime
      common/preion/ipre
      common/nhion/np1hs
      common/inidi/initdi
      common/debug/ideb
      integer nionstage
      common/num_ion_stages/nionstage(14)
      common/line_em/cinx(14,26,401),taulinex(14,26,401),wlix(14,26,401)
     &     ,ilabcfx(14,26,401)
      data nionstage/1,2,6,7,8,10,3,3,4,14,16,7,3,26/
      common/kmax/kmax(14,27)
      common/hgsrec/algs_h1
      DIMENSION GS(50),Z(50),KS(50),P(500),ENL(500)
      DIMENSION INDI(50),INDIC(10,50),JMAX(50),JIO(50)
      DIMENSION ALV(100),EX(100)
      DIMENSION kqq(100),EMCHECK(NL),WLF(300)
      DIMENSION WLR(8,5),FRR(8,5),XO(8)
      DIMENSION WLDR(100),threshfb(2),fbh(ne1:ne3)
      DIMENSION IOION(100)
      save ar,BD,cr,dr,er,IOION,WLDR,ndie
c     recombination emission from O I-V
      DATA WLR/8*0.,
     &        673., 617., 518., 515., 485., 539., 430., 0.,
     &        435., 396., 345., 328., 321., 374., 306., 304., 
     &        279.8, 238.5, 6*0., 
     &        629.7,248.5,172.2,220.3,215.2,192.9,1216.,0./
      DATA FRR/8*0.,
     &            .026,.053,.014,.020,.038,.30,.035,0.,
     &            .014,.074,.013,.023,.032,.27,.069,.040,
     &            .22,.27,6*0.,
     &            .11,.0057,.023,.041,.11,.16,.37,0./
C     GS = (STATISTICAL WEIGHT OF GROUND STATE)/2.
      DATA GS/2.,0.5,2.,2.,.5,2.,.5,.5,2.,.5,2.,.83333,6.,1.25,
     &        .44444,1.5,6.,.44444,1.5,6.,.5,2.,.5,2.,26*0./
C     CHARGE OF THE RECOMBINING ION
      DATA Z/1.,1.,2.,6.,7.,8.,5.,3.,4.,5.,6.,1.,2.,1.,2.,3.,4.,
     &       1.,2.,3.,4.,5.,6.,7.,26*0./
C     IF KS=1 THEN USE K-SHELL CROSSECTION
      DATA KS/1,1,1,-1,1,1,-1,-1,-1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,
     &        -1,-1,-1,1,1,26*-1/
      DATA JIO/2,7,14,47*0/
      DATA PI/3.1415926E0/,ELCH/1.60219E-12/,AMU/1.660531E-24/
      if(te<250.)  ideb=1

      if(initdi.eq.0) then
      INITDI=1
      OPEN(39,FILE='diel.dat')
      do k=1,100
         read(39,*,end=88)IOION(k),WLDR(k),ar(k),BD(k),cr(k),dr(k),er(k)
      enddo
88    ndie=k-1 
      rewind (39)
c      close (39)
      endif
c      stop
      emtot=0.
      bah=0.
      texc=0.
      TEV=TE/1.1609E4
      T4=TE/1.E4
      DO  J=JMIN,JJ
         EM(I,J)=0.d0
         EMC(I,J)=0.d0
      ENDDO
C      ALV(I) = RECOMBINATION COEFFICIENT TO GROUND STATE OF ION I
      DO K=1,100
         EX(K)=0.5
         ALV(K)=ALGS(K)*1.E13*T4**.5
      ENDDO
      ALV(1)=1.54
      EX(1)=0.51
C     HE I CONSISTENT WITH RECOMB()
      ALV(2)=1.59
      EX(2)=0.488
c     do He II separately below
      ALV(3)=0.
      EX(3)=0.503
      ALV(4)=8.98
      EX(4)=0.50
      ALV(7)=3.27
      EX(7)=0.50
      ALV(8)=1.30
      EX(8)=0.50
      ALV(9)=3.61
      EX(9)=0.50
      ALV(10)=28.1
      EX(10)=0.50
      ALV(11)=65.0
      EX(11)=0.50
      ALV(13)=13.35
      EX(13)=0.50
c!!   Same al as in RECOMB()
      ALV(14)=1.22
      EX(14)=0.70
      ALV(42)=2.46E-2
c!!   guesses for 42 ad 43 = Fe I and II
      ALV(42)=0.4*1.42
      ALV(43)=0.4*10.2
      ALV(50)=5.61E-4
      ALV(51)=4.50E-2
      ALV(53)=2.83E-3

      DO K=1,100
         ALV(K)=ALV(K)*1.E-13/T4**EX(K)
      enddo
      imin=11111111
c!!!
      kqq(12)=1
      kqq(14)=1
      kqq(25)=1
      kqq(39)=1
c!! stryk denna
      em(i,78)=1.d-80

c     excited continua of O III - O IV
      DO K=1,3
        IF(K.EQ.1) THEN
C          O III N=3
           REV=16.3*ELCH*ABn(5)*xion(i,5,4)*1.11E-12/T4**.769
           EN=16.3
        ELSEIF(K.EQ.2) THEN
C          O IV N=4
           REV=14.1*ELCH*ABn(5)*xion(i,5,5)*1.23E-12/T4**.843
           EN=14.10
        ELSEIF(K.EQ.3) THEN
C          O IV N=3
           REV=27.45*ELCH*ABn(5)*xion(i,5,5)*1.53E-12/T4**.769
           EN=27.45
        ENDIF
           REV=REV*XEL/(4.*PI)
           DO J=JMIN,JJ
              IF(EN.GT.E(J-1).AND.EN.LT.E(J)) THEN
                EM(I,J-1)=REV/(E(J)-E(J-1))+EM(I,J-1)
                demq=REV/(E(J)-E(J-1))
                if(kradius> 511.and.rev > 1.e-30) then
                   write(6,9121)j,rev,em(i,j-1)
 9121              format('O III exi ',i5,1pe12.3,10e12.3)
                endif
                   if(j.eq.8) then                                            
                      rev8=rev
 922                  format(' em, oc ',i5,1pe12.3,10e12.3)
                   endif
                     if(em(i,j-1)< 0.) then
                        write(6,*)' OI-III con',j,k,rev,e(j),em(i,j-1)
                     endif                     
                GOTO 752
              ENDIF
            ENDDO
752       CONTINUE
      ENDDO

c     new from pwn
      DO iel=1,14
         do ion=1,nionstage(iel)
            do k=1,401
               wl0=wlix(iel,ion,k)
               IF(wl0.gt.0.) then
                  EN=12398.54/WL0
                  IF(EN > E(JMIN)) then
                     DO J=JMIN,JJ
                        IF(en > e(j-1) .and. en < e(j)) then
                           EMJ=DEL(I)*CINx(iel,ion,k)/(4.*PI*(E(J)-E(J-1)))
                           EM(I,J-1)=EM(I,J-1)+EMJ
                           if(kradius>511.and.emj>1.e-25) then
                              write(6,9251)j,iel,ion,k,wl0,cinx(iel,ion,k),em(i,j-1)
 9251                         format('iel,ion,k,wl0,cinx ',4i5,f10.3,1pe12.3,10e12.3)
                           endif
                        endif
                     enddo
                  endif
               endif
            enddo
         enddo
      enddo

      DO  KA=1,24
         DO N=1,10
            INDIC(N,KA)=-1
         enddo
         INDI(KA)=-1
      enddo
      EMFF=0.
      
      DO K1=1,3
         if(k1.eq.1) then
            chi = ek(1,1)
         elseif(k1.eq.2.or.k1.eq.3) then
            chi = ek(2,k1-1)
         endif
         E2=CHI+1.5*TE/1.16E4
         DO J=JIO(K1),JJ
            IF(E2.lT.E(J)) then
               JMAX(K1)=J-1
               GOTO 136
            endif
         enddo
 136     IF(JMAX(K1).EQ.JIO(K1)) JMAX(K1)=JIO(K1)+1
      enddo

      DO K=1,NL
         EMCHECK(K)=0.
      ENDDO

C     EXCITED O I EMISSION
      CALL ATDATO
      NP1O=14
      DO J=JMIN,JJ
        DO N=4,9
c     CALL RECR(2,N,J,CHI,AZG,KS,TE,RE)
           CALL RECR(2,N,J,CHI,TE,RE)
            DRECO=ABN(5)*BOL(NP1O)*RE*XEL/(4.*PI)
            EM(I,J)=EM(I,J)+DRECO
            if(kradius> 511.and.re>1.d-30) then
               write(6,9124)j,bol(np1o),re,dreco,em(i,j)
 9124          format('O I exi ',i5,1pe12.3,10e12.3)
            endif
          EMC(I,J)=EMC(I,J)+DRECO
        ENDDO
      ENDDO
      
C     *****************************************************************
C     LINE RADIATION
C     *****************************************************************
c add x-ray lines (not used anymore. only for mewe list)
c
c note that bin em(i,j) is between j and j+1  (see main)

      etotxray= 0.
      do j=jmin,jj
         em(i,j) = em(i,j) + fxlin(j)/((E(J+1)-E(J))*4.*PI)
         dfxl =  fxlin(j)/((E(J+1)-E(J))*4.*PI)
         dem=4.*pi*(e(j+1)-e(j))*em(i,j)
         etotxray= etotxray + fxlin(j)
c     if(fxlin(j) > 0. .and. kradius>50) then
         if(j.eq.8.or.j==32.or.j==107) then
            wlbin1=12398./e(j)
            wlbin3=12398./e(j+1)
         endif
      enddo
      if(kradius > 500) then
         write(6,*)'kradius ',kradius,i
         write(6,9381)kradius,dreco8,re8,emj8,rev8,revo8
         do j=jmin,jj
            if(isnan(em(i,j))) then
               write(6,9462)j,e1(j),em(i,j)
 9462          format('em ',i5,1pe12.3,10e12.3)
            endif
           if (isnan(em(i,j))) then
              write(6,9382)kradius,j,te,dreco8,re8,emj8,rev8,revo8,em(i,8)
 9381         format('NaN dreco8,re8,emj8,rev8,revo8',i5,1pe12.3,10e12.3)
 9382         format('NaN dreco8,re8,emj8,rev8,revo8',2i5,1pe12.3,10e12.3)              
           endif
        enddo
      endif

                                       
      RETURN
      END

      subroutine romb(fx,aa,bb,err,fi,fail,k)
      implicit real*8(a-h,o-z)
      external fx
      dimension tr(100),w(100)
      a=aa
      b=bb
      dx=(b-a)/2.
      fi1=0.5*(fx(a)+fx(b))
      fi2=fx(a+dx)
      tr(1)=dx*(fi1+fi2)
      n=1
      k=1
      w(1)=4.
1     fi=tr(1)
      w(k+1)=4.*w(k)
      n=2*n
      k=k+1
      tdx=dx
      dx=dx/2.
      x=a+dx
      do i=1,n
         fi2=fi2+fx(x)
         x=x+tdx
      enddo
      tr(k)=dx*(fi1+fi2)
      kk=k-1
      do  j=1,kk
         m=k-j
         tr(m)=(w(j)*tr(m+1)-tr(m))/(w(j)-1.)
      enddo
      fail=abs(tr(1)-fi)/abs(tr(1))
      if(fail-err)6,6,4
4     if(k-100)1,6,6
6     fi=tr(1)
      return
      end



