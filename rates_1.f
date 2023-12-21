      subroutine energy_eq(dt,coolf,xe,dens,teold,tenew)      
      implicit real*8(a-h,o-z)
      real*8 kb
      data kb/1.3806e-16/
c
c from Kozma & CF 1998. neglect dx_e/dt and adiabatic cooling terms
c
      tenew = teold + 2.*coolf*xe*dens*dt/(3*kb*(1.+xe))
      dte = 2.*coolf*xe*dens*dt/(3*kb*(1.+xe))
      return
      end

      subroutine elecdens(ik,xe,ze)
      implicit real*8(a-h,o-z)
      include 'PARAM'
      common/ionx/xion(md,14,27)
      common/abl/abn(15)
      dimension nionel(14)
      data nionel/1,2,6,7,8,10,11,12,13,14,16,18,20,26/

      xeo=xe
      xe = 0.
      ze = 0.
      do iel=3,14
         do i=2,nionel(iel)+1      
            xe = xe + abn(iel)*(i-1.)*xion(ik,iel,i)
            if(iel.ge.3) then
               ze = ze + abn(iel)*(i-1.)*xion(ik,iel,i)
            endif            
         enddo
      enddo

      return
      end

      subroutine timestep(iss,k,ns,lambda,xtot,vel,dens,dxtotsh,
     &     dtold,dtnew)
      implicit real*8(a-h,o-z)

      include 'PARAM'

      common/heitau/tauheiold,tauhei,dtauhei
      COMMON/DIF/EM(MD,NE1:NE2),TAU(MD,NE1:NE2),TAUTOT(MD,NE1:NE2),
     &           EMC(MD,NE1:NE2)

      tauhi = tautot(k,3)

      tauhei = tautot(k,104)

      tauheii = tautot(k,143)

      return
      end


      subroutine ionization_aug(ik,dt,denel,denh,dent)
c     implicit real*8(a-h,o-z)
      implicit none
      integer nl,nlp1,nispec,iion,ipre,ideb,itime,nz,nion,nshell
      PARAMETER (NL=340,NLP1=NL+1)
      real*8 abn,xion,te_old,xion_old,time,dt,denel,denh,dent
c     real*8 aa(nlp1,nlp1),aaold(nlp1,nlp1),ction,alrec,di,eps,geion
      real*8 aa(nlp1,nlp1),ction,alrec,di,eps,geion
      parameter (nispec=14,iion=26)
      include 'PARAM'
      common/abl/abn(15)
      common/ionx/xion(md,14,27)
      common/ionxold/te_old(md),xion_old(md,14,27)
      common/preion/ipre
      common/timecheck/time,itime
      common/debug/ideb
      parameter(nz=30,nion=27,nshell=10)
      integer iel,ion,shell,epsi,fr(10)
      integer ns,kmax,n,ik,i,ielnat,ierr,il,init,init_aug,iout,ipr,iz,izmax
      integer j,k,na,nionel,nrc,nsi,init_augfrac
      real*8 fr_aug,eion,en_aug,en_augi,eioni,x,xo,zion,ztot,ztot2,collion
      real*8 xold,simul,dt_old
      common/augerfrac/eion(nz,nion,nshell),en_aug(nz,nion,nshell),
     &     fr_aug(nz,nion,nshell,10),kmax(nz,nion,nshell),
     &     init_augfrac
      common/shells/ns(30,27)
      common/rec_coll/alrec(30,30),collion(30,30),ction(14,27)
      common/raug/zion(30,27,7),geion(30,27,7)
      dimension nionel(14)
c      data nionel/1,2,6,7,8,8,2,3,4,14,10,6,2,15/
      data nionel/1,2,6,7,8,10,11,12,13,14,16,18,20,26/
      dimension ztot(30)
c     dimension AA(nlp1,nlp1),X(nl),XO(NL),xold(nl),aaold(nlp1,nlp1)
      dimension X(nl),XO(NL),xold(nl)

c     NOTE!
c     iel = real atomic number 1-26 for alrec and collion
      

c  dn(i)/dt = (n(i) - nold(i))/dt = 
c             a_(i+1) n(i+1) -(a(i) + c(i) + p(i)) n(i) + (p(i-1)+c(i-1)) n(i-1)

c  n(1) + n(2) +.....+ n(nionel(iel)+1) = 1.


c    n(iel) = 2

c    n(1)/dt - a_(2) n(2) + (c(1) + p(1)) n(1) =  nold(1))/dt

c    n(2)/dt - a_(3) n(3) + (c(2) + p(2)) n(2) - (c(1) + p(1)) n(1) =  nold(2))/dt

c    n(1) + n(2)  + n(3) = 1.

c non-hydrogenic electron density

      do iel=3,14
         izmax=nionel(iel)
         ipr=0
c NOTE! here iel = CF numbering 1-14(Fe)
c     do i=1,nionel(iel)+2
c            do j=1,nionel(iel)+2           
         do i=1,nlp1
            do j=1,nlp1
               aa(i,j) = 0.            
            enddo
         enddo

         do i=1,nionel(iel)+1

            xold(i)=xion(ik,iel,i)

            do j=1,nionel(iel)
            
               aa(i,j) = 0.d0
            
            enddo
         enddo

c     rates from ion to ion+ion+j, includig auger fractions from Kaastra & Mewe
c NOTE! zion(iel,ion,i) here iel=1-14                  
         do ion=1,izmax
            ztot(ion)=0.
            ztot2=0.
            do i=1,ns(iel,ion)
               ztot2= zion(iel,ion,i) + ztot2
               ztot(ion)=ztot2
            enddo
         enddo
         ztot(izmax+1)=0.

         do i=1,nionel(iel)+1
c     neutral has no losses to lower stages!
            if(i.eq.1) then
c      auger, C-T and coll. ion. losses to higher stages
               aa(i,i) = 1./dt + ztot(i) + 
     &              ction(iel,i)*dent + collion(izmax,i)*denel
c     recombination from higher stages.. not needed to do separately                              
               aa(i,i+1) = -alrec(izmax,i+1)*denel

            else
c     now the higher stages
               
c     Auger
               do k=1,10
c     lower ionization stage added to i            
                  il=i-k
                  if(il>0 ) then
c     go through all the shells of il               
                     do nsi=1,ns(iel,il)
                        if(kmax(iel,il,nsi)>0) then
c     Note, iel = CF number. fr_aug in CF numbering in aug_fr()
                           aa(i,il)=-fr_aug(iel,il,nsi,k)*zion(iel,il,nsi) +
     &                          aa(i,il)
                        endif
                     enddo
                  endif
               enddo
c add coll. ion. and CT              
               aa(i,i-1) = -(ction(iel,i-1)*dent + collion(izmax,i-1)*denel) + aa(i,i-1)
c total losses from level i               
               aa(i,i) = 1./dt + alrec(izmax,i)*denel + ztot(i) + ction(iel,i)*dent + 
     &              collion(izmax,i)*denel
               if(iel==-5) then
                  write(6,9281)i,izmax,dt,denel,alrec(izmax,i),ztot(i),ction(iel,i),dent,collion(izmax,i)
 9281             format('aa(i,i) ',2i4,1pe12.3,10e12.3)
               endif
c     recombination from higher stages               
               aa(i,i+1) = -alrec(izmax,i+1)*denel

            endif

c            aa(i,nionel(iel)) = 1./dt
            
c number cons.

c     *            aa(nionel(iel)+1,i) = 1.
            aa(1,i) = 1. 

c RHS of stat. eqns.

            aa(i,nionel(iel)+2)= xold(i)/dt

         enddo

c ionized stage

c     *         aa(nionel(iel)+1,nionel(iel)+1) = 1.
         aa(1,nionel(iel)+1) = 1. 

c RHS of cons.


c     *         aa(nionel(iel)+1,nionel(iel)+2) = 1.
         aa(1,nionel(iel)+2) = 1. 
c number of eqns
         na=nionel(iel)+1
c number of columns in matrix incl RHS
         nrc=na+1

         eps=1.d-30

         iout=6
         di=simul(na,aa,x,eps,1,nrc)
         ierr=0
         
         do i=1,nionel(iel)+1
            if(x(i).lt.0.) then
               x(i)=xold(i)/2.
            endif
            xion(ik,iel,i) = x(i)
            if(xion(ik,iel,i).gt.1.01) then
               ierr=1
            endif
         enddo
      enddo
      return
      end



      
      subroutine ionization_h_he(ik,dt,dent,denel,
     &     photok,xel,ze)
      implicit real*8(a-h,o-z)
      PARAMETER (NL=340,NLP1=NL+1)
      include 'PARAM'
      common/ionx/xion(md,14,27)
      common/ionxold/te_old(md),xion_old(md,14,27)
      common/abl/abn(15)
      common/timecheck/time,itime
      COMMON/RECAL/RECO(NL)
      COMMON/HPOP/XN(6,NL),XN1,XN2,XN3
      common/rec_coll/alrec(30,30),collion(30,30),ction(14,27)
      common/debug/ideb
      dimension nionel(14)
      data nionel/1,2,6,7,8,10,11,12,13,14,16,18,20,26/
      dimension photok(14,26)
      dimension AA(nlp1,nlp1),X(nl)
      dimension dx(nl)

      do i=1,6
         x(i)=1.e-30
      enddo
      xion(ik,1,1) = x(1)
      
      xion(ik,1,2) = x(2) 
      
      xion(ik,2,1) = x(3) 
      
      xion(ik,2,2) = x(4) 
      
      xion(ik,2,3) = x(5)

      return

      end

      subroutine ctionrate(ik,te,ction)
      implicit real*8(a-h,o-z)
      PARAMETER (NL=340,NLP1=NL+1)
      include 'PARAM'
      common/ionx/xion(md,14,27)
      common/abl/abn(15)
      common/CTIon/ Ctionp(7,4,30)
      dimension nionel(14)
c      data nionel/1,2,6,7,8,8,2,3,4,14,10,6,2,15/
      data nionel/1,2,6,7,8,10,11,12,13,14,16,18,20,26/

      dimension ction(14,26)

      t4 = te/1.e4

      do iel=1,14

         do i=1,nionel(iel)

            ction(iel,i) = 0.
            
         enddo
         
      enddo

      xhii = xion(ik,1,2)*abn(1)
c!!
c      xhii=0.

c O I + H II FS 71

      ction(5,1) = 1.04E-9  * xhii

c Na I + H II

c charge transfer ionization of Na I with H II from 3S at 8000K
c rate from Natta & Giovardi ApJ 356, 646 -90

      ction(7,1) = 4.d-11*xhii
 

C     CT IONIZATION OF FE I BY H II
C     FROM WATSON IN LES HOUCHES 1974 ATOMIC AND MOLECULAR... P 257
C!    OBS! RATE IS AT 300 K
c!!  fudge factor skip!!
      fac=1.

      ction(14,1) = fac*7.4E-9*xhii

C     CT IONIZATION OF FE II BY H II
C     FROM D A NEUFELD AND A. DALGARNO, PR A 35, 3142 (1987)
      ction(14,2) = 2.3E-9*EXP(-3.0054/T4)*xhii

c include CT reactions from oak ridge data base
c this includes the reactions above

      do iel=3,14
         do ion=1,nionel(iel)

            if(ion.le.3)  then

               call convrtel(iel,nelem)

               if(CTIonp(1,ion,nelem).gt.0.) then
                  
                  ction(iel,ion) = HCTIon(ion,nelem,te)*xhii

               else
                  
                  ction(iel,ion) = 0.

               endif

            else

               ction(iel,ion)=0.

            endif
         enddo
      enddo

      return
      end

      subroutine convrtel(iel,nelem)

c convert from elem numbering to atomic number

c eg iel=5 and nelem=8 for O

      implicit real*8(a-h,o-z)
      if(iel.ge.3.and.iel.le.5) then
         nelem=iel+3
      elseif(iel.ge.6.and.iel.le.10) then
         nelem=iel+4
      elseif(iel.eq.11) then
         nelem=iel+5
      elseif(iel.eq.12) then
         nelem=iel+6
      elseif(iel.eq.13) then
         nelem=iel+7
      elseif(iel.eq.14) then
         nelem=iel+12
      endif
      return
      end


      SUBROUTINE RATE(te,XEL,photok,photolm)
C     ****************************************************************
C     *****
C     RATE CALCULATES PHOTOIONIZATION AND PHOTOHEATING RATES
C     OBS! FLUX(J) = MEAN INTENSITY (ERG/CM**2 EV)
C
C
C     *****
C     ****************************************************************
C
C      1 = H I      2 = HE I     3 = HE II    4 = O VI     5 = O VII
C      6 = O VIII   7 = O V      8 = C III    9 = C IV    10 = C V
C     11 = C VI    12 = C I     13 = C II    14 = O I     15 = O II
C     16 = O III   17 = O IV    18 = N I     19 = N II    20 = N III
C     21 = N IV    22 = N V     23 = N VI    24 = N VII   25 = SI I
C     26 = SI II   27 = SI III  28 = SI IV   29 = SI V    30 = SI VI
C     31 = SI VII  32 = SI VIII 33 = SI IX   34 = SI X    35 = SI XI
C     36 = SI XII  37 = SI XIII 38 = SI XIV  39 = MG I    40 = MG II
C     41 = MG III  42 = FE I    43 = FE II   44 = FE III  45 = FE IV
C     46 = AL I    47 = AL II   48 = AL III  49 = AL IV   50 = CA I
C     51 = CA II   52 = CA III  53 = NA I    54 = NA II   55 = NA III
C     56 = S I     57 = S II    58 = S III   59 = S IV    60 = S V
C     61 = S VI    62 = S VII   63 = S VIII  64 = S IX    65 = S X
C     66 = S XI    67 = S XII   68 = S XIII  69 = S XIV   70 = S XV
C     71 = S XVI   72 = NE I    73 = NE II   74 = NE III  75 = NE IV 
C     76 = NE V    77 = NE VI   78 = NE VII  79 = NE VIII 80 = NE IX 
C     81 = NE X    82 = AR I    83 = AR II   84 = FE V    85 = FE VI
C     86 = FE VII  87 = FE VIII 88 = FE IX   89 = FE X    90 = FE XI
C     91 = FE XII  92 = FE XIII 93 = FE XIV  94 = FE XV   95 = 
C     96 = AR III  97 = AR IV   98 = AR V    99 = AR VI  100 = AR VII
C     
C
      
      implicit real*8(a-h,o-z)
      save
c      parameter (md=350,mdp1=md+1)
      include 'param'
c      parameter (ne1=-100,ne2=250,ne3=ne2+1)
      parameter (nl=340,nlp1=nl+1)
      character*8 lab(200)
      common/fre/nint,jmin,jj
      common/ind/ik
      common/phq/ze(14,27),ge(14,27),zk(14,27)
      common/int/fl(2,ne1:ne2),si(14,27,ne1:ne2),e1(ne1:ne3),e(ne1:ne3)
      common/sik/sk(14,27,ne1:ne2)
      common/tres/ el(14,27),ek(14,27)
      common/spect/tel,fd(md,ne1:ne2),f0(ne1:ne2),ipar
      common/dtau/flux(ne1:ne2)
      common/chec/gec(14,27)
      common/hpop/xn(6,nl),xn1,xn2,xn3
      common/heiion/zheico,zheila,zotshi,zbalm,zhbalm
      common/abl/abn(15)
      common/ionx/xion(md,14,27)
      parameter(nlp=30000)
      common/linepoint/nlines,jline(nlp),iion(nlp)
      common/lineopac/totopl(nlp)
      common/linephoto/phline(nlp)
      common/fracphot/frphotline(14,27,nlp)
      common/timecheck/time,itime
      common/preion/ipre
      common/ionabun/iab(100)
      COMMON/RECAL/RECO(NL)
      common/debug/ideb
      dimension jlowl(14,27),jupl(14,27),jlowk(14,27),jupk(14,27),
     &                                    zet(14,27),zkt(14,27)
c    ions with ionization potentials less than 4 rydberg
c    and having l- or m-shells 
      dimension ionl4(33)
      dimension iaug(14,27)
      dimension photok(14,26),photolm(14,26)
       data ionl4/8, 12, 13, 14, 15, 18, 19, 20, 25, 26, 27, 28, 
     &     39, 40, 42, 43, 44, 46, 47, 48, 50, 51, 52, 53, 
     &     54, 56, 57, 58, 59, 72, 73, 82, 83/

      dimension nionel(14)
      data nionel/1,2,6,7,8,10,11,12,13,14,16,18,20,26/
      dimension zeold(14,27),gemax(14,27),flm(14,27)
      dimension jgm(14,27),lm(14,27)
      data pi/3.1415926e0/,elch/1.60219e-12/,amu/1.660531e-24/
      ikold=0

      call rateaug(1,1)
      call rateaug(2,2)
      call rateaug(6,6)
      call rateaug(7,7)
      call rateaug(8,8)
      call rateaug(10,10)
      call rateaug(12,3)
      call rateaug(14,14)
      call rateaug(16,16)
      call rateaug(18,7)
      call rateaug(26,26)
     
      zkhi=zk(1,1)

      do iel=1,14
         do iz=1,nionel(iel)
            zk(iel,iz) = 1.d-38
            zeold(iel,iz) =  ze(iel,iz) 
            ze(iel,iz) = 1.d-38
            gec(iel,iz) = 1.d-38
            ge(iel,iz) = 1.d-38
         enddo
      enddo

      IF(IK.NE.IKOLD) THEN
         IKOLD=IK
C
C     DETERMINE THE LOWER AND UPPER LIMIT IN ENERGY FOR THE INTEGRATION
C           OF THE PHOTOIONIZATION RATES
C
         do iel=1,14
            do iz=1,nionel(iel)
               
               ZKT(iel,iz)=0.
               ZET(iel,iz)=0.
               DO J=JMIN,JJ
                  ZA=4.*PI*FL(2,J)*SI(iel,iz,J)*(E(J+1)-E(J))/
     &                 (ELCH*E1(J))
                  ZB=4.*PI*FL(2,J)*SK(iel,iz,J)*(E(J+1)-E(J))/
     &                 (ELCH*E1(J))
                  ZKT(iel,iz)=ZKT(iel,iz)+ZB
                  ZET(iel,iz)=ZET(iel,iz)+ZA
               ENDDO
            ENDDO
         enddo


         do iel=1,14
            do iz=1,nionel(iel)

               ZA=0.
               ZB=0.
               ZABT=ZET(iel,iz)
               IF(ZABT.LE.0.) GOTO 23
               DO J=JMIN,JJ
                  ZA=4.*PI*FL(2,J)*SI(iel,iz,J)*(E(J+1)-E(J))/
     &                 (ELCH*E1(J))+ZA
                  IF(ZA.EQ.0.) JLOWL(iel,iz)=J
                  IF(ABS(ZA-ZABT)/ZABT.LT.1.E-5) GOTO 22
               ENDDO
 22            JUPL(iel,iz)=J
 23            CONTINUE
               ZABT=ZKT(iel,iz)
               IF(ZABT.LE.0.) GOTO 25
               DO J=JMIN,JJ
                  ZB=4.*PI*FL(2,J)*SK(IEL,IZ,J)*(E(J+1)-E(J))/
     &                 (ELCH*E1(J))+ZB
                  IF(ZB.EQ.0.) JLOWK(IEL,IZ)=J
                  IF(ABS(ZB-ZABT)/ZABT.LT.1.E-5) GOTO 24
               ENDDO
 24            JUPK(IEL,IZ)=J
 25            CONTINUE
            ENDDO
         ENDDO
      ENDIF

      GAA=0.0E0

      dphmax=0.

      do iel=1,14
         do iz=1,nionel(iel)
            GAA=0.0E0
            ZA=0.0E0
            ZB=0.E0
c!!!!!
            jlowl(iel,iz)=jmin
            jupl(iel,iz)=jj
            jlowk(iel,iz)=jmin
            jupk(iel,iz)=jj
            gemax(iel,iz)=0.
            DO J=JLOWL(IEL,IZ),JUPL(IEL,IZ)
C     ONLY GROUND STATES. 5 EV IS THE LOWEST IONIZATION POTENTIAL = NA
               IF(E(J).GT.5.0) THEN
                  ZA=4.*PI*FL(2,J)*SI(IEL,IZ,J)*(E(J+1)-E(J))/
     &                 (ELCH*E1(J))
                  ZE(IEL,IZ)=ZE(IEL,IZ)+ZA
                  if((iel.eq.-12.and.iz.eq.1).and.itime.ge.3500.and.
     &                 za.gt.0.) then
c                     write(6,928)iel,iz,j,e1(j),fl(2,j),si(iel,iz,j),
c     &                    sk(iel,iz,j),za,ze(iel,iz)
 928                 format(' Ar ze ',3i5,1pe12.3,10e12.3)
                  endif
                  IF(E1(J).GT.EL(IEL,IZ)) THEN
                     GAA=ELCH*ZA*(E1(J)-EL(IEL,IZ))
                  ELSE
                     GAA=0.0E0
                  ENDIF
                  GE(IEL,IZ)=GE(IEL,IZ)+GAA
                  if(iel.le.2) then
                     if(gaa.ge.gemax(iel,iz)) then
                        gemax(iel,iz)=gaa
                        jgm(iel,iz)=j
                        flm(iel,iz)=fl(2,j)
                     endif
                  endif
                   
               ENDIF
            ENDDO


c Auger electron energy

c?? !! check if this is ok 

            if(iz.le.nionel(iel)-3) then
               iaug(iel,iz) = 1
            else
               iaug(iel,iz) = 0
            endif

            if(iaug(iel,iz).eq.1) then
               enaug=0.
            elseif(iaug(iel,iz).eq.0) then
               enaug=ek(iel,iz)
            endif

            gemax(iel,iz)=0.

            DO J=JLOWK(IEL,IZ),JUPK(IEL,IZ)
C     ONLY GROUND STATES. 5 EV IS THE LOWEST IONIZATION POTENTIAL = NA
               IF(E(J).GT.5.0) THEN
                  if(iel.eq.1.and.iz.eq.1) then
C     COMPTON IONIZATION
C          SEE HALPERN AND GRINDLAY AP J 242:1042
C
                     COSTH=1.-511.E3*13.6/E1(J)**2
                     IF(COSTH.GT.-1.) THEN
                        sigmak=SK(1,1,J)+0.665E-24*(0.5+3.*COSTH/8.+
     &                       COSTH**3/8.)
                     else
                        sigmak=SK(1,1,J)
                     ENDIF
                  else
                     sigmak=SK(iel,iz,J)
                  endif
                  
                  ZB=4.*PI*FL(2,J)*sigmak*(E(J+1)-E(J))/
     &                 (ELCH*E1(J))
                  ZK(IEL,IZ)=ZK(IEL,IZ)+ZB
                  if (isnan(zk(iel,iz))) then
                     write(6,9169)iel,iz,j,fl(2,j),sigmak,zk(iel,iz)

 9169                format('iel,iz,j,fl(2,j),sigmak,zk(iel,iz) ',3i5,1pe12.3,10e12.3)
                     stop
                  endif

                  IF(E1(J).GT.EK(IEL,IZ)) THEN
                     GAA=ELCH*ZB*(E1(J)-enaug)
                  ELSE
                     GAA=0.0E0
                  ENDIF
                  GE(IEL,IZ)=GE(IEL,IZ)+GAA
                  if(iel.le.2) then
                     if(gaa.ge.gemax(iel,iz)) then
                        gemax(iel,iz)=gaa
                        jgm(iel,iz)=j
                        flm(iel,iz)=fl(2,j)
                     endif
                  endif
                  
               ENDIF
            ENDDO

c     include photoionization due to continous absorption of line
c     photons from rlossd
c     
c     
c     check if there are lines with energy higher than the 
c     ionization potential. i.e. treshold less than 4 Rydberg
            ilines=0
            phlines=0.
            phhlines=0.
c            do il=1,33
c               if(k.eq.ionl4(il)) then
c                  ilines=1
c               endif
c            enddo

c!!! include all ions for line photoionization

            ilines=1

            gemax(iel,iz)=0.

c     now add line photoionization
            if(ilines.eq.1.or.iel.le.2) then
               do l=1,nlines
                  if(jline(l).gt.jlowk(iel,iz).or.jline(l).gt.
     &                 jlowl(iel,iz)) then
c     first photoionization
                     dphl=phline(l)*(si(iel,iz,jline(l))+
     &                    sk(iel,iz,jline(l)))
                     if (isnan(dphl)) then
                        write(6,9178)iel,iz,l,jline(l),phline(l)
 9178                   format('iel,iz,l,jline(l),phline(l) ',4i5,1pe12.3)
                     endif
                     phlines=phlines+dphl

                     dphhll = 0.
                     dphhlk = 0.
                     if(e1(jline(l)).gt.el(iel,iz)) then
                        dphhll=elch*(e1(jline(l))-el(iel,iz))*phline(l)
     &                       *si(iel,iz,jline(l))
                     endif
c  now photo-heating
                     if(e1(jline(l)).gt.ek(iel,iz)) then
                        dphhlk=elch*(e1(jline(l))-enaug)*phline(l)
     &                       *sk(iel,iz,jline(l))
                     endif
                     phhlines=phhlines+dphhll+dphhlk
                     dphl=abn(iel)*xion(ik,iel,iz)*(dphhll+dphhlk)
                     if(dphl.gt.dphmax.and.iel.eq.1) then
                        ielmax=iel
                        ionmax=iz
                        dphmax=dphl
                        llmax=l
                        dphmaxl=dphhll
                        dphmaxk=dphhlk
                        e1max=e1(jline(l))
                     endif
                     if(dphl.gt.1.e-3*ze(iel,iz).or.dphl.gt.
     &                    1.e-10*zk(iel,iz)) then
c     if(dphl.gt.1.e-3*ze(k).or.dphl.gt.1.e-3*zk(k)) then
c     fraction of photoionization due to lines
                        if(iel.ge.3.and.ze(iel,iz).gt.1.e-20) then
                           frphotline(iel,iz,l)=dphl/ze(iel,iz)
                        endif
                        if(iel.le.2.and.zk(iel,iz).gt.1.e-20) then
                           frphotline(iel,iz,l)=dphl/zk(iel,iz)
                        endif
     
 923                    format(' line phot ',4i5,1pe12.4,20e12.4)
                     endif
                  endif

                  if(iel.le.2) then
                     if(phhlines.ge.gemax(iel,iz)) then
                        gemax(iel,iz)=phhlines
                        lm(iel,iz)=l
                     endif
                  endif

               enddo
               if(iel.ge.3) then
                  ze(iel,iz)=ze(iel,iz)+phlines
                  if(phlines.lt.0.) then
                     write(6,*)'phlines lt 0 ',k,l,phlines
                  endif
               elseif(iel.le.2) then
                  zk(iel,iz)=zk(iel,iz)+phlines
                  if (isnan(zk(iel,iz))) then
                     write(6,9179)iel,iz,phlines,zk(iel,iz)
 9179                format('iel,iz,phlines ',2i5,1pe12.3,10e12.3)
                     endif
               endif
c     if(k.eq.1.and.phlines.gt.0..and.itime.ge.3100) then
c     write(6,917)l,e1(j),phlines,ze(k)
c 917           format(' phlines(1) ',i5,1pe12.3,10e12.3)
c     endif
               
               ge(iel,iz)=ge(iel,iz)+phhlines

            endif

         enddo
      enddo



C     HE I PHOTOIONIZATION OF H I 584 AND 504 A. Check if inices ok!
      ZHEILA=4.*PI*FL(2,17)*SK(1,1,17)*(E(18)-E(17))/
     &     (ELCH*E1(17)*ZK(1,1))
      ZHEICO=4.*PI*FL(2,22)*SK(1,1,22)*(E(23)-E(22))/
     &     (ELCH*E1(22)*ZK(1,1))
      ZOTSHI=4.*PI*FL(2,2)*SK(1,1,2)*(E(3)-E(2))/
     &     (ELCH*E1(2)*ZK(1,1))
      IF(XN(5,1).GT.0.) ZHBALM=XN(5,2)*ZBALM/(XN(5,1)*ZK(1,1))

c
c 1 = H, 2 = He, 3 = C, 4 = N, 5 = O, 6 = Ne, 7 = Na, 8 = Mg, 
c 9 = Al, 10 = Si, 11 = S, 12 = Ar, 13 = Ca, 14 = Fe

      do iel=1,14
         do iz=1,nionel(iel)
            photok(iel,iz)=0.
            photolm(iel,iz)=0.
         enddo
      enddo
      do iel=1,14
         do iz=1,nionel(iel)
            photok(iel,iz)=zk(iel,iz)
            photolm(iel,iz)=ze(iel,iz)
         enddo
      enddo

      RETURN
      END

      SUBROUTINE RECR(ION,N,J,CHI,TE,RE)      
C
C     RECOMBINATION (USING THE MILNE RELATION)
C     AVERAGE THE EMISSION OVER THE INTERVAL SO THAT THE TOTAL EMISSION
C     IS CORRECT
C
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NL=340,NLP1=NL+1)
      include 'param'
c      PARAMETER (NE1=-100,NE2=250,NE3=NE2+1)
      COMMON/INT/FL(2,NE1:NE2),SI(14,27,NE1:NE2),E1(NE1:NE3),E(NE1:NE3)
      COMMON/A16/SIG(NL),GA(NL),BP(NL)
      COMMON/SIK/SK(14,27,NE1:NE2)
      COMMON/A26/SIGO(20,NE1:NE2),SIGCA(10,NE1:NE2),SIGHE(16,NE1:NE2),
     &     SIGFE(NL,NE1:NE2),sigfei(nl,ne1:ne2),SIGH(50,NE1:NE2)
      COMMON/NLEV/e000,NION,NH,NP1H,NHMAX,NHMIN
      COMMON/A14/C(NL,NL),CI(NL),G(NL),EN(NL),A(NL,NL),
     &     WL(NL,NL),DCDT(NL,NL),DCIDT(NL)
      COMMON/mdIND/kradius
      E00=E000-EN(N)
      CHI=E00
      IF(E(J).LT.E00) GOTO 250
      F2=(CHI-E1(J))*1.1609E4/TE
      DELT=E(J+1)-E(J)
      F3=DELT*1.1609E4/(2.*TE)
      F2AB=ABS(F2)
      F3AB=ABS(F3)
      IF(F3AB.GT.75.) GOTO 250
      IF(F2AB.GT.75.) GOTO 250
            IF(N.LE.6) THEN
              GAUNT=gbf(N,e1(J))
            ELSE
              GAUNT=1.
            ENDIF
            IF(ION.EQ.5) SIGM=1.e-18*SIG(N)*GAUNT*(E00/E1(J))**GA(N)
            IF(ION.EQ.2) SIGM=1.E-18*SIGOX(N,E1(J))
C           AVERAGE OVER THE ENERGY INTERVAL
            RE=G(N)*1.31159E-4*E1(J)**3.*(EXP(F2+F3)-EXP(F2-F3))
     &           *SIGM/(G(NH+1)*2.*F3*TE*SQRT(TE))
            if(kradius > 511) then
               write(6,9278)ion,n,j,e1(j),te,sigm,re
 9278          format('ion,n,j,te,sigm,re ',3i5,1pe12.3,10e12.3)
            endif
            GOTO 275
  250 RE=0.0E0
  275 CONTINUE
      if(re.gt.0.and.kradius>511) then
         write(6,922)ion,n,j,te,e00,sigm,gaunt,f2,f3,re
 0922    format('RE<0 ',3i5,1pe15.6,10e15.6)
      endif
      RETURN
      END

c$$$      SUBROUTINE RECEX(Z,J,N,TE,RECO)
c$$$      IMPLICIT REAL*8(A-H,O-Z)
c$$$c
c$$$c NOT USED! REPLACED BY RECR above       
c$$$C      ***********************************************************
c$$$C      *****
c$$$C     THIS ROUTINE CALCULATES THE REC. EMISSION DUE TO EXCITED STATES
c$$$C     USING A HYDROGENIC APPROXIMATION AND UNIT GAUNT FACTOR.
c$$$C     Z = CHARGE
c$$$C     J = ENERGY BIN
c$$$C     N = LEVEL
c$$$C      *****
c$$$C      ***********************************************************
c$$$      include 'param'
c$$$c      PARAMETER (NE1=-100,NE2=250,NE3=NE2+1)
c$$$      COMMON/INT/FL(2,NE1:NE2),SI(14,27,NE1:NE2),E1(NE1:NE3),E(NE1:NE3)
c$$$      COMMON/mdIND/kradius
c$$$      GAUNT=gbf(n,e1(j))
c$$$      RN=DBLE(N)
c$$$      RIZ=13.597*Z**2./RN**2.
c$$$      TEV=TE/1.1609E4
c$$$      EX=ABS((RIZ-E1(J))/TEV)
c$$$      RECO=0.
c$$$      IF(EX.LT.100.) THEN
c$$$            RECO=4.*3.14159*4.1478E-19*Z**4.*GAUNT*
c$$$     &                  EXP((RIZ-E1(J))/TEV)/(RN**3.*TE**1.5)
c$$$      ENDIF
c$$$      RETURN
c$$$      END

      SUBROUTINE XIONS(X,CHI)
C
C     CALCULATE IONIZATION POTENTIALS AND TOTAL FRACTIONS OF ALL IONS
C     FOR CALCULATION OF RECOMBINATION RADIATION
C
C     ENUMERATION OF THE IONS:
C
C      1 = H I      2 = HE I     3 = HE II    4 = O VI     5 = O VII
C      6 = O VIII   7 = O V      8 = C III    9 = C IV    10 = C V
C     11 = C VI    12 = C I     13 = C II    14 = O I     15 = O II
C     16 = O III   17 = O IV    18 = N I     19 = N II    20 = N III
C     21 = N IV    22 = N V     23 = N VI    24 = N VII   25 = SI I
C     26 = SI II   27 = SI III  28 = SI IV   29 = SI V    30 = SI VI
C     31 = SI VII  32 = SI VIII 33 = SI IX   34 = SI X    35 = SI XI
C     36 = SI XII  37 = SI XIII 38 = SI XIV  39 = MG I    40 = MG II
C     41 = MG III  42 = FE I    43 = FE II   44 = FE III  45 = FE IV
C     46 = AL I    47 = AL II   48 = AL III  49 = AL IV   50 = CA I
C     51 = CA II   52 = CA III  53 = NA I    54 = NA II   55 = NA III
C     56 = S I     57 = S II    58 = S III   59 = S IV    60 = S V
C     61 = S VI    62 = S VII   63 = S VIII  64 = S IX    65 = S X
C     66 = S XI    67 = S XII   68 = S XIII  69 = S XIV   70 = S XV
C     71 = S XVI   72 = NE I    73 = NE II   74 = NE III  75 = NE IV 
C     76 = NE V    77 = NE VI   78 = NE VII  79 = NE VIII 80 = NE IX 
C     81 = NE X    82 = AR I    83 = AR II   84 = FE V    85 = FE VI
C     86 = FE VII  87 = FE VIII 88 = FE IX   89 = FE X    90 = FE XI
C     91 = FE XII  92 = FE XIII 93 = FE XIV  94 = FE XV   
C     
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/TRES/EL(14,27),EK(14,27)
      COMMON/ABU/XA(2,100)
      common/abl/abn(15)
      include 'PARAM'
      common/ind/ik
      common/ionx/xion(md,14,27)
      DIMENSION X(100),CHI(100)
      CHI(1)=13.6
      CHI(2)=24.6
      CHI(3)=54.4
      CHI(5)=EK(5,7)
      CHI(6)=EK(5,8)
      CHI(10)=EK(3,5)
      CHI(11)=EK(3,6)
      CHI(23)=EK(4,6)
      CHI(24)=EK(4,7)
      chi(80)=ek(6,9)
      chi(81)=ek(6,10)

      X(1)=ABn(1)*XION(IK,1,2)
      X(2)=ABn(2)*XION(IK,2,2)
      X(3)=ABn(2)*XION(IK,2,3)
c O
      X(4)=ABN(5)*XION(IK,5,7)
      X(5)=ABN(5)*XION(IK,5,8)
      X(6)=ABN(5)*Xion(ik,5,9)
      X(7)=ABN(5)*XION(IK,5,6)
      X(14)=ABN(5)*XION(IK,5,2)
      X(15)=ABN(5)*XION(IK,5,3)
      X(16)=ABN(5)*XION(IK,5,4)
      X(17)=ABN(5)*XION(IK,5,5)
c C
      X(8)=ABN(3)*XION(IK,3,4)
      X(9)=ABN(3)*XION(IK,3,5)
      X(10)=ABN(3)*XION(IK,3,6)
      X(11)=ABN(3)*XION(IK,3,7)
      X(12)=ABN(3)*XION(IK,3,2)
      X(13)=ABN(3)*XION(IK,3,3)
c N
      X(18)=ABN(4)*XION(IK,4,2)
      X(19)=ABN(4)*XION(IK,4,3)
      X(20)=ABN(4)*XION(IK,4,4)
      X(21)=ABN(4)*XION(IK,4,5)
      X(22)=ABN(4)*XION(IK,4,6)
      X(23)=ABN(4)*XION(IK,4,7)
      X(24)=ABN(4)*XION(IK,4,8)
c Si
      do i=1,13
         x(i+24)=abn(10)*xion(ik,10,i+1)
      enddo
      x(38)=abn(10)*xion(ik,10,15)
c Mg
      X(39)=ABn(8)*XION(IK,8,2)
      X(40)=ABn(8)*XION(IK,8,3)
      X(41)=ABn(8)*XION(IK,8,4)
c Fe
      X(42)=ABn(14)*XION(IK,14,2)
      X(43)=ABn(14)*XION(IK,14,3)
      X(44)=ABn(14)*XION(IK,14,4)
      X(45)=ABn(14)*xion(ik,14,5)

      do i=5,15
         x(i+79)=abn(14)*xion(ik,14,i+1)
      enddo

c Al
      do i=1,4
         x(i+45)=abn(9)*xion(ik,9,i+1)
      enddo

c Ca
      X(50)=ABn(13)*XION(IK,13,2)
      X(51)=ABn(13)*XION(IK,13,3)
      X(52)=ABn(13)*XION(IK,13,4)
c Na
      X(53)=ABn(7)*XION(IK,7,2)
      X(54)=ABn(7)*XION(IK,7,3)
      X(55)=ABn(7)*XION(IK,7,4)
c S
      do i=1,16
         x(i+55)=abn(11)*xion(ik,11,i+1)
      enddo

c Ne

      do i=1,10
         x(i+71)=abn(6)*xion(ik,6,i+1)
      enddo

c Ar
      x(82)=abn(12)*xion(ik,12,2)
      x(83)=abn(12)*xion(ik,12,3)

      do i=3,7
         x(i+93)=abn(12)*xion(ik,12,i+1)
      enddo

c O III 1D
      x(95)=0.

      DO K=1,100
            IF(X(K).LT.0.) X(K)=0.
      ENDDO

      RETURN
      END

      double precision function gbf(n,e)
C
C     Calculation of bound-free gaunt factor from Mihalas, D., 
C     Ap.J. 149:169 (1967)
C     e in eV
C
      implicit real*8 (a-h,o-z)
      parameter(nl=6)
      dimension a0(nl),a1(nl),a2(nl),a3(nl),am1(nl),am2(nl),am3(nl)
      dimension wl0(nl)
      data wl0/915.3291d0,3704.9034d0,8504.7783d0,15560.594d0,
     &         25260.706d0,38194.187d0/
      data a0/1.2302628d0,1.1595421d0,1.1450949d0,1.1306695d0,
     &        1.1190904d0,1.1168376d0/
      data a1/-2.9094219d-3,-2.0735860d-3,-1.9366592d-3,-1.3482273d-3,
     &        -1.0401085d-3,-8.9466573d-4/
      data a2/7.3993579d-6,2.7033384d-6,2.3572356d-6,-4.6949424d-6,
     &        -6.9943488d-6,-8.8393113e-6/
      data a3/-8.7356966d-9,0.d0,0.d0,2.3548636d-8,2.8496742d-8,
     &        3.4696768d-8/
      data am1/-5.5759888d0,-1.2709045d0,-0.55936432d0,-0.31190730d0,
     &         -0.16051018d0,-0.13075417d0/
      data am2/12.803223d0,2.1325684d0,0.52471924d0,0.19683564d0,
     &         5.5545091d-2,4.1921183d-2/
      data am3/0.d0,-2.0244141d0,-0.23387146d0,-5.4418565d-2,
     &         -8.9182854d-3,-5.5303574d-3/
      wlmy=1.239854/e
      gbf=0.
      if(wlmy.lt.wl0(n)/1.d4) then
            x=1./wlmy
            gbf=a0(n)+a1(n)*x+a2(n)*x*x+a3(n)*x*x*x+am1(n)/x+
     &                                am2(n)/(x*x)+am3(n)/(x*x*x)
      endif
c mihalas approximation only good for > 50 AA.

      if(e.gt.250.) then
         gbf=1.
      endif

      return
      end


      SUBROUTINE NONTH(X,FH,FIH,FEH,FIHE,FEHE)
C     NON-THERMAL ELECTRON DEPOSITION FROM SHULL AND VAN STEENBERG (1985)
      IMPLICIT REAL*8(A-H,O-Z)
      IF(X.LT.0.95) THEN
       FH=0.9971*(1.-(1-X**.2663)**1.3163)
       FIH=0.3908*(1.-X**0.4092)**1.7592
       FIHE=0.0554*(1.-X**0.4614)**1.666
       FEH=0.4766*(1.-X**0.2735)**1.5221
       FEHE=0.0246*(1.-X**0.4049)**1.6594
      ELSE
       FH=1.
       FIH=0.
       FIHE=0.
       FEH=0.
       FEHE=0.
      ENDIF
      RETURN
      END      

c23456789012345678901234567890123456789012345678901234567890123456789012
