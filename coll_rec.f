C************************************************************           
C************************************************************           
      double precision function arf(a,b,c,d,pot,tev) 
C************************************************************           
C Collisional ionization according to Arnaud & Rothenflug.
C************************************************************           
      IMPLICIT REAL*8(A-H,O-Z) 
C
      y=0.0d0
      x=pot/tev
      if(x.gt.1.d2) goto 100
      f1=arf1(x)
      f2=arf2(x)
      y1=a*(1.d0-x*f1)+b*(1.d0+x-x*(2.d0+x)*f1)+c*f1+d*x*f2
      y=(6.69d-7*y1*dexp(-x))/(x*tev**1.5d0)
  100 continue
      arf=y
C
      return
      end
C************************************************************           
C************************************************************           
      double precision function arauli(iz,tev) 
C************************************************************           
C Auto-ionization according to Arnaud & Rothenflug (Li-sequence).
C************************************************************           
      IMPLICIT REAL*8(A-H,O-Z) 
C
      q=0.0d0
      z=dble(iz)
      potea=13.6d0*((z-0.835d0)**2-0.25d0*(z-1.62d0)**2)
      y=potea/tev
      if(y.gt.1.d2) goto 100
       zeff=z-0.43d0
       b=1.d0/(1.d0+2.d-4*z**3)
       f1=arf1(y)
       gy=2.22d0*f1+0.67d0*(1.d0-y*f1)+0.49d0*y*f1+1.2d0*y*(1.d0-
     & y*f1)
       q=(1.92d-7*b*gy*dexp(-y))/(zeff**2*(dsqrt(tev)))
       if(iz.eq.6) q=0.6d0*q
       if(iz.eq.7) q=0.8d0*q
       if(iz.eq.8) q=1.25d0*q
  100 continue
      arauli=q
C
      return
      end
C************************************************************           
C************************************************************           
      double precision function arau(iz,iseq,tev) 
C************************************************************           
C Auto-ionization according to Arnaud & Rothenflug (Except Li-
C sequence).
C************************************************************           
      IMPLICIT REAL*8(A-H,O-Z) 
C
      q=0.0d0
      z=dble(iz)
      if(iseq.eq.11.and.iz.le.16) potea=26.d0*(z-10.d0)
      if(iseq.eq.11.and.iz.gt.16) potea=11.d0*(z-10.d0)**1.5
      if(iseq.eq.12) potea=10.3d0*(z-10.d0)**1.52
      if(iseq.eq.13) potea=18.d0*(z-11.d0)**1.33
      if(iseq.eq.14) potea=18.4d0*(z-12.d0)**1.36
      if(iseq.eq.15) potea=23.7d0*(z-13.d0)**1.29
      if(iseq.eq.16) potea=40.1d0*(z-14.d0)**1.1
      if(iseq.eq.19) potea=29.d0
      if(iseq.eq.20) potea=25.d0
      if(iseq.eq.21) potea=73.d0
      if(iseq.eq.22) potea=60.d0
      y=potea/tev
      if(y.gt.1.d2) goto 100
       f1=arf1(y)
       if(iseq.eq.11) then
        if(iz.le.16) then
         a=2.8d-17/(z-10.d0)**0.7
         y1=1.d0-y*f1
        else
         a=1.3d-14/(z-10.d0)**3.73
         y1=1.d0-0.5d0*(y-y**2+f1*y**3)
        endif
       elseif(iseq.ge.12.and.iseq.le.16) then
        a=4.0d-13/(potea*z**2)
        y1=1.d0-0.5d0*(y-y**2+f1*y**3)
       elseif(iseq.ge.19) then
        if(iseq.eq.19) then
         b=1.12d0
         a=9.8d-17
        elseif(iseq.eq.20) then
         b=1.12d0
         a=6.0d-17
        elseif(iseq.eq.21) then
         b=1.0d0
         a=5.0d-18
        elseif(iseq.eq.22) then
         b=1.0d0
         a=1.8d-17
        endif
        y1=1.d0+b*f1
       endif
       q=(6.69d7*a*potea*y1*dexp(-y))/dsqrt(tev)
  100 continue
      arau=q
C
      return
      end
C************************************************************           
C************************************************************           
      double precision function arf1(x) 
C
      IMPLICIT REAL*8(A-H,O-Z) 
C
      if(x.le.0.02d0) then
       f1=exp(x)*(x-0.5772d0-dlog(x))
      elseif(x.gt.0.02d0.and.x.lt.10.d0) then
       a1=dlog((x+1.d0)/x)
       if(x.lt.1.5d0) then 
        a=-0.5d0
       else
        a=0.5d0
       endif
       a2=(0.36d0+0.03d0*(x+0.01d0)**a)/(x+1.d0)**2
       f1=a1-a2
      elseif(x.ge.10.d0) then
       f1=(1.d0/x)*(1.d0-1.d0/x+2.d0/x**2-6.d0/x**3+24.d0/x**4)
      endif
      arf1=f1
C
      return
      end
C************************************************************           
C************************************************************           
      double precision function arf2(x) 
C
      IMPLICIT REAL*8(A-H,O-Z) 
C
      dimension p(14),q(15)
      DATA p/1.0d0,2.1658d2,2.0336d4,1.0911d6,3.7114d7,8.3963d8,
     &1.2889d10,1.3449d11,9.4002d11,4.2571d12,1.1743d13,1.7549d13,
     &1.0806d13,4.9776d11/ 
      DATA q/1.0d0,2.1958d2,2.0984d4,1.1517d6,4.0349d7,9.4900d8,
     &1.5345d10,1.7182d11,1.3249d12,6.9071d12,2.3531d13,4.9432d13,
     &5.7760d13,3.0225d13,3.3641d12/
      sump=0.0d0
      sumq=0.0d0
       do 100 i=1,14
       y=dble(1-i)
       sump=sump+p(i)*x**y
       sumq=sumq+q(i)*x**y
  100  continue
      sumq=sumq+q(15)/x**1.4d1
      arf2=sump/(sumq*x**2)
C
      return
      end

C************************************************************           
C************************************************************           

      double precision function  rec_colli(ik,te,dens,xel)

C************************************************************           
C                                                                       
C  This function calculates all collisional rates and uses
C  these to calculate steady-state populations.
C************************************************************
C************************************************************        
C     Enumeration of species:
C
C     1 = H   2 = He  3 = C   4 = N   5 = O   6 = Ne  7 = Na
C     8 = Mg  9 = Al 10 = Si 11 = S  12 = Ar 13 = Ca 14 = Fe
C    15 = Ni
C
C     Number of ionization stages included for each specy:
C
C     H : 2 stages     He : 3 stages      C : 7 stages
C     N : 8 stages      O : 9 stages     Ne :11 stages
C    Na :12 stages     Mg :13 stages     Al :14 stages
C    Si :15 stages      S :17 stages     Ar :19 stages
C    Ca :21 stages     Fe :27 stages     Ni :29 stages
C************************************************************           
C                                                                       
      IMPLICIT REAL*8(A-H,O-Z)                                          
      character*120 chr
      parameter (mh=2,mhe=3,mc=7,mn=8,mo=9,mne=11,mna=12,mmg=13,
     &mal=14,msi=15,ms=17,mar=19,mca=21,mfe=27,mni=29,
     &mmat=mh+mhe+mc+mn+mo+mne+mna+mmg+mal+msi+ms+mar+mca+mfe+mni,
     &ispecy=15)
C                                                                       
c      COMMON/MSTUFF/XX(mmat+1),AM(mmat+1,mmat+1),EPS,NERR 
      COMMON/CONSTS/PI,PISQRT,CLIGHT 
      COMMON/SPOT/OTSP(16)
      COMMON/NUTSK/YI(ispecy+1),YEX1(ispecy+1),FI1(ispecy+1)
      COMMON/REC/AL2(16)  
c      COMMON/COLLIS/COLION(10)  
      COMMON/ABC/AL(16)
                                       
      common/abc1/alc(mc),aln(mn),alo(mo),alne(mne),alna(mna),
     &almg(mmg),alal(mal),alsi(msi),alsu(ms),alar(mar),alca(mca),
     &alfe(mfe),alni(mni)
      common/abc2/coh(1),cohe(2),coc(mc-1),con(mn-1),coo(mo-1),
     &cone(mne-1),cona(mna-1),comg(mmg-1),coal(mal-1),cosi(msi-1),
     &cosu(ms-1),coar(mar-1),coca(mca-1),cofe(mfe-1),coni(mni-1)
      common/abc3/chtrac(mc),chtran(mn),chtrao(mo)
c      include 'PARAM'
      include "parameters.h"
      common/ionx/xion(md,14,27)
      common/initr/initrec
      common/abl/abn(15)
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
c     2302 next
c      common/rec_coll/alrec(30,30),collion(30,30),cthr(14,27)
      common/rec_coll/alrec(30,30),collion(30,30),ction(14,27)
      common/gsreco/al_gs(26,26,10)
      dimension aldfe(27),cofeea(26),cofedir(26),alsi_nahar(15)
      common/rr_diel_badn/dielbadn(30,30),alrecbadn(30,30),altotbadn(30,30)
      common/hgsrec/algs_h1
c     2302 next
c      dimension cther(14,26)      
      dimension cthr(14,26),cther(14,26)
      dimension tlr(100),recc(10,100),recn(10,100),reco(10,100),recsi(10,100)
      save tlr,recc,recn,reco,recsi
      dimension nion(14)
      data nion/1,2,6,7,8,10,11,12,13,14,16,18,20,26/
c      write(6,*)' in coll rec '
c      write(6,*)'al_gs(14,1,1),al_gs(16,1,1) collrec ',al_gs(14,1,1),al_gs(16,1,1)
      tq=1.e6           
c
      TELOG=DLOG10(TE)                                                  
      TEV=TE/1.1603D4
      T4=TE/1.D4                                                        
C
      DO I=1,mc-1     
         COC(I)=0.0D0
         CHTRAC(I)=0.0D0
         ALC(I)=0.0D0
      enddo
      ALC(mc)=0.0d0
      CHTRAC(mc)=0.0d0
      DO I=1,mn-1     
         CON(I)=0.0D0
         CHTRAN(I)=0.0D0
         ALN(I)=0.0D0
      enddo
      ALN(mn)=0.0d0
      CHTRAN(mn)=0.0d0
      DO  I=1,mo-1     
         COO(I)=0.0D0
         CHTRAO(I)=0.0D0
         ALO(I)=0.0D0
      enddo
      ALO(mo)=0.0d0
      CHTRAO(mo)=0.0d0
      DO I=1,mne-1     
         CONE(I)=0.0D0
         ALNE(I)=0.0D0
      enddo
      ALNE(mne)=0.0d0
      DO I=1,mna-1     
         CONA(I)=0.0D0
         ALNA(I)=0.0D0
      enddo
      ALNA(mna)=0.0d0
      DO I=1,mmg-1     
         COMG(I)=0.0D0
         ALMG(I)=0.0D0
      enddo
      ALMG(mmg)=0.0d0
      DO  I=1,mal-1     
         COAL(I)=0.0D0
         ALAL(I)=0.0D0
      enddo
      ALAL(mal)=0.0d0
      DO I=1,msi-1     
         COSI(I)=0.0D0
         ALSI(I)=0.0D0
      enddo
      ALSI(msi)=0.0d0
      DO I=1,ms-1     
         COSU(I)=0.0D0
         ALSU(I)=0.0D0
      enddo
      ALSU(ms)=0.0d0
      do  I=1,mar-1     
         COAR(I)=0.0D0
         ALAR(I)=0.0D0
      enddo
      ALAR(mar)=0.0d0
      DO I=1,mca-1     
         COCA(I)=0.0D0
         ALCA(I)=0.0D0
      enddo
      ALCA(mca)=0.0d0
      DO I=1,mfe-1     
         COFE(I)=0.0D0
         ALFE(I)=0.0D0
      enddo
      ALFE(mfe)=0.0d0
      DO I=1,mni-1     
         CONI(I)=0.0D0
         ALNI(I)=0.0D0
      enddo
      ALNI(mni)=0.0d0
      DO I=1,ispecy+1
         OTSP(I)=2.0D0
         YEX1(I)=0.0D0
         AL2(I)=0.0D0
         AL(I)=0.0D0
         FI1(I)=0.0D0
      enddo
c or HI, He I, He II use OTS!      
      OTSP(1)=0.
      OTSP(2)=0.
      OTSP(3)=0.
      YI(1)=1.5789D5/TE                                                 
      YI(2)=YI(1)                                                       
      YI(3)=4.*YI(1)                                                    
      YI(4)=36.*YI(1)                                                   
      YI(5)=49.*YI(1)                                                   
      YI(6)=64.*YI(1)                                                   
      YI(7)=100.*YI(1)                                                  
      YI(8)=121.*YI(1)                                                  
      YI(9)=144.*YI(1)                                                  
      YI(10)=169.*YI(1)                                                  
      YI(11)=196.*YI(1)                                                  
      YI(12)=256.*YI(1)                                                  
      YI(13)=324.*YI(1)                                                  
      YI(14)=400.*YI(1)                                                  
      YI(15)=676.*YI(1)                                                  
      YI(16)=784.*YI(1)                                                  

C     ***************************************************************   
C     *****                                                             
C     RECOMBINATION COEFF. FOR HYDROGEN LIKE IONS AND HELIUM            
C     SEATON (1959) AND GOULD AND TAKHUR (1966)                         
C     RATES ARE ADJUSTED TO AGREE WITH OSTERBROCK (74) FOR T=1E4        
C     (On-the-spot if OTSP.LE.1.)                                       
C     FI is the asymptotic TOTAL recombination term, accurate 
C     to better than 3 % for T*Z**2 < 1E6. 
C     FI2 is the term for recombination to n=1 only. 
C     *****                                                             
C     ***************************************************************   
C                                                                       
      DO 922 K1=1,16                                                     
      Y=YI(K1)                                                          
      FI=(0.4288+0.5*DLOG(Y)+0.469/(Y**.33333))*0.99D0                  
C     DEXP. INTEGRAL APPROX. FOR RECOMB. DIRECTLY TO N=1                
C  Y < 1                                                                
      if(y.gt.1.) goto 926                                              
      FI1(K1)=-.5772-DLOG(Y)+Y-Y**2/4.+Y**3/18.-Y**4/96.+               
     &Y**5/600.                                                         
      FI2=Y*DEXP(Y)*FI1(K1)                                             
      goto 920                                                          
C  1 < Y < 20                                                           
  926 FI2=(Y**2.+2.3347*Y+.2506)/(Y**2.+3.3306*Y+1.6815)                
      IF(Y.GT.20.) GOTO 3114                                            
      FI3=-.17501+.54157/Y-2.4153/(Y*Y)+11.428/Y**3.-36.749/Y**4.       
     &+70.31/Y**5.-75.169/Y**6.+40.428/Y**7.-8.2368/Y**8.               
      GOTO 3112                                                         
C  20 < Y         Seaton 1959 Eq 26. FI3 = x S(1)                                                      
 3114 FI3=-0.1728*(1.-8./(3.*Y)+70./(9.*Y*Y)-800./(27.*Y**3.)           
     &+11440./(81.*Y**4.))                                              
      FI4=-0.0496*(1.-1./Y)                                             
      FI3=FI3+FI4                                                       
 3112 CONTINUE                                                          
      YEX1(K1)=FI2                                                      
      FI2=FI2+FI3                                                       
      IF(Y.GT.100.) GOTO 907                                            
      FI1(K1)=FI2                                                       
      GOTO 910                                                          
  907 FI1(K1)=0.0                                                       
  910 CONTINUE                                                          
  920 IF(K1.EQ.2) GOTO 921                                              
      AL2(K1)=1.38D-16*DSQRT(TE)*Y*(FI-FI2)                             
      IF(K1.EQ.1) AL2(1)=AL2(1)*0.994D0*(TE/1.D4)**(0.03354)            
      IF(OTSP(K1).LE.1.) GOTO 815                                       
      FI2=0.0D0                                                         
      AL(K1)=1.34D-16*DSQRT(TE)*Y*(FI-FI2)                              
      GOTO 8915                                                         
  815 AL(K1)=AL2(K1)                                                    
 8915 CONTINUE                                                          
      GOTO 922                                                          
C     HELIUM I                                                          
  921 AL2(2)=1.38D-16*DSQRT(TE)*Y*(FI-FI2)*(TE/1.D4)**0.1375            
      AL(2)=AL2(2)
      algs_h1=1.33D-16*Y*FI2*DSQRT(TE)      
      IF(OTSP(2).LE.1.) GOTO 923                                        
      AL(2)=AL2(2)+1.33D-16*Y*FI2*DSQRT(TE)      
 923  CONTINUE

C                                                                       
C     DIELECTRONIC RECOMBINATION (ALDROV. AND PEQU.)                    
C                                                                       
      ALDHE=0.                                                          
      IF(TE.LT.5000.) GOTO 945                                          
      ALA=1.+0.3*DEXP(-9.4D4/TE)                                        
      ALDHE=1.9D-3*DEXP(-4.7D5/TE)*ALA/(TE*DSQRT(TE))                   
      GOTO 946                                                          
  945 ALDHE=0.0                                                         
  946 CONTINUE                                                          
  922 CONTINUE       

      alchh=al2(1)
      alchhe=al2(2)
      alchheii=al2(3)

      alrec(1,2)=al(1)
      alrec(2,2)=al(2)
      alrec(2,3)=al(3)


C     ***************************************************************   
C     *****                                                             
C     COLLISIONAL IONIZATION RATES INCLUDING AUTOIONIZATION
C     (Arnaud & Rothenflug)
C     *****                                                             
C     ***************************************************************   

C -Carbon-
      coc(1)=arf(6.d0,-16.d0,12.d0,-15.1d0,11.3d0,tev)
      coc(1)=coc(1)+arf(24.3d0,-7.8d0,2.5d0,-24.d0,16.6d0,tev)
      coc(2)=arf(16.d0,-9.d0,2.5d0,-10.5d0,24.4d0,tev)
      coc(2)=coc(2)+arf(23.7d0,-7.6d0,2.5d0,-21.7d0,30.9d0,tev)
      coc(3)=arf(23.2d0,-7.4d0,2.5d0,-19.4d0,47.9d0,tev)
      coc(3)=coc(3)+arf(20.d0,-5.6d0,4.1d0,-18.d0,325.d0,tev)
      coc(4)=arf(8.2d0,-2.7d0,1.4d0,-6.6d0,64.5d0,tev)
      coarn=coc(4)
      coc(4)=coc(4)+arf(20.d0,-5.6d0,4.1d0,-18.d0,343.d0,tev)
      coc(5)=arf(20.4d0,-6.1d0,4.5d0,-18.d0,392.d0,tev)
      coc(6)=arf(12.2d0,-3.9d0,1.9d0,-10.3d0,490.d0,tev)
C -Nitrogen-
      con(1)=arf(19.5d0,-30.5d0,15.d0,-29.d0,14.5d0,tev)
      con(1)=con(1)+arf(19.d0,-4.5d0,2.8d0,-20.2d0,20.3d0,tev)
      con(2)=arf(21.d0,-9.d0,5.3d0,-22.5d0,29.6d0,tev)
      con(2)=con(2)+arf(18.5d0,-4.3d0,2.8d0,-18.d0,36.7d0,tev)
      con(3)=arf(16.d0,-7.5d0,2.3d0,-10.d0,47.4d0,tev)
      con(3)=con(3)+arf(18.1d0,-4.d0,2.8d0,-15.8d0,55.8d0,tev)
      con(4)=arf(17.6d0,-3.8d0,2.8d0,-13.6d0,77.5d0,tev)
      con(4)=con(4)+arf(20.5d0,-5.8d0,4.1d0,-18.d0,471.d0,tev)
      con(5)=arf(10.5d0,-3.3d0,1.4d0,-7.7d0,97.9d0,tev)
      con(5)=con(5)+arf(20.5d0,-5.8d0,4.1d0,-18.d0,493.d0,tev)
      con(6)=arf(20.8d0,-6.3d0,4.4d0,-18.2d0,552.d0,tev)
      con(7)=arf(12.3d0,-4.d0,1.9d0,-10.3d0,667.d0,tev)
C -Oxygen-
      coo(1)=arf(9.5d0,-17.5d0,12.5d0,-19.5d0,13.6d0,tev)
      coo(1)=coo(1)+arf(18.2d0,-4.d0,2.8d0,-20.2d0,28.5d0,tev)
      coo(2)=arf(25.d0,-8.d0,8.4d0,-29.5d0,35.1d0,tev)
      coo(2)=coo(2)+arf(17.8d0,-3.8d0,2.9d0,-18.1d0,42.6d0,tev)
      coo(3)=arf(25.d0,-7.d0,5.d0,-18.d0,54.9d0,tev)
      coo(3)=coo(3)+arf(17.3d0,-3.5d0,2.9d0,-16.1d0,63.8d0,tev)
      coo(4)=arf(15.d0,-5.d0,2.2d0,-10.5d0,77.4d0,tev)
      coo(4)=coo(4)+arf(16.8d0,-3.3d0,2.8d0,-14.1d0,87.6d0,tev)
      coo(5)=arf(16.4d0,-3.d0,2.9d0,-12.d0,114.d0,tev)
      coo(5)=coo(5)+arf(20.8d0,-6.d0,4.1d0,-18.d0,644.d0,tev)
      coo(6)=arf(10.4d0,-3.3d0,1.4d0,-7.4d0,138.d0,tev)
      coo(6)=coo(6)+arf(20.8d0,-6.d0,4.1d0,-18.d0,670.d0,tev)
      coo(7)=arf(21.2d0,-6.5d0,4.3d0,-18.4d0,739.d0,tev)
      coo(8)=arf(12.3d0,-4.d0,1.9d0,-10.3d0,871.d0,tev)
C -Neon-
      cone(1)=arf(40.d0,-42.d0,18.d0,-56.d0,21.6d0,tev)
      cone(1)=cone(1)+arf(19.d0,-4.9d0,2.8d0,-22.d0,48.5d0,tev)
      cone(2)=arf(37.d0,-33.d0,15.5d0,-46.d0,41.1d0,tev)
      cone(2)=cone(2)+arf(18.6d0,-4.6d0,2.8d0,-20.2d0,66.4d0,tev)
      cone(3)=arf(33.d0,-17.5d0,11.2d0,-33.d0,63.5d0,tev)
      cone(3)=cone(3)+arf(18.2d0,-4.4d0,2.8d0,-18.4d0,86.2d0,tev)
      cone(4)=arf(34.d0,-10.d0,7.5d0,-25.d0,97.1d0,tev)
      cone(4)=cone(4)+arf(17.8d0,-4.d0,2.8d0,-16.7d0,108.d0,tev)
      cone(5)=arf(25.5d0,-8.5d0,4.5d0,-16.8d0,126.d0,tev)
      cone(5)=cone(5)+arf(17.4d0,-3.8d0,2.8d0,-14.9d0,139.d0,tev)
      cone(6)=arf(14.5d0,-4.6d0,1.9d0,-8.5d0,158.d0,tev)
      cone(6)=cone(6)+arf(16.9d0,-3.4d0,2.8d0,-13.2d0,172.d0,tev)
      cone(7)=arf(16.5d0,-3.1d0,2.8d0,-11.4d0,207.d0,tev)
      cone(7)=cone(7)+arf(21.5d0,-6.4d0,4.1d0,-18.d0,644.d0,tev)
      cone(8)=arf(10.1d0,-3.1d0,1.4d0,-7.1d0,239.d0,tev)
      cone(8)=cone(8)+arf(21.5d0,-6.4d0,4.1d0,-18.d0,1107.d0,tev)
      cone(9)=arf(21.9d0,-6.8d0,4.2d0,-18.7d0,1196.d0,tev)
      cone(10)=arf(12.5d0,-4.1d0,1.9d0,-10.4d0,1362.d0,tev)
C -Sodium-
      cona(1)=arf(16.d0,-1.d0,0.2d0,-13.5d0,5.1d0,tev)
      cona(1)=cona(1)+arf(63.9d0,-27.d0,33.d0,-80.d0,34.d0,tev)
      cona(2)=arf(40.d0,-28.d0,19.4d0,-60.d0,47.3d0,tev)
      cona(2)=cona(2)+arf(19.2d0,-5.3d0,2.8d0,-21.2d0,80.1d0,tev)
      cona(3)=arf(50.1d0,-20.2d0,14.8d0,-41.7d0,71.7d0,tev)
      cona(3)=cona(3)+arf(18.8d0,-5.d0,2.8d0,-19.6d0,102.d0,tev)
      cona(4)=arf(43.3d0,-16.3d0,10.7d0,-33.4d0,141.d0,tev)
      cona(4)=cona(4)+arf(18.4d0,-4.7d0,2.8d0,-18.d0,126.d0,tev)
      cona(5)=arf(35.1d0,-12.4d0,7.2d0,-25.1d0,138.d0,tev)
      cona(5)=cona(5)+arf(18.d0,-4.3d0,2.8d0,-16.3d0,151.d0,tev)
      cona(6)=arf(25.5d0,-8.5d0,4.2d0,-16.8d0,172.d0,tev)
      cona(6)=cona(6)+arf(17.6d0,-4.d0,2.8d0,-14.7d0,186.d0,tev)
      cona(7)=arf(14.5d0,-4.6d0,1.8d0,-8.5d0,208.d0,tev)
      cona(7)=cona(7)+arf(17.2d0,-3.7d0,2.8d0,-13.1d0,224.d0,tev)
      cona(8)=arf(16.8d0,-3.4d0,2.8d0,-11.4d0,264.d0,tev)
      cona(8)=cona(8)+arf(21.7d0,-6.5d0,4.1d0,-18.d0,1328.d0,tev)
      cona(9)=arf(10.d0,-3.d0,1.4d0,-6.9d0,300.d0,tev)
      cona(9)=cona(9)+arf(21.7d0,-6.5d0,4.1d0,-18.d0,1366.d0,tev)
      cona(10)=arf(22.2d0,-7.d0,4.2d0,-18.8d0,1465.d0,tev)
      cona(11)=arf(12.5d0,-4.1d0,1.9d0,-10.4d0,1649.d0,tev)
C -Magnesium-
      comg(1)=arf(18.d0,-1.d0,0.6d0,-4.d0,7.6d0,tev)
      comg(1)=comg(1)+arf(37.7d0,-30.d0,24.8d0,-62.d0,54.d0,tev)
      comg(1)=comg(1)+arf(17.6d0,-5.2d0,3.3d0,-19.d0,92.2d0,tev)
      comg(2)=arf(9.d0,-3.6d0,0.3d0,-5.4d0,15.d0,tev)
      comg(2)=comg(2)+arf(37.7d0,-30.d0,24.8d0,-62.d0,65.d0,tev)
      comg(2)=comg(2)+arf(17.6d0,-5.2d0,3.3d0,-19.d0,104.5d0,tev)
      comg(3)=arf(55.5d0,-24.1d0,18.7d0,-65.d0,80.1d0,tev)
      comg(3)=comg(3)+arf(19.3d0,-5.6d0,2.8d0,-20.5d0,119.d0,tev)
      comg(4)=arf(50.1d0,-20.2d0,14.2d0,-41.7d0,109.d0,tev)
      comg(4)=comg(4)+arf(19.d0,-5.3d0,2.8d0,-19.d0,144.d0,tev)
      comg(5)=arf(43.3d0,-16.3d0,10.3d0,-33.4d0,141.d0,tev)
      comg(5)=comg(5)+arf(18.6d0,-4.9d0,2.8d0,-17.5d0,172.d0,tev)
      comg(6)=arf(35.1d0,-12.4d0,6.9d0,-25.1d0,187.d0,tev)
      comg(6)=comg(6)+arf(18.2d0,-4.6d0,2.8d0,-16.d0,201.d0,tev)
      comg(7)=arf(25.5d0,-8.5d0,4.1d0,-16.8d0,225.d0,tev)
      comg(7)=comg(7)+arf(18.d0,-4.3d0,2.8d0,-14.5d0,241.d0,tev)
      comg(8)=arf(14.5d0,-4.6d0,1.8d0,-8.5d0,266.d0,tev)
      comg(8)=comg(8)+arf(17.5d0,-4.d0,2.8d0,-13.d0,283.d0,tev)
      comg(9)=arf(17.1d0,-3.6d0,2.7d0,-11.5d0,328.d0,tev)
      comg(9)=comg(9)+arf(22.d0,-6.7d0,4.1d0,-18.d0,1611.d0,tev)
      comg(10)=arf(10.d0,-3.d0,1.4d0,-6.8d0,367.d0,tev)
      comg(10)=comg(10)+arf(22.d0,-6.7d0,4.1d0,-18.d0,1653.d0,tev)
      comg(11)=arf(22.4d0,-7.1d0,4.1d0,-18.9d0,1762.d0,tev)
      comg(12)=arf(12.6d0,-4.2d0,1.9d0,-10.4d0,1963.d0,tev)
C -Aluminium-
      coal(1)=arf(47.d0,-26.d0,0.6d0,-39.d0,6.d0,tev)
      coal(1)=coal(1)+arf(55.1d0,-37.2d0,1.4d0,-41.d0,10.6d0,tev)
      coal(2)=arf(17.d0,-6.d0,1.d0,-8.d0,18.8d0,tev)
      coal(2)=coal(2)+arf(31.3d0,-22.7d0,21.d0,-44.1d0,90.d0,tev)
      coal(2)=coal(2)+arf(12.1d0,-3.5d0,3.3d0,-13.1d0,131.d0,tev)
      coal(3)=arf(6.3d0,-2.4d0,0.5d0,-4.1d0,28.4d0,tev)
      coal(3)=coal(3)+arf(31.3d0,-22.7d0,21.d0,-44.1d0,103.d0,tev)
      coal(3)=coal(3)+arf(12.1d0,-3.5d0,3.3d0,-13.1d0,145.6d0,tev)
      coal(4)=arf(72.d0,-24.1d0,18.d0,-50.d0,120.d0,tev)
      coal(4)=coal(4)+arf(19.5d0,-5.9d0,2.8d0,-19.8d0,164.d0,tev)
      coal(5)=arf(60.8d0,-20.2d0,13.7d0,-41.7d0,154.d0,tev)
      coal(5)=coal(5)+arf(19.1d0,-5.5d0,2.8d0,-18.4d0,194.d0,tev)
      coal(6)=arf(49.5d0,-16.3d0,9.9d0,-33.4d0,190.d0,tev)
      coal(6)=coal(6)+arf(18.9d0,-5.2d0,2.8d0,-17.1d0,225.d0,tev)
      coal(7)=arf(38.3d0,-12.4d0,6.7d0,-25.1d0,241.d0,tev)
      coal(7)=coal(7)+arf(18.4d0,-4.8d0,2.8d0,-15.7d0,258.d0,tev)
      coal(8)=arf(27.d0,-8.5d0,3.9d0,-16.8d0,285.d0,tev)
      coal(8)=coal(8)+arf(18.2d0,-4.5d0,2.8d0,-14.3d0,302.d0,tev)
      coal(9)=arf(14.d0,-4.6d0,1.7d0,-8.5d0,330.d0,tev)
      coal(9)=coal(9)+arf(17.9d0,-4.1d0,2.8d0,-13.d0,350.d0,tev)
      coal(10)=arf(17.4d0,-3.8d0,2.7d0,-11.6d0,399.d0,tev)
      coal(10)=coal(10)+arf(22.2d0,-6.8d0,4.1d0,-18.d0,1921.d0,tev)
      coal(11)=arf(9.9d0,-3.d0,1.4d0,-6.7d0,442.d0,tev)
      coal(11)=coal(11)+arf(22.2d0,-6.8d0,4.1d0,-18.d0,1967.d0,tev)
      coal(12)=arf(22.7d0,-7.2d0,4.1d0,-19.d0,2086.d0,tev)
      coal(13)=arf(12.6d0,-4.2d0,1.9d0,-10.4d0,2304.d0,tev)

C -Silicon-
      cosi(1)=arf(74.5d0,-49.4d0,1.3d0,-54.6d0,8.1d0,tev)
      cosi(1)=cosi(1)+arf(53.8d0,-35.8d0,1.4d0,-40.7d0,13.5d0,tev)
      cosi(2)=arf(50.4d0,-33.4d0,0.6d0,-36.9d0,16.3d0,tev)
      cosi(2)=cosi(2)+arf(55.1d0,-37.2d0,1.4d0,-41.d0,22.9d0,tev)
      cosi(3)=arf(19.8d0,-5.7d0,1.3d0,-11.9d0,33.5d0,tev)
      cosi(3)=cosi(3)+arf(66.7d0,-24.8d0,18.7d0,-65.d0,133.d0,tev)
      cosi(3)=cosi(3)+arf(22.d0,-7.2d0,3.3d0,-20.9d0,176.6d0,tev)
      cosi(4)=arf(9.d0,-3.d0,0.6d0,-5.8d0,45.1d0,tev)
      cosi(4)=cosi(4)+arf(66.7d0,-24.8d0,18.7d0,-65.d0,148.d0,tev)
      cosi(4)=cosi(4)+arf(22.d0,-7.2d0,3.3d0,-20.9d0,193.5d0,tev)
      cosi(5)=arf(72.d0,-24.1d0,17.4d0,-50.d0,167.d0,tev)
      cosi(5)=cosi(5)+arf(19.6d0,-6.2d0,2.8d0,-19.d0,217.d0,tev)
      cosi(6)=arf(60.8d0,-20.2d0,13.2d0,-41.7d0,205.d0,tev)
      cosi(6)=cosi(6)+arf(19.3d0,-5.8d0,2.8d0,-17.8d0,250.d0,tev)
      cosi(7)=arf(49.5d0,-16.3d0,9.6d0,-33.4d0,246.d0,tev)
      cosi(7)=cosi(7)+arf(19.d0,-5.4d0,2.8d0,-16.6d0,285.d0,tev)
      cosi(8)=arf(38.3d0,-12.4d0,6.4d0,-25.1d0,303.d0,tev)
      cosi(8)=cosi(8)+arf(18.6d0,-5.1d0,2.8d0,-15.4d0,321.d0,tev)
      cosi(9)=arf(27.d0,-8.5d0,3.8d0,-16.8d0,351.d0,tev)
      cosi(9)=cosi(9)+arf(18.3d0,-4.7d0,2.8d0,-14.1d0,371.d0,tev)
      cosi(10)=arf(14.d0,-4.6d0,1.6d0,-8.5d0,401.d0,tev)
      cosi(10)=cosi(10)+arf(18.d0,-4.3d0,2.8d0,-12.9d0,423.d0,tev)
      cosi(11)=arf(17.7d0,-4.d0,2.7d0,-11.7d0,476.d0,tev)
      cosi(11)=cosi(11)+arf(22.4d0,-6.9d0,4.1d0,-18.d0,2259.d0,tev)
      cosi(12)=arf(9.8d0,-2.9d0,1.4d0,-6.6d0,523.d0,tev)
      cosi(12)=cosi(12)+arf(22.4d0,-6.9d0,4.1d0,-18.d0,2309.d0,tev)
      cosi(13)=arf(22.9d0,-7.3d0,4.d0,-19.1d0,2438.d0,tev)
      cosi(14)=arf(12.7d0,-4.3d0,1.9d0,-10.4d0,2673.d0,tev)
C -Sulphur-
      cosu(1)=arf(6.d0,-22.d0,20.d0,-20.d0,10.4d0,tev)
      cosu(1)=cosu(1)+arf(51.3d0,-33.2d0,1.4d0,-40.2d0,20.2d0,tev)
      cosu(2)=arf(98.7d0,-65.4d0,1.9d0,-72.3d0,23.4d0,tev)
      cosu(2)=cosu(2)+arf(52.5d0,-34.5d0,1.4d0,-40.5d0,30.7d0,tev)
      cosu(3)=arf(74.5d0,-49.4d0,1.3d0,-54.6d0,35.d0,tev)
      cosu(3)=cosu(3)+arf(53.8d0,-35.8d0,1.4d0,-40.7d0,43.8d0,tev)
      cosu(4)=arf(50.4d0,-33.4d0,0.6d0,-36.9d0,47.3d0,tev)
      cosu(4)=cosu(4)+arf(55.1d0,-37.2d0,1.4d0,-41.d0,57.6d0,tev)
      cosu(5)=arf(19.8d0,-5.7d0,1.6d0,-11.9d0,72.7d0,tev)
      cosu(5)=cosu(5)+arf(73.2d0,-27.d0,15.8d0,-61.1d0,239.d0,tev)
      cosu(5)=cosu(5)+arf(23.1d0,-8.d0,3.3d0,-19.5d0,288.2d0,tev)
      cosu(6)=arf(9.d0,-2.8d0,0.7d0,-5.4d0,88.1d0,tev)
      cosu(6)=cosu(6)+arf(73.2d0,-27.d0,15.8d0,-61.1d0,257.d0,tev)
      cosu(6)=cosu(6)+arf(23.1d0,-8.d0,3.3d0,-19.5d0,309.7d0,tev)
      cosu(7)=arf(72.d0,-24.1d0,14.2d0,-50.d0,281.d0,tev)
      cosu(7)=cosu(7)+arf(19.6d0,-6.8d0,2.8d0,-17.5d0,343.d0,tev)
      cosu(8)=arf(60.8d0,-20.2d0,10.9d0,-41.7d0,328.d0,tev)
      cosu(8)=cosu(8)+arf(19.3d0,-6.3d0,2.8d0,-16.6d0,384.d0,tev)
      cosu(9)=arf(49.5d0,-16.3d0,8.d0,-33.4d0,379.d0,tev)
      cosu(9)=cosu(9)+arf(19.1d0,-5.9d0,2.8d0,-15.6d0,426.d0,tev)
      cosu(10)=arf(38.3d0,-12.4d0,5.5d0,-25.1d0,447.d0,tev)
      cosu(10)=cosu(10)+arf(18.8d0,-5.5d0,2.8d0,-14.7d0,469.d0,tev)
      cosu(11)=arf(27.d0,-8.5d0,3.3d0,-16.8d0,505.d0,tev)
      cosu(11)=cosu(11)+arf(18.6d0,-5.1d0,2.8d0,-13.7d0,528.d0,tev)
      cosu(12)=arf(14.d0,-4.6d0,1.5d0,-8.5d0,564.d0,tev)
      cosu(12)=cosu(12)+arf(18.3d0,-4.7d0,2.8d0,-12.8d0,589.d0,tev)
      cosu(13)=arf(18.1d0,-4.4d0,2.7d0,-11.8d0,652.d0,tev)
      cosu(13)=cosu(13)+arf(22.8d0,-7.1d0,4.1d0,-18.d0,3017.d0,tev)
      cosu(14)=arf(9.7d0,-2.8d0,1.4d0,-6.4d0,707.d0,tev)
      cosu(14)=cosu(14)+arf(22.8d0,-7.1d0,4.1d0,-18.d0,3075.d0,tev)
      cosu(15)=arf(23.3d0,-7.6d0,4.d0,-19.3d0,3224.d0,tev)
      cosu(16)=arf(12.8d0,-4.3d0,1.9d0,-10.5d0,3493.d0,tev)
C -Argon-
      coar(1)=arf(171.d0,-78.d0,3.8d0,-169.d0,15.8d0,tev)
      coar(1)=coar(1)+arf(48.7d0,-30.5d0,1.4d0,-39.7d0,29.2d0,tev)
      coar(2)=arf(147.d0,-97.4d0,3.2d0,-107.7d0,27.6d0,tev)
      coar(2)=coar(2)+arf(50.d0,-31.8d0,1.4d0,-40.d0,41.7d0,tev)
      coar(3)=arf(122.8d0,-81.4d0,2.6d0,-90.d0,40.9d0,tev)
      coar(3)=coar(3)+arf(51.3d0,-33.2d0,1.4d0,-40.2d0,55.5d0,tev)
      coar(4)=arf(98.7d0,-65.4d0,1.9d0,-72.3d0,59.7d0,tev)
      coar(4)=coar(4)+arf(52.5d0,-34.5d0,1.4d0,-40.5d0,70.4d0,tev)
      coar(5)=arf(74.5d0,-49.4d0,1.3d0,-54.6d0,75.2d0,tev)
      coar(5)=coar(5)+arf(53.8d0,-35.8d0,1.4d0,-40.7d0,87.6d0,tev)
      coar(6)=arf(50.4d0,-33.4d0,0.6d0,-36.9d0,91.2d0,tev)
      coar(6)=coar(6)+arf(55.1d0,-37.2d0,1.4d0,-41.d0,105.d0,tev)
      coar(7)=arf(19.8d0,-5.7d0,1.9d0,-11.9d0,125.d0,tev)
      coar(7)=coar(7)+arf(74.8d0,-27.d0,14.1d0,-58.6d0,373.d0,tev)
      coar(7)=coar(7)+arf(23.4d0,-8.3d0,3.3d0,-18.5d0,427.d0,tev)
      coar(8)=arf(9.d0,-2.7d0,0.8d0,-5.4d0,143.d0,tev)
      coar(8)=coar(8)+arf(74.8d0,-27.d0,14.1d0,-58.6d0,394.d0,tev)
      coar(8)=coar(8)+arf(23.4d0,-8.3d0,3.3d0,-18.5d0,453.1d0,tev)
      coar(9)=arf(72.d0,-24.1d0,11.9d0,-50.d0,423.d0,tev)
      coar(9)=coar(9)+arf(19.6d0,-7.3d0,2.8d0,-16.d0,498.d0,tev)
      coar(10)=arf(60.8d0,-20.2d0,9.3d0,-41.7d0,479.d0,tev)
      coar(10)=coar(10)+arf(19.4d0,-6.8d0,2.8d0,-15.3d0,545.d0,tev)
      coar(11)=arf(49.5d0,-16.3d0,6.9d0,-33.4d0,539.d0,tev)
      coar(11)=coar(11)+arf(19.2d0,-6.4d0,2.8d0,-14.6d0,594.d0,tev)
      coar(12)=arf(38.3d0,-12.4d0,4.8d0,-25.1d0,618.d0,tev)
      coar(12)=coar(12)+arf(18.9d0,-5.9d0,2.8d0,-14.d0,644.d0,tev)
      coar(13)=arf(27.d0,-8.5d0,3.d0,-16.8d0,686.d0,tev)
      coar(13)=coar(13)+arf(18.7d0,-5.4d0,2.8d0,-13.3d0,713.d0,tev)
      coar(14)=arf(14.d0,-4.6d0,1.4d0,-8.5d0,755.d0,tev)
      coar(14)=coar(14)+arf(18.5d0,-5.d0,2.8d0,-12.6d0,784.d0,tev)
      coar(15)=arf(18.4d0,-4.6d0,2.7d0,-12.d0,855.d0,tev)
      coar(15)=coar(15)+arf(23.1d0,-7.3d0,4.1d0,-18.d0,3855.d0,tev)
      coar(16)=arf(9.6d0,-2.8d0,1.4d0,-6.2d0,918.d0,tev)
      coar(16)=coar(16)+arf(23.1d0,-7.3d0,4.1d0,-18.d0,3951.d0,tev)
      coar(17)=arf(23.7d0,-7.8d0,3.9d0,-19.5d0,4121.d0,tev)
      coar(18)=arf(12.8d0,-4.4d0,1.9d0,-10.5d0,4426.d0,tev)
C -Calcium-
      coca(1)=arf(2.5d0,-2.5d0,8.d0,-5.5d0,6.1d0,tev)
      coca(1)=coca(1)+arf(74.3d0,-24.2d0,7.d0,-63.d0,28.d0,tev)
      coca(1)=coca(1)+arf(17.6d0,-3.8d0,1.9d0,-13.8d0,40.3d0,tev)
      coca(2)=arf(7.9d0,-2.d0,0.2d0,-6.d0,11.9d0,tev)
      coca(2)=coca(2)+arf(74.3d0,-24.2d0,7.d0,-68.d0,37.d0,tev)
      coca(2)=coca(2)+arf(17.6d0,-3.8d0,1.9d0,-13.8d0,45.2d0,tev)
      coca(3)=arf(74.3d0,-24.3d0,7.d0,-63.d0,51.2d0,tev)
      coca(3)=coca(3)+arf(17.6d0,-3.8d0,1.9d0,-13.8d0,70.1d0,tev)
      coca(4)=arf(55.8d0,-15.8d0,6.4d0,-44.5d0,67.3d0,tev)
      coca(4)=coca(4)+arf(16.2d0,-3.2d0,1.8d0,-11.6d0,86.4d0,tev)
      coca(5)=arf(47.1d0,-14.5d0,4.8d0,-35.5d0,84.5d0,tev)
      coca(5)=coca(5)+arf(18.9d0,-5.1d0,1.6d0,-13.2d0,104.d0,tev)
      coca(6)=arf(40.9d0,-13.6d0,3.4d0,-30.1d0,109.d0,tev)
      coca(6)=coca(6)+arf(20.4d0,-6.3d0,2.1d0,-13.8d0,123.d0,tev)
      coca(7)=arf(22.9d0,-7.4d0,2.8d0,-15.9d0,128.d0,tev)
      coca(7)=coca(7)+arf(21.9d0,-7.7d0,1.9d0,-14.9d0,144.d0,tev)
      coca(8)=arf(11.1d0,-3.4d0,1.3d0,-7.3d0,148.d0,tev)
      coca(8)=coca(8)+arf(22.7d0,-8.6d0,1.9d0,-15.5d0,165.d0,tev)
      coca(9)=arf(19.8d0,-5.7d0,1.8d0,-11.9d0,189.d0,tev)
      coca(9)=coca(9)+arf(76.1d0,-27.d0,12.8d0,-56.6d0,534.d0,tev)
      coca(9)=coca(9)+arf(23.5d0,-8.4d0,3.3d0,-17.8d0,593.1d0,tev)
      coca(10)=arf(9.d0,-2.6d0,0.9d0,-5.4d0,211.d0,tev)
      coca(10)=coca(10)+arf(76.1d0,-27.d0,12.8d0,-56.6d0,559.d0,tev)
      coca(10)=coca(10)+arf(23.5d0,-8.4d0,3.3d0,-17.8d0,623.7d0,tev)
      coca(11)=arf(72.d0,-24.1d0,10.3d0,-50.d0,592.d0,tev)
      coca(11)=coca(11)+arf(19.5d0,-7.8d0,2.8d0,-14.5d0,680.d0,tev)
      coca(12)=arf(60.8d0,-20.2d0,8.1d0,-41.7d0,657.d0,tev)
      coca(12)=coca(12)+arf(19.4d0,-7.3d0,2.8d0,-14.1d0,734.d0,tev)
      coca(13)=arf(49.5d0,-16.3d0,6.1d0,-33.4d0,727.d0,tev)
      coca(13)=coca(13)+arf(19.2d0,-6.8d0,2.8d0,-13.7d0,790.d0,tev)
      coca(14)=arf(38.3d0,-12.4d0,4.3d0,-25.1d0,818.d0,tev)
      coca(14)=coca(14)+arf(19.d0,-6.3d0,2.8d0,-13.2d0,847.d0,tev)
      coca(15)=arf(27.d0,-8.5d0,2.7d0,-16.8d0,894.d0,tev)
      coca(15)=coca(15)+arf(18.9d0,-5.8d0,2.8d0,-12.8d0,925.d0,tev)
      coca(16)=arf(14.d0,-4.6d0,1.3d0,-8.5d0,974.d0,tev)
      coca(16)=coca(16)+arf(18.7d0,-5.3d0,2.8d0,-12.4d0,1006.d0,tev)
      coca(17)=arf(18.6d0,-4.6d0,2.7d0,-12.1d0,1087.d0,tev)
      coca(17)=coca(17)+arf(23.4d0,-7.4d0,4.1d0,-18.d0,4865.d0,tev)
      coca(18)=arf(9.5d0,-2.7d0,1.4d0,-6.1d0,1157.d0,tev)
      coca(18)=coca(18)+arf(23.4d0,-7.4d0,4.1d0,-18.d0,4939.d0,tev)
      coca(19)=arf(24.0d0,-7.9d0,3.9d0,-19.6d0,5129.d0,tev)
      coca(20)=arf(12.9d0,-4.4d0,1.9d0,-10.5d0,5470.d0,tev)

            call coion_fit(1,1,tq,cov)

C -Iron-
      call collionfeea(te,cofeea)
      call collionfe(te,cofedir)
      do k=1,26
         cofe(k)=cofedir(k)+cofeea(k)
      enddo
      call coion_fit(1,1,tq,cov)
C -Nickel-
      coni(1)=arf(2.5d0,-0.8d0,0.2d0,-1.2d0,8.7d0,tev)
      coni(1)=coni(1)+arf(12.6d0,-4.d0,0.4d0,-6.1d0,10.d0,tev)
      coni(1)=coni(1)+arf(69.9d0,-23.7d0,9.5d0,-51.7d0,73.d0,tev)
      coni(2)=arf(32.d0,-10.d0,1.d0,-15.4d0,18.2d0,tev)
      coni(2)=coni(2)+arf(69.9d0,-23.7d0,9.5d0,-51.7d0,97.d0,tev)
      coni(2)=coni(2)+arf(19.2d0,-5.7d0,2.3d0,-12.7d0,142.d0,tev)
      coni(3)=arf(44.4d0,-14.1d0,1.4d0,-21.5d0,35.2d0,tev)
      coni(3)=coni(3)+arf(69.9d0,-23.7d0,9.5d0,-51.7d0,122.d0,tev)
      coni(3)=coni(3)+arf(19.2d0,-5.7d0,2.3d0,-12.7d0,166.d0,tev)
      coni(4)=arf(50.3d0,-16.d0,1.6d0,-24.3d0,54.9d0,tev)
      coni(4)=coni(4)+arf(69.9d0,-23.7d0,9.5d0,-51.7d0,146.d0,tev)
      coni(4)=coni(4)+arf(19.2d0,-5.7d0,2.3d0,-12.7d0,190.d0,tev)
      coni(5)=arf(49.9d0,-15.9d0,1.6d0,-24.1d0,75.5d0,tev)
      coni(5)=coni(5)+arf(69.9d0,-23.7d0,9.5d0,-51.7d0,171.d0,tev)
      coni(5)=coni(5)+arf(19.2d0,-5.7d0,2.3d0,-12.7d0,215.d0,tev)
      coni(6)=arf(50.8d0,-16.1d0,1.6d0,-24.6d0,108.d0,tev)
      coni(6)=coni(6)+arf(69.9d0,-23.7d0,9.5d0,-51.7d0,196.d0,tev)
      coni(6)=coni(6)+arf(19.2d0,-5.7d0,2.3d0,-12.7d0,239.d0,tev)
      coni(7)=arf(43.2d0,-13.7d0,1.3d0,-20.9d0,133.d0,tev)
      coni(7)=coni(7)+arf(69.9d0,-23.7d0,9.5d0,-51.7d0,221.d0,tev)
      coni(7)=coni(7)+arf(19.2d0,-5.7d0,2.3d0,-12.7d0,264.d0,tev)
      coni(8)=arf(34.5d0,-10.9d0,1.1d0,-16.7d0,162.d0,tev)
      coni(8)=coni(8)+arf(69.9d0,-23.7d0,9.5d0,-51.7d0,246.d0,tev)
      coni(8)=coni(8)+arf(19.2d0,-5.7d0,2.3d0,-12.7d0,288.d0,tev)
      coni(9)=arf(24.1d0,-7.7d0,0.7d0,-11.7d0,193.d0,tev)
      coni(9)=coni(9)+arf(69.9d0,-23.7d0,9.5d0,-51.7d0,271.d0,tev)
      coni(9)=coni(9)+arf(19.2d0,-5.7d0,2.3d0,-12.7d0,331.d0,tev)
      coni(10)=arf(12.5d0,-4.d0,0.4d0,-6.d0,225.d0,tev)
      coni(10)=coni(10)+arf(69.9d0,-23.7d0,9.5d0,-51.7d0,296.d0,tev)
      coni(10)=coni(10)+arf(19.2d0,-5.7d0,2.3d0,-12.7d0,338.d0,tev)
      coni(11)=arf(69.9d0,-23.7d0,9.5d0,-51.7d0,321.d0,tev)
      coni(11)=coni(11)+arf(19.2d0,-5.7d0,2.3d0,-12.7d0,363.d0,tev)
      coni(12)=arf(57.7d0,-18.6d0,7.8d0,-40.3d0,352.d0,tev)
      coni(12)=coni(12)+arf(21.d0,-7.1d0,2.3d0,-14.1d0,393.d0,tev)
      coni(13)=arf(45.6d0,-13.9d0,6.2d0,-30.d0,384.d0,tev)
      coni(13)=coni(13)+arf(22.8d0,-8.4d0,2.3d0,-15.4d0,423.d0,tev)
      coni(14)=arf(33.4d0,-9.7d0,4.6d0,-20.8d0,430.d0,tev)
      coni(14)=coni(14)+arf(24.6d0,-9.8d0,2.3d0,-16.8d0,458.d0,tev)
      coni(15)=arf(21.3d0,-5.9d0,3.d0,-12.6d0,464.d0,tev)
      coni(15)=coni(15)+arf(26.4d0,-11.2d0,2.3d0,-18.1d0,494.d0,tev)
      coni(16)=arf(9.1d0,-2.6d0,1.4d0,-5.6d0,499.d0,tev)
      coni(16)=coni(16)+arf(28.2d0,-12.5d0,2.3d0,-19.4d0,531.d0,tev)
      coni(17)=arf(19.8d0,-5.7d0,2.2d0,-11.9d0,571.d0,tev)
      coni(17)=coni(17)+arf(79.6d0,-27.d0,10.1d0,-51.9d0,1458.d0,tev)
      coni(17)=coni(17)+arf(23.5d0,-7.8d0,3.3d0,-16.2d0,1529.1d0,tev)
      coni(18)=arf(9.d0,-2.6d0,1.d0,-5.4d0,608.d0,tev)
      coni(18)=coni(18)+arf(79.6d0,-27.d0,10.1d0,-51.9d0,1500.d0,tev)
      coni(18)=coni(18)+arf(23.5d0,-7.8d0,3.3d0,-16.2d0,1578.d0,tev)
      coni(19)=arf(82.2d0,-26.4d0,5.9d0,-49.8d0,1546.d0,tev)
      coni(19)=coni(19)+arf(19.4d0,-5.7d0,2.8d0,-12.3d0,1694.d0,tev)
      coni(20)=arf(68.5d0,-22.d0,4.9d0,-41.5d0,1648.d0,tev)
      coni(20)=coni(20)+arf(19.4d0,-5.7d0,2.8d0,-12.3d0,1775.d0,tev)
      coni(21)=arf(54.8d0,-17.6d0,3.9d0,-33.2d0,1756.d0,tev)
      coni(21)=coni(21)+arf(19.4d0,-5.7d0,2.8d0,-12.3d0,1858.d0,tev)
      coni(22)=arf(41.1d0,-13.2d0,2.9d0,-24.9d0,1894.d0,tev)
      coni(22)=coni(22)+arf(19.4d0,-5.7d0,2.8d0,-12.3d0,1938.d0,tev)
      coni(23)=arf(27.4d0,-8.8d0,2.d0,-16.6d0,2011.d0,tev)
      coni(23)=coni(23)+arf(19.4d0,-5.7d0,2.8d0,-12.3d0,2056.d0,tev)
      coni(24)=arf(13.7d0,-4.4d0,1.d0,-8.3d0,2131.d0,tev)
      coni(24)=coni(24)+arf(19.4d0,-5.7d0,2.8d0,-12.3d0,2178.d0,tev)
      coni(25)=arf(19.3d0,-5.4d0,2.7d0,-12.3d0,2295.d0,tev)
      coni(25)=coni(25)+arf(24.3d0,-8.d0,4.1d0,-18.d0,9914.d0,tev)
      coni(26)=arf(9.2d0,-2.6d0,1.4d0,-5.7d0,2399.d0,tev)
      coni(26)=coni(26)+arf(24.3d0,-8.d0,4.1d0,-18.d0,10020.d0,tev)
      coni(27)=arf(25.d0,-8.5d0,3.7d0,-20.1d0,10280.d0,tev)
      coni(28)=arf(13.d0,-4.5d0,1.9d0,-10.6d0,10790.d0,tev)
c     direct collisional ionization from Voronov 1997, ADNDT, 65, 1 (via verner)
      call coion_fit(1,1,tq,cov)
      do iel=1,14
         do ion=1,nion(iel)            
            in=nion(iel)-ion+1
            call coion_fit(nion(iel),in,te,cov)
c            call coion_fit(1,1,te,cov)
c            stop
            if(nion(iel)==1) then
               carn=ch
               coh(ion)=cov
            elseif(nion(iel)==2) then
               if(ion==1) then
                  carn=che0
               elseif(ion==2) then
                  carn=che
               endif
               cohe(ion)=cov
            elseif(nion(iel)==6) then
               coc(ion)=cov
            elseif(nion(iel)==7) then
               con(ion)=cov
            elseif(nion(iel)==8) then
               carn=coo(ion)
               coo(ion)=cov
            elseif(nion(iel)==10) then
               cone(ion)=cov
            elseif(nion(iel)==11) then
               cona(ion)=cov
            elseif(nion(iel)==12) then
               comg(ion)=cov
            elseif(nion(iel)==13) then
               coal(ion)=cov
            elseif(nion(iel)==14) then
               cosi(ion)=cov
            elseif(nion(iel)==16) then
               cosu(ion)=cov               
            elseif(nion(iel)==18) then
               coar(ion)=cov
            elseif(nion(iel)==20) then
               coca(ion)=cov
            elseif(nion(iel)==26) then
               carn=cofe(ion)
               cofe(ion)=cov
            endif

         enddo
      enddo

C
C - Autoionization: Lithium sequence -
      coc(4)=coc(4)+arauli(6,tev)
      con(5)=con(5)+arauli(7,tev)
      coo(6)=coo(6)+arauli(8,tev)
      cone(8)=cone(8)+arauli(10,tev)
      cona(9)=cona(9)+arauli(11,tev)
      comg(10)=comg(10)+arauli(12,tev)
      coal(11)=coal(11)+arauli(13,tev)
      cosi(12)=cosi(12)+arauli(14,tev)
      cosu(14)=cosu(14)+arauli(16,tev)
      coar(16)=coar(16)+arauli(18,tev)
      coca(18)=coca(18)+arauli(20,tev)
      coni(26)=coni(26)+arauli(28,tev)
C - Autoionization: Sodium sequence -
      cona(1)=cona(1)+arau(11,11,tev)
      comg(2)=comg(2)+arau(12,11,tev)
      coal(3)=coal(3)+arau(13,11,tev)
      cosi(4)=cosi(4)+arau(14,11,tev)
      cosu(6)=cosu(6)+arau(16,11,tev)
      coar(8)=coar(8)+arau(18,11,tev)
      coca(10)=coca(10)+arau(20,11,tev)
      coni(18)=coni(18)+arau(28,11,tev)
C - Autoionization: Magnesium sequence -
      comg(1)=comg(1)+arau(12,12,tev)
      coal(2)=coal(2)+arau(13,12,tev)
      cosi(3)=cosi(3)+arau(14,12,tev)
      cosu(5)=cosu(5)+arau(16,12,tev)
      coar(7)=coar(7)+arau(18,12,tev)
      coca(9)=coca(9)+arau(20,12,tev)
      coni(17)=coni(17)+arau(28,12,tev)
C - Autoionization: Aluminium sequence -
      coal(1)=coal(1)+arau(13,13,tev)
      cosi(2)=cosi(2)+arau(14,13,tev)
      cosu(4)=cosu(4)+arau(16,13,tev)
      coar(6)=coar(6)+arau(18,13,tev)
      coca(8)=coca(8)+arau(20,13,tev)
      coni(16)=coni(16)+arau(28,13,tev)
C - Autoionization: Silicon sequence -
      cosi(1)=cosi(1)+arau(14,14,tev)
      cosu(3)=cosu(3)+arau(16,14,tev)
      coar(5)=coar(5)+arau(18,14,tev)
      coca(7)=coca(7)+arau(20,14,tev)
      coni(15)=coni(15)+arau(28,14,tev)
C - Autoionization: Phosphorus sequence -
      cosu(2)=cosu(2)+arau(16,15,tev)
      coar(4)=coar(4)+arau(18,15,tev)
      coca(6)=coca(6)+arau(20,15,tev)
      coni(14)=coni(14)+arau(28,15,tev)
C - Autoionization: Sulphur sequence -
      cosu(1)=cosu(1)+arau(16,16,tev)
      coar(3)=coar(3)+arau(18,16,tev)
      coca(5)=coca(5)+arau(20,16,tev)
      coni(13)=coni(13)+arau(28,16,tev)
C - Autoionization: Potassium sequence -
      coca(2)=coca(2)+arau(20,19,tev)
C - Autoionization: Calcium sequence -
      coca(1)=coca(1)+arau(20,20,tev)

      do i=1,2
         collion(2,i)=cohe(i)
      enddo

      do i=1,6
         collion(6,i)=coc(i)
      enddo

      do i=1,7
         collion(7,i)=con(i)
      enddo


      do i=1,8
         collion(8,i)=coo(i)
      enddo

      do i=1,10
         collion(10,i)=cone(i)
      enddo

      do i=1,11
         collion(11,i)=cona(i)
      enddo
      do i=1,12
         collion(12,i)=comg(i)
      enddo

      do i=1,13
         collion(13,i)=coal(i)
      enddo

      do i=1,14
         collion(14,i)=cosi(i)
      enddo

      do i=1,16
         collion(16,i)=cosu(i)
      enddo

      do i=1,18
         collion(18,i)=coar(i)
      enddo

      do i=1,20
         collion(20,i)=coca(i)
      enddo

      do i=1,26
         collion(26,i)=cofe(i)
      enddo

C *********************************************************************   
C   RECOMBINATION RATES FOR METALS.
C                                                                       
C   For the dielectronic recombination rates, use Shull & Steenberg's   
C   fits (Ap.J. Suppl. 48, 95, 1982) to the results of Jacobs et al..
C   (Correct rates listed in Landini & Monsignori Fossi.) 
C   Correction factors according to Bely-Dubau et al. (1984) are;
C   Carbon = 1., Nitrogen = 1., Oxygen = 0.81, Neon = 0.7, 
C   Sodium = 1(?), Magnesium = 0.56, Aluminium = 1(?), Silicon = 0.6,
C   Sulphur = 0.6, Argon = 0.6, Calcium = 0.64, Iron = 0.59, 
C   Nickel = 0.6.
C   Low temperature rates from Nussbaumer and Storey's calcs.                  
C   (Check update by Romanik (Ap.J 330,1022, 1988).)
C                                                                       
C   Direct rad.recomb. for C IV-V, N V-VI, O VI-VIII and Ne VIII-IX from
C   the compilation by Arnaud & Rothenflug, A&A Suppl. 60, 425, 1985.   
C   (Also check Landini and Monsignori Fossi, 1990?, Astr.Ap.Suppl.)
C     *********************************************************************
c call Badnell et al rates for the Na, Mg, Al, Si, P isoelectronic ions
      call rr_diel_badnell(te)      

c skip the rec rates for C,N,O below and use rates from Nahar below instead

C     ***************************************************************   
C     *****                                                             
C     OXYGEN 
C     *****                                                             
C     ***************************************************************   



      BDCORR=0.81D0                                                     
      SUPR=1.                                                           
      IF(DENS.GT.1.D9) SUPR=0.26                                     
c
      ALO(9)=AL(6)                                                      
c
      CALL DIEL(6.23D-2,7.01D6,0.304D0,1.47D6,TE,ALD)                  
      ALD=ALD*SUPR*BDCORR                                              
      ALO(8)=7.19D-11/T4**.834+ALD                                      
c
cf!!
      bdcorr=1.
      CALL DIEL(1.06D-1,6.25D6,0.34D0,1.12D6,TE,ALD)                    
      ALD=ALD*SUPR*BDCORR                                              
      ALO(7)=2.44D-11/T4**.774+ALD                                      
c
      CALL DIEL(4.13D-3,1.25D5,0.162D0,2.27D5,TE,ALD)                   
      CALL DIELB(-2.8425D0,.2283D0,4.0407D1,-3.4956D0,1.7558D0,T4,ALDB) 
c     if(telog.lt.5.d0) then
c      badnel=0.8838d0
c     else if(telog.ge.5.d0.and.telog.lt.7.d0) then
c      badnel=-69.517d0+31.494d0*telog-4.5871d0*telog**2+.22337d0*
c    & telog**3
c     else if(telog.ge.7.d0) then
c      badnel=2.476d0
c     endif
c     ALD=badnel*ALD
      ALD=(ALD+ALDB)*SUPR                                              
      ALO(6)=1.59D-11/T4**.759+ALD                                      
C                                                                       
      CALL DIEL(1.84D-2,2.12D5,0.1D0,2.83D5,TE,ALD)                     
      CALL DIELB(6.1D-3,0.2269D0,3.2149D1,1.9939D0,-6.46D-2,T4,ALDB)
      ALD=(ALD*BDCORR+ALDB)*SUPR   
C
C     C-T BD 80
C
C*********

      CT=((.16+.096*T4**1.14)*1.E-9*xhi+6.5E-10*xhei)/xel
                                            
      ALO(5)=9.6D-12/T4**.670+ALD                             
c
      CALL DIEL(1.48D-2,2.41D5,0.305d0,2.83D5,TE,ALD)                     
      CALL DIELB(0.0D0,21.88D0,1.6273D1,-.7020D0,1.1899D0,T4,ALDB)      
      ALD=(ALD*BDCORR+ALDB)*SUPR    

C
C     C-T (BHD 80)
C
      CT=(.45+8.18*T4**.47)*1.E-9*xhi/xel
      CTHE=1.E-9*xhei/xel
                                   
      ALO(4)=5.1D-12/T4**.666+ALD
c
      CALL DIEL(5.07D-3,1.98D5,0.181D0,3.35D5,TE,ALD)                   
      CALL DIELB(-3.6D-3,.7519D0,1.5252D0,-8.38D-2,.2769D0,T4,ALDB)     
      ALD=(ALD*BDCORR+ALDB)*SUPR  
C
C     C-T (BHD 80)
C

      CT=(.28+.49*T4**.61)*1.E-9*xhi/xel
      CTHE=T4*2.E-10*xhei/xel
      ALO(3)=2.0D-12/T4**.646+ALD

      alchoii=alo(3)
c
      CALL DIEL(1.11D-3,1.75D5,0.0925D0,1.45D5,TE,ALD)                  
      CALL DIELB(0.0D0,2.38D-2,6.59D-2,3.49D-2,.5334D0,T4,ALDB)         
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       

C
C     C-T FS 71
C

      cth=(.66+.25*T4**.5)*1.E-9*xhi/xel

 
      ALO(2)=3.1D-13/T4**.678 + ALD
      alcho=alo(2)

c
      ALO(1)=0.0d0


C     ***************************************************************   
C     *****                                                             
C     CARBON
C     *****                                                             
C     ***************************************************************   
      IF(DENS.GT.1.D9) SUPR=0.17                                     
c
      ALC(1)=0.0D0                                                      
c
      CALL DIEL(2.54D-3,1.57D5,0.0442D0,3.74D5,TE,ALD)                  
      CALL DIELB(1.08D-2,-.1075D0,.281D0,-1.93D-2,-.1127D0,T4,ALDB)     
      ALD=(ALD+ALDB)*SUPR                                              
      ALC(2)=4.7D-13/T4**.624+ALD                                       
      alchc=alc(2)


c
      CALL DIEL(6.15D-3,1.41D5,0.0588D0,1.41D5,TE,ALD)                  
      CALL DIELB(1.8267D0,4.1012D0,4.8443D0,.2261D0,.596D0,T4,ALDB)     
      ALD=(ALD+ALDB)*SUPR   

C
C     C-T BHD 80
C
      CT=1.18E-12*xhi/xel
                                           
      ALC(3)=2.3D-12/T4**.645+ALD + ct
      alchcii=alc(3)
c
      CALL DIEL(1.62D-3,8.19D4,0.343D0,1.59D5,TE,ALD)             
c     if(telog.lt.4.8d0) then
c      badnel=0.852d0
c     else if(telog.ge.4.8d0.and.telog.lt.6.8d0) then
c      badnel=-76.482d0+36.941d0*telog-5.7961d0*telog**2+.30345d0*
c    & telog**3
c     else if(telog.ge.6.8d0) then
c      badnel=2.119d0
c     endif
c     ALD=badnel*ALD
      CALL DIELB(2.3196D0,1.0733D1,6.883D0,-.1824D0,.4101D0,T4,ALDB)    
      ALD=(ALD+ALDB)*SUPR    

C
C     C-T BHD 80
C
      CT=(1.49+2.09*T4**.39)*1.E-9*XHi/xel
C
                                          
      ALC(4)=4.90D-12/T4**.803+ALD + ct
c
      CALL DIEL(4.78D-2,3.44D6,0.362D0,5.87D5,TE,ALD)                   
      ALD=ALD*SUPR                

C
C     C-T (B D 80)
C
      CT=(-.05+.81*T4**1.353)*1.E-9*XHi/xEL

      ct=max(0.d0,ct)
                                     
      ALC(5)=9.16D-12/T4**.791+ALD + ct
c
      CALL DIEL(3.22D-2,4.06D6,0.315D0,8.31D5,TE,ALD)                   
      ALD=ALD*SUPR                                                     
      ALC(6)=1.7D-11/T4**.721+ALD                                       
c
      ALC(7)=AL(4)                                                      


C     ***************************************************************   
C     *****                                                             
C     NITROGEN 
C     *****                                                             
C     ***************************************************************   
      ALN(1)=0.                                                         
c
      IF(DENS.GT.1.D9) SUPR=0.20                                     
      CALL DIEL(2.98D-3,2.20D5,0.0D0,0.0D0,TE,ALD)                      
      CALL DIELB(0.D0,.631D0,.199D0,-1.97D-2,.4398D0,T4,ALDB)           
      ALD=(ALD+ALDB)*SUPR   

C
C     C-T BD 79
C
      CT=1.0E-12*xhi/xel
                                           
      ALN(2)=4.1D-13/T4**.608+ALD + ct

      alchn = aln(2)
c
      CALL DIEL(7.41D-3,2.01D5,0.0764D0,7.37D4,TE,ALD)                  
      CALL DIELB(3.2D-2,-.6624D0,4.3191D0,3.D-4,.5946D0,T4,ALDB)        
      ALD=(ALD+ALDB)*SUPR    

C
C     C-T (BHD 1980 )
C
      CT=(.57+.29*T4**.46)*1.E-9*xhi/xel
                                          
      ALN(3)=2.2D-12/T4**.639+ALD+ct 
      alchnii = aln(3)
c
      CALL DIEL(1.13D-2,1.72D5,0.164D0,2.25D5,TE,ALD)                   
      CALL DIELB(-.8806D0,1.124D1,3.071D1,-1.1721D0,.6127D0,T4,ALDB)    
      ALD=(ALD+ALDB)*SUPR 

C
C     C-T BHD 80
C
      CT=(-.82+3.75*T4**.67)*1.E-9*XHi/xel

      ct = max(ct,0.)
                                             
      ALN(4)=5.0D-12/T4**.676+ALD + ct               
c
      CALL DIEL(2.62D-3,1.02D5,0.243D0,1.25D5,TE,ALD)         
c     if(telog.lt.5.d0) then
c      badnel=1.261d0
c     else if(telog.ge.5.d0.and.telog.lt.7.d0) then
c      badnel=-68.517d0+31.725d0*telog-4.7333d0*telog**2+.23588d0*
c    & telog**3
c     else if(telog.ge.7.d0) then
c      badnel=2.533d0
c     endif
c     ALD=badnel*ALD
      CALL DIELB(.4134D0,-4.6319D0,2.5917D1,-2.229D0,.236D0,T4,ALDB)    
      ALD=(ALD+ALDB)*SUPR      

C
C     CT (BUTLER AND DALG. 1980 )
C
      CT=(-.0036+.164*T4**1.455)*1.E-9*xhi/xel
      
      ct = max(ct,0.)
                                        
      ALN(5)=9.40D-12/T4**.765+ALD + ct                
c
      CALL DIEL(7.50D-2,4.75D6,0.35D0,8.35D5,TE,ALD)                    
      ALD=ALD*SUPR                                                     
      ALN(6)=1.57D-11/T4**.791+ALD                                      
c
      CALL DIEL(4.61D-2,5.44D6,0.309D0,1.14D6,TE,ALD)                   
      ALD=ALD*SUPR                                                     
      ALN(7)=2.9D-11/T4**.750+ALD                                       
c
      ALN(8)=AL(5)                                                      


      if(initrec.eq.0) then

c read recomb data from Nahar for C, N, O
         
c C
         open(71,file='./ATDAT/trrc.tbl.cions.dat',status='old')
         do i=1,33
            read(71,933)chr
 933        format(a)
         enddo

         do i=1,81
            read(71,*)tlr(i),(recc(k,i),k=1,6)
            do k=1,6
               recc(k,i)=log10(recc(k,i))
            enddo
         enddo
         close(71)

c N
         open(71,file='./ATDAT/trrc.tbl.nions.dat',status='old')
         do i=1,4
            read(71,933)chr
         enddo

         do i=1,81
            read(71,*)tlr(i),(recn(k,i),k=1,7)
            do k=1,7
               recn(k,i)=log10(recn(k,i))
            enddo
         enddo
         close(71)

c O
         open(71,file='./ATDAT/trrc.tbl.oions.dat',status='old')
         do i=1,24
            read(71,933)chr
         enddo

         do i=1,81
            read(71,*)tlr(i),(reco(k,i),k=1,8)
            do k=1,8
               reco(k,i)=log10(reco(k,i))
 2          enddo
         enddo
         close(71)
         
c rec from Nahar for Si II-> I, III--> II, S III-->II, IV --> III. Only Si III--> II used. Others from BAdnell + below
         open(71,file='./ATDAT/rec_si_2_3_s_3_4_nahar.dat',status='old')
         do i=1,81
            read(71,*)temp,(recsi(k,i),k=1,4)
            tlr(i)=log10(temp)
            do k=1,2
               recsi(k,i)=log10(recsi(k,i))
            enddo
         enddo
         close(71)
         
         initrec=1

      endif
      do i=1,81
         telog=log10(te)
         if(telog.ge.tlr(81)) then
            irec=80
            fac= (telog-tlr(irec))/(tlr(irec+1)-tlr(irec))
            goto 33
         elseif(telog.ge.tlr(i).and.telog.lt.tlr(i+1)) then
            irec=i
            fac= (telog-tlr(i))/(tlr(i+1)-tlr(i))
           goto 33 
         endif
      enddo
 33   continue
      do k=1,6
         alcl=recc(k,irec)+(recc(k,irec+1)-recc(k,irec))*fac
         alc(k+1)=10.**alcl
      enddo

      do k=1,7
         alnl=recn(k,irec)+(recn(k,irec+1)-recn(k,irec))*fac
         aln(k+1)=10.**alnl
      enddo

      do k=1,8
         alol=reco(k,irec)+(reco(k,irec+1)-reco(k,irec))*fac
         alo(k+1)=10.**alol
      enddo

      do k=1,2
         alsil=recsi(k,irec)+(recsi(k,irec+1)-recsi(k,irec))*fac
         alsi_nahar(k+1)=10.**alsil
      enddo

c
      ALC(1)=0.0D0                                                      

C
C     C-T BHD 80
C
      CT=1.18E-12*xhi/xel
                                           
c      ALC(3)=alc(3) + ct

C
C     C-T BHD 80
C
      CT=(1.49+2.09*T4**.39)*1.E-9*XHi/xel
                                          
c      ALC(4)=alc(4) + ct

C
C     C-T (B D 80)
C
      CT=(-.05+.81*T4**1.353)*1.E-9*XHi/xEL

      ct=max(0.d0,ct)
                                     
c      ALC(5)=alc(5) + ct

C
C     C-T BD 79
C
      CT=1.0E-12*xhi/xel
                                           
c      ALN(2)=aln(2) + ct

C
C     C-T (BHD 1980 )
C
      CT=(.57+.29*T4**.46)*1.E-9*xhi/xel
                                          
c      ALN(3)=aln(3)+ct 

C
C     C-T BHD 80
C
      CT=(-.82+3.75*T4**.67)*1.E-9*XHi/xel

      ct = max(ct,0.)
                                             
c      ALN(4)=aln(4) + ct               

C
C     CT (BUTLER AND DALG. 1980 )
C
      CT=(-.0036+.164*T4**1.455)*1.E-9*xhi/xel
      
      ct = max(ct,0.)
c      ALN(5)=aln(5) + ct                



C
C     C-T BD 80
C


      CT=((.16+.096*T4**1.14)*1.E-9*xhi+6.5E-10*xhei)/xel
                                            
c      ALO(5)=alo(5)+CT                             

C
C     C-T (BHD 80)
C
      CT=(.45+8.18*T4**.47)*1.E-9*xhi/xel
      CTHE=1.E-9*xhei/xel
                                   
c      ALO(4)=alo(4)+CT+CTHE                             

C
C     C-T (BHD 80)
C

      CT=(.28+.49*T4**.61)*1.E-9*xhi/xel
      CTHE=T4*2.E-10*xhei/xel
c      ALO(3)=alo(3)+CT+CTHE

C
C     C-T FS 71
C

      cth=(.66+.25*T4**.5)*1.E-9*xhi/xel
 
c      ALO(2)=alo(2) + cth                             

c
      ALO(1)=0.0d0


C     ***************************************************************   
C     *****                                                             
C     NEON      
C     *****                                                             
C     ***************************************************************   
      ALNE(1)=0.                                                       
c

cf!!
      bdcorr=1.
      SUPR=1.D0
      IF(DENS.GT.1.D9) SUPR=0.20                                     
C
      CALL DIEL(9.77D-4,3.11D5,7.30D-2,2.06D5,TE,ALD)                
      ALDB=0.0D0
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALNE(2)=2.2D-13/T4**.759+ALD 
      
      alchne = alne(2)
c
      CALL DIEL(2.65D-3,2.84D5,0.242D0,3.07D5,TE,ALD)                  
      CALL DIELB(.0129D0,-.1779D0,.9353D0,-.0682D0,.4516D0,T4,ALDB)
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALNE(3)=1.5D-12/T4**.693+ALD                             
      alchneii = alne(3)
C
      CALL DIEL(3.69D-3,2.24D5,1.01D0,2.94D5,TE,ALD)                   
      CALL DIELB(3.6781D0,14.1481D0,17.1175D0,-.5017D0,.2313D0,T4,      
     &ALDB)
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALNE(4)=4.4D-12/T4**.675+ALD                             
C
      CALL DIEL(1.18D-2,2.70D5,.391D0,5.50D5,TE,ALD)                   
      CALL DIELB(-.0254D0,5.5365D0,1.70727D1,-0.7225D0,.1702D0,T4,    
     &ALDB)
      ALD=(ALD*BDCORR+ALDB)*SUPR        
c!! PL
      chtrne=0.
      ALNE(5)=9.1D-12/T4**.668+ALD+CHTRNE                             
C
      CALL DIEL(2.44D-2,3.09D5,2.52D0,9.91D5,TE,ALD)                    
      CALL DIELB(-.0141D0,33.8479D0,43.1608D0,-1.6072D0,.1942D0,T4,    
     &ALDB)
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALNE(6)=1.50D-11/T4**.684+ALD                                   
C
      CALL DIEL(3.02D-2,2.83D5,0.445D0,1.73D6,TE,ALD)                   
      CALL DIELB(19.928D0,235.0536D0,152.5096D0,9.1413D0,0.1282D0,
     &T4,ALDB)
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALNE(7)=2.30D-11/T4**.704+ALD                                   
C
      CALL DIEL(6.10D-3,1.68D5,0.254D0,6.13D5,TE,ALD)                   
c     if(telog.lt.5.d0) then
c      badnel=1.152d0
c     else if(telog.ge.5.d0.and.telog.lt.7.6d0) then
c      badnel=-54.407d0+23.205d0*telog-3.1176d0*telog**2+.13976d0*
c    & telog**3
c     else if(telog.ge.7.6d0) then
c      badnel=3.238d0
c     endif
c     ALD=badnel*ALD
      CALL DIELB(5.4751D0,203.9751D0,86.9016D0,-7.4568D0,2.5145D0,
     &T4,ALDB)
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALNE(8)=3.46D-11/T4**.742+ALD                                   
C
      CALL DIEL(2.52D-1,1.40D7,0.304D0,1.80D6,TE,ALD)                   
      ALD=ALD*BDCORR*SUPR                                       
      ALNE(9)=5.70D-11/T4**.803+ALD                                   
C
      CALL DIEL(1.02D-1,1.10D7,0.296D0,2.240D6,TE,ALD)                   
cf!!
      BDCORR=0.7D0
      ALD=ALD*BDCORR*SUPR                                       
      ALNE(10)=8.60D-11/T4**.769+ALD                           
C
      ALNE(11)=AL(7)                           
C     ***************************************************************   
C     *****                                                             
C     SODIUM       
C     *****                                                             
C     ***************************************************************   
      ALNA(1)=0.0d0
c
      BDCORR=1.D0
      SUPR=1.D0
      ALDB=0.0D0
      IF(DENS.GT.1.D9) SUPR=0.20                                     
C
      CALL DIEL(2.07d-3,4.566D5,2.999D-2,1.932d6,TE,ALD)    
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALNA(2)=3.467D-13/T4**.846+ALD                         
      alchna = alna(2)

C                                                                       
      CALL DIEL(2.237D-3,3.819D5,.1536D0,3.944D5,TE,ALD)                  
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALNA(3)=8.775d-13/T4**.7464+ALD   
      alchnaii = alna(3)
C
      CALL DIEL(4.529D-3,3.259D5,.3923D0,4.918D5,TE,ALD)            
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALNA(4)=3.399D-12/T4**.7054+ALD
C
      CALL DIEL(6.571D-3,2.711D5,.9028D0,5.476D5,TE,ALD)   
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALNA(5)=7.849D-12/T4**.6952+ALD

C
      CALL DIEL(2.087D-2,3.6D5,.3705D0,7.315D5,TE,ALD)        
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALNA(6)=1.447D-11/T4**.6814+ALD                   
C
      CALL DIEL(2.976D-2,3.463D5,1.175D0,8.552d5,TE,ALD) 
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALNA(7)=2.191D-11/T4**.6875+ALD                 
C
      CALL DIEL(3.54d-2,3.097D5,0.322D0,8.129D5,TE,ALD)            
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALNA(8)=3.253D-11/T4**.7075+ALD                             
C
      CALL DIEL(7.352d-3,1.883d5,0.2842D0,9.716d5,TE,ALD)       
      ALD=ALD*BDCORR*SUPR                                       
      ALNA(9)=4.03D-11/T4**.7873+ALD                
C
      CALL DIEL(2.52D-1,1.40D7,0.3095D0,2.18D6,TE,ALD)      
      ALD=ALD*BDCORR*SUPR                                       
      ALNA(10)=4.90D-11/T4**.753+ALD           
C
      CALL DIEL(1.21D-1,1.285D7,0.2935D0,2.631D6,TE,ALD)      
      ALD=ALD*BDCORR*SUPR                                       
      ALNA(11)=8.11D-11/T4**.816+ALD           
C
      ALNA(12)=AL(8)                           
C     ***************************************************************   
C     *****                                                             
C     MAGNESIUM
C     *****                                                             
C     ***************************************************************   
      ALMG(1)=0.                                                       
c
cf!!
      bdcorr=1.
      SUPR=1.D0
      IF(DENS.GT.1.D9) SUPR=0.20                                     
      CALL DIEL(4.49D-4,5.01D4,2.1D-2,2.81D4,TE,ALD)                
      CALL DIELB(1.2044D0,-4.6836D0,7.662D0,-.593D0,1.626D0,T4,ALDB)
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALMG(2)=1.4D-13/T4**.855+ALD                             

      alchmg = almg(2)

C                                                                       
      ALDB=0.0D0
      CALL DIEL(1.95D-3,6.06D5,0.074D0,1.44D6,TE,ALD)                  
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALMG(3)=8.80D-13/T4**.838+ALD                             
      alchmgii = almg(3)
C
      CALL DIEL(5.12D-3,4.69D5,0.323D0,7.55D5,TE,ALD)                   
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALMG(4)=3.50D-12/T4**.734+ALD
C
      CALL DIEL(7.74D-3,3.74D5,.636D0,7.88D5,TE,ALD)                   
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALMG(5)=7.7D-12/T4**.718+ALD
C
      CALL DIEL(1.17D-2,3.28D5,0.807D0,1.02D6,TE,ALD)                    
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALMG(6)=1.40D-11/T4**.716+ALD                                   
C
      CALL DIEL(3.69D-2,4.80D5,0.351D0,9.73D5,TE,ALD)                   
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALMG(7)=2.30D-11/T4**.695+ALD                                   
C
      CALL DIEL(3.63D-2,3.88D5,0.548D0,7.38D5,TE,ALD)                   
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALMG(8)=3.20D-11/T4**.691+ALD                                   
C
      CALL DIEL(4.15D-2,3.39D5,0.233D0,3.82D5,TE,ALD)                   
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALMG(9)=4.6D-11/T4**.711+ALD                                   
C
      CALL DIEL(8.86D-3,2.11D5,0.318D0,1.54D6,TE,ALD)                   
c     if(telog.lt.5.d0) then
c      badnel=0.7197d0
c     else if(telog.ge.5.d0.and.telog.lt.7.6d0) then
c      badnel=-96.484d0+42.163d0*telog-5.9487d0*telog**2+.28085d0*
c    & telog**3
c     else if(telog.ge.7.6d0) then
c      badnel=3.6456d0
c     endif
      ALD=(ALD*BDCORR+ALDB)*SUPR 
      ALMG(10)=5.8D-11/T4**.804+ALD                                   
C
      CALL DIEL(2.52D-1,1.40D7,0.315D0,2.64D6,TE,ALD)                   
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALMG(11)=6.8D-11/T4**.765+ALD                                   
C
      CALL DIEL(1.44D-1,1.50D7,0.291D0,3.09D6,TE,ALD)                   
cf!!
      BDCORR=0.56D0
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALMG(12)=1.10D-10/T4**.829+ALD                                   
C
      ALMG(13)=AL(9)
C     ***************************************************************   
C     *****                                                             
C     ALUMINIUM      
C     *****                                                             
C     ***************************************************************   
      ALAL(1)=0.                                                       
c
      BDCORR=1.D0
      SUPR=1.D0
      IF(DENS.GT.1.D9) SUPR=0.20                                     
      CALL DIEL(2.547D-3,7.311D4,0.4018D0,5.794D4,TE,ALD)                
      CALL DIELB(0.0219D0,-0.4528D0,2.5427D0,-.1678D0,0.2276D0,T4,
     &ALDB)
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALAL(2)=3.98D-13/T4**.8019+ALD                             

      alchal = alal(2)
C                                                                       
      CALL DIEL(1.503D-3,6.621D4,0.06283D0,3.638D4,TE,ALD)                  
      CALL DIELB(0.7086D0,-3.1083D0,7.0422D0,0.5998D0,.4194D0,T4,
     &ALDB)
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALAL(3)=7.197D-13/T4**.7697+ALD                             
      alchalii = alal(3)
C
      ALDB=0.0D0
      CALL DIEL(3.254D-3,7.977D5,0.1825D0,1.072D6,TE,ALD)                   
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALAL(4)=2.20D-12/T4**.8295+ALD                             
C
      CALL DIEL(6.735D-3,7.312D5,5.683D-16,8.689D-13,TE,ALD)                   
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALAL(5)=6.481D-12/T4**.7345+ALD
C
      CALL DIEL(1.14D-2,4.295D5,1.07D0,9.009D5,TE,ALD)                    
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALAL(6)=1.275D-11/T4**.717+ALD                                   
C
      CALL DIEL(1.707D-2,3.689D5,1.232D0,1.396D6,TE,ALD)                   
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALAL(7)=2.049D-11/T4**.709+ALD                                   
C
      CALL DIEL(3.398D-2,4.191D5,.8399D0,1.433D6,TE,ALD)                   
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALAL(8)=3.145D-11/T4**.6915+ALD                                   
C
      CALL DIEL(3.928D-2,3.753D5,0.8177D0,1.257D6,TE,ALD)                   
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALAL(9)=4.308D-11/T4**.697+ALD                                   
C
      CALL DIEL(5.064D-2,3.627D5,0.2657D0,6.541D5,TE,ALD)                   
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALAL(10)=5.951D-11/T4**.7125+ALD                                   
C
      CALL DIEL(1.106D-2,2.301D5,0.672D0,2.46D6,TE,ALD)                   
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALAL(11)=8.343D-11/T4**.8291+ALD                                   
C
      CALL DIEL(2.871D-1,1.622D7,0.3105D0,3.083D6,TE,ALD)                   
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALAL(12)=9.40D-11/T4**.778+ALD                                   
C
      CALL DIEL(1.65D-1,1.728D7,0.2885D0,3.577D6,TE,ALD)                   
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALAL(13)=1.60D-10/T4**.843+ALD                                   
C
      ALAL(14)=AL(10)
C     ***************************************************************   
C     *****                                                             
C     SILICON      
C     *****                                                             
C     ***************************************************************   
      ALSI(1)=0.                                                       
c
c     f!!

      bdcorr=1.
      SUPR=1.D0
      IF(DENS.GT.1.D9) SUPR=0.20                                     
      CALL DIEL(1.100D-3,7.7D4,0.0D0,0.0D0,TE,ALD)                
      CALL DIELB(-0.0219D0,0.4364D0,0.0684D0,-.0032D0,0.1342D0,T4,
     &ALDB)
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALSI(2)=5.90D-13/T4**.601+ALD
c use this instead!      
      call rec_fit_bad(14,2,te,alsi(2))
      call rec_rr_dr(14,2,te,alsi(2))
      
c     replace with total rate from rr_diel_badnell(te)
      alsi(2)=altotbadn(14,2)
      
      alchsi = alsi(2)                          
C                                                                       
      CALL DIEL(5.87D-3,9.63D4,0.753D0,6.46D4,TE,ALD)                  
      CALL DIELB(3.2163D0,-12.0571D0,16.2118D0,-0.5886D0,0.5613D0,T4,
     &ALDB)
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALSI(3)=1.0D-12/T4**.786+ALD
c     use Nahar above instead!
      alsi(3)=alsi_nahar(3)
c     replace with total rate from rr_diel_badnell(te)
      alsi(3)=altotbadn(14,3)
      alchsiii = alsi(3)
      call rec_rr_dr(14,3,te,alsi(3))
C
      CALL DIEL(5.03D-3,8.75D4,0.188D0,4.71D4,TE,ALD)                   
      CALL DIELB(0.1203D0,-2.69D0,19.1943D0,-.1479D0,0.1118D0,T4,
     &ALDB)
      ALD=(ALD*BDCORR+ALDB)*SUPR              
      ALSI(4)=3.70D-12/T4**.693+ALD
      call rec_rr_dr(14,4,te,alsi(4))
c     replace with total rate from rr_diel_badnell(te)
      alsi(4)=altotbadn(14,4)
C
      ALDB=0.0D0
      CALL DIEL(5.43D-3,1.05D6,0.45D0,7.98D5,TE,ALD)                   
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALSI(5)=5.50D-12/T4**.821+ALD
      call rec_rr_dr(14,5,te,alsi(5))
C
      CALL DIEL(8.86D-3,1.14D6,0.0D0,0.0D0,TE,ALD)
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALSI(6)=1.2D-11/T4**.735+ALD
      call rec_rr_dr(14,6,te,alsi(6))
C
      CALL DIEL(1.68D-2,4.85D5,1.8D0,1.03D6,TE,ALD)                   
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALSI(7)=2.11D-11/T4**.716+ALD
      call rec_rr_dr(14,7,te,alsi(7))

      CALL DIEL(2.49D-2,4.15D5,1.88D0,1.91D6,TE,ALD)                   
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALSI(8)=3.00D-11/T4**.702+ALD                                   
C
      CALL DIEL(3.13D-2,3.66D5,2.01D0,2.11D6,TE,ALD)                   
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALSI(9)=4.30D-11/T4**.688+ALD                                   
C
      CALL DIEL(4.25D-2,3.63D5,1.22D0,2.14D6,TE,ALD)                   
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALSI(10)=5.80D-11/T4**.703+ALD                                   
C
      CALL DIEL(6.18D-2,3.88D5,0.303D0,1.12D6,TE,ALD)                   
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALSI(11)=7.70D-11/T4**.714+ALD                                   
C
      CALL DIEL(1.38D-2,2.51D5,1.42D0,3.93D6,TE,ALD)                   
c     if(telog.lt.5.2d0) then
c      badnel=1.424d0
c     else if(telog.ge.5.2d0.and.telog.lt.7.8d0) then
c      t5log=telog-5.d0
c      badnel=1.2821d0+.29422d0*t5log-2.4287d0*t5log**2-1.7778d0*
c    & t5log**3+.3314d0*t5log**4
c     else if(telog.ge.7.8d0) then
c      badnel=2.4893d0
c     endif
      ALD=(ALD*BDCORR+ALDB)*SUPR 
      ALSI(12)=1.20D-10/T4**.855+ALD                                   
C
      CALL DIEL(3.27D-1,1.88D7,0.306D0,3.60D6,TE,ALD)                   
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALSI(13)=1.40D-10/T4**.823+ALD
c     fit to nahar ok to 1e5 K
      alsi(13)=1.05e-10/t4**0.58

C
      CALL DIEL(1.89D-1,1.99D7,0.286D0,4.14D6,TE,ALD)                   
cf!!
      BDCORR=0.6D0
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALSI(14)=2.20D-10/T4**.858+ALD
c     fit to nahar      
      ALSI(14)=1.40D-10/T4**.6096
C
      ALSI(15)=AL(11)
c     fit to nahar
      ALSI(14)=1.80D-10/T4**.57                                  
C     ***************************************************************   
C     *****                                                             
C     SULPHUR
C     *****                                                             
C     ***************************************************************   
      ALSU(1)=0.                                                       
C
cf!!
      BDCORR=1.
      SUPR=1.D0
      IF(DENS.GT.1.D9) SUPR=0.20                                     
      ALDB=0.0D0
C
      CALL DIEL(1.62D-3,1.25D5,0.0D0,0.0D0,TE,ALD)                
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALSU(2)=4.1D-13/T4**.630+ALD
      call rec_fit_bad(16,2,te,alsu(2))

      
c     replace with total rate from rr_diel_badnell(te)
      alsu(2)=altotbadn(16,2)
      alchs = alsu(2)

C                                                                       
      CALL DIEL(1.09D-2,1.92D5,0.012D0,1.80D4,TE,ALD)                  
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALSU(3)=1.8D-12/T4**.686+ALD                             
      alchsii = alsu(3)
      call rec_fit_bad(16,3,te,alsu(3))
      
c     replace with total rate from rr_diel_badnell(te)
      alsu(3)=altotbadn(16,3)
C
      CALL DIEL(3.35D-2,1.89D5,0.0659D0,1.59D5,TE,ALD)              
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALSU(4)=2.7D-12/T4**.745+ALD
      call rec_fit_bad(16,4,te,alsu(4))
      call rec_rr_dr(16,4,te,alsu(4))

c     replace with total rate from rr_diel_badnell(te)
      alsu(4)=altotbadn(16,4)      
C
      CALL DIEL(3.14D-2,1.68D5,0.0689D0,8.04D4,TE,ALD)              
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALSU(5)=5.7D-12/T4**.755+ALD
      call rec_rr_dr(16,5,te,alsu(5))

c     replace with total rate from rr_diel_badnell(te)
      alsu(5)=altotbadn(16,5)      
      
C
      CALL DIEL(1.27D-2,1.38D5,0.187D0,1.71D5,TE,ALD)               
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALSU(6)=1.2D-11/T4**.701+ALD
      call rec_rr_dr(16,6,te,alsu(6))

c     replace with total rate from rr_diel_badnell(te)
      alsu(6)=altotbadn(16,6)      
      
c
      CALL DIEL(1.47D-2,1.80D6,0.129D0,1.75D6,TE,ALD)                   
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALSU(7)=1.70D-11/T4**.849+ALD
      call rec_rr_dr(16,7,te,alsu(7))
C
      CALL DIEL(1.34D-2,6.90D5,1.04D0,2.15D6,TE,ALD)                   
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALSU(8)=2.70D-11/T4**.733+ALD
      call rec_rr_dr(16,8,te,alsu(8))
C
      CALL DIEL(2.38D-2,5.84D5,1.12D0,2.59D6,TE,ALD)                   
      ALD=ALD*BDCORR*SUPR                                       
      ALSU(9)=4.00D-11/T4**.696+ALD
      call rec_rr_dr(16,9,te,alsu(9))
C
      CALL DIEL(3.19D-2,5.17D5,1.40D0,2.91D6,TE,ALD)                   
      ALD=ALD*BDCORR*SUPR
      ALSU(10)=5.50D-11/T4**.711+ALD                                   
C
      CALL DIEL(7.13D-2,6.66D5,1.00D0,2.32D6,TE,ALD)                   
      ALD=ALD*BDCORR*SUPR
      ALSU(11)=7.40D-11/T4**.716+ALD                                   
C
      CALL DIEL(8.0D-2,6.0D5,0.555D0,2.41D6,TE,ALD)                   
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALSU(12)=9.20D-11/T4**.714+ALD                                   
C
      CALL DIEL(7.96D-2,5.09D5,1.63D0,6.37D6,TE,ALD)                   
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALSU(13)=1.40D-10/T4**.755+ALD                                   
C
      CALL DIEL(1.34D-2,2.91D5,0.304D0,1.04D6,TE,ALD)                   
c     if(telog.lt.5.2d0) then
c      badnel=2.029D0
c     else if(telog.ge.5.2d0.and.telog.lt.8.2d0) then
c      t5log=telog-5.d0
c      badnel=2.5237d0-3.6277d0*t5log+6.1681d0*t5log**2-2.0387d0*
c    & t5log**3+.1995d0*t5log**4
c     else if(telog.ge.8.2d0) then
c      badnel=8.19D0
c     endif
      ALD=(ALD*BDCORR+ALDB)*SUPR 
      ALSU(14)=1.70D-10/T4**.832+ALD                                   
C
      CALL DIEL(4.02D-1,2.41D7,0.298D0,4.67D6,TE,ALD)                   
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALSU(15)=2.00D-10/T4**.806+ALD                                   
C
      CALL DIEL(2.41D-1,2.54D7,0.281D0,5.3D6,TE,ALD)                   
cf!!
      BDCORR=0.6D0
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALSU(16)=2.90D-10/T4**.84+ALD                                   
C
      ALSU(17)=AL(12)
C     ***************************************************************   
C     *****                                                             
C     ARGON
C     *****                                                             
C     ***************************************************************   
      ALAR(1)=0.                                                       
C
cf!!
      BDCORR=1.
      SUPR=1.D0
      IF(DENS.GT.1.D9) SUPR=0.20                                     
      ALDB=0.0D0
C
      CALL DIEL(1.0D-3,3.20D5,5.0D-3,3.10D5,TE,ALD)                
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALAR(2)=3.77D-13/T4**.651+ALD
      alchar = alar(2)

C                                                                       
      CALL DIEL(1.10D-2,2.90D5,0.045D0,5.50D5,TE,ALD)                  
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALAR(3)=1.95D-12/T4**.752+ALD                             
      alcharii = alar(3)
C
      CALL DIEL(3.40D-2,2.39D5,0.057D0,6.00D5,TE,ALD)              
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALAR(4)=3.23D-12/T4**.869+ALD
      
c     replace with total rate from rr_diel_badnell(te)
      alar(4)=altotbadn(18,4)            
C
      CALL DIEL(6.85D-2,2.56D5,0.087D0,3.81D5,TE,ALD)              
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALAR(5)=6.03D-12/T4**.812+ALD
      
c     replace with total rate from rr_diel_badnell(te)
      alar(5)=altotbadn(18,5)      

c     replace with total rate from rr_diel_badnell(te)
      alar(5)=altotbadn(18,5)                  
C
      CALL DIEL(9.00D-2,2.50D5,0.0769D0,3.30D5,TE,ALD)               
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALAR(6)=9.12D-12/T4**.811+ALD
c use this instead!
      call rec_fit_bad(18,6,te,alar(6))
      call rec_rr_dr(18,6,te,alar(6))

c     replace with total rate from rr_diel_badnell(te)
      alar(6)=altotbadn(18,6)                  
c
      CALL DIEL(6.35D-2,2.10D5,0.14D0,2.15D5,TE,ALD)                   
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALAR(7)=1.58D-11/T4**.793+ALD
      call rec_rr_dr(18,7,te,alar(7))
c     replace with total rate from rr_diel_badnell(te)
      alar(7)=altotbadn(18,7)                  
C
      CALL DIEL(2.60D-2,1.80D5,0.12D0,2.15D5,TE,ALD)                   
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALAR(8)=2.69D-11/T4**.744+ALD
      call rec_rr_dr(18,8,te,alar(8))
C
      CALL DIEL(1.70D-2,2.70D6,0.10D0,3.30D6,TE,ALD)                   
      ALD=ALD*BDCORR*SUPR                                       
      ALAR(9)=3.55D-11/T4**.910+ALD
      call rec_rr_dr(18,9,te,alar(9))
C
      CALL DIEL(2.10D-2,8.30D5,1.92D0,3.50D6,TE,ALD)                   
      ALD=ALD*BDCORR*SUPR
      ALAR(10)=4.90D-11/T4**.801+ALD
      call rec_rr_dr(18,10,te,alar(10))
C
      CALL DIEL(3.50D-2,6.95D5,1.66D0,3.60D6,TE,ALD)                   
      ALD=ALD*BDCORR*SUPR
      ALAR(11)=6.92D-11/T4**.811+ALD
      call rec_rr_dr(18,11,te,alar(11))
C
      CALL DIEL(4.30D-2,6.05D5,1.67D0,3.80D6,TE,ALD)                   
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALAR(12)=9.55D-11/T4**.793+ALD                                   
C
      CALL DIEL(7.13D-2,6.68D5,1.400d0,2.90D6,TE,ALD)                   
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALAR(13)=1.23D-10/T4**.702+ALD                                   
C
      CALL DIEL(9.60D-2,6.50D5,1.31D0,3.60D6,TE,ALD)                   
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALAR(14)=1.58D-10/T4**.79+ALD                                   
C
      CALL DIEL(8.50D-2,5.30D5,1.02D0,2.80D6,TE,ALD)                   
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALAR(15)=2.14D-10/T4**.774+ALD


      CALL DIEL(1.70D-2,3.55D5,0.245D0,1.10D6,TE,ALD)                   
c     if(telog.lt.5.2d0) then
c      badnel=3.428D0
c     else if(telog.ge.5.2d0.and.telog.lt.8.2d0) then
c      t5log=telog-5.d0
c      badnel=4.8511d0-8.9917d0*t5log+9.9552d0*t5log**2-2.8855d0*
c    & t5log**3+.25046d0*t5log**4
c     else if(telog.ge.8.2d0) then
c      badnel=9.73D0
c     endif
      ALD=(ALD*BDCORR+ALDB)*SUPR 
      ALAR(16)=2.63D-10/T4**.907+ALD                                   
C
      CALL DIEL(4.76D-1,3.01D7,0.294D0,6.05D6,TE,ALD)                   
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALAR(17)=3.10D-10/T4**.819+ALD                                   
C
      CALL DIEL(2.97D-1,3.13D7,0.277D0,6.54D6,TE,ALD)                   
cf!!
      BDCORR=0.6D0
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALAR(18)=4.40D-10/T4**.854+ALD                                   
C
      ALAR(19)=AL(12)

      
C     ***************************************************************   
C     *****                                                             
C     CALCIUM
C     *****                                                             
C     ***************************************************************   
      ALCA(1)=0.                                                       
C
cf!!
      BDCORR=1.
      SUPR=1.D0
      IF(DENS.GT.1.D9) SUPR=0.20                                     
      ALDB=0.0D0
C
      CALL DIEL(3.28D-4,3.46D4,0.0907D0,1.64D4,TE,ALD)                
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALCA(2)=1.12D-13/T4**.900+ALD                             
      alchca = alca(2)
C                                                                       
      CALL DIEL(5.84D-2,3.85D5,0.11D0,2.45D5,TE,ALD)                  
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALCA(3)=6.78D-13/T4**.800+ALD                             
      alchcaii = alca(3)
C
      CALL DIEL(1.12D-1,4.08D5,0.0174D0,4.27D5,TE,ALD)              
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALCA(4)=3.96D-12/T4**.700+ALD
C
      CALL DIEL(1.32D-1,3.82D5,0.132D0,6.92D5,TE,ALD)              
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALCA(5)=7.08D-12/T4**.780+ALD
C
      CALL DIEL(1.33D-1,3.53D5,0.114D0,8.78D5,TE,ALD)               
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALCA(6)=1.07D-11/T4**.840+ALD
c     replace with total rate from rr_diel_badnell(te)
      alca(6)=altotbadn(20,6)                  
C
      CALL DIEL(1.26D-1,3.19D5,0.162D0,7.43D5,TE,ALD)                   
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALCA(7)=1.80D-11/T4**.82+ALD
c     replace with total rate from rr_diel_badnell(te)
      alca(7)=altotbadn(20,7)                        
C
      CALL DIEL(1.39D-1,3.22D5,0.0878D0,6.99D5,TE,ALD)                   
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALCA(8)=2.40D-11/T4**.820+ALD
c     replace with total rate from rr_diel_badnell(te)
      alca(8)=altotbadn(20,8)                        
C
      CALL DIEL(9.55D-2,2.47D5,0.263D0,4.43D5,TE,ALD)                   
      ALD=ALD*BDCORR*SUPR                                       
      ALCA(9)=3.76D-11/T4**.810+ALD
c     replace with total rate from rr_diel_badnell(te)
      alca(9)=altotbadn(20,9)                        
C
      CALL DIEL(4.02D-2,2.29D5,0.0627D0,2.81D5,TE,ALD)                   
      ALD=ALD*BDCORR*SUPR
      ALCA(10)=5.04D-11/T4**.780+ALD
c     replace with total rate from rr_diel_badnell(te)
      alca(10)=altotbadn(20,10)                        
C
      CALL DIEL(4.19D-2,3.73D6,0.0616D0,5.84D6,TE,ALD)                   
      ALD=ALD*BDCORR*SUPR
      ALCA(11)=6.46D-11/T4**.900+ALD                                   
C
      CALL DIEL(2.57D-2,9.26D5,2.77D0,4.89D6,TE,ALD)                   
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALCA(12)=8.51D-11/T4**.820+ALD                                   
C
      CALL DIEL(4.45D-2,7.96D5,2.23D0,4.62D6,TE,ALD)                   
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALCA(13)=1.18D-10/T4**.810+ALD                                   
C
      CALL DIEL(5.48D-2,6.90D5,2.00D0,4.52D6,TE,ALD)                   
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALCA(14)=1.58D-10/T4**.80+ALD                                   
C
      CALL DIEL(7.13D-2,6.70D5,1.82D0,3.20D6,TE,ALD)                   
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALCA(15)=2.04D-10/T4**.730+ALD                                   
C
      CALL DIEL(9.03D-2,4.72D5,0.424D0,1.37D6,TE,ALD)                   
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALCA(16)=2.60D-10/T4**.800+ALD                                   
C
      CALL DIEL(1.10D-1,5.67D7,0.243D0,4.41D6,TE,ALD)                   
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALCA(17)=3.24D-10/T4**.780+ALD                                   
C
      CALL DIEL(2.05D-2,4.21D5,0.185D0,2.27D6,TE,ALD)                   
c     if(telog.lt.5.2d0) then
c      badnel=4.047D0
c     else if(telog.ge.5.2d0.and.telog.lt.8.2d0) then
c      t5log=telog-5.d0
c      badnel=5.9299d0-11.571d0*t5log+11.271d0*t5log**2-2.4872d0*
c    & t5log**3+.10303d0*t5log**4
c     else if(telog.ge.8.2d0) then
c      badnel=13.623D0
c     endif
      ALD=(ALD*BDCORR+ALDB)*SUPR 
      ALCA(18)=3.81D-10/T4**.850+ALD                                   
C
      CALL DIEL(5.49D-1,3.65D7,0.292D0,7.25D6,TE,ALD)                   
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALCA(19)=4.60D-10/T4**.833+ALD                                   
C
      CALL DIEL(3.55D-1,3.78D7,0.275D0,7.68D6,TE,ALD)                   
cf!!
      BDCORR=0.64D0
      ALD=(ALD*BDCORR+ALDB)*SUPR                                       
      ALCA(20)=5.60D-10/T4**.839+ALD                                   
C
      ALCA(21)=AL(13)
C     ***************************************************************   
C     *****                                                             
C     IRON
C     *****                                                             
C     ***************************************************************   
      ALFE(1)=0.                                                       
      ALFE(2)=1.42D-13/T4**.891
      ALFE(3)=1.02D-12/T4**.843
      ALFE(4)=3.32D-12/T4**.746
      ALFE(5)=7.80D-12/T4**.682
      ALFE(6)=1.51D-11/T4**.699
      ALFE(7)=2.62D-11/T4**.728
      ALFE(8)=4.12D-11/T4**.759
      ALFE(9)=6.05D-11/T4**.79
      ALFE(10)=8.13D-11/T4**.81
      ALFE(11)=1.09D-10/T4**.829
      ALFE(12)=1.33D-10/T4**.828
      ALFE(13)=1.64D-10/T4**.834
      ALFE(14)=2.00D-10/T4**.836

c total rate from dielectronic!

      do k=1,5

         alfe(k+1) = 0.

      enddo

      alfe(14) = 0.



      call radreclfe(1.46d-10,0.597d0,5.22d-2,te,alfe(15))
      call radreclfe(1.68d-10,0.602d0,5.07d-2,te,alfe(16))
      call radreclfe(1.91d-10,0.601d0,5.10d-2,te,alfe(17))
      call radreclfe(2.24d-10,0.579d0,5.49d-2,te,alfe(18))
      call radreclfe(2.59d-10,0.567d0,5.65d-2,te,alfe(19))
      call radreclfe(2.96d-10,0.557d0,5.79d-2,te,alfe(20))
      call radreclfe(3.16d-10,0.534d0,6.02d-2,te,alfe(21))
      call radreclfe(3.49d-10,0.521d0,6.22d-2,te,alfe(22))
      call radreclfe(3.91d-10,0.523d0,6.15d-2,te,alfe(23))
      call radreclfe(4.33d-10,0.531d0,5.77d-2,te,alfe(24))
      call radreclfe(4.77d-10,0.537d0,5.49d-2,te,alfe(25))
      call radreclfe(5.84d-10,0.546d0,4.02d-2,te,alfe(26))
      ALFE(27)=AL(14)
      call dielrecfe(te,aldfe)
      do k=2,26
         alfe(k)=alfe(k)+aldfe(k)
      enddo
      do k=2,27
      enddo
 9239 format('k,aldfe,alffe',i5,1pe12.4,10e12.4)

      alchfe = alfe(2)
      alchfeii = alfe(3)

C     ***************************************************************   
C     *****                                                             
C     NICKEL
C     *****                                                             
C     ***************************************************************   
      ALNI(1)=0.                                                       
C
cf!!
      BDCORR=1.
      SUPR=1.D0
      IF(DENS.GT.1.D9) SUPR=0.20                                     
      ALDB=0.0D0
C
      CALL DIEL(1.41D-3,9.82D4,0.469D0,1.01D5,TE,ALD)
      ALD=BDCORR*SUPR*ALD
      ALNI(2)=3.60D-13/T4**.700+ALD                             
C                                                                      
      CALL DIEL(5.20D-3,2.01D5,0.357D0,1.91D5,TE,ALD)
      ALD=BDCORR*SUPR*ALD
      ALNI(3)=1.0D-12/T4**.700+ALD
C
      CALL DIEL(1.38D-2,3.05D5,0.281D0,2.32D5,TE,ALD)
      ALD=BDCORR*SUPR*ALD
      ALNI(4)=1.4D-12/T4**.700+ALD
C
      CALL DIEL(2.30D-2,4.2D5,0.128D0,3.18D5,TE,ALD)
      ALD=BDCORR*SUPR*ALD
      ALNI(5)=1.60D-12/T4**.700+ALD
C
      CALL DIEL(4.19D-2,5.56D5,0.0417D0,4.55D5,TE,ALD)
      ALD=BDCORR*SUPR*ALD
      ALNI(6)=3.85D-12/T4**.746+ALD
C
      CALL DIEL(6.83D-2,6.72D5,0.0558D0,5.51D5,TE,ALD)
      ALD=BDCORR*SUPR*ALD
      ALNI(7)=9.05D-12/T4**.682+ALD
C
      CALL DIEL(1.22D-1,7.93D5,0.0346D0,5.28D5,TE,ALD)
      ALD=BDCORR*SUPR*ALD
      ALNI(8)=1.75D-11/T4**.699+ALD
C
      CALL DIEL(3.00D-1,9.00D5,0.0D0,0.0D0,TE,ALD)
      ALD=BDCORR*SUPR*ALD
      ALNI(9)=3.04D-11/T4**.728+ALD
C
      CALL DIEL(1.50D-1,1.00D6,1.90D0,5.50D5,TE,ALD)
      ALD=BDCORR*SUPR*ALD
      ALNI(10)=8.91D-11/T4**.759+ALD
C
      CALL DIEL(6.97D-1,7.81D5,0.277D0,8.87D5,TE,ALD)
      ALD=BDCORR*SUPR*ALD
      ALNI(11)=1.19D-10/T4**.790+ALD
C
      CALL DIEL(7.09D-1,7.64D5,0.135D0,1.80D6,TE,ALD)
      ALD=BDCORR*SUPR*ALD
      ALNI(12)=1.50D-10/T4**.810+ALD
C
      CALL DIEL(6.44D-1,7.44D5,0.134D0,1.25D6,TE,ALD)
      ALD=BDCORR*SUPR*ALD
      ALNI(13)=1.91D-10/T4**.829+ALD
C
      CALL DIEL(5.25D-1,6.65D5,0.192D0,1.89D6,TE,ALD)
      ALD=BDCORR*SUPR*ALD
      ALNI(14)=2.29D-10/T4**.828+ALD
C
      CALL DIEL(4.46D-1,5.97D5,0.332D0,8.84D5,TE,ALD)
      ALD=BDCORR*SUPR*ALD
      ALNI(15)=2.63D-10/T4**.834+ALD
C
      CALL DIEL(3.63D-1,5.24D5,0.337D0,1.29D6,TE,ALD)
      ALD=BDCORR*SUPR*ALD
      ALNI(16)=3.16D-10/T4**.836+ALD
C
      CALL DIEL(3.02D-1,4.96D5,0.121D0,6.24D5,TE,ALD)
      ALD=BDCORR*SUPR*ALD
      ALNI(17)=3.63D-10/T4**.84+ALD
C
      CALL DIEL(1.02D-1,4.46D5,0.0514D0,1.59D6,TE,ALD)
      ALD=BDCORR*SUPR*ALD
      ALNI(18)=4.03D-10/T4**.846+ALD
C
      CALL DIEL(2.70D-1,8.49D6,0.183D0,8.01D6,TE,ALD)
      ALD=BDCORR*SUPR*ALD
      ALNI(19)=4.73D-10/T4**.850+ALD
C
      CALL DIEL(4.67D-2,1.36D6,7.56D0,9.32D6,TE,ALD)
      ALD=BDCORR*SUPR*ALD
      ALNI(20)=5.25D-10/T4**.836+ALD
C
      CALL DIEL(8.35D-2,1.23D6,4.55D0,9.45D6,TE,ALD)
      ALD=BDCORR*SUPR*ALD
      ALNI(21)=5.75D-10/T4**.824+ALD
C
      CALL DIEL(9.96D-2,1.06D6,4.87D0,9.45D6,TE,ALD)
      ALD=BDCORR*SUPR*ALD
      ALNI(22)=6.38D-10/T4**.816+ALD
C
      CALL DIEL(1.99D-1,1.25D6,2.19D0,8.01D6,TE,ALD)
      ALD=BDCORR*SUPR*ALD
      ALNI(23)=7.08D-10/T4**.811+ALD
C
      CALL DIEL(2.40D-1,1.23D6,1.15D0,7.57D6,TE,ALD)
      ALD=BDCORR*SUPR*ALD
      ALNI(24)=7.94D-10/T4**.808+ALD
C
      CALL DIEL(1.15D-1,3.32D5,1.23D0,2.64D6,TE,ALD)
      ALD=BDCORR*SUPR*ALD
      ALNI(25)=8.71D-10/T4**.800+ALD
C
      CALL DIEL(3.16D-2,6.45D5,0.132D0,1.93D6,TE,ALD)
c     if(telog.lt.5.6d0) then
c      badnel=5.547D0
c     else if(telog.ge.5.6d0.and.telog.lt.8.6d0) then
c      t5log=telog-5.d0
c      badnel=20.01d0-38.556d0*t5log+27.173d0*t5log**2-5.2864d0*
c    & t5log**3+.23104d0*t5log**4
c     else if(telog.ge.8.6d0) then
c      badnel=25.532D0
c     endif
      ALD=bdcorr*SUPR*ALD
      ALNI(26)=8.91D-10/T4**.718+ALD
C
      CALL DIEL(8.03D-1,6.65D7,0.289D0,1.19D7,TE,ALD)
      ALD=BDCORR*SUPR*ALD
      ALNI(27)=1.40D-9/T4**.842+ALD
C
      CALL DIEL(5.75D-1,6.81D7,0.286D0,9.08D6,TE,ALD)
cf!!
      BDCORR=0.59D0
      ALD=BDCORR*SUPR*ALD
      ALNI(28)=1.90D-9/T4**.868+ALD
C
      ALNI(29)=AL(15)

      do i=1,6
         alrec(6,i+1)=alc(i+1)
      enddo

      do i=1,7
         alrec(7,i+1)=aln(i+1)
      enddo


      do i=1,8
         alrec(8,i+1)=alo(i+1)
      enddo

      do i=1,10
         alrec(10,i+1)=alne(i+1)
      enddo

      do i=1,11
         alrec(11,i+1)=alna(i+1)
      enddo


      do i=1,12
         alrec(12,i+1)=almg(i+1)
      enddo

      do i=1,13
         alrec(13,i+1)=alal(i+1)
      enddo

      do i=1,14
         alrec(14,i+1)=alsi(i+1)
      enddo

      do i=1,16
         alrec(16,i+1)=alsu(i+1)
      enddo

      do i=1,18
         alrec(18,i+1)=alar(i+1)
      enddo

      do i=1,20
         alrec(20,i+1)=alca(i+1)
      enddo

      do i=1,26
         alrec(26,i+1)=alfe(i+1)
      enddo

c replace the fits below with the ones from dielbadnell.f !!      
      


c add charge transfer

      xhi=xion(ik,1,1)*abn(1)
      xhei=xion(ik,2,1)*abn(2)
c!!!
c      xhi=0.
c      xhei=0.

      do iel=3,14
         do ion=1,nion(iel)

            if(ion.le.5)  then
               
c CT with H I

               if(ion.ge.2.and.ion.le.5) then

                  call convrtel(iel,nelem)

                  cthr(iel,ion) = HCTRecom(ion,nelem,te)*xhi/xel

               else
                  cthr(iel,ion) =0.
               endif

c CT with He I

               if(ion.ge.2.and.ion.le.5) then

                  cther(iel,ion) = hectrec(ion,iel,te)*xhei/xel
               else
                  cther(iel,ion) = 0.
               endif            
               
               alrold= alrec(iel,ion)
               alrec(iel,ion) = alrec(iel,ion) + cthr(iel,ion) + 
     &              cther(iel,ion)

            endif
         enddo
      enddo

      alchfe = alfe(2)
      alchfeii = alfe(3)

      alchfe = alrec(14,2)
      alchfeii = alrec(14,3)

      RETURN                                                            
      END                                                               

      subroutine dielrecfe(te,alfe)
      implicit real*8(a-h,o-z)
      character*80 ch
      character*100 dum
      character*3 str
      character*4 str4
      common/initferates/initferec,initfecolli,initfecolliea
      dimension tefel(100),recfenahar(5,100)
      dimension c14(10),e14(10)
      data c14/3.55E-4, 2.40E-3, 7.83E-3, 1.10E-2, 3.30E-2, 1.45E-1,
     & 8.50E-2, 2.59E-2, 8.93E-3, 9.80E-3/
      data e14/2.19E-2, 1.79E-1, 7.53E-1, 2.21E0, 9.57E0, 3.09E1, 
     &  6.37E1, 2.19E2, 1.50E3, 7.86E3/
      save edr,cdr,tefel,recfenahar
      real*8 aq(15),bq(15),t0(15),t1(15),cq(15),t2(15)
      real*8 ed(15,9),cd(15,9)
      real*8 alrr(27),ald(27),alfe(27),aldfe(27)
      dimension edr(4,26),cdr(4,26),cfe(26,6),efe(26,6)
c     fit params for Fe+6 to Fe+14 (from cloudy)
      real*8 cfe7(6),cfe8(6),cfe9(6),cfe10(6),cfe11(6),cfe12(6),cfe13(6),
     &     cfe14(6),cfe15(6)
      real*8 csh(26,13),esh(26,13),tminsh(26),tmaxsh(26)
c data from Cloudy c17.03 from Gu. rec of Fe VII->Fe VI to Fe XV->Fe XIV
      data cfe7/2.50507e-11, 5.60226e-11, 1.85001e-10, 3.57495e-9, 1.66321e-7,0./
      data cfe8/ 9.19610e-11, 2.92460e-10, 1.02120e-9, 1.14852e-8, 3.25418e-7, 0./
      data cfe9/9.02625e-11, 6.22962e-10, 5.77545e-9, 1.78847e-8, 3.40610e-7, 0./ 
      data cfe10/9.04286e-12, 9.68148e-10, 4.83636e-9, 2.48159e-8, 3.96815e-7, 0./ 
      data cfe11/6.77873e-10, 1.47252e-9, 5.31202e-9, 2.54793e-8, 3.47407e-7, 0./
      data cfe12/1.29742e-9, 4.10172e-9, 1.23605e-8, 2.33615e-8, 2.97261e-7, 0./
      data cfe13/8.78027e-10,2.31680e-9,3.49333e-9,1.16927e-8,8.18537e-8,1.54740e-7/
      data cfe14/2.23178e-10, 1.87313e-9, 2.86171e-9, 1.38575e-8, 1.17803e-7, 1.06251e-7/
      data cfe15/2.17263e-10, 7.35929e-10, 2.81276e-9, 1.32411e-8, 1.15761e-7, 4.80389e-8/

      real*8 efe7(6),efe8(6),efe9(6),efe10(6),efe11(6),efe12(6),efe13(6),
     &     efe14(6),efe15(6)
      data efe7/8.30501e-2, 8.52897e-1, 3.40225e0, 2.23053e1, 6.80367e1, 0./
      data efe8/1.44392e-1, 9.23999e-1, 5.45498e0, 2.04301e1, 7.06112e1, 0. /
      data efe9/ 5.79132e-2, 1.27852e0, 3.22439e0, 1.79602e1, 6.96277e1, 0./
      data efe10/ 1.02421e-1, 1.79393e0, 4.83226e0, 1.91117e1, 6.80858e1, 0./
      data efe11/ 1.24630e-1, 6.86045e-1, 3.09611e0, 1.44023e1, 6.42820e1, 0./
      data efe12/ 1.34459e-1, 6.63028e-1, 2.61753e0, 1.30392e1, 6.10222e1, 0./
      data efe13/ 7.79748e-2, 5.35522e-1, 1.88407e0, 8.38459e0, 3.38613e1, 7.89706e1/
      data efe14/ 8.83019e-2, 6.12756e-1, 2.36035e0, 9.61736e0, 3.64467e1, 8.72406e1/
      data efe15/ 1.51322e-1, 5.63155e-1, 2.57013e0, 9.08166e0, 3.69528e1, 1.08067e2/
      save aq,bq,t0,t1,cq,t2,cfe,efe
      save ed,cd,csh,esh,tminsh,tmaxsh
      dimension abadr(27),bbadr(27),cbadr(27),t0badr(27),t1badr(27),t2badr(27)
      dimension cbadi(7),ebadi(7),cbad(27,7),ebad(27,7),nkb(27)
      save abadr,bbadr,cbadr,t0badr,t1badr,t2badr,cbad,ebad,nkb

      if(initferec.eq.0) then

      do j=1,6
         cfe(7,j)=cfe7(j)
         cfe(8,j)=cfe8(j)
         cfe(9,j)=cfe9(j)
         cfe(10,j)=cfe10(j)
         cfe(11,j)=cfe11(j)
         cfe(12,j)=cfe12(j)
         cfe(13,j)=cfe13(j)
         cfe(14,j)=cfe14(j)
         cfe(15,j)=cfe15(j)
         efe(7,j)=efe7(j)
         efe(8,j)=efe8(j)
         efe(9,j)=efe9(j)
         efe(10,j)=efe10(j)
         efe(11,j)=efe11(j)
         efe(12,j)=efe12(j)
         efe(13,j)=efe13(j)
         efe(14,j)=efe14(j)
         efe(15,j)=efe15(j)
      enddo

         do k=1,15
            do i=1,9
               cd(k,i)=0.
               ed(k,i)=0.
            enddo
         enddo
c      stop
c diel (?) from unknown source
         open(23,file='./ATDAT/recrate.dat',status='old')
         do i=1,25
            read(23,*)(edr(k,i),k=1,4),(cdr(k,i),k=1,4)
         enddo

         close (23)


c Fe II - Fe VI recomb. from Nahar. Fe VI -> Fe V from Nahar 1999

         open(23,file='./ATDAT/feii_vi_rec_nahar.dat',status='old')

         Do i=1,36
            read(23,91)ch
 91         format(a)
         enddo

         do i=1,61
            read(23,*)tefel(i),(recfenahar(k,i),k=1,5)
         enddo
         
         close (23)

c     Fe XXVII - Fe XIII + Fe IX from Badnells compilations
         open(23,file='./ATDAT/Fe_rec_rates_Badnell_DR.dat',status='old')
         do i=1,27
            do k=1,7
               cbad(i,k)=0.
               ebad(i,k)=0.
            enddo
         enddo
         do i=1,15
            if(i <=2) then
               nk=3
            elseif(i <=4) then
               nk=5
            elseif(i <=9) then
               nk=6
            elseif(i <=15) then
               nk=7
            endif            
            read(23,*)iz,nel,mm,iw,(cbadi(k),k=1,nk)
            read(23,*)iz,nel,mm,iw,(ebadi(k),k=1,nk)
            io=iz-nel+1
            nkb(io)=nk
            do k=1,nk
               cbad(io,k)=cbadi(k)
               ebad(io,k)=ebadi(k)
            enddo
         enddo

         close(23)
         
c RR rates from Badnells compil
         open(23,file='./ATDAT/Fe_rec_rates_Badnell_RR.dat',status='old')
         do i=1,27
            abadr(i)=0.
            bbadr(i)=0.
            t0badr(i)=0.
            t1badr(i)=0.
            cbadr(i)=0.
            t2badr(i)=0.
         enddo

         do i=1,15
            if(i<=14) then
               read(23,*)iz,nel,mm,iw,ai,bi,t0i,t1i
               ci=0.
               t2i=0.
            elseif(i==15) then
               read(23,*)iz,nel,mm,iw,ai,bi,t0i,t1i,ci,t2i
            endif
            io=iz-nel+1
            abadr(io)=ai
            bbadr(io)=bi
            t0badr(io)=t0i
            t1badr(io)=t1i
            cbadr(io)=ci
            t2badr(io)=t2i
         enddo

         close(23)

c Fe VIII-VII and IX->VIII from Schmidt+ 2008 
         open(23,file='./ATDAT/rec_fe_viii_ix.dat',status='old')
         do i=1,9
            read(23,*)iq,cd(8,i),ed(8,i),cd(9,i),ed(9,i)
         enddo
         close (23)

c Fe X-IX and XI -> X Lei.. 2009         
         open(23,file='./ATDAT/rec_fe_x_xi.dat',status='old')
         read(23,*)dum
         read(23,*)dum
         do i=1,9
            read(23,*)iq,cd(10,i),ed(10,i),cd(11,i),ed(11,i)
         enddo
         close (23)

c     Fe XII-XI Novotny Lei.. 2012
         open(23,file='./ATDAT/rec_fe_xii.dat',status='old')
         read(23,*)dum
         do i=1,7
            read(23,*)iq,cd(12,i),ed(12,i)
         enddo
         close (23)

c     Fe XIII-XII Hahn+ 2014
         open(23,file='./ATDAT/rec_fe_xiii.dat',status='old')
         read(23,*)dum
         do i=1,8
            read(23,*)iq,ed(13,i),cd(13,i)
         enddo
         close (23)                  
         
c radiative rec fits from Badnell (2006)         
         open(23,file='./ATDAT/rad_rec_fe_badnell.dat',status='old')
         read(23,*)dum
         do i=1,8
c     Badnell uses # of el in M shell
            ion=15-i
            read(23,*)iq,aq(ion),bq(ion),t0(ion),t1(ion),cq(ion),t2(ion)

         enddo
         close (23)


c     Fe X-XXII Schippers 2010 astroph 1002.3678.pdf
         do i=8,23
            do j=1,13
               csh(i,j)=0.
               esh(i,j)=0.
            enddo
         enddo

         open(23,file='./ATDAT/schippers_2010_dr.dat',status='old')
         do i=1,5
            read(23,*)dum
         enddo
         do j=1,6
            read(23,*)str,(csh(i,j),i=8,11),(csh(i,j),i=14,16)
 920        format(a3,1pe9.3,10e9.3)
         enddo
         read(23,*)str,(csh(i,7),i=8,11),(csh(i,7),i=14,15)
         read(23,*)str,(csh(i,8),i=8,11)
         read(23,*)str,(csh(i,9),i=8,11)
         do j=10,13
            read(23,*)str,csh(8,j)
         enddo
         
         do j=1,6
            read(23,*)str,(esh(i,j),i=8,11),(esh(i,j),i=14,16)
         enddo
         read(23,*)str,(esh(i,7),i=8,11),(esh(i,7),i=14,15)
         read(23,*)str,(esh(i,8),i=8,11)
         read(23,*)str,(esh(i,9),i=8,11)
         do j=10,13
            read(23,*)str,esh(8,j)
         enddo
         read(23,*)str4,(tminsh(i),i=8,11),(tminsh(i),i=14,16)
         read(23,*)str4,(tmaxsh(i),i=8,11),(tmaxsh(i),i=14,16)

         do i=1,2
            read(23,*)dum
         enddo
         do j=1,2
            read(23,*)str,(csh(i,j),i=17,23)
         enddo
         do j=3,7
            read(23,*)str,(csh(i,j),i=17,22)
         enddo
         read(23,*)str,(csh(i,8),i=17,21)
         read(23,*)str,(csh(i,9),i=17,18)
         read(23,*)str,csh(17,10)

         do j=1,2
            read(23,*)str,(esh(i,j),i=17,23)
         enddo
         do j=3,7
            read(23,*)str,(esh(i,j),i=17,22)
         enddo
         read(23,*)str,(esh(i,8),i=17,21)
         read(23,*)str,(esh(i,9),i=17,18)
         read(23,*)str,esh(17,10)
         
         read(23,*)str4,(tminsh(i),i=17,23)
         read(23,*)str4,(tmaxsh(i),i=17,23)
         
         close (23)                                   

         initferec=1
      endif

      do k=1,26
         alfe(k)=0.d0
         aldfe(k)=0.d0
         alrr(k)=0.d0
         ald(k)=0.d0
      enddo
      

c     diel of Fe VII - Fe XVI from Cloudy

      do i=7,15
         ald(i)=0.
         do j=1,6
            ald(i)=ald(i) + cfe(i,j)*exp(-1.1609e4*efe(i,j)/te)/(te/1.1609e4)**1.5
         enddo
c!! neglect rad recomb!
         alfe(i)=ald(i)
 928     format('diel from Cl ',i5,1pe12.3)
      enddo

c Fe II -> I -- Fe VI -> Fe V rec. Nahar

      tel = log10(te)

      do i=1,61

         if(tel.gt.tefel(i).and.tel.le.tefel(i+1)) then

            it = i

            goto 11

         elseif(tel.le.tefel(1)) then

            it = 0

            goto 11

         elseif(tel.gt.tefel(61)) then

            it = 61

            goto 11

         endif
         

      enddo

 11   continue

      do k=1,5

         if(it.eq.0) then

            aldfe(k+1) = recfenahar(k,1)

         elseif(it.eq.61) then

            aldfe(k+1) = recfenahar(k,61)

         else

            alil=log10(recfenahar(k,i))

            alilp1=log10(recfenahar(k,i+1))

            all = alil + (tel-tefel(i))*(alilp1 - alil)/
     &           (tefel(i+1) - tefel(i))

            aldfe(k+1) = 10.**all
c!!   neglect rad rem. OR is it included (as is usual for Nahar)
            alfe(k+1)=aldfe(k+1)

         endif
      enddo

c     radiative and dielectronic from Badnell's compil. for Fe IX and Fe XIII-XXVII

      do io=9,27
         if(cbad(io,1)>0.) then
            aldfe(io)=0.
            do k=1,nkb(io)
               ald(io) = ald(io) + cbad(io,k)*exp(-ebad(io,k)/te)/te**1.5
            enddo
            if(io==9) then
               b=-bbadr(io) + cbadr(io)*exp(-t2badr(io)/te)
            else
               b=-bbadr(io)
            endif
            alrr(io)=abadr(io)/( (te/t0badr(io))**0.5*
     &           (1.+(te/t0badr(io))**0.5)**(1.-b)*
     &           (1.+(te/t1badr(io))**0.5)**(1.+b))
         endif
         alfe(io)=alrr(io)+ald(io)
      enddo
      
      
c     Fe VIII-XI, XIV-XXIII from Schippers 2010         
      do k=8,23
         if(k.le.11.or.k.ge.14) then
            if(te > tminsh(k).and.te < tmaxsh(k)) then
               ald(k)=0.
               do i=1,13
                  ald(k)=ald(k)+csh(k,i)*exp(-esh(k,i)/te)/te**1.5
 912              format('k,i,cd,ed,ald',2i4,1pe12.3,10e12.3)
               enddo
               alfe(k)=alrr(k)+ald(k)
            endif
         endif
      enddo

c     Fe VIII-> VII and IX -> VIII
      

      do k=8,14
c radiative            
         bp=bq(k)+cq(k)*exp(-t2(k)/te)
         alrr(k)=aq(k)/((te/t0(k))**0.5*(1+(te/t0(k))**0.5)**(1-bp)*
     &        (1.+(te/t0(k))**0.5)*(1.+bp))
c     dielectronic 
         ald(k)=0.
         do i=1,9
            ald(k)=ald(k)+cd(k,i)*exp(-ed(k,i)/te)/te**1.5
         enddo
         alfe(k)=alrr(k)+ald(k)
      enddo
      

c Fe XIV-> XIII rate from Schmidt et al 2006 astroph-0603340
c     this icludes both radiative and dielectronic

c note that e14 is in eV not in K as the other fits      

      tev = te/1.1609e4
      ald(14) = 0.
      do i=1,10
         ald(14) = ald(14) + c14(i)*exp(-e14(i)/tev)/te**1.5
      enddo
      alfe(14)=ald(14)+alrr(14)

      return
      end

      

      subroutine collionfe(te,cofe)
      implicit real*8(a-h,o-z)
      common/initferates/initferec,initfecolli,initfecolliea
      save rf,pot,a,b,c,d
      dimension pot(4,26),a(4,26),b(4,26),c(4,26),d(4,26),cofe(26)
      dimension kmax(26)
      data kmax/8*3,6*2,2*3,8*2,2*1/
      if(initfecolli.eq.0) then
         open(24,file='./ATDAT/colljonfe.dat',status='old')
         do i=1,26
            do k=1,kmax(i)
               read(24,*)rf,pot(k,i),a(k,i),b(k,i),c(k,i),d(k,i)
            enddo
         enddo
         initfecolli=1
      endif
      tev=te/1.1609e4
      do i=1,26
         cofe(i)=0.
         do k=1,kmax(i)
            x=pot(k,i)/tev
            if(x.lt.1.d2) then
               f1=arf1(x)
               f2=arf2(x)
               y1=a(k,i)*(1.d0-x*f1)+b(k,i)*(1.d0+x-x*(2.d0+x)*f1)+
     &            c(k,i)*f1+d(k,i)*x*f2
               cofe(i)=cofe(i)+6.69d-7*y1*dexp(-x)/(x*tev**1.5d0)
            endif               
         enddo
      enddo
      return
      end


      subroutine collionfeea(te,cofeea)
      implicit real*8(a-h,o-z)
      common/initferates/initferec,initfecolli,initfecolliea
      save pot,a,b,c,d,e
      dimension pot(26),a(26),b(26),c(26),d(26),e(26),cofeea(26)
      if(initfecolliea.eq.0) then
         open(25,file='./ATDAT/colljoneafe.dat',status='old')
         do i=1,26
               read(25,*)if,pot(i),a(i),b(i),c(i),d(i),e(i)
         enddo
         initfecolliea=1
      endif
      tev=te/1.1609e4
      do i=1,26
         cofeea(i)=0.
         x=pot(i)/tev
         if(x.lt.1.d2.and.x.gt.0.) then
            f1=arf1(x)
            y1=a(i)+b(i)*(1.d0-x*f1)+c(i)*(1.d0-x*(1.d0-x*f1))+
     &            d(i)*(1.d0-0.5d0*(x-x*x+x*x*x*f1))+e(i)*f1
            cofeea(i)=6.69d-9*y1*dexp(-x)/sqrt(tev)
         endif               
      enddo
      return
      end

      real*8 function HECTRecom(ion,nelem,te)
      implicit real*8(a-h,o-z)
*     ion is stage of ionization, 2 for the ion going to the atom
*     nelem is atomic number of element, 2 up to 30
      integer ion , nelem
      real*8 te
      common/CTHERecomb/ CTHERecomb(6,4,30)

**     local variables
      real*8 tused 
      integer ipIon
*
      ipIon = ion - 1
*

      if( ipIon.gt.4 ) then
*       use statistical charge transfer for ion &gt; 4
        HECTRecom = 1.92e-9 * ipIon
        return
      endif
*
*     Make sure te is between temp. boundaries; set constant outside of range
      tused = max( te,CTHERecomb(5,ipIon,nelem) )
      tused = min( tused , CTHERecomb(6,ipIon,nelem) )
      tused = tused * 1e-4
*
*     the interpolation equation
      HECTRecom = CTHERecomb(1,ipIon,nelem)* 1e-9 * 
     1 (tused**CTHERecomb(2,ipIon,nelem)) *
     2 (1. + 
     3 CTHERecomb(3,ipIon,nelem) * 
     4 exp(CTHERecomb(4,ipIon,nelem)*tused) )
c      write(6,*)' hectrecom ',CTHERecomb(1,ipIon,nelem),hectrecom
*
      end


      subroutine radreclfe(a,al,be,te,alfe)
      implicit real*8(a-h,o-z)
      alfe=a*(te/1.d4)**(-al-be*log10(te/1.e4))
      return
      end


Code by Jim Kingdon, in collaboration with G.J. Ferland
      real*8 function HCTRecom(ion,nelem,te)
      implicit real*8(a-h,o-z)
*     ion is stage of ionization, 2 for the ion going to the atom
*     nelem is atomic number of element, 2 up to 30
*     Example:  O+ + H =&gt; O + H+ is HCTRecom(2,8)
      integer ion , nelem
      real*8 te
      common/CTRecomb/ CTRecomb(6,4,30)
*
*     local variables
      real tused 
      integer ipIon
*
      ipIon = ion - 1
*
c      write(6,*)' hctr ',ion,ipion,nelem

      if( ipIon.gt.4 ) then
*       use statistical charge transfer for ion &gt; 4
        HCTRecom = 1.92e-9 * ipIon
        return
      endif
*
*     Make sure te is between temp. boundaries; set constant outside of range
      tused = max( te,CTRecomb(5,ipIon,nelem) )
      tused = min( tused , CTRecomb(6,ipIon,nelem) )
      tused = tused * 1e-4
*
*     the interpolation equation
      HCTRecom = CTRecomb(1,ipIon,nelem)* 1e-9 * 
     1 (tused**CTRecomb(2,ipIon,nelem)) *
     2 (1. + 
     3 CTRecomb(3,ipIon,nelem) * exp(CTRecomb(4,ipIon,nelem)*tused) )
*
c      write(6,*)' hctrecom ',CTRecomb(1,ipIon,nelem),hctrecom
      end

******************************************************************************
      real*8 function HCTIon(ion,nelem,te)
*     ion is stage of ionization, 1 for atom
*     nelem is atomic number of element, 2 up to 30
*     Example:  O + H+ =&gt; O+ + H is HCTIon(1,8)
      implicit real*8(a-h,o-z)
      real*8 te
      integer ion , nelem
      common/CTIon/ Ctionp(7,4,30)
*
*     local variables
      real*8 tused 
      integer ipIon
*
      ipIon = ion
*
*     Make sure te is between temp. boundaries; set constant outside of range

      tused = max( te,Ctionp(5,ipIon,nelem) )
      tused = min( tused , Ctionp(6,ipIon,nelem) )
      tused = tused * 1e-4
*
*     the interpolation equation
      HCTIon = Ctionp(1,ipIon,nelem)* 1e-9 * 
     1 (tused**Ctionp(2,ipIon,nelem)) *
     2 (1. + 
     3 Ctionp(3,ipIon,nelem) * exp(Ctionp(4,ipIon,nelem)*tused) ) *
     4 exp(-Ctionp(7,ipIon,nelem)/tused)
*
      end


      real*8 function hectrec(ion,iel,te)

c CT recomb. with He I from Oak Ridge base

      implicit real*8(a-h,o-z)
      integer ion , nelem
      real*8 te
      common/CTHERecomb/ CTHERecomb(6,4,30)

      integer ipIon
*
      hectrec=0.

c      write(6,*)'ion, iel,te a ',ion, iel, te

      if(iel.eq.3) then

         if(ion.eq.4) then

            ionm1 = ion-1

            hectrec = HECTRecom(ion,6,te)

         elseif(ion.eq.5) then

c C V + He I -> C IV + He II

            if(te.le.1.d3) then
               telow=1.d3
               hectrec = fkingdon(3.12d-7,-7.37d-2,3.50d1,2.40d0,
     &              1.d+3,1d+4,telow)
            elseif(te.gt.1.d3.and.te.lt.1.d4) then
               hectrec = fkingdon(3.12d-7,-7.37d-2,3.50d1,2.40d0,
     &              1.d+3,1d+4,te)
            elseif(te.ge.1.d4.and.te.lt.3.5d5) then
               hectrec = fkingdon(1.49d-5,2.73d0,5.93d0,-8.74d-2,1d+4,
     &              3.5d+5,te)
            elseif(te.gt.3.5d5.and.te.lt.1.d7) then
               hectrec = fkingdon(5.80d-2,0.73d0,-0.86d0,-9.60d-3,
     &              3.5d+5,1d+7,te)
            elseif(te.ge.1.d7) then
               tehigh=1.d7
               hectrec = fkingdon(1.96d-2,0.43d0,-1.09d0,-1.23d-2,8d+4,
     &              1d+7,tehigh)
            endif

         endif

c         write(6,*)'c v rec te,hectrec ',ion,te,hectrec

      elseif(iel.eq.4) then

         if(ion.eq.3) then

c N III + He I -> N II + He II

            if(te.le.1.d3) then
               telow=1.d3
               hectrec = fkingdon(4.84d-1,0.92d0,2.37d0,-1.02d1,1.d+3,
     &              4.d+4,telow)
            elseif(te.gt.1.d3.and.te.lt.4.d4) then
               hectrec = fkingdon(4.84d-1,0.92d0,2.37d0,-1.02d1,1.d+3,
     &              4.d+4,te)
            elseif(te.ge.4.d4.and.te.lt.1.d7) then
               hectrec = fkingdon(3.17d0,0.20d0,-0.72d0,-4.81d-2,4.d+4,
     &              1.d+7,te)
            elseif(te.ge.1.d7) then
               tehigh=1.d7
               hectrec = fkingdon(3.17d0,0.20d0,-0.72d0,-4.81d-2,4.d+4,
     &              1d+7,tehigh)
            endif

         elseif(ion.eq.4) then

            hectrec = HECTRecom(ion,7,te)


         elseif(ion.eq.5) then

c N V + He I -> N IV + He II

            if(te.le.1.d3) then
               telow=1.d3
               hectrec = fkingdon(1.26d-2,1.55d0,1.12d1,-7.82d0,1d+3,
     &              9d+4,telow)
            elseif(te.gt.1.d3.and.te.lt.9.d4) then
               hectrec = fkingdon(1.26d-2,1.55d0,1.12d1,-7.82d0,1d+3,
     &              9d+4,te)
            elseif(te.ge.9.d4.and.te.lt.1.d7) then
               hectrec = fkingdon(3.75d-1,0.54d0,-0.82d0,-2.07d-2,9.d+4,
     &              1.d+7,te)
            elseif(te.ge.1.d7) then
               tehigh=1.d7
               hectrec = fkingdon(3.75d-1,0.54d0,-0.82d0,-2.07d-2,9.d+4,
     &              1d+7,tehigh)
            endif

         endif

      elseif(iel.eq.5) then

         if(ion.eq.3) then

c O III + He I -> O II + He II

            if(te.le.1.d3) then
               telow=1.d3
               hectrec = fkingdon(7.10d-3,2.60d0,8.99d0,-0.78d0,1d+3,
     &              5d4,telow)
            elseif(te.gt.1.d3.and.te.lt.5.d4) then
               hectrec = fkingdon(7.10d-3,2.60d0,8.99d0,-0.78d0,1d+3,
     &              5d4,te)
            elseif(te.ge.5.d4.and.te.lt.1.d7) then
               hectrec = fkingdon(6.21d-1,0.53d0,-0.66d0,-2.22d-2,5d+4,
     &              1.d+7,te)
            elseif(te.ge.1.d7) then
               tehigh=1.d7
               hectrec = fkingdon(6.21d-1,0.53d0,-0.66d0,-2.22d-2,5d+4,
     &              1.d+7,tehigh)
            endif

         elseif(ion.eq.4) then
            
            hectrec = HECTRecom(ion,8,te)

         elseif(ion.eq.5) then
            
            hectrec = HECTRecom(ion,8,te)

         endif

c         write(6,*)'o rec te,hectrec ',ion,te,hectrec


      elseif(iel.eq.6) then

         if(ion.eq.3) then

c Ne III + He I -> Ne II + He II

            if(te.le.1.d3) then
               telow=1.d3
               hectrec = fkingdon(1.00D-5,0.d0,0.d0,0.d0,1d+3,5d+3,
     &              telow)
            elseif(te.gt.1.d3.and.te.lt.5.d3) then
               hectrec = fkingdon(1.00D-5,0.d0,0.d0,0.d0,1d+3,5d+3,te)
            elseif(te.gt.5.d3.and.te.lt.8.d3) then
               hectrec = fkingdon(8.48d-3,3.35d0,-1.92d0,-1.50d0,5d+3,
     &              8d+3,te)
            elseif(te.gt.8.d3.and.te.lt.1.d7) then
               hectrec = fkingdon(2.52d-2,0.14d0,-1.99d0,-0.91d0,8d+3,
     &              1d+7,te)
            elseif(te.ge.1.d7) then
               tehigh=1.d7
               hectrec = fkingdon(2.52d-2,0.14d0,-1.99d0,-0.91d0,8d+3,
     &              1d+7,tehigh)
            endif

         elseif(ion.eq.4) then

c Ne IV + He I -> Ne III + He II

            if(te.le.1.d3) then
               telow=1.d3
               hectrec = fkingdon(1.00D-5,0.d0,0.d0,0.d0,1d+3,2.5d+4,
     &              telow)
            elseif(te.gt.1.d3.and.te.lt.2.5d4) then
               hectrec = fkingdon(1.00D-5,0.d0,0.d0,0.d0,1d+3,2.5d+4,te)
            elseif(te.gt.2.5d4.and.te.lt.9.5d4) then
               hectrec = fkingdon(1.34d-4,2.33d0,-2.55d0,-0.37d0,2.5d+4,
     &              9.5d+4,te)
            elseif(te.gt.9.5d4.and.te.lt.1.d7) then
               hectrec = fkingdon(1.00d-1,0.24d0,-1.09d0,-2.47d-2,
     &              9.5d+4,1d+7,te)
            elseif(te.ge.1.d7) then
               tehigh=1.d7
               hectrec = fkingdon(1.00d-1,0.24d0,-1.09d0,-2.47d-2,
     &              9.5d+4,1d+7,tehigh)
            endif

         elseif(ion.eq.5) then

c Ne V + He I -> Ne IV + He II

            if(te.le.1.d3) then
               telow=1.d3
               hectrec = fkingdon(1.77d0,0.14d0,4.88d-2,-3.35d0,1d+3,
     &              5d+4,telow)
            elseif(te.gt.1.d3.and.te.lt.5.d4) then
               hectrec = fkingdon(1.77d0,0.14d0,4.88d-2,-3.35d0,
     &              1d+3,5d+4,te)
            elseif(te.gt.5.d4.and.te.lt.1.d7) then
               hectrec = fkingdon(2.67d-1,0.54d0,0.91d0,-1.88d-2,5d+4,
     &              1d+7,te)
            elseif(te.ge.1.d7) then
               tehigh=1.d7
               hectrec = fkingdon(2.67d-1,0.54d0,0.91d0,-1.88d-2,5d+4,
     &              1d+7,tehigh)
            endif

         endif

      elseif(iel.eq.7) then

         if(ion.eq.3) then

c Na III + He I -> Na II + He II

            if(te.le.1.d3) then
               telow=1.d3
               hectrec = fkingdon(5.98d-6,-0.29d0,1.34d-2,4.58d0,1d+3,
     &              2.d+4,telow)
            elseif(te.gt.1.d3.and.te.lt.2.d4) then
               hectrec = fkingdon(5.98d-6,-0.29d0,1.34d-2,4.58d0,1d+3,
     &              2.d+4,te)
            elseif(te.gt.2.d4.and.te.lt.1.d7) then
               hectrec = fkingdon(1.30d-2,0.17d0,-1.38d0,-0.20d0,2d+4,
     &              1d+7,te)
            elseif(te.ge.1.d7) then
               tehigh=1.d7
               hectrec = fkingdon(1.30d-2,0.17d0,-1.38d0,-0.20d0,2d+4,
     &              1d+7,tehigh)
            endif

         elseif(ion.eq.4) then

c Na IV + He I -> Na III + He II

            if(te.le.1.d3) then
               telow=1.d3
               hectrec = fkingdon(1.00D-5,0.d0,0.d0,0.d0,1d+3,5d+3,
     &              telow)
            elseif(te.gt.1.d3.and.te.lt.5d3) then
               hectrec = fkingdon(1.00D-5,0.d0,0.d0,0.d0,1d+3,5d+3,te)
            elseif(te.ge.5.d3.and.te.lt.2.d4) then
               hectrec = fkingdon(3.64d-6,8.00d0,1.44d3,-2.96d0,5d+3,
     &              2.d+4,te)
            elseif(te.gt.2.d4.and.te.lt.1.d7) then
               hectrec = fkingdon(6.45d-2,0.20d0,-0.93d0,-3.83d-2,2d+4,
     &              1d+7,te)
            elseif(te.ge.1.d7) then
               tehigh=1.d7
               hectrec = fkingdon(6.45d-2,0.20d0,-0.93d0,-3.83d-2,2d+4,
     &              1d+7,tehigh)
            endif
            
         elseif(ion.eq.5) then
            
            hectrec = HECTRecom(ion,11,te)
            
         endif

      elseif(iel.eq.8) then

         if(ion.eq.4) then

c Mg IV + He I -> Mg III + He II

            if(te.le.1.d3) then
               telow=1.d3
               hectrec = fkingdon(1.88d-2,0.13d0,0.83d0,-4.94d0,1d+3,
     &              7d+3,telow)
            elseif(te.gt.1.d3.and.te.lt.7.e3) then
               hectrec = fkingdon(1.88d-2,0.13d0,0.83d0,-4.94d0,1d+3,
     &              7d+3,te)
            elseif(te.gt.7.d3.and.te.lt.1.d7) then
               hectrec = fkingdon(3.41d-2,0.15d0,-0.45d0,-4.83d-2,7d+3,
     &              1d+7,te)
            elseif(te.ge.1.d7) then
               tehigh=1.d7
               hectrec = fkingdon(3.41d-2,0.15d0,-0.45d0,-4.83d-2,7d+3,
     &              1d+7,tehigh)
            endif
            
         elseif(ion.eq.5) then
            
            hectrec = HECTRecom(ion,11,te)
            
         endif

      elseif(iel.eq.9) then

c Al

         if(ion.eq.4.or.ion.eq.5) then

            hectrec = HECTRecom(ion,13,te)

         endif

      elseif(iel.eq.10) then

c Si

        if(ion.eq.4.or.ion.eq.5) then

            hectrec = HECTRecom(ion,14,te)

         endif

      elseif(iel.eq.11) then

c S

         if(ion.eq.4.or.ion.eq.5) then

            hectrec = HECTRecom(ion,16,te)

         endif

      elseif(iel.eq.12) then

c Ar
         if(ion.ge.3.and.ion.le.5) then

            hectrec = HECTRecom(ion,18,te)

         endif

      elseif(iel.eq.13) then

         if(ion.eq.4) then

c Ca IV + He I -> Ca III + He II

            if(te.le.1.d3) then
               telow=1.d3
               hectrec = fkingdon(1.00D-5,0.d0,0.d0,0.d0,1d+3,3d+4,
     &              telow)
            elseif(te.gt.1.d3.and.te.lt.3.d4) then
               hectrec = fkingdon(1.00D-5,0.d0,0.d0,0.d0,1d+3,3d+4,te)
            elseif(te.gt.3.d4.and.te.lt.8.d4) then
               hectrec = fkingdon(2.61d-8,4.69d0,1.52d6,-4.65d0,3d+4,
     &              8d+4,te)
            elseif(te.gt.8.d4.and.te.lt.1.d7) then
               hectrec = fkingdon(1.96d-2,0.43d0,-1.09d0,-1.23d-2,8d+4,
     &              1d+7,te)
            elseif(te.ge.1.d7) then
               tehigh=1.d7
               hectrec = fkingdon(1.96d-2,0.43d0,-1.09d0,-1.23d-2,8d+4,
     &              1d+7,tehigh)
            endif


         elseif(ion.eq.5) then
            
            hectrec = HECTRecom(ion,20,te)
            
         endif

      elseif(iel.eq.14) then

         if(ion.eq.4.or.ion.eq.5) then

            hectrec = HECTRecom(ion,26,te)

         endif

      endif

      return
      end



      real*8 function fkingdon(a,b,c,d,tlow,thigh,te)
      implicit real*8(a-h,o-z)
      t4=te/1.e4
c      write(6,9)a,b,c,d,tlow,thigh,te
      if(te.ge.tlow.and.te.le.thigh) then

         rate = a*t4**b*(1 + c*exp(d*t4))

      else
         write(6,*)' te out of range ',a,b,c,d,tlow,thigh,te
         rate=0.
c         stop
      endif
      fkingdon = 1.e-9*rate
c      write(6,9)a,b,c,d,tlow,thigh,te
 9    format(1pe12.3,10e12.3)
      return
      end






********************************************************************************
      block data ctdata
*
      implicit real*8(a-h,o-z)
      real*8 CTIonp
*     second dimension is ionization stage,
*     1=+0 for parent, etc
*     third dimension is atomic number of atom
      common/CTIon/ Ctionp(7,4,30)
      real*8 CTRecomb
*     second dimension is ionization stage,
*     1=+1 for parent, etc
*     third dimension is atomic number of atom
      common/CTRecomb/ CTRecomb(6,4,30)
*
*     local variables
      integer i
*
*     digital form of the fits to the charge transfer
*     ionization rate coefficients 
*
*     Note: First parameter is in units of 1e-9!
*     Note: Seventh parameter is in units of 1e4 K
*     ionization
      data (Ctionp(i,1,3),i=1,7)/2.84e-3,1.99,375.54,-54.07,1e2,1e4,0.0/
      data (Ctionp(i,2,3),i=1,7)/7*0./
      data (Ctionp(i,3,3),i=1,7)/7*0./
      data (Ctionp(i,1,4),i=1,7)/7*0./
      data (Ctionp(i,2,4),i=1,7)/7*0./
      data (Ctionp(i,3,4),i=1,7)/7*0./
      data (Ctionp(i,1,5),i=1,7)/7*0./
      data (Ctionp(i,2,5),i=1,7)/7*0./
      data (Ctionp(i,3,5),i=1,7)/7*0./
      data (Ctionp(i,1,6),i=1,7)/1.07e-6,3.15,176.43,-4.29,1e3,1e5,0.0/
      data (Ctionp(i,2,6),i=1,7)/7*0./
      data (Ctionp(i,3,6),i=1,7)/7*0./
      data (Ctionp(i,1,7),i=1,7)/4.55e-3,-0.29,-0.92,-8.38,1e2,5e4,
     &     1.086/
      data (Ctionp(i,2,7),i=1,7)/7*0./
      data (Ctionp(i,3,7),i=1,7)/7*0./
      data (Ctionp(i,1,8),i=1,7)/7.40e-2,0.47,24.37,-0.74,1e1,1e4,
     &     0.023/
      data (Ctionp(i,2,8),i=1,7)/7*0./
      data (Ctionp(i,3,8),i=1,7)/7*0./
      data (Ctionp(i,1,9),i=1,7)/7*0./
      data (Ctionp(i,2,9),i=1,7)/7*0./
      data (Ctionp(i,3,9),i=1,7)/7*0./
      data (Ctionp(i,1,10),i=1,7)/7*0./
      data (Ctionp(i,2,10),i=1,7)/7*0./
      data (Ctionp(i,3,10),i=1,7)/7*0./
      data (Ctionp(i,1,11),i=1,7)/3.34e-6,9.31,2632.31,-3.04,1e3,2e4,
     &     0.0/
      data (Ctionp(i,2,11),i=1,7)/7*0./
      data (Ctionp(i,3,11),i=1,7)/7*0./
      data (Ctionp(i,1,12),i=1,7)/9.76e-3,3.14,55.54,-1.12,5e3,3e4,0.0/
      data (Ctionp(i,2,12),i=1,7)/7.60e-5,0.00,-1.97,-4.32,1e4,3e5,
     &     1.670/
      data (Ctionp(i,3,12),i=1,7)/7*0./
      data (Ctionp(i,1,13),i=1,7)/7*0./
      data (Ctionp(i,2,13),i=1,7)/7*0./
      data (Ctionp(i,3,13),i=1,7)/7*0./
      data (Ctionp(i,1,14),i=1,7)/0.92,1.15,0.80,-0.24,1e3,2e5,0.0/
      data (Ctionp(i,2,14),i=1,7)/2.26,7.36e-2,-0.43,-0.11,2e3,1e5,
     1 3.031/
      data (Ctionp(i,3,14),i=1,7)/7*0./
      data (Ctionp(i,1,15),i=1,7)/7*0./
      data (Ctionp(i,2,15),i=1,7)/7*0./
      data (Ctionp(i,3,15),i=1,7)/7*0./
      data (Ctionp(i,1,16),i=1,7)/1.00e-5,0.00,0.00,0.00,1e3,1e4,0.0/
      data (Ctionp(i,2,16),i=1,7)/7*0./
      data (Ctionp(i,3,16),i=1,7)/7*0./
      data (Ctionp(i,1,17),i=1,7)/7*0./
      data (Ctionp(i,2,17),i=1,7)/7*0./
      data (Ctionp(i,3,17),i=1,7)/7*0./
      data (Ctionp(i,1,18),i=1,7)/7*0./
      data (Ctionp(i,2,18),i=1,7)/7*0./
      data (Ctionp(i,3,18),i=1,7)/7*0./
      data (Ctionp(i,1,19),i=1,7)/7*0./
      data (Ctionp(i,2,19),i=1,7)/7*0./
      data (Ctionp(i,3,19),i=1,7)/7*0./
      data (Ctionp(i,1,20),i=1,7)/7*0./
      data (Ctionp(i,2,20),i=1,7)/7*0./
      data (Ctionp(i,3,20),i=1,7)/7*0./
      data (Ctionp(i,1,21),i=1,7)/7*0./
      data (Ctionp(i,2,21),i=1,7)/7*0./
      data (Ctionp(i,3,21),i=1,7)/7*0./
      data (Ctionp(i,1,22),i=1,7)/7*0./
      data (Ctionp(i,2,22),i=1,7)/7*0./
      data (Ctionp(i,3,22),i=1,7)/7*0./
      data (Ctionp(i,1,23),i=1,7)/7*0./
      data (Ctionp(i,2,23),i=1,7)/7*0./
      data (Ctionp(i,3,23),i=1,7)/7*0./
      data (Ctionp(i,1,24),i=1,7)/7*0./
      data (Ctionp(i,2,24),i=1,7)/4.39,0.61,-0.89,-3.56,1e3,3e4,3.349/
      data (Ctionp(i,3,24),i=1,7)/7*0./
      data (Ctionp(i,1,25),i=1,7)/7*0./
      data (Ctionp(i,2,25),i=1,7)/2.83e-1,6.80e-3,6.44e-2,-9.70,1e3,3e4,
     1 2.368/
      data (Ctionp(i,3,25),i=1,7)/7*0./
      data (Ctionp(i,1,26),i=1,7)/7*0./
      data (Ctionp(i,2,26),i=1,7)/2.10,7.72e-2,-0.41,-7.31,1e4,1e5,
     &     3.005/
      data (Ctionp(i,3,26),i=1,7)/7*0./
      data (Ctionp(i,1,27),i=1,7)/7*0./
      data (Ctionp(i,2,27),i=1,7)/1.20e-2,3.49,24.41,-1.26,1e3,3e4,
     &     4.044/
      data (Ctionp(i,3,27),i=1,7)/7*0./
      data (Ctionp(i,1,28),i=1,7)/7*0./
      data (Ctionp(i,2,28),i=1,7)/7*0./
      data (Ctionp(i,3,28),i=1,7)/7*0./
      data (Ctionp(i,1,29),i=1,7)/7*0./
      data (Ctionp(i,2,29),i=1,7)/7*0./
      data (Ctionp(i,3,29),i=1,7)/7*0./
      data (Ctionp(i,1,30),i=1,7)/7*0./
      data (Ctionp(i,2,30),i=1,7)/7*0./
      data (Ctionp(i,3,30),i=1,7)/7*0./
*
*     digital form of the fits to the charge transfer
*     recombination rate coefficients (total)
*
*     Note: First parameter is in units of 1e-9!
*     recombination
      data (CTRecomb(i,1,2),i=1,6)/7.47e-6,2.06,9.93,-3.89,6e3,1e5/
      data (CTRecomb(i,2,2),i=1,6)/1.00e-5,0.,0.,0.,1e3,1e7/
      data (CTRecomb(i,1,3),i=1,6)/6*0./
      data (CTRecomb(i,2,3),i=1,6)/1.26,0.96,3.02,-0.65,1e3,3e4/
      data (CTRecomb(i,3,3),i=1,6)/1.00e-5,0.,0.,0.,2e3,5e4/
      data (CTRecomb(i,1,4),i=1,6)/6*0./
      data (CTRecomb(i,2,4),i=1,6)/1.00e-5,0.,0.,0.,2e3,5e4/
      data (CTRecomb(i,3,4),i=1,6)/1.00e-5,0.,0.,0.,2e3,5e4/
      data (CTRecomb(i,4,4),i=1,6)/5.17,0.82,-0.69,-1.12,2e3,5e4/
      data (CTRecomb(i,1,5),i=1,6)/6*0./
      data (CTRecomb(i,2,5),i=1,6)/2.00e-2,0.,0.,0.,1e3,1e9/
      data (CTRecomb(i,3,5),i=1,6)/1.00e-5,0.,0.,0.,2e3,5e4/
      data (CTRecomb(i,4,5),i=1,6)/2.74,0.93,-0.61,-1.13,2e3,5e4/
      data (CTRecomb(i,1,6),i=1,6)/4.88e-7,3.25,-1.12,-0.21,5.5e3,1e5/
      data (CTRecomb(i,2,6),i=1,6)/1.67e-4,2.79,304.72,-4.07,5e3,5e4/
      data (CTRecomb(i,3,6),i=1,6)/3.25,0.21,0.19,-3.29,1e3,1e5/
      data (CTRecomb(i,4,6),i=1,6)/332.46,-0.11,-9.95e-1,-1.58e-3,1e1,
     1 1e5/
      data (CTRecomb(i,1,7),i=1,6)/1.01e-3,-0.29,-0.92,-8.38,1e2,5e4/
      data (CTRecomb(i,2,7),i=1,6)/3.05e-1,0.60,2.65,-0.93,1e3,1e5/
      data (CTRecomb(i,3,7),i=1,6)/4.54,0.57,-0.65,-0.89,1e1,1e5/
      data (CTRecomb(i,4,7),i=1,6)/2.95,0.55,-0.39,-1.07,1e3,1e6/
      data (CTRecomb(i,1,8),i=1,6)/1.04,3.15e-2,-0.61,-9.73,1e1,1e4/
      data (CTRecomb(i,2,8),i=1,6)/1.04,0.27,2.02,-5.92,1e2,1e5/
      data (CTRecomb(i,3,8),i=1,6)/3.98,0.26,0.56,-2.62,1e3,5e4/
      data (CTRecomb(i,4,8),i=1,6)/2.52e-1,0.63,2.08,-4.16,1e3,3e4/
      data (CTRecomb(i,1,9),i=1,6)/6*0./
      data (CTRecomb(i,2,9),i=1,6)/1.00e-5,0.,0.,0.,2e3,5e4/
      data (CTRecomb(i,3,9),i=1,6)/9.86,0.29,-0.21,-1.15,2e3,5e4/
      data (CTRecomb(i,4,9),i=1,6)/7.15e-1,1.21,-0.70,-0.85,2e3,5e4/
      data (CTRecomb(i,1,10),i=1,6)/6*0./
      data (CTRecomb(i,2,10),i=1,6)/1.00e-5,0.,0.,0.,5e3,5e4/
      data (CTRecomb(i,3,10),i=1,6)/14.73,4.52e-2,-0.84,-0.31,5e3,5e4/
      data (CTRecomb(i,4,10),i=1,6)/6.47,0.54,3.59,-5.22,1e3,3e4/
      data (CTRecomb(i,1,11),i=1,6)/6*0./
      data (CTRecomb(i,2,11),i=1,6)/1.00e-5,0.,0.,0.,2e3,5e4/
      data (CTRecomb(i,3,11),i=1,6)/1.33,1.15,1.20,-0.32,2e3,5e4/
      data (CTRecomb(i,4,11),i=1,6)/1.01e-1,1.34,10.05,-6.41,2e3,5e4/
      data (CTRecomb(i,1,12),i=1,6)/6*0./
      data (CTRecomb(i,2,12),i=1,6)/8.58e-5,2.49e-3,2.93e-2,-4.33,1e3,
     1 3e4/
      data (CTRecomb(i,3,12),i=1,6)/6.49,0.53,2.82,-7.63,1e3,3e4/
      data (CTRecomb(i,4,12),i=1,6)/6.36,0.55,3.86,-5.19,1e3,3e4/
      data (CTRecomb(i,1,13),i=1,6)/6*0./
      data (CTRecomb(i,2,13),i=1,6)/1.00e-5,0.,0.,0.,1e3,3e4/
      data (CTRecomb(i,3,13),i=1,6)/7.11e-5,4.12,1.72e4,-22.24,1e3,3e4/
      data (CTRecomb(i,4,13),i=1,6)/7.52e-1,0.77,6.24,-5.67,1e3,3e4/
      data (CTRecomb(i,1,14),i=1,6)/6*0./
      data (CTRecomb(i,2,14),i=1,6)/6.77,7.36e-2,-0.43,-0.11,5e2,1e5/
      data (CTRecomb(i,3,14),i=1,6)/4.90e-1,-8.74e-2,-0.36,-0.79,1e3,
     1 3e4/
      data (CTRecomb(i,4,14),i=1,6)/7.58,0.37,1.06,-4.09,1e3,5e4/
      data (CTRecomb(i,1,15),i=1,6)/6*0./
      data (CTRecomb(i,2,15),i=1,6)/1.74e-4,3.84,36.06,-0.97,1e3,3e4/
      data (CTRecomb(i,3,15),i=1,6)/9.46e-2,-5.58e-2,0.77,-6.43,1e3,3e4/
      data (CTRecomb(i,4,15),i=1,6)/5.37,0.47,2.21,-8.52,1e3,3e4/
      data (CTRecomb(i,1,16),i=1,6)/3.82e-7,11.10,2.57e4,-8.22,1e3,1e4/
      data (CTRecomb(i,2,16),i=1,6)/1.00e-5,0.,0.,0.,1e3,3e4/
      data (CTRecomb(i,3,16),i=1,6)/2.29,4.02e-2,1.59,-6.06,1e3,3e4/
      data (CTRecomb(i,4,16),i=1,6)/6.44,0.13,2.69,-5.69,1e3,3e4/
      data (CTRecomb(i,1,17),i=1,6)/6*0./
      data (CTRecomb(i,2,17),i=1,6)/1.00e-5,0.,0.,0.,1e3,3e4/
      data (CTRecomb(i,3,17),i=1,6)/1.88,0.32,1.77,-5.70,1e3,3e4/
      data (CTRecomb(i,4,17),i=1,6)/7.27,0.29,1.04,-10.14,1e3,3e4/
      data (CTRecomb(i,1,18),i=1,6)/6*0./
      data (CTRecomb(i,2,18),i=1,6)/1.00e-5,0.,0.,0.,1e3,3e4/
      data (CTRecomb(i,3,18),i=1,6)/4.57,0.27,-0.18,-1.57,1e3,3e4/
      data (CTRecomb(i,4,18),i=1,6)/6.37,0.85,10.21,-6.22,1e3,3e4/
      data (CTRecomb(i,1,19),i=1,6)/6*0./
      data (CTRecomb(i,2,19),i=1,6)/1.00e-5,0.,0.,0.,1e3,3e4/
      data (CTRecomb(i,3,19),i=1,6)/4.76,0.44,-0.56,-0.88,1e3,3e4/
      data (CTRecomb(i,4,19),i=1,6)/1.00e-5,0.,0.,0.,1e3,3e4/
      data (CTRecomb(i,1,20),i=1,6)/6*0./
      data (CTRecomb(i,2,20),i=1,6)/0.,0.,0.,0.,1e1,1e9/
      data (CTRecomb(i,3,20),i=1,6)/3.17e-2,2.12,12.06,-0.40,1e3,3e4/
      data (CTRecomb(i,4,20),i=1,6)/2.68,0.69,-0.68,-4.47,1e3,3e4/
      data (CTRecomb(i,1,21),i=1,6)/6*0./
      data (CTRecomb(i,2,21),i=1,6)/0.,0.,0.,0.,1e1,1e9/
      data (CTRecomb(i,3,21),i=1,6)/7.22e-3,2.34,411.50,-13.24,1e3,3e4/
      data (CTRecomb(i,4,21),i=1,6)/1.20e-1,1.48,4.00,-9.33,1e3,3e4/
      data (CTRecomb(i,1,22),i=1,6)/6*0./
      data (CTRecomb(i,2,22),i=1,6)/0.,0.,0.,0.,1e1,1e9/
      data (CTRecomb(i,3,22),i=1,6)/6.34e-1,6.87e-3,0.18,-8.04,1e3,3e4/
      data (CTRecomb(i,4,22),i=1,6)/4.37e-3,1.25,40.02,-8.05,1e3,3e4/
      data (CTRecomb(i,1,23),i=1,6)/6*0./
      data (CTRecomb(i,2,23),i=1,6)/1.00e-5,0.,0.,0.,1e3,3e4/
      data (CTRecomb(i,3,23),i=1,6)/5.12,-2.18e-2,-0.24,-0.83,1e3,3e4/
      data (CTRecomb(i,4,23),i=1,6)/1.96e-1,-8.53e-3,0.28,-6.46,1e3,3e4/
      data (CTRecomb(i,1,24),i=1,6)/6*0./
      data (CTRecomb(i,2,24),i=1,6)/5.27e-1,0.61,-0.89,-3.56,1e3,3e4/
      data (CTRecomb(i,3,24),i=1,6)/10.90,0.24,0.26,-11.94,1e3,3e4/
      data (CTRecomb(i,4,24),i=1,6)/1.18,0.20,0.77,-7.09,1e3,3e4/
      data (CTRecomb(i,1,25),i=1,6)/6*0./
      data (CTRecomb(i,2,25),i=1,6)/1.65e-1,6.80e-3,6.44e-2,-9.70,1e3,
     1 3e4/
      data (CTRecomb(i,3,25),i=1,6)/14.20,0.34,-0.41,-1.19,1e3,3e4/
      data (CTRecomb(i,4,25),i=1,6)/4.43e-1,0.91,10.76,-7.49,1e3,3e4/
      data (CTRecomb(i,1,26),i=1,6)/6*0./
      data (CTRecomb(i,2,26),i=1,6)/1.26,7.72e-2,-0.41,-7.31,1e3,1e5/
      data (CTRecomb(i,3,26),i=1,6)/3.42,0.51,-2.06,-8.99,1e3,1e5/
      data (CTRecomb(i,4,26),i=1,6)/14.60,3.57e-2,-0.92,-0.37,1e3,3e4/
      data (CTRecomb(i,1,27),i=1,6)/6*0./
      data (CTRecomb(i,2,27),i=1,6)/5.30,0.24,-0.91,-0.47,1e3,3e4/
      data (CTRecomb(i,3,27),i=1,6)/3.26,0.87,2.85,-9.23,1e3,3e4/
      data (CTRecomb(i,4,27),i=1,6)/1.03,0.58,-0.89,-0.66,1e3,3e4/
      data (CTRecomb(i,1,28),i=1,6)/6*0./
      data (CTRecomb(i,2,28),i=1,6)/1.05,1.28,6.54,-1.81,1e3,1e5/
      data (CTRecomb(i,3,28),i=1,6)/9.73,0.35,0.90,-5.33,1e3,3e4/
      data (CTRecomb(i,4,28),i=1,6)/6.14,0.25,-0.91,-0.42,1e3,3e4/
      data (CTRecomb(i,1,29),i=1,6)/6*0./
      data (CTRecomb(i,2,29),i=1,6)/1.47e-3,3.51,23.91,-0.93,1e3,3e4/
      data (CTRecomb(i,3,29),i=1,6)/9.26,0.37,0.40,-10.73,1e3,3e4/
      data (CTRecomb(i,4,29),i=1,6)/11.59,0.20,0.80,-6.62,1e3,3e4/
      data (CTRecomb(i,1,30),i=1,6)/6*0./
      data (CTRecomb(i,2,30),i=1,6)/1.00e-5,0.,0.,0.,1e3,3e4/
      data (CTRecomb(i,3,30),i=1,6)/6.96e-4,4.24,26.06,-1.24,1e3,3e4/
      data (CTRecomb(i,4,30),i=1,6)/1.33e-2,1.56,-0.92,-1.20,1e3,3e4/
*
      end

      block data cthedata
      implicit real*8(a-h,o-z)
*
      real*8 CTHEIon
*     second dimension is ionization stage,
*     1=+0 for parent, etc
*     third dimension is atomic number of atom
      common/CTHEIon/ CTHEIon(7,4,30)
      real*8 CTHERecomb
*     second dimension is ionization stage,
*     1=+1 for parent, etc
*     third dimension is atomic number of atom
      common/CTHERecomb/ CTHERecomb(6,4,30)
*
*     local variables
      integer i
*
*     digital form of the fits to the charge transfer
*     ionization rate coefficients 
*
*     Note: First parameter is in units of 1e-9!
*     Note: Seventh parameter is in units of 1e4 K
*     ionization

      data (CTHERecomb(i,3,6),i=1,6)/1.12,0.42,-0.69,-0.34,1e+3,1e+7/  
c      data (CTHERecomb(i,4,6),i=1,6)/3.12e-7,-7.37e-2,3.50e1,2.40,
c     &     1e+3,1e+4/
c      data (CTHERecomb(i,4,6),i=1,6)/1.49e-5,2.73,5.93,-8.74e-2,1e+4,
c     &     3.5e+5/
c      data (CTHERecomb(i,4,6),i=1,6)/5.80e-2,0.73,-0.86,-9.60e-3,
c     &     3.5e+5,1e+7/
c      data (CTHERecomb(i,2,7),i=1,6)/4.84e-1,0.92,2.37,-1.02e1,1e+3,
c     &     4e+4/
c      data (CTHERecomb(i,2,7),i=1,6)/3.17,0.20,-0.72,-4.81e-2,4e+4,
c     &     1e+7/
      data (CTHERecomb(i,3,7),i=1,6)/2.05,0.23,-0.72,-0.19,1e+3,1e+7/
c      data (CTHERecomb(i,4,7),i=1,6)/1.26e-2,1.55,1.12e1,-7.82,1e+3,
c     &     9e+4/
c      data (CTHERecomb(i,4,7),i=1,6)/3.75e-1,0.54,-0.82,-2.07e-2,
c     &     9e+4,1e+7/
c      data (CTHERecomb(i,2,8),i=1,6)/7.10e-3,2.60,8.99,-0.78,1e+3,5e+4/
c      data (CTHERecomb(i,2,8),i=1,6)/6.21e-1,0.53,-0.66,-2.22e-2,5e+4,
c     &     1e+7/
      data (CTHERecomb(i,3,8),i=1,6)/1.12,0.42,-0.71,-1.98e-2,1e+3,1e+7/
      data (CTHERecomb(i,4,8),i=1,6)/9.97e-1,0.40,-0.46,-0.35,1e+3,1e+7/
c      data (CTHERecomb(i,2,10),i=1,6)/1.00e-5,0.00,0.00,0.00,1e+3,5e+3/
c      data (CTHERecomb(i,2,10),i=1,6)/8.48e-3,3.35,-1.92,-1.50,5e+3,
c     &     8e+3/
c      data (CTHERecomb(i,2,10),i=1,6)/2.52e-2,0.14,-1.99,-0.91,8e+3,
c     &     1e+7/
c      data (CTHERecomb(i,3,10),i=1,6)/1.00e-5,0.00, 0.00,0.00,1e+3,
c     &     2.5e+4/
c      data (CTHERecomb(i,3,10),i=1,6)/1.34e-4,2.33,-2.55,-0.37,2.5e+4,
c     &     9.5e+4/
c      data (CTHERecomb(i,3,10),i=1,6)/1.00e-1,0.24,-1.09,-2.47e-2,
c     &     9.5e+4,1e+7/
c      data (CTHERecomb(i,4,10),i=1,6)/1.77,0.14, 4.88e-2,-3.35,1e+3,
c     &     5e+4/
c      data (CTHERecomb(i,4,10),i=1,6)/2.67e-1,0.54, 0.91,-1.88e-2 5e+4,
c     &     1e+7/
c      data (CTHERecomb(i,2,11),i=1,6)/5.98e-6,-0.29, 1.34e-2,4.58,1e+3,
c     &     2e+4/
c      data (CTHERecomb(i,2,11),i=1,6)/1.30e-2,0.17,-1.38,-0.20,2e+4,
c     &     1e+7/
c      data (CTHERecomb(i,3,11),i=1,6)/1.00e-5,0.00,0.00,0.00,1e+3,5e+3/
c      data (CTHERecomb(i,3,11),i=1,6)/3.64e-6,8.00,1.44e3,-2.96,5e+3,
c     &     2e+4/
c      data (CTHERecomb(i,3,11),i=1,6)/6.45e-2,0.20,-0.93,-3.83e-2,2e+4,
c     &     1e+7/
      data (CTHERecomb(i,4,11),i=1,6)/3.97,0.13,-0.59,-0.13,1e+3,1e+7/
c      data (CTHERecomb(i,3,12),i=1,6)/1.88e-2,0.13,0.83,-4.94,1e+3,7e+3/
c      data (CTHERecomb(i,3,12),i=1,6)/3.41e-2,0.15,-0.45,-4.83e-2,7e+3,
c     &     1e+7/
      data (CTHERecomb(i,4,12),i=1,6)/1.37,0.21,-0.59,-5.94e-2,1e+3,
     &     1e+7/
      data (CTHERecomb(i,3,13),i=1,6)/1.53e-6,1.60e-2,-0.25,-9.71,1e+3,
     &     1e+7/
      data (CTHERecomb(i,4,13),i=1,6)/3.16,0.27,-0.64,-2.47e-2,1e+3,
     &     1e+7/
      data (CTHERecomb(i,3,14),i=1,6)/1.03,0.60,-0.61,-1.42,1e+2,1e+6/
      data (CTHERecomb(i,4,14),i=1,6)/5.75e-1,0.93,1.33,-0.29,1e+3,5e+5/
      data (CTHERecomb(i,3,16),i=1,6)/3.58,7.77e-3,-0.94,-0.30,1e+3,
     &     3.1e+4/
      data (CTHERecomb(i,4,16),i=1,6)/7.44e-4,0.34,3.74,-5.18,1e+3,
     &     3.1e+4/
      data (CTHERecomb(i,2,18),i=1,6)/1.30e-1,0.00,0.00,0.00,1e+3,
     &     3.1e+4/
      data (CTHERecomb(i,3,18),i=1,6)/1.00e-5,0.00,0.00,0.00,1e+3,
     &     3.1e+4/
      data (CTHERecomb(i,4,18),i=1,6)/9.72e-1,-0.14,1.05,-4.80,1e+3,
     &     1e+7/

c      data (CTHERecomb(i,3,20),i=1,6)/1.00e-5,0.00,0.00,0.00,1e+3,3e+4/
c      data (CTHERecomb(i,3,20),i=1,6)/2.61e-8,4.69,1.52e6,-4.65,3e+4,
c     &     8e+4/
c      data (CTHERecomb(i,3,20),i=1,6)/1.96e-2,0.43,-1.09,-1.23e-2,8e+4,
c     &     1e+7/
      data (CTHERecomb(i,4,20),i=1,6)/5.21,0.24,-0.51,-6.43e-2,1e+3,
     &     1e+7/
      data (CTHERecomb(i,3,26),i=1,6)/3.05e-2,0.00, 0.00,0.00,1e+3,1e+7/
      data (CTHERecomb(i,4,26),i=1,6)/1.38,0.43,1.12,-0.16,1e+3,1e+7/
      data (CTHERecomb(i,3,27),i=1,6)/7.10e-1,0.19,0.63,-4.91e-3,1e+3,
     &     1e+7/
c      data (CTHERecomb(i,4,27),i=1,6)/1.75,0.31, 1.25,1.14e2,1e+3,1e+6/
c      data (CTHERecomb(i,4,27),i=1,6)/7.14e3,-1.18,-1.25,-4.83e-3,1e+6,
c     &     1e+7/
c      data (CTHERecomb(i,3,28),i=1,6)/1.59,0.26,-0.63,-0.19,1e+3,1e+6/
c      data (CTHERecomb(i,3,28),i=1,6)/4.05e3,-1.16,-1.25,-5.34e-3,1e+6,
c     &     1e+7/
c      data (CTHERecomb(i,4,28),i=1,6)/2.96,0.29,-0.25,-0.13,1e+3,1e+7/

c!!! ion     data (CTHEIon(i,1,14),i=1,6)/1.30,0.00,0.00,0.00,1e+1,1e+4/

      end

      subroutine rec_fit_bad(iel,ion,te,alrr_di)
      implicit real*8 (a-h,o-z)
      common/recbad/ c(26,26,8),e(26,26,8),rr(26,26,6)
      aldi=0.
      do i=1,8
         aldi=aldi + c(iel,ion,i)*exp(-e(iel,ion,i)/te)/te**1.5
      enddo
c      write(6,91)(c(iel,ion,i),i=1,8)
c      write(6,91)(e(iel,ion,i),i=1,8)

      a=rr(iel,ion,1)
      b=rr(iel,ion,2)
      t0=rr(iel,ion,3)
      t1=rr(iel,ion,4)
      c0=rr(iel,ion,5)
      t2=rr(iel,ion,6)
c      write(6,91)a,b,t0,t1,c0,t2
 91   format(1pe12.4,10e12.4)
c      stop
      b=b+c0*exp(-t2/te)

      alrr=a*sqrt(t0/te)/((1+sqrt(te/t0))**(1.-b)*
     &     (1.+sqrt(te/t1))**(1.+b))

      alrr_di=alrr+aldi

      return
      end
      

      block data rec_fit
      implicit real*8 (a-h,o-z)
      common/recbad/ c(26,26,8),e(26,26,8),rr(26,26,6)
c Al-like from Abdel-Naby (2012) 

      data (c(14,2,i),i=1,8)/3.463d-8,1.9347d-7,1.652d-7,6.273d-6,5.465d-5,4.103d-3,2*0.d0/
      data (e(14,2,i),i=1,8)/2.411d+1,1.289d+2,4.201d+2,1.009d+4,5.042d+4,1.291d+5,2*0.d0/

      data (c(16,4,i),i=1,8)/5.817d-7,1.391d-6,1.123d-5,1.521d-4,1.875d-3,2.097d-2,2*0.d0/
      data (e(16,4,i),i=1,8)/3.628d+2,1.058d+3,7.160d+3,3.260d+4,1.235d+5,2.070d+5,2*0.d0/

      data (c(18,6,i),i=1,8)/1.136d-5,1.302d-4,8.017d-4,4.318d-2,7.582d-3,7.582d-3,2*0.d0/
      data (e(18,6,i),i=1,8)/7.953d+2,7.596d+3,3.505d+4,2.299d+5,4.576d+5,4.576d+5,2*0.d0/
      
      data (rr(14,2,i),i=1,6)/3.501d-11,0.6344d0,1.410d+1,4.343d+7,0.2248d0,5.831d+4/
      data (rr(16,4,i),i=1,6)/2.664d-10,0.6896d0,2.107d+1,2.028d+7,0.0840d0,6.752d+4/
      data (rr(18,6,i),i=1,6)/4.137d-10,0.7064d0,6.853d+1,2.108d+7,0.0316d0,1.020d+5/


      data (c(20, 8,i),i=1,8)/2.978d-4,7.585d-4,5.973d-3,8.950d-2,2.076d-3,2.078d-3,0.d0,0.d0/
      data (c(26,14,i),i=1,8)/2.305d-3,1.072d-2,3.512d-2,2.105d-1,3.622d-2,0.0d0,0.d0,0.d0/
      data (e(20, 8,i),i=1,8)/3.540d+3,1.530d+4,7.536d+4,3.313d+5,8.562d+5,8.562d+5,0.d0,0.d0/
      data (e(26,14,i),i=1,8)/4.747d+3,1.877d+4,1.190d+5,5.090d+5,1.595d+6,0.d0,0.d0,0.d0/
      data (rr(20,8,i),i=1,6)/4.279d-10,0.6860d0,2.289d+2,2.459d+7,0.0262d0,9.582d+4/
      data (rr(26,14,i),i=1,6)/5.370d-10,0.6338d0,1.531d+3,4.632d+7,0.0278d0,9.070d+4/

c  Phosphoros-like from Bleda+ (2022)          
      data (c(16,2,i),i=1,8)/7.300d-8,2.577d-7,4.961d-8,9.520d-7,9.586d-7,6.849d-4,6.539d-4,0.d0/
      data (c(18,4,i),i=1,8)/1.688d-6,4.430d-6,2.567d-5,1.028d-4,1.413d-3,4.267d-2,0.d0,0.d0/
      data (e(16,2,i),i=1,8)/5.077d+2,6.007d+2,2.342d+3,7.269d+3,2.190d+4,1.483d+5,1.906d+5,0.d0/
      data (e(18,4,i),i=1,8)/1.013d+3,1.620d+3,8.418d+3,2.451d+4,1.525d+5,3.233d+5,0.d0,0.0d0/
      data (c(20,6,i),i=1,8)/9.288d-6,5.118d-5,4.784d-4,1.906d-3,1.046d-1,3.473d-3,0.d0,0.d0/
      data (c(26,12,i),i=1,8)/1.791d-3,5.627d-3,2.131d-2,2.097d-1,1.936d-1,0.d0,0.d0,0.d0/
      data (e(20, 6,i),i=1,8)/9.012d+2,3.985d+3,1.738d+4,8.190d+4,4.066d+5,1.009d+6,0.d0,0.d0/
      data (e(26,12,i),i=1,8)/1.463d+3,1.289d+4,7.234d+4,4.470d+5,1.042d+6,0.d0,0.d0,0.d0/
      data (rr(16,2,i),i=1,6)/7.438d-11,0.6530d0,3.138d0,1.944d+6,0.1730,1.675d+4/
      data (rr(18,4,i),i=1,6)/3.464d-10,0.6893d0,1.393d+1,1.187d+7,0.1555,4.430d+5/     
      data (rr(20,6,i),i=1,6)/4.817d-10,0.6852d0,5.279d+1,2.706d+7,0.0239,4.049d+5/
      data (rr(26,12,i),i=1,6)/3.456d-10,0.5956d0,1.795d+3,4.097d+7,0.0405,6.769d+4/

c     Si-like from Kaur (2018)

      data (c(16,3,i),i=1,8)/3.040d-07,4.393d-07,1.609d-06,4.980d-06,3.457d-05,8.617d-03,9.284d-04,0.d0/
      data (e(16,3,i),i=1,8)/5.016d+01,3.266d+02,3.102d+03,1.210d+04,4.969d+04,2.010d+05,2.575d+05,0.d0/
      data (c(18,5,i),i=1,8)/1.590d-05,1.636d-05,7.566d-05,3.805d-04,5.247d-03,3.272d-02,1.060d-04,0.d0/
      data (e(18,5,i),i=1,8)/2.879d+02,1.717d+03,9.917d+03,5.769d+04,2.178d+05,3.191d+05,1.250d+06,0.d0/
      data (c(20,7,i),i=1,8)/4.836E-04,3.208E-04,9.281E-04,5.307E-02,2.175E-02,0.d0,0.d0,0.d0/
      data (e(20,7,i),i=1,8)/3.497E+02,2.664E+03,3.433E+04,3.149E+05,6.358E+05,0.d0,0.d0,0.d0/
      data (c(26,13,i),i=1,8)/4.469E-03,8.538E-03,1.741E-02,1.630E-01,8.680E-02,0.d0,0.d0,0.d0/
      data (e(26,13,i),i=1,8)/2.462E+03,1.261E+04,9.330E+04,4.887E+05,1.312E+06,0.d0,0.d0,0.d0/

      data (rr(16,3,i),i=1,6)/2.478d-11,0.4642,3.294d+02,2.166d+07,0.3351d0,7.630d+05/
      data (rr(18,5,i),i=1,6)/3.939d-10,0.6607,3.207d+01,3.043d+07,0.0761d0,6.360d+05/
      data (rr(20,7,i),i=1,6)/1.427E-08,0.7285,3.790E-01,3.977E+07,0.d0,0.d0/
      data (rr(26,13,i),i=1,6)/1.984E-09,0.7101,1.158E+02,4.400E+07,0.d0,0.d0/



      end

            
      
