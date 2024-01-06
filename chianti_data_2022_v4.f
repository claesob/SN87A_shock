      subroutine popchianti_new(iel,ion,te,rltot)
      IMPLICIT none
      include 'PARAM'
      integer nl,nlp1,nfel,nlp,ik,IAGN,ionq,iel,ion,idep
      integer ieli,il,init,initca,initch,initfe,initfei,inith,inithe,ionnr
      integer iout,ipopch,ir,itcon,iu,n,nionel,nmaxh,np1,nq,nrc,nvar
      integer i,ierr,ilabcfx,im1,imiss,j,k,ioni,ipr,nups,nchi22
      integer NION,NP1H,NMAX,NMIN,itime,ideb
      real*8 C3,C33,dens,del,DRQ,EM,ESC
      real*8 FRH,SR,WOBS,e00,CS,CI,G,E,abxa,bb
      real*8 bol,bolca,bolfe,bolfei,bolh,bolhe,c1,cistand,csaha,cstand
      real*8 e_ion,estand,gstand,astand,rs,upsil,weh,wlstand,xelec
      real*8 told,toldca,toldfe,toldfei,toldh,toldhe,toldo,ts,zi
      real*8 DCDT,DCIDT,RECNET,RECNDT,RECCO,RTE
      real*8  PH,DRDT,RECT,PHET,RECO,XN,XN1,XN2,XN3
      real*8  xion,abn,time,dr_dv,c,be,xi,aa,betot
      real*8  elch,amu,te,rltot,a,cinx,den,denel
      real*8  tev,wl,de,xinc

      PARAMETER (NL=340,NLP1=NL+1)
      PARAMETER(NFEL=3000)
      parameter (nlp=30000)
      common/abl/abn(15)
      common/ionx/xion(md,14,27)
      common/ichia/ipopch
      COMMON/NION/IONQ
      COMMON/INI/INITH,INIT,INITCA,INITHE,INITFE,initfei
      COMMON/IND/IK
      COMMON/ELEC/DEL(MD)
      COMMON/NLEV/e00,NION,N,NP1H,NMAXh,NMIN
      COMMON/BOLD/TOLD,TOLDH,TOLDCA,TOLDO,TOLDHE,TOLDFE,toldfei,BB(NL),
     &       BOLH(NL),BOL(NL),BOLCA(NL),BOLHE(NL),BOLFE(NL),bolfei(nl)
      COMMON/RECAL/RECO(NL)
      COMMON/HPOP/XN(6,NL),XN1,XN2,XN3
      COMMON/A14/Cstand(NL,NL),CIstand(NL),Gstand(NL),Estand(NL),
     &     Astand(NL,NL),WLstand(NL,NL),DCDT(NL,NL),DCIDT(NL)
      parameter(nchi22=20)
      common/atdatch/e_ion(nchi22),e(nchi22,nl),g(nchi22,nl),wl(nchi22,nl,nl),a(nchi22,nl,nl),
     &     c(nchi22,nl,nl),ci(nchi22,nl),nmax(nchi22),ionnr(14,26)
      COMMON/A17/RECNET(NL),RECNDT(NL),RECCO(NL),RTE(NL),
     &     PH(NL),DRDT(NL),RECT(6,NL),PHET(NL)
      COMMON/EQUIH/FRH(2,401),WEH(401),SR(NFEL+nlp),WOBS(NL,NL)
      common/timecheck/time,itime
      common/initchianti/initch
      PARAMETER (nups=65)
      common/chianti_2022_cs/upsil(nchi22,nups,nups)
      
      DIMENSION XINC(NL),XI(NL)
c     &     A(NLP1,NLP1),BS(NL)
      DATA CSAHA/2.0708d-16/
      dimension nionel(14)
      data nionel/1,2,6,7,8,10,11,12,13,14,16,18,20,26/
      ionq=ion
      abxa = abn(iel)*xion(ik,iel,ion)
      if(initch==1) then
         te=1.e6
         do ir=1,20
            if(ir.eq.1) then
               ieli=3
               ioni=2
            elseif(ir==2) then
               ieli=5
               ioni=2
            elseif(ir==3) then
               ieli=8
               ioni=2
            elseif(ir==4) then
               ieli=14
               ioni=3
            elseif(ir==5) then
               ieli=4
               ioni=2
            elseif(ir==6) then
               ieli=14
               ioni=4
            elseif(ir==7) then
               ieli=14
               ioni=5
            elseif(ir==8) then
               ieli=14
               ioni=6
            elseif(ir==9) then
               ieli=14
               ioni=7
            elseif(ir==10) then
               ieli=14
               ioni=8
            elseif(ir==11) then
               ieli=14
               ioni=9
            elseif(ir==12) then
               ieli=4
               ioni=3
            elseif(ir==13) then
               ieli=4
               ioni=4
            elseif(ir==14) then
               ieli=3
               ioni=3
            elseif(ir==15) then
               ieli=11
               ioni=4
            elseif(ir==16) then
               ieli=11
               ioni=5
            elseif(ir==17) then
               ieli=13
               ioni=5
            elseif(ir==18) then
               ieli=13
               ioni=2                                             
            elseif(ir==19) then
               ieli=14
               ioni=1
            elseif(ir==20) then
               ieli=14
               ioni=2                                            
            endif

            call ATDAT_chianti_2023(initch,ieli,ioni,te)
         enddo
         initch=0
      endif
      n=nmax(ionnr(iel,ion))
      ipopch=0
      if(itime.le.3.or.abxa.gt.1.e-8) then
         ipopch=1
         ieli=iel
         ioni=ion

         call el_coll_rate_chianti(ionnr(ieli,ioni),ieli,ioni,
     &        nmax(ionnr(ieli,ioni)),te)

         NMIN=n
         NQ=n
         NP1=N+1
         NP1H=NP1
         NVAR=NP1
         NRC=NP1+1
         nmaxh=nmin
         xelec=del(ik)

c!!   no recombination or photoionization to these

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

c     put the data from chianti into the usual arrays
         
         do il=1,nmax(ionnr(iel,ion))
            estand(il)=e(ionnr(iel,ion),il)
            gstand(il)=g(ionnr(iel,ion),il)
         enddo
         do iu=1,nmax(ionnr(iel,ion))
            do il=1,nmax(ionnr(iel,ion))
               wlstand(iu,il)=wl(ionnr(iel,ion),iu,il)
               astand(iu,il)=a(ionnr(iel,ion),iu,il)
               cstand(iu,il)=c(ionnr(iel,ion),iu,il)
            enddo
         enddo

         CALL MULTISIMPq(iel,ion,ZI,Te,XI,rltot)

         BB(NP1)=XI(NP1)
         DENEL=DEL(IK)*DE(RS)
         DO J=1,N
            BB(J)=XI(J)/(G(ionnr(iel,ion),J)*BB(NP1)*DENEL*C1*
     &           EXP((E00-E(ionnr(iel,ion),J))/TEV)
     &           /TS**1.5)
         ENDDO
         ITCON=1
         
      else


         do k=1,401
            weh(k)=1.d-40
         enddo

         rltot = 1.d-40

      endif


      
      end
      
      SUBROUTINE ATDAT_chianti_2022(init,iel,ion,te)
      IMPLICIT REAL*8(A-H,O-Z)
      character*43 conf
      character*5 levlab
      character*5 lorb
      integer twosp1
      PARAMETER (NL=340,NLP1=NL+1)
      parameter(nchi22=20)
      real*8 jlev,eicm(nchi22,nl)
      common/atdatch/e_ion(nchi22),e(nchi22,nl),g(nchi22,nl),wl(nchi22,nl,nl),a(nchi22,nl,nl),
     &     c(nchi22,nl,nl),ci(nchi22,nl),nmax(nchi22),ionnr(14,26)
      COMMON/NBA/NBACK
      COMMON/NLEV/e00,NION,NLEV,NP1H,NMAXQ,NMIN
      common/timecheck/time,itime
      common/colls/u
      integer n1,i2sp1
      inew=0
      
      if(iel.eq.3.and.ion.eq.2) then
         open(86,file='./ATDAT/c_2.elvlc',status='old')
         open(87,file='./ATDAT/c_2.wgfa',status='old')
         ionnr(iel,ion)=1                 
         n=100
         nmax(ionnr(iel,ion))=53
         inew=1
      elseif(iel.eq.4.and.ion.eq.2) then
         open(86,file='./ATDAT/n_2.elvlc',status='old')
         open(87,file='./ATDAT/n_2.wgfa',status='old')
         ionnr(iel,ion)=5
         n=100
         nmax(ionnr(iel,ion))=12
         inew=1
      elseif(iel.eq.5.and.ion.eq.2) then
         open(86,file='./ATDAT/o_2_rev.elvlc',status='old')
         open(87,file='./ATDAT/o_2_rev.wgfa',status='old')
         ionnr(iel,ion)=2
         n=100
         nmax(ionnr(iel,ion))=35
         inew=1
      elseif(iel.eq.8.and.ion.eq.2) then
         open(86,file='./ATDAT/mg_2.elvlc',status='old')
         open(87,file='./ATDAT/mg_2.wgfa',status='old')
         ionnr(iel,ion)=3               
         n=100
         nmax(ionnr(iel,ion))=41
         inew=1
      elseif(iel.eq.11.and.ion.eq.4) then
         open(86,file='./ATDAT/s_4.elvlc',status='old')
         open(87,file='./ATDAT/s_4.wgfa',status='old')
         ionnr(iel,ion)=15
         n=100
         nmax(ionnr(iel,ion))=52
         inew=1
         e_ion(ionnr(iel,ion))=47.222
      elseif(iel.eq.11.and.ion.eq.5) then
         open(86,file='./ATDAT/s_5.elvlc',status='old')
         open(87,file='./ATDAT/s_5.wgfa',status='old')
         ionnr(iel,ion)=16
         n=100
         nmax(ionnr(iel,ion))=52
         inew=1
         e_ion(ionnr(iel,ion))=72.4945         
      elseif(iel.eq.13.and.ion.eq.2) then
         open(86,file='./ATDAT/ca_2.elvlc',status='old')
         open(87,file='./ATDAT/ca_2.wgfa',status='old')
         ionnr(iel,ion)=18             
         n=41
         nmax(ionnr(iel,ion))=41
         inew=1
         e_ion(ionnr(iel,ion))=11.871719 

      elseif(iel.eq.13.and.ion.eq.5) then
         open(86,file='./ATDAT/ca_5.elvlc',status='old')
         open(87,file='./ATDAT/ca_5.wgfa',status='old')
         ionnr(iel,ion)=17             
         n=100
c Truncated at 61 (ie 5p 1P. Extend?
         nmax(ionnr(iel,ion))=5
         inew=1
         e_ion(ionnr(iel,ion))=84.34

      elseif(iel.eq.14.and.ion.eq.3) then
         open(86,file='./ATDAT/fe_3.elvlc',status='old')
         open(87,file='./ATDAT/fe_3.wgfa',status='old')
         ionnr(iel,ion)=4               
         n=100
         nmax(ionnr(iel,ion))=55
         inew=1
      elseif(iel.eq.14.and.ion.eq.4) then
         open(86,file='./ATDAT/fe_4.elvlc',status='old')
         open(87,file='./ATDAT/fe_4.wgfa',status='old')
         ionnr(iel,ion)=6             
         n=100
         nmax(ionnr(iel,ion))=37
         inew=1
      elseif(iel.eq.14.and.ion.eq.5) then
         open(86,file='./ATDAT/fe_5.elvlc',status='old')
         open(87,file='./ATDAT/fe_5.wgfa',status='old')
         ionnr(iel,ion)=7             
         n=100
         nmax(ionnr(iel,ion))=34
         inew=1
      elseif(iel.eq.14.and.ion.eq.6) then
         open(86,file='./ATDAT/fe_6.elvlc',status='old')
         open(87,file='./ATDAT/fe_6.wgfa',status='old')
         ionnr(iel,ion)=8           
         n=100
         nmax(ionnr(iel,ion))=55
         inew=1
      elseif(iel.eq.14.and.ion.eq.7) then
         open(86,file='./ATDAT/fe_7.elvlc',status='old')
         open(87,file='./ATDAT/fe_7.wgfa',status='old')
         ionnr(iel,ion)=9             
         n=100
c Trncated at 17 since obs and theoretical energies not monotoneous!                  
         nmax(ionnr(iel,ion))=17
         inew=1
      elseif(iel.eq.14.and.ion.eq.8) then
         open(86,file='./ATDAT/fe_8.elvlc',status='old')
         open(87,file='./ATDAT/fe_8.wgfa',status='old')
         ionnr(iel,ion)=10             
         n=100
c Trncated at 18 since obs and theoretical energies not monotoneous!         
         nmax(ionnr(iel,ion))=18
         inew=1
      elseif(iel.eq.14.and.ion.eq.9) then
         open(86,file='./ATDAT/fe_9.elvlc',status='old')
         open(87,file='./ATDAT/fe_9.wgfa',status='old')
         ionnr(iel,ion)=11             
         n=100
c Truncated at 17 inlucing only  obs energies
         nmax(ionnr(iel,ion))=17
         inew=1
      elseif(iel.eq.4.and.ion.eq.3) then
         open(86,file='./ATDAT/n_3.elvlc',status='old')
         open(87,file='./ATDAT/n_3.wgfa',status='old')
         ionnr(iel,ion)=12             
         n=100
c Truncated at 53 (ie 3d 4D. Extend?
         nmax(ionnr(iel,ion))=53
         inew=1
      elseif(iel.eq.4.and.ion.eq.4) then
         open(86,file='./ATDAT/n_4.elvlc',status='old')
         open(87,file='./ATDAT/n_4.wgfa',status='old')
         ionnr(iel,ion)=13             
         n=100
c Truncated at 53 (ie 4 1F. Extend?
         nmax(ionnr(iel,ion))=60
         inew=1
      elseif(iel.eq.3.and.ion.eq.3) then
         open(86,file='./ATDAT/c_3.elvlc',status='old')
         open(87,file='./ATDAT/c_3.wgfa',status='old')
         ionnr(iel,ion)=14             
         n=100
c Truncated at 61 (ie 5p 1P. Extend?
         nmax(ionnr(iel,ion))=61
         inew=1
         
      endif
      
      if(inew==1) then
c         write(6,*)'iel,ion,nmax',iel,ion,nmax(ionnr(iel,ion))
         eijcm=0.      
         do i=1,nmax(ionnr(iel,ion))
            eijcm_old=eijcm
c     read(86,9646)iq,ch29,jlev,wn
            if(iel==8.and.ion==2) then
               read(86,99,err=33)n1,conf,jlev,eijcm,eijcmth
 99            format(i7,a43,f7.1,2f15.3)
            else
               read(86,991,err=33)n1,conf,levlab,twosp1,lorb,jlev,eijcm,
     &              eijcmth
 991           format(i7,a30,a5,i5,a5,f5.1,2f15.3)
            endif

            eijcm_old=eijcm
            if(n1.le.nmax(ionnr(iel,ion))) then
               if(eijcm.gt.0..and.n1.ne.0) then
                  eicm(ionnr(iel,ion),i)=eijcm
               elseif(eijcm.le.0..and.n1.ne.0) then
                  eicm(ionnr(iel,ion),i)=eijcmth
c     split degenerate levels by arbtrary 1 cm-1
                  if(eijcm.eq.eijcm_old.and.n1.gt.1) then
                     eicm(ionnr(iel,ion),i)=eicm(ionnr(iel,ion),i-1)+1.
                  endif

                  if(eijcm.eq.0..and.n1.gt.1) then
                     eicm(ionnr(iel,ion),i)=eicm(ionnr(iel,ion),i-1)+1.
                  endif
               endif
c     energy in eV
               e(ionnr(iel,ion),i)=eicm(ionnr(iel,ion),i)/8065.54445
c     stat weight
               g(ionnr(iel,ion),i)=2.*jlev+1.

c     if(iel.eq.8.and.ion.eq.2) then
c     wlq=12398.54/E(I)
               wlq=12398.41/e(ionnr(iel,ion),i)
               wn=1.e8/wlq
c     endif

            endif
         enddo
 33      continue
         close(86)
         do i=1,100000
            read(87,*,err=11,end=11)il,iu,wlq,gf,a21

            if(iu.le.nmax(ionnr(iel,ion))) then
               a(ionnr(iel,ion),iu,il)=a21
            endif
            
c     if(iel.eq.5.and.ion==2) then
            if(iu.le.nmax(ionnr(iel,ion))) then
               de=e(ionnr(iel,ion),iu)-e(ionnr(iel,ion),il)
               wll=12398.41/de
            endif
c     endif
         enddo

 11      continue
         close(87)
         init=1
         call upsilon_chianti(init,ionnr(iel,ion),iel,ion,
     &        nmax(ionnr(iel,ion)),te)
      endif
c      nelv=nmax(ionnr(iel,ion))
      return
      end

      SUBROUTINE ATDAT_chianti_2023(init,iel,ion,te)
      IMPLICIT REAL*8(A-H,O-Z)
      character*1 dum
      character*43 conf
      character*30 conf30
      character*5 levlab
      character*5 lorb
      integer twosp1,cfform
      character*6 lang
      character*120 clab

      PARAMETER (NL=340,NLP1=NL+1)
      parameter(nchi22=20)
      PARAMETER (nups=65)
      real*8 jlev,eicm(nchi22,nl),omega,glev
      common/atdatch/e_ion(nchi22),e(nchi22,nl),g(nchi22,nl),wl(nchi22,nl,nl),a(nchi22,nl,nl),
     &     c(nchi22,nl,nl),ci(nchi22,nl),nmax(nchi22),ionnr(14,26)
      COMMON/NBA/NBACK
      COMMON/NLEV/e00,NION,NLEV,NP1H,NMAXQ,NMIN
      common/timecheck/time,itime
      common/colls/u
      common/chianti_2022_cs/upsil(nchi22,nups,nups)
      integer n1,i2jp1(nl),i2sp1
      integer gu,gl
      inew=0
      cfform=0
      if(iel.eq.3.and.ion.eq.2) then
         open(86,file='./ATDAT/c_2.elvlc',status='old')
         open(87,file='./ATDAT/c_2.wgfa',status='old')
         ionnr(iel,ion)=1                 
         n=100
         nmax(ionnr(iel,ion))=53
         inew=1
      elseif(iel.eq.4.and.ion.eq.2) then
         open(86,file='./ATDAT/n_2.elvlc',status='old')
         open(87,file='./ATDAT/n_2.wgfa',status='old')
         ionnr(iel,ion)=5
         n=100
         nmax(ionnr(iel,ion))=12
         inew=1
      elseif(iel.eq.5.and.ion.eq.2) then
         open(86,file='./ATDAT/o_2_rev.elvlc',status='old')
         open(87,file='./ATDAT/o_2_rev.wgfa',status='old')
         ionnr(iel,ion)=2
         n=100
         nmax(ionnr(iel,ion))=35
         inew=1
      elseif(iel.eq.8.and.ion.eq.2) then
         open(86,file='./ATDAT/mg_2.elvlc',status='old')
         open(87,file='./ATDAT/mg_2.wgfa',status='old')
         ionnr(iel,ion)=3               
         n=100
         nmax(ionnr(iel,ion))=41
         inew=1
      elseif(iel.eq.11.and.ion.eq.4) then
         open(86,file='./ATDAT/s_4.elvlc',status='old')
         open(87,file='./ATDAT/s_4.wgfa',status='old')
         ionnr(iel,ion)=15
         n=100
         nmax(ionnr(iel,ion))=52
         inew=1
         e_ion(ionnr(iel,ion))=47.222
      elseif(iel.eq.11.and.ion.eq.5) then
         open(86,file='./ATDAT/s_5.elvlc',status='old')
         open(87,file='./ATDAT/s_5.wgfa',status='old')
         ionnr(iel,ion)=16
         n=100
         nmax(ionnr(iel,ion))=52
         inew=1
         e_ion(ionnr(iel,ion))=72.4945         
      elseif(iel.eq.13.and.ion.eq.2) then
         open(86,file='./ATDAT/ca_2.elvlc',status='old')
         open(87,file='./ATDAT/ca_2.wgfa',status='old')
         ionnr(iel,ion)=18             
         n=41
         nmax(ionnr(iel,ion))=41
         inew=1
         e_ion(ionnr(iel,ion))=11.871719 

      elseif(iel.eq.13.and.ion.eq.5) then
         open(86,file='./ATDAT/ca_5.elvlc',status='old')
         open(87,file='./ATDAT/ca_5.wgfa',status='old')
         ionnr(iel,ion)=17             
         n=100
c Truncated at 61 (ie 5p 1P. Extend?
         nmax(ionnr(iel,ion))=5
         inew=1
         e_ion(ionnr(iel,ion))=84.34
      elseif(iel.eq.14.and.ion.eq.1) then
         open(86,file='./ATDAT/fe_1.elvlc',status='old')
         open(87,file='./ATDAT/fe_1.wgfa',status='old')
         ionnr(iel,ion)=19               
         n=121
         nmax(ionnr(iel,ion))=121
c     Obs! Truncated at 65 to fit in ups array. Should b 121!!
         nmax(ionnr(iel,ion))=65
         inew=1         
         cfform=1 
      elseif(iel.eq.14.and.ion.eq.2) then
         open(86,file='./ATDAT/fe_2.elvlc',status='old')
         open(87,file='./ATDAT/fe_2.wgfa',status='old')
         ionnr(iel,ion)=20               
         n=191
         nmax(ionnr(iel,ion))=191
c     Obs! Truncated at 65 to fit in ups array. Should b 191!!
         nmax(ionnr(iel,ion))=65
         inew=1         
         cfform=1 
      elseif(iel.eq.14.and.ion.eq.3) then
         open(86,file='./ATDAT/fe_3.elvlc',status='old')
         open(87,file='./ATDAT/fe_3.wgfa',status='old')
         ionnr(iel,ion)=4               
         n=100
         nmax(ionnr(iel,ion))=55
         inew=1
      elseif(iel.eq.14.and.ion.eq.4) then
         open(86,file='./ATDAT/fe_4.elvlc',status='old')
         open(87,file='./ATDAT/fe_4.wgfa',status='old')
         ionnr(iel,ion)=6             
         n=100
         nmax(ionnr(iel,ion))=37
         inew=1
      elseif(iel.eq.14.and.ion.eq.5) then
         open(86,file='./ATDAT/fe_5.elvlc',status='old')
         open(87,file='./ATDAT/fe_5.wgfa',status='old')
         ionnr(iel,ion)=7             
         n=100
         nmax(ionnr(iel,ion))=34
         inew=1
      elseif(iel.eq.14.and.ion.eq.6) then
         open(86,file='./ATDAT/fe_6.elvlc',status='old')
         open(87,file='./ATDAT/fe_6.wgfa',status='old')
         ionnr(iel,ion)=8           
         n=100
         nmax(ionnr(iel,ion))=55
         inew=1
      elseif(iel.eq.14.and.ion.eq.7) then
         open(86,file='./ATDAT/fe_7.elvlc',status='old')
         open(87,file='./ATDAT/fe_7.wgfa',status='old')
         ionnr(iel,ion)=9             
         n=100
c Trncated at 17 since obs and theoretical energies not monotoneous!                  
         nmax(ionnr(iel,ion))=17
         inew=1
      elseif(iel.eq.14.and.ion.eq.8) then
         open(86,file='./ATDAT/fe_8.elvlc',status='old')
         open(87,file='./ATDAT/fe_8.wgfa',status='old')
         ionnr(iel,ion)=10             
         n=100
c Trncated at 18 since obs and theoretical energies not monotoneous!         
         nmax(ionnr(iel,ion))=18
         inew=1
      elseif(iel.eq.14.and.ion.eq.9) then
         open(86,file='./ATDAT/fe_9.elvlc',status='old')
         open(87,file='./ATDAT/fe_9.wgfa',status='old')
         ionnr(iel,ion)=11             
         n=100
c Truncated at 17 inlucing only  obs energies
         nmax(ionnr(iel,ion))=17
         inew=1
      elseif(iel.eq.4.and.ion.eq.3) then
         open(86,file='./ATDAT/n_3.elvlc',status='old')
         open(87,file='./ATDAT/n_3.wgfa',status='old')
         ionnr(iel,ion)=12             
         n=100
c Truncated at 53 (ie 3d 4D. Extend?
         nmax(ionnr(iel,ion))=53
         inew=1
      elseif(iel.eq.4.and.ion.eq.4) then
         open(86,file='./ATDAT/n_4.elvlc',status='old')
         open(87,file='./ATDAT/n_4.wgfa',status='old')
         ionnr(iel,ion)=13             
         n=100
c Truncated at 53 (ie 4 1F. Extend?
         nmax(ionnr(iel,ion))=60
         inew=1
      elseif(iel.eq.3.and.ion.eq.3) then
         open(86,file='./ATDAT/c_3.elvlc',status='old')
         open(87,file='./ATDAT/c_3.wgfa',status='old')
         ionnr(iel,ion)=14             
         n=100
c Truncated at 61 (ie 5p 1P. Extend?
         nmax(ionnr(iel,ion))=61
         inew=1
         
      endif
      
      if(inew==1) then
         if(cfform==0) then
            eijcm=0.      
            do i=1,nmax(ionnr(iel,ion))
               eijcm_old=eijcm
               if(iel==8.and.ion==2) then
                  read(86,99,err=33)n1,conf,jlev,eijcm,eijcmth
 99               format(i7,a43,f7.1,2f15.3)
               else
                  read(86,991,err=33)n1,conf,levlab,twosp1,lorb,jlev,eijcm,
     &                 eijcmth
 991              format(i7,a30,a5,i5,a5,f5.1,2f15.3)
               endif

               eijcm_old=eijcm
               if(n1.le.nmax(ionnr(iel,ion))) then
                  if(eijcm.gt.0..and.n1.ne.0) then
                     eicm(ionnr(iel,ion),i)=eijcm
                  elseif(eijcm.le.0..and.n1.ne.0) then
                     eicm(ionnr(iel,ion),i)=eijcmth
c     split degenerate levels by arbtrary 1 cm-1
                     if(eijcm.eq.eijcm_old.and.n1.gt.1) then
                        eicm(ionnr(iel,ion),i)=eicm(ionnr(iel,ion),i-1)+1.
                     endif

                     if(eijcm.eq.0..and.n1.gt.1) then
                        eicm(ionnr(iel,ion),i)=eicm(ionnr(iel,ion),i-1)+1.
                     endif
                  endif
c     energy in eV
                  e(ionnr(iel,ion),i)=eicm(ionnr(iel,ion),i)/8065.54445
c     stat weight
                  g(ionnr(iel,ion),i)=2.*jlev+1.

c     if(iel.eq.8.and.ion.eq.2) then
c     wlq=12398.54/E(I)
                  wlq=12398.41/e(ionnr(iel,ion),i)
                  wn=1.e8/wlq
c     endif

               endif
            enddo
 33         continue
         elseif(cfform==1) then
            if(iel==14.and.ion==1) then
               ntxt=2
            elseif(iel==14.and.ion==2) then
               ntxt=10
            endif
            do i=1,ntxt
               read(86,*,err=34)dum
            enddo
            do i=1,nmax(ionnr(iel,ion))
               read(86,*,err=34)n1,eijcm,glev,conf
               if(n1.le.nmax(ionnr(iel,ion))) then
                  eicm(ionnr(iel,ion),i)=eijcm
               endif
c     energy in eV
               e(ionnr(iel,ion),i)=eicm(ionnr(iel,ion),i)/8065.54445
c     stat weight
               g(ionnr(iel,ion),i)=glev
c     wlq=12398.54/E(I)
               wlq=12398.41/e(ionnr(iel,ion),i)
               wn=1.e8/wlq
            enddo
 34       continue
          endif         
         close(86)
         if(cfform==0) then
            do i=1,100000
               read(87,*,err=11,end=11)il,iu,wlq,gf,a21
               if(iu.le.nmax(ionnr(iel,ion))) then
                  a(ionnr(iel,ion),iu,il)=a21
               endif            
c     if(iel.eq.5.and.ion==2) then
               if(iu.le.nmax(ionnr(iel,ion))) then
                  de=e(ionnr(iel,ion),iu)-e(ionnr(iel,ion),il)
                  wll=12398.41/de
               endif
            enddo
         elseif(cfform==1) then
            nmax_lev=nmax(ionnr(iel,ion))
            do il=1,nmax_lev-1
               do iu=il+1,nmax_lev
                  upsil(ionnr(iel,ion),il,iu)=0.
               enddo
            enddo
            read(87,*)dum
            read(87,*)dum           
            do i=1,100000
               if(iel==14.and.ion==1) then
                  read(87,*,err=11,end=11)iu,il,wlq,a21,omega
               elseif(iel==14.and.ion==2) then
                  read(87,*,err=11,end=11)iu,il,dum,gu,dum,gl,wlq,a21,omega
               endif
               if(iu.le.nmax(ionnr(iel,ion))) then
                  a(ionnr(iel,ion),iu,il)=a21
c???  upsil(ionnr(iel,ion),il,iu)=omega
                  upsil(ionnr(iel,ion),il,iu)=omega
               endif            
               if(iu.le.nmax(ionnr(iel,ion))) then
                  de=e(ionnr(iel,ion),iu)-e(ionnr(iel,ion),il)
                  wll=12398.41/de
               endif               
            enddo
         endif
 11      continue
         close(87)
         init=1
         if(cfform==0) then
            call upsilon_chianti(init,ionnr(iel,ion),iel,ion,
     &           nmax(ionnr(iel,ion)),te)
         endif
      endif
      nelv=nmax(ionnr(iel,ion))
      return
      end
      
      subroutine el_coll_rate_chianti(ionnri,iel,ion,nmaxi,te)
      implicit none
      integer nl,nlp1,nups,ionnr,ionnri,nmaxi,nchi22
      PARAMETER (NL=340,NLP1=NL+1)
      real*8 upsil,c,ci,g,e,a,wl,dcdt,dcidt,te,tex,e_ion
      integer i,iu,il,iel,ion,nmax,init
c      COMMON/A14/C(NL,NL),CI(NL),G(NL),E(NL),A(NL,NL),
c     &     WL(NL,NL),DCDT(NL,NL),DCIDT(NL)
      parameter(nchi22=20)
      common/atdatch/e_ion(nchi22),e(nchi22,nl),g(nchi22,nl),wl(nchi22,nl,nl),a(nchi22,nl,nl),
     &     c(nchi22,nl,nl),ci(nchi22,nl),nmax(nchi22),ionnr(14,26)
      PARAMETER (nups=65)
      common/chianti_2022_cs/upsil(nchi22,nups,nups)

      init=0
      if(ionnri.ne.19.and.ionnri.ne.20) then
         call upsilon_chianti(init,ionnri,iel,ion,nmaxi,te)
      endif

      do iu=2,nmaxi
         do il=1,iu-1
            tex=1.602e-12*(e(ionnri,iu)-e(ionnri,il))/1.38e-16
c same as 8.63e-6 * .....            
            c(ionnri,il,iu)=2.1716e-8*(13.6068*1.602e-12/(1.38e-16*te))**0.5*
     &           upsil(ionnri,il,iu)*exp(-tex/te)/g(ionnri,il)
            c(ionnri,iu,il)=2.1716e-8*(13.6068*1.602e-12/(1.38e-16*te))**0.5*
     &           upsil(ionnri,il,iu)/g(ionnri,iu)
         enddo
      enddo

      return
      end      
      
      subroutine upsilon_chianti(init,ionnr,iel,ion,nmax_lev,te)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NMAX=500,nups=65)
      parameter(nchi22=20,nknot=20)
      REAL*8 yp1,ypn,x,ya,y2(nknot),ybis(nchi22,nups,nups,nknot)
      real*8 delta_e(nchi22,nups,nups),t_high(nchi22,nups,nups),
     &     c_val(nchi22,nups,nups)
      real*8 t_knot(nchi22,nups,nups,nknot),scups_knot(nchi22,nups,nups,nknot)
      real*8 tn(nknot),scup(nknot),t_knoti(nknot),scups_knoti(nknot)
      common/chianti_2022_cs/upsil(nchi22,nups,nups)
      INTEGER i,k,nti,t_typei
      integer nt(nchi22,nups,nups),t_type(nchi22,nups,nups),ilmax(nchi22),
     &     iumax(nchi22),iumaxil(nchi22,1000)
      save ybis,ilmax,iumax,t_type,delta_e,nt,t_high,c_val,t_knot,
     &     scups_knot,iumaxil
      if(init==1) then
         yp1=1.e31
         ypn=1.e31
         if(iel==3.and.ion==2) then
            open(11,file='./ATDAT/c_2.scups',status='old')
         elseif(iel==4.and.ion==2) then
            open(11,file='./ATDAT/n_2.scups',status='old')
         elseif(iel==5.and.ion==2) then
            open(11,file='./ATDAT/o_2_rev.upsdata',status='old')
         elseif(iel==8.and.ion==2) then
            open(11,file='./ATDAT/mg_2.scups',status='old')
         elseif(iel==14.and.ion==3) then
            open(11,file='./ATDAT/fe_3.scups',status='old')
         elseif(iel==14.and.ion==4) then
            open(11,file='./ATDAT/fe_4.scups',status='old')
         elseif(iel==14.and.ion==5) then
            open(11,file='./ATDAT/fe_5.scups',status='old')
         elseif(iel==14.and.ion==6) then
            open(11,file='./ATDAT/fe_6.scups',status='old')
         elseif(iel==14.and.ion==7) then
            open(11,file='./ATDAT/fe_7.scups',status='old')
         elseif(iel==14.and.ion==8) then
            open(11,file='./ATDAT/fe_8.scups',status='old')
         elseif(iel==14.and.ion==9) then
            open(11,file='./ATDAT/fe_9.scups',status='old')
         elseif(iel==4.and.ion==3) then
            open(11,file='./ATDAT/n_3.scups',status='old')
         elseif(iel==4.and.ion==4) then
            open(11,file='./ATDAT/n_4.scups',status='old')
         elseif(iel==3.and.ion==3) then
            open(11,file='./ATDAT/c_3.scups',status='old')
         elseif(iel==11.and.ion==4) then
            open(11,file='./ATDAT/s_4.scups',status='old')
         elseif(iel==11.and.ion==5) then
            open(11,file='./ATDAT/s_5.scups',status='old')
         elseif(iel==13.and.ion==2) then
            open(11,file='./ATDAT/ca_2.scups',status='old')
         elseif(iel==13.and.ion==5) then
            open(11,file='./ATDAT/ca_5.scups',status='old')
         endif
         ilm=0
         ium=0
         do i=1,20000
            read(11,*,end=22,err=22)il,iu,delta_ei,gfi,t_highi,nti,
     &           t_typei,c_vali
            read(11,*,err=22,end=22)(t_knoti(k),k=1,nti)
            read(11,*,err=22,end=22)(scups_knoti(k),k=1,nti)
            if(iu.le.nmax_lev) then
               nt(ionnr,il,iu)=nti
               do k=1,nt(ionnr,il,iu)
                  t_knot(ionnr,il,iu,k)=t_knoti(k)
                  scups_knot(ionnr,il,iu,k)=scups_knoti(k)
               enddo
               delta_e(ionnr,il,iu)=delta_ei
               t_type(ionnr,il,iu)=t_typei
               t_high(ionnr,il,iu)=t_highi
               c_val(ionnr,il,iu)=c_vali         

               do k=1,nt(ionnr,il,iu)
                  tn(k)=t_knot(ionnr,il,iu,k)
                  scup(k)=scups_knot(ionnr,il,iu,k)
               enddo
               do k=1,nknot
                  y2(k)=0.
               enddo
               call splinec2(tn,scup,nt(ionnr,il,iu),yp1,ypn,y2)
               do k=1,nknot
                  ybis(ionnr,il,iu,k)=y2(k)
               enddo

               ilm=max(ilm,il)
               ium=max(ium,iu)
               iumaxil(ionnr,il)=iu
            endif
         enddo
 22      continue
         ilmax(ionnr)=ilm
         iumax(ionnr)=ium         

      endif      

      do il=1,nmax_lev-1
         do iu=il+1,nmax_lev
            upsil(ionnr,il,iu)=0.
         enddo
      enddo

      do il=1,ilmax(ionnr)
         do iu=il+1,iumax(ionnr)            
            if(iu.le.iumaxil(ionnr,il)) then
               if(ybis(ionnr,il,iu,2).ne.0.) then
                  t_e_sc=abs(te/1.57888e5*delta_e(ionnr,il,iu))
c     if((t_e_sc.lt.t_high(iel,ion,il,iu)) .or.
c     &           (t_high(iel,ion,il,iu) < 0.)) then
                  if(t_type(ionnr,il,iu)==1.or.
     &                 t_type(ionnr,il,iu)==4) then
                     x=1.-log(c_val(ionnr,il,iu))
     &                    /log(t_e_sc+c_val(ionnr,il,iu))
                  elseif(t_type(ionnr,il,iu)==2.or.
     &                    t_type(ionnr,il,iu)==3) then
                     x=t_e_sc/(t_e_sc+c_val(ionnr,il,iu))
                  endif
                  do k=1,nt(ionnr,il,iu)
                     tn(k)=t_knot(ionnr,il,iu,k)
                     scup(k)=scups_knot(ionnr,il,iu,k)
                     y2(k)=ybis(ionnr,il,iu,k)
                  enddo
                  call splint2(tn,scup,y2,
     &                 nt(ionnr,il,iu),x,ya)
                  if(t_type(ionnr,il,iu)==1) then
                     ups=ya*log(t_e_sc+2.718281)
                  elseif(t_type(ionnr,il,iu)==2) then
                     ups=ya
                  elseif(t_type(ionnr,il,iu)==3) then
                     ups=ya/log(t_e_sc+1.)
                  elseif(t_type(ionnr,il,iu)==4) then
                     ups=ya*log(t_e_sc+c_val(ionnr,il,iu))
                  endif
                  r=1
                  
                  upsil(ionnr,il,iu)=ups
               else
                  upsil(ionnr,il,iu)=0.
               endif
            else
               ups=0.
            endif
         
         enddo
      enddo
      return
      end


      SUBROUTINE splint2(xa,ya,y2a,n,x,y)
      implicit real*8(a-h,o-z)
      INTEGER n
      REAL*8 x,y,xa(n),y2a(n),ya(n)
      INTEGER k,khi,klo
      REAL*8 a,b,h
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
c      if (h.eq.0.) pause 'bad xa input in splint'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**
     *2)/6.
      return
      END

      SUBROUTINE splinec2(x,y,n,yp1,ypn,y2)
      implicit real*8(a-h,o-z)
      INTEGER n,NMAX
      REAL*8 yp1,ypn,x(n),y(n),y2(n)
      PARAMETER (NMAX=500)
      INTEGER i,k
      REAL*8 p,qn,sig,un,u(NMAX)
      if (yp1.gt..99e30) then
        y2(1)=0.
        u(1)=0.
      else
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+
     *1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*
     *u(i-1))/p
11    continue
      if (ypn.gt..99e30) then
        qn=0.
        un=0.
      else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      return
      END
