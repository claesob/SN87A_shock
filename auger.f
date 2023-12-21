      SUBROUTINE AUGION(iel,IZMAX,xel,XA)
      IMPLICIT REAL*8(A-H,O-Z)
c     PARAMETER (MD=300,MDP1=MD+1)
      include "parameters.h"

      parameter (mh=2,mhe=3,mc=7,mn=8,mo=9,mne=11,mna=12,mmg=13,
     &     mal=14,msi=15,ms=17,mar=19,mca=21,mfe=27,mni=29,
     &     mmat=mh+mhe+mc+mn+mo+mne+mna+mmg+mal+msi+ms+mar+mca+mfe+mni,
     &     ispecy=15)
      common/abc1/alc(mc),aln(mn),alo(mo),alne(mne),alna(mna),
     &     almg(mmg),alal(mal),alsi(msi),alsu(ms),alar(mar),alca(mca),
     &     alfe(mfe),alni(mni)
      common/abc2/coh(1),cohe(2),coc(mc-1),con(mn-1),coo(mo-1),
     &     cone(mne-1),cona(mna-1),comg(mmg-1),coal(mal-1),cosi(msi-1),
     &     cosu(ms-1),coar(mar-1),coca(mca-1),cofe(mfe-1),coni(mni-1)      
      COMMON/SECX/CSEC(20),DCSDX(20),CISEC(20),DCISDX
      COMMON/PHY/DEN(MD)
      COMMON/IND/IK
      COMMON/RSI/ZS(18,5)
      common/raug/zion(30,27,7),geion(30,27)
      COMMON/AUG/AUG
      parameter(nz=30,nion=27,nshell=10)
      integer iel,ion,shell,epsi,fr(10)
      integer ns,kmax,n
      real*8 fr_aug,eion,en_aug,en_augi,eioni
      common/augerfrac/eion(nz,nion,nshell),en_aug(nz,nion,nshell),
     &     fr_aug(nz,nion,nshell,10),kmax(nz,nion,nshell),
     &     init_augfrac
      common/shells/ns(30,27)
      common/rec_coll/alrec(30,30),collion(30,30),ction(14,27)
      real*8 x(30),dx(30),aa(31,31),ztot(30)
      DIMENSION XA(30),ZSA(30,18,5),phrate(30,30)
c       dimension nionel(14)
c       data nionel/1,2,6,7,8,10,11,12,13,14,16,18,20,26/

C     IONIZATION BALANCE EQUATIONS TAKING INTO ACCOUNT AUGER
C      IONIZATION AND COLLISIONAL IONIZATION 
C     XR(I)=N(I+1)/N(I)
c      CALL COSI

c     set up ionization balance

      nionel=iel

      do i=1,izmax
         do j=1,izmax+1
            aa(i,j)=0.
         enddo
      enddo
c$$$      do ion=1,izmax
c$$$         if(iel==8) then
c$$$            alrec(ion)=alo(ion)
c$$$         elseif(iel==10) then
c$$$            alrec(ion)=alne(ion)
c$$$         elseif(iel==12) then
c$$$            alrec(ion)=almg(ion)
c$$$         elseif(iel==14) then
c$$$            alrec(ion)=alsi(ion)
c$$$         elseif(iel==16) then
c$$$            alrec(ion)=alsu(ion)
c$$$         elseif(iel==18) then
c$$$            alrec(ion)=alar(ion)
c$$$         endif
c         write(6,91)iel,ion,alrec(ion),(zion(iel,ion,is),is=1,5)
c 91      format('Auger  ',2i5,1pe11.3,10e11.3)
c      enddo

c rates from ion to ion+ion+j, includig auger fractions from Kaastra & Mewe
      do ion=1,izmax
         ztot(ion)=0.
         ztot2=0.
         do i=1,ns(iel,ion)
c$$$            do j=1,kmax(iel,ion,i)
c$$$               phrate(ion,ion+j)=zsa(iel,ion,i)*fr_aug(iel,ion,i,j)
c$$$               ztot(ion)= phrate(ion,ion+j) + ztot(ion)
c$$$               write(6,921)ion,i,j,fr_aug(iel,ion,i,j),zsa(iel,ion,i),
c$$$     &              phrate(ion,ion+j),ztot(ion)
c$$$ 921           format('ztot ',3i4,1pe12.3,10e12.3)
c$$$            enddo
            ztot2= zion(iel,ion,i) + ztot2
            ztot(ion)=ztot2
c            write(6,922)ion,i,zion(iel,ion,i),ztot2
 922        format('ztot2 ',2i4,1pe12.3,10e12.3)
         enddo
      enddo

      do i=1,izmax+1
c         write(6,9287)i,alrec(iel,i)
 9287    format(' alpha ',i5,1pe12.3)
      enddo

      do i=1,izmax

c     right hand side
         aa(i,izmax+1)=0.

         do k=1,10
c lower ionization stage added to i            
            il=i-k
            if(il>0 ) then
c go through all the shells of il               
               do nsi=1,ns(iel,il)
                  if(kmax(iel,il,nsi)>0) then                     
c                     write(6,*)'il,nsi,kmax(iel,il,nsi) ',
c     &                    il,nsi,kmax(iel,il,nsi)
                     aa(i,il)=fr_aug(iel,il,nsi,k)*zion(iel,il,nsi) +
     &                    aa(i,il)
c                     write(6,927)i,il,k,fr_aug(iel,il,nsi,k),
c     &                    zion(iel,il,nsi),aa(i,il)
 927                 format('i,il,k,fr_aug,zion,aa ',
     &                    3i4,1pe11.3,10e11.3)
                  endif
               enddo
            endif
         enddo

c     recombination from higher stages
         aa(i,i+1)=alrec(iel,i+1)*xel*den(ik)
c     recombination to lower stages and auger losses to higher stages 
         aa(i,i)=-alrec(iel,i)*xel*den(ik)-ztot(i)
      enddo

c     replace highest stage with number cons
      do i=1,izmax
         aa(izmax,i)=1.d0
      enddo

      aa(izmax,izmax+1)=1.


      x1_x2=-aa(1,2)/aa(1,1)

      x3_x2=(x1_x2*aa(2,1)-aa(2,2))/aa(2,3)

      x2=1./(x1_x2+1.d0+x3_x2)

      x1=x1_x2*x2

      x3=x3_x2*x2

      write(6,*)'x1/x2,x3/x2 ',x1/x2,x3/x2
      write(6,*)'x1,x2,x3 ',x1,x2,x3


      do i=1,izmax
c         write(6,9781)i,(aa(i,j),j=1,izmax+1)
 9781    format('aa aug ion ',i5,1pe12.3,20e12.3)
      enddo

      di=simul_ion(izmax,aa,x,eps,1,nrc)

      dxmax =0.

      do i=1,izmax
         xa(i)=x(i)
         write(6,9128)iel,i,ztot(i),alrec(iel,i+1),x(i)
 9128    format('balance ',2i5,1pe12.3,10e12.3)         
      enddo
      RETURN
      END      


      subroutine auger(izmax,deel,xa)
      implicit real*8(a-h,o-z)
c     parameter (md=300,mdp1=md+1)
      include "parameters.h"
c      include 'PARAM'
      parameter (mh=2,mhe=3,mc=7,mn=8,mo=9,mne=11,mna=12,mmg=13,
     &     mal=14,msi=15,ms=17,mar=19,mca=21,mfe=27,mni=29,
     &     mmat=mh+mhe+mc+mn+mo+mne+mna+mmg+mal+msi+ms+mar+mca+mfe+mni,
     &     ispecy=15)
      common/abc1/alc(mc),aln(mn),alo(mo),alne(mne),alna(mna),
     &     almg(mmg),alal(mal),alsi(msi),alsu(ms),alar(mar),alca(mca),
     &     alfe(mfe),alni(mni)
      common/abc2/coh(1),cohe(2),coc(mc-1),con(mn-1),coo(mo-1),
     &     cone(mne-1),cona(mna-1),comg(mmg-1),coal(mal-1),cosi(msi-1),
     &     cosu(ms-1),coar(mar-1),coca(mca-1),cofe(mfe-1),coni(mni-1)      

      common/phy/den(md)
      common/ind/ik
      common/rar/zar(18,5)
      common/augalu/alaug(27),coaug(26)
      common/aug/aug
      dimension xa(30),zsa(18,5)
c     ionization balance equations taking into account auger
c      ionization and collisional ionization
c     xr(i)=n(i+1)/n(i)
      do  iz=1,izmax
         do  mi=1,5
            zsa(iz,mi)=zar(iz,mi)
         enddo
      enddo
      do iz=1,izmax
         ze=0.
         do mi=1,5
            ze=ze+zar(iz,mi)
         enddo
         if(iz.le.1) then
            zeff=ze
            xa(1)=(coaug(1)+ze/(den(ik)*deel))/alaug(iz+1)
         else
            xs=1.
            ia1=1
            ia2=iz-1
            if(iz.gt.6) then
               ia1=iz-5
            endif
c     ia=i (weisheit)
c     i=z-i (weisheit)
            s=0.
            do ia=ia2,ia1,-1
               i=iz-ia
               if(xa(ia).le.0.) xa(ia)=1.e-20
               xs=xs/xa(ia)
               n=18+1-ia
               call psi(ia,i,n,zsa,ps)
               s=s+xs*pse
            enddo
            zeff=ze+zeff/xa(iz-1)-s
            if(aug.lt.0.) zeff=ze
            xa(iz)=(coaug(iz)+zeff/(deel*den(ik)))/alaug(iz+1)
         endif
c         write(6,9129)iz,xa(iz),coaug(iz),ze,zeff,deel,
c     &        den(ik),alaug(iz+1),xa(iz)
 9129    format(' araug ',i5,1pe12.3,10e12.3)
      enddo
      return
      end

      SUBROUTINE RATEAUG(iel,izmax)
C     IEL IN CF (3 FFOR C)
c gsv2 in CF system 3 for C      
      IMPLICIT REAL*8(A-H,O-Z)
c      PARAMETER (NE1=-200,NE2=130,NE3=NE2+1)
      include 'parameters.h'
      COMMON/FRE/NINT,JMIN,JJ
      COMMON/RSI/ZS(18,5)
      common/raug/zion(30,27,7),geion(30,27)
      COMMON/IND/IK
      COMMON/DTAU/FLUX(NE1:NE2)
      COMMON/INT/FL(2,NE1:NE2),SI(14,27,NE1:NE2),E1(NE1:NE3),E(NE1:NE3)
c      COMMON/CSI/GS(14,5,NE1:NE2)
c      COMMON/CSI_vern/GSV(14,5,NE1:NE2)
      COMMON/CSI_vern2/GSV2(14,30,7,NE1:NE2)
      common/ethresh/ethresh(30,30,7)
      common/shells/ns(30,27)
      DATA PI/3.1415926E0/,ELCH/1.60219E-12/,AMU/1.660531E-24/
c     correspondance between Z and element enumeration
      integer elcf(30),elcfi
      data elcf/1,2,0,0,0,3,4,5,0,6,7,8,9,10,0,11,0,12,0,13,
     &     0,0,0,0,0,14,0,15,0,0/
      elcfi=elcf(iel)
c      write(6,*)' i rateaug ',elcfi,izmax,jmin,jj,(ns(elcfi,iz),iz=1,izmax)
      ZA=0.
      ga=0.
      DO IZ=1,IZMAX
         geion(elcfi,IZ)=0.
         DO  IS=1,ns(elcfi,iz)
            Zion(elcfi,IZ,IS)=0.
            DO J=JMIN,JJ
c               ZA=4.*PI*FLUX(J)*GSv2(elcfi,IZ,IS,J)*(E(J+1)-E(J))/
c     &              (ELCH*E1(J))
               ZA=4.*PI*FL(2,J)*GSv2(elcfi,IZ,IS,J)*(E(J+1)-E(J))/
     &              (ELCH*E1(J))         
               ZION(elcfi,IZ,IS)=ZION(elcfi,IZ,IS)+ZA
               GA=4.*PI*FL(2,J)*GSv2(elcfi,IZ,IS,J)*(1. - ethresh(elcfi,IZ,IS)/e1(j))* (E(J+1)-E(J))         
               GeION(elcfi,iz)=GeION(elcfi,iz)+GA
c               if(iel==-12.and.j>0) then
c               write(6,91)elcfi,iz,is,j,e1(j),gsv2(elcfi,iz,is,j),fl(2,j),
c     &              zion(elcfi,iz,is),geion(elcfi,iz)
c               endif
            enddo
         enddo
c         write(6,9)iel,iz,(zion(elcfi,iz,is),is=1,5)
 9       format('Ion rate ',2i5,1pe11.3,10e11.3)
 91      format('iel,iz,is,j,e1,gsv2,zion,geion',4i5,1pe11.3,10e11.3)
      enddo
      RETURN
      END
      
