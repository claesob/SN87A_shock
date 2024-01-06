      SUBROUTINE EMISS_cont(te,xel,dens,cool)
C   ******************************************************
C     em(j) = emissivity in erg/(s cm3 eV sr)
C     cool = 4*pi*energy integration of em(j) (= total cooling) 
C     coolco = continuum cooling only
C     (em (j) and cool both include Xe*dens**2)
C
C     Continuum emission from Gronenschild & Mewe (A&A Suppl.
C     32, 283 (1978)). Note: Only emission to ground state. See 
c     Mewe et al 1986 for more accurate routine.
c
C   ******************************************************
      IMPLICIT REAL*8(A-H,O-Z)                                          
      parameter (mh=2,mhe=3,mc=7,mn=8,mo=9,mne=11,mna=12,mmg=13,
     &mal=14,msi=15,ms=17,mar=19,mca=21,mfe=27,mni=29,
     &mmat=mh+mhe+mc+mn+mo+mne+mna+mmg+mal+msi+ms+mar+mca+mfe+mni,
     &ispecy=15)
C                                                                       
      include 'param'
c      PARAMETER (NE1=-100,NE2=250,NE3=NE2+1)
      common/abl/abn(15)
      common/ionx/xion(md,14,27)
      COMMON/IND/Ik
      COMMON/INT/FL(2,NE1:NE2),SI(14,27,NE1:NE2),E1(NE1:NE3),E(NE1:NE3)
      COMMON/DIF/EM(MD,NE1:NE2),TAU(MD,NE1:NE2),TAUTOT(MD,NE1:NE2),
     &           EMC(MD,NE1:NE2)
      COMMON/FRE/ninq,JMIN,JJ                                                
      common/opac/cappa(ne1:ne2)
      common/timecheck/time,itime
      COMMON/SIKq/sikh(mh-1,ne1:ne2),sikhe(mhe-1,ne1:ne2),sikc
     &(mc-1,-190:901),sikn(mn-1,ne1:ne2),siko(mo-1,ne1:ne2),
     &sikne(mne-1,ne1:ne2),sikna(mna-1,ne1:ne2),sikmg(mmg-1,
     &ne1:ne2),sikal(mal-1,ne1:ne2),siksi(msi-1,ne1:ne2),siks
     &(ms-1,ne1:ne2),sikar(mar-1,ne1:ne2),sikca(mca-1,ne1:ne2),
     &sikfe(mfe-1,ne1:ne2),sikni(mni-1,ne1:ne2) 
C
      dimension zh(2),zhe(3),zc(7),zn(8),zo(9),zne(11),zna(12),
     &zmg(13),zal(14),zsi(15),zs(17),zar(19),zca(21),zfe(27),
     &zni(29),zz(29)
      dimension poth(1),pothe(2),potc(6),potn(7),poto(8),potne(10),
     &potna(11),potmg(12),potal(13),potsi(14),pots(16),potar(18),
     &potca(20),potfe(26),potni(28),potz(29)
      dimension bfh(1),bfhe(2),bfc(6),bfn(7),bfo(8),bfne(10),
     &bfna(11),bfmg(12),bfal(13),bfsi(14),bfs(16),bfar(18),bfca
     &(20),bffe(26),bfni(28),bfz(29)
      dimension ff(ispecy),free(ne1:ne2),ff1(ispecy,mni),
     &ff2(ispecy,mni)
      dimension fb(ispecy),fbound(ne1:ne2),fb1(ispecy,mni),
     &fb2(ispecy,mni)
      dimension itrap(mni),exche(14),f2h2(ispecy),f1h2(ispecy),
     &f2he2(ispecy),f1he2(ispecy),ftwo(ne1:ne2),ffisp(ispecy)
      dimension fx_ind(ispecy,mni),fx_el(ispecy)
      dimension aeff(15,mni)

C  -- Effective ion charges (see Gronenschild & Mewe (A&A), 1978).
      data zh/0.,1./,zhe/0.,1.34,2./
      data zc/0.,1.82,2.68,3.75,4.36,5.37,6./
      data zn/0.,2.07,2.95,3.73,4.77,5.37,6.37,7./
      data zo/0.,2.00,3.21,4.02,4.77,5.79,6.37,7.37,8./
      data zne/0.,2.52,3.48,4.32,5.34,6.09,6.82,7.80,8.38,9.38,
     &10./
      data zna/0.,1.84,3.73,4.59,5.40,6.37,7.11,7.82,8.81,9.39,
     &10.38,11./
      data zmg/0.,2.24,3.15,4.85,5.66,6.44,7.42,8.13,8.85,9.82,
     &10.39,11.38,12./
      data zal/0.,1.99,3.53,4.34,5.94,6.73,7.48,8.42,9.16,9.85,
     &10.83,11.40,12.38,13./
      data zsi/0.,2.32,3.28,4.71,5.46,7.01,7.76,8.51,9.44,10.16,
     &10.86,11.83,12.40,13.39,14./
      data zs/0.,2.62,3.94,4.81,5.59,6.94,7.64,9.09,9.82,10.56,
     &11.47,12.19,12.88,13.85,14.42,15.40,16./
      data zar/0.,3.23,4.27,5.20,6.29,7.05,7.77,9.10,9.73,11.15,
     &11.87,12.59,13.48,14.20,14.90,15.86,16.43,17.41,18./
      data zca/0.,2.68,3.74,5.82,6.67,7.48,8.49,9.20,9.90,11.18,
     &11.82,13.20,13.90,14.62,15.51,16.22,16.93,17.88,18.45,19.42,
     &20./
      data zfe/0.,3.05,4.37,4.51,6.02,7.05,8.09,9.10,10.00,12.47,
     &13.17,13.85,14.80,15.46,16.11,17.39,18.01,19.29,19.99,20.69,
     &21.57,22.29,23.00,23.95,24.52,25.48,26./
      data zni/0.,3.20,3.47,4.83,6.03,7.07,8.45,9.38,10.35,11.30,
     &12.20,14.57,15.26,15.94,16.87,17.52,18.17,19.44,20.06,21.32,
     &22.02,22.73,23.60,24.32,25.04,25.98,26.56,27.49,28./
C
C - Lowest ionization potentials from ground state. 
      data poth/13.598/,pothe/24.587,54.416/
      data potc/11.3,24.4,47.9,64.5,392.,490./
      data potn/14.5,29.6,47.4,77.5,97.9,552.,667./
      data poto/13.6,35.1,54.9,77.4,114.,138.,739.,871./
      data potne/21.6,41.1,63.5,97.1,126.,158.,207.,239.,1196.,
     &1362./
      data potna/5.1,47.3,71.7,99.,138.,172.,208.,264.,300.,1465.,
     &1649./
      data potmg/7.6,15.,80.1,109.,141.,187.,225.,266.,328.,367.,
     &1762.,1963./
      data potal/6.,18.8,28.4,120.,154.,190.,241.,285.,330.,399.,
     &442.,2086.,2304./
      data potsi/8.1,16.3,33.5,45.1,167.,205.,246.,303.,351.,401.,
     &476.,523.,2438.,2673./
      data pots/10.4,23.4,35.,47.3,72.7,88.1,281.,328.,379.,447.,
     &505.,564.,652.,707.,3224.,3493./
      data potar/15.8,27.6,40.9,59.7,75.2,91.2,125.,143.,423.,
     &479.,539.,618.,686.,755.,855.,918.,4121.,4426./
      data potca/6.1,11.9,51.2,67.3,84.5,109.,128.,148.,189.,211.,
     &592.,657.,727.,818.,894.,974.,1087.,1157.,5129.,5470./
      data potfe/7.9,16.2,30.7,54.8,75.,99.,125.,151.,235.,262.,
     &290.,331.,361.,392.,457.,490.,1265.,1358.,1456.,1582.,
     &1689.,1799.,1950.,2045.,8828.,9278./
      data potni/8.7,18.2,35.2,54.9,75.5,108.,133.,162.,193.,225.,
     &321.,352.,384.,430.,464.,499.,571.,608.,1546.,1648.,1756.,
     &1894.,2011.,2131.,2295.,2399.,10280.,10790./

C -- 'number of empty states in valence shell'.

      data bfh/2./,bfhe/1.,2./,bfc/2.5,3.,3.5,4.,1.,2./,
     &bfn/2.,2.5,3.,3.5,4.,1.,2./,bfo/1.5,2.,2.5,3.,3.5,4.,1.,2./
      data bfne/0.5,1.,1.5,2.,2.5,3.,3.5,4.,1.,2./
      data bfna/6.,0.5,1.,1.5,2.,2.5,3.,3.5,4.,1.,2./
      data bfmg/5.67,6.,0.5,1.,1.5,2.,2.5,3.,3.5,4.,1.,2./
      data bfal/5.33,5.67,6.,0.5,1.,1.5,2.,2.5,3.,3.5,4.,1.,2./
      data bfsi/5.,5.33,5.67,6.,0.5,1.,1.5,2.,2.5,3.,3.5,4.,1.,2./
      data bfs/4.33,4.67,5.,5.33,5.67,6.,0.5,1.,1.5,2.,2.5,3.,3.5,
     &4.,1.,2./
      data bfar/3.67,4.,4.33,4.67,5.,5.33,5.67,6.,0.5,1.,1.5,2.,2.5,
     &3.,3.5,4.,1.,2./
      data bfca/7.75,8.,3.67,4.,4.33,4.67,5.,5.33,5.67,6.,0.5,1.,1.5,
     &2.,2.5,3.,3.5,4.,1.,2./
      data bffe/7.75,8.,1.67,2.,2.33,2.67,3.,3.33,3.67,4.,4.33,4.67,
     &5.,5.33,5.67,6.,0.5,1.,1.5,2.,2.5,3.,3.5,4.,1.,2./
      data bfni/8.,0.67,1.,1.33,1.67,2.,2.33,2.67,3.,3.33,3.67,4.,
     &4.33,4.67,5.,5.33,5.67,6.,0.5,1.,1.5,2.,2.5,3.,3.5,4.,1.,2./
C 
C -- Two-photon excitation energies for He-like transitions.
      data exche/20.62,304.,425.,569.,910.,1120.,1350.,1590.,1840.,
     &2430.,3120.,3890.,6525.,7700./

      DATA PI/3.1415926d0/
C
      do i=1,ispecy
         do j=1,mni
            fx_ind(i,j)=0.
         enddo
      enddo
      xelde=xel*dens**2
      tev=te/1.1603d4

c      write(6,*)' te,xel,dens cont ',ik,te,xel,dens

      gh=0.055d0
      ghe=0.045d0
      fh=0.416d0
      fhe=0.832d0
      do j=jmin,jj
         free(j)=0.0d0
         fbound(j)=0.0d0
         ftwo(j)=0.0d0
c         em(ik,j)=0.0d0
         cappa(j)=0.0d0
      enddo

      do isp=1,ispecy
         f2h2(isp)=0.0d0
         f2he2(isp)=0.0d0
         do i=1,mni
            fb1(isp,i)=0.0d0
            ff1(isp,i)=0.0d0
            fb2(isp,i)=0.0d0
            ff2(isp,i)=0.0d0
         enddo
      enddo
      cool=0.0d0
      coolco=0.0d0
      coolff=0.0d0
C
C   ******************
C   Continuum emission
C   ******************
C
C  =>   Energy loop   =>
C  *********************
      do isp=1,15
        ffisp(isp)=0.
        do i=1,mni
           aeff(isp,i)=0.           
        enddo
      enddo

      do j=jmin,jj
         ej=e(j)
         ej1=e1(j)
c!! skip Ni
         do isp=1,ispecy-1
            f1h2(isp)=f2h2(isp)
            f2h2(isp)=0.0d0
            f1he2(isp)=f2he2(isp)
            f2he2(isp)=0.0d0
            do  i=1,mni
               ff1(isp,i)=ff2(isp,i)
               fb1(isp,i)=fb2(isp,i)
               ff2(isp,i)=0.0d0
               fb2(isp,i)=0.0d0
            enddo
            fb(isp)=0.0d0
            ff(isp)=0.0d0
         enddo
         if((ej/tev).gt.1.d2) then
            etexp=0.0d0
         else
            etexp=dexp(-ej/tev)
         endif
C
C    - free-free  and  free-bound -
C
         cff=(1.6545d-23*xel*dens**2)/dsqrt(te)
         if (isnan(cff)) then
            write(6,9291)xel,dens,te
 9291    format('fb xel,dens,te ',1pe12.3,10e12.3)
               endif

c         do isp=1,ispecy
c!!! skip H and He (done in emiss)

         do isp=3,ispecy
            if(isp.eq.9) goto 102
            if(abn(isp).lt.1.d-10) goto 102
            if(isp.eq.1) mion=mh
            if(isp.eq.2) mion=mhe
            if(isp.eq.3) mion=mc
            if(isp.eq.4) mion=mn
            if(isp.eq.5) mion=mo
            if(isp.eq.6) mion=mne
            if(isp.eq.7) mion=mna
            if(isp.eq.8) mion=mmg
            if(isp.eq.9) mion=mal
            if(isp.eq.10) mion=msi
            if(isp.eq.11) mion=ms
            if(isp.eq.12) mion=mar
            if(isp.eq.13) mion=mca
            if(isp.eq.14) mion=mfe
            if(isp.eq.15) mion=mni
            do i=2,mion
               abxa=abn(isp)*xion(ik,isp,i)

               itrap(i)=0
               if(abxa.lt.1.d-10) goto 103
               if(isp.eq.1) zz(i)=zh(i)
               if(isp.eq.2) zz(i)=zhe(i)
               if(isp.eq.3) zz(i)=zc(i)
               if(isp.eq.4) zz(i)=zn(i)
               if(isp.eq.5) zz(i)=zo(i)
               if(isp.eq.6) zz(i)=zne(i)
               if(isp.eq.7) zz(i)=zna(i)
               if(isp.eq.8) zz(i)=zmg(i)
               if(isp.eq.9) zz(i)=zal(i)
               if(isp.eq.10) zz(i)=zsi(i)
               if(isp.eq.11) zz(i)=zs(i)
               if(isp.eq.12) zz(i)=zar(i)
               if(isp.eq.13) zz(i)=zca(i)
               if(isp.eq.14) zz(i)=zfe(i)
               if(isp.eq.15) zz(i)=zni(i)
               if(isp.eq.1) potz(i)=poth(i-1)
               if(isp.eq.2) potz(i)=pothe(i-1)
               if(isp.eq.3) potz(i)=potc(i-1)
               if(isp.eq.4) potz(i)=potn(i-1)
               if(isp.eq.5) potz(i)=poto(i-1)
               if(isp.eq.6) potz(i)=potne(i-1)
               if(isp.eq.7) potz(i)=potna(i-1)
               if(isp.eq.8) potz(i)=potmg(i-1)
               if(isp.eq.9) potz(i)=potal(i-1)
               if(isp.eq.10) potz(i)=potsi(i-1)
               if(isp.eq.11) potz(i)=pots(i-1)
               if(isp.eq.12) potz(i)=potar(i-1)
               if(isp.eq.13) potz(i)=potca(i-1)
               if(isp.eq.14) potz(i)=potfe(i-1)
               if(isp.eq.15) potz(i)=potni(i-1)
               if(isp.eq.1) bfz(i)=bfh(i-1)
               if(isp.eq.2) bfz(i)=bfhe(i-1)
               if(isp.eq.3) bfz(i)=bfc(i-1)
               if(isp.eq.4) bfz(i)=bfn(i-1)
               if(isp.eq.5) bfz(i)=bfo(i-1)
               if(isp.eq.6) bfz(i)=bfne(i-1)
               if(isp.eq.7) bfz(i)=bfna(i-1)
               if(isp.eq.8) bfz(i)=bfmg(i-1)
               if(isp.eq.9) bfz(i)=bfal(i-1)
               if(isp.eq.10) bfz(i)=bfsi(i-1)
               if(isp.eq.11) bfz(i)=bfs(i-1)
               if(isp.eq.12) bfz(i)=bfar(i-1)
               if(isp.eq.13) bfz(i)=bfca(i-1)
               if(isp.eq.14) bfz(i)=bffe(i-1)
               if(isp.eq.15) bfz(i)=bfni(i-1)
               if(etexp.ne.0.0d0) then
                  ff2(isp,i)=cff*etexp*abxa*zz(i)**2*
     &                 ffgau(ej,tev,zz(i))
c               write(6,*)'isp,i,etexp,cff,abxa ',isp,i,etexp,cff,abxa
c               write(6,*)'ej,tev,zz(i),ff2(isp,i) ',ej,tev,zz(i),
c     &              ff2(isp,i)


               endif

               if(potz(i).le.ej) then
                  potexp=(ej-potz(i))/tev
                  if(potexp.gt.1.d2) then
                     fb2(isp,i)=0.0d0
                  else
                     fb2(isp,i)=cff*(poth(1)/tev)*abxa*bfz(i)*
     &                    (potz(i)/poth(1))**2*dexp(-potexp)
                  endif
               else
                  fb2(isp,i)=0.0d0
               endif

               if(potz(i).le.ej) then

                  if(potz(i).gt.e(j-1).and.potz(i).le.e(j)) then

                     itrap(i)=1

                     potexp=(e(j)-potz(i))/tev

                     if(potexp.gt.1.d2) then
                        fac1 = 1.
                     else
                        fac1 = (1. - exp(-(e(j)-potz(i))/tev))
                     endif

                  else

                     itrap(i)=0

                     potexp=(e(j-1)-potz(i))/tev

                     if(potexp.gt.1.d2) then

                        fac1 = 0.
                     else

                     fac1= exp(-(e(j-1)-potz(i))/tev) - 
     &                    exp(-(e(j)-potz(i))/tev)

                     endif

                  endif

                  fac=tev*fac1
                     
                  fb2(isp,i)=cff*(poth(1)/tev)*abxa*bfz(i)*
     &                 (potz(i)/poth(1))**2*fac/(e(j)-e(j-1))
                  if (isnan(fb2(isp,i))) then
                     write(6,9182)isp,i,j,cff,tev,poth(1),abxa,bfz(i),potz(i),fac,e(j),e(j-1)
 9182                format('fb2 a ',3i5,1pe12.3,10e12.3)
                  endif


               else
                  fb2(isp,i)=0.0d0
               endif

c note that fb includes all ions of a given element

               fb(isp)=fb(isp) + fb2(isp,i)
               if (isnan(fb2(isp,i))) then
                  write(6,9282)isp,i,j,fb2(isp,i)
 9282             format('fb2 ',3i5,1pe12.3,10e12.3)
               endif

               daeff=fb2(isp,i)*(e(j)-e(j-1))/(1.602e-12*e1(j)*
     &              xel*abxa*dens**2)

               aeff(isp,i)=aeff(isp,i)+daeff

               if((isp.eq.8.or.isp.eq.14).and.fb2(isp,i).gt.
     &              1.d-50.and.j.ge.-33.and.j.le.-30) then
                  aleffo=fb(isp)*(e(j)-e(j-1))/(1.602e-12*e1(j)*
     &                 xel*abxa*dens**2)
                  
                  
c                  write(6,9387)j,itrap(i),e1(j),fb(isp),
c     &                 fac1,fac,fb2(isp,i),aleffo,aeff(12,i)
 9387             format('aleffq ',2i5,1pe12.3,10e12.3)
               endif


c               if(te.lt.1.e5.and.isp.eq.12.and.j.eq.-28) then
c                  write(6,928)j,itrap(i),ej,fb1(isp,i),
c     &                 fb2(isp,i),potz(i)
c 928              format('fb1 ',2i5,1pe12.3,10e12.3)
c               endif

 103           continue
            enddo
            if(j.ne.jmin) then
               do i=2,mion
                  ff(isp)=ff(isp)+0.5d0*(ff2(isp,i)+ff1(isp,i))
               enddo
               free(j)=free(j)+ff(isp)
               fbound(j)=fbound(j)+fb(isp)
               if (isnan(fb(isp))) then
                  write(6,9281)isp,j,fb(isp),fbound(j)
 9281             format('fbound ',2i5,1pe12.3,10e12.3)

               endif
               ffisp(isp)=ffisp(isp)+ff(isp)*(e(j)-e(j-1))/xelde

               if(j.eq.-36.or.j.eq.-28) then
c                  write(6,9992)j,isp,i,abxa,xelde,pi,fb(isp),
c     &                 ff(isp)*(e(j)-e(j-1))/xelde,
c     &                 fbound(j),
c     &                 fbound(j)/(4.d0*pi*xelde)
 9992             format('emconbfff ',3i5,1pe12.3,10e12.3)
               endif

            endif
 102        continue
            if(j.eq.jj) then
c               write(6,9277)isp,te,(aeff(isp,i),i=1,mion)
 9277          format(' aleff ',i5,1pe12.3,30e12.3)
            endif
         enddo




C     
C    - two-photon -
C
            ctwo=(2.727d-15*ej*xel*dens**2)/dsqrt(te)
            do isp=1,ispecy
               if(abn(isp).lt.1.d-10) goto 5001
               if(isp.eq.1) mion=mh
               if(isp.eq.2) mion=mhe
               if(isp.eq.3) mion=mc
               if(isp.eq.4) mion=mn
               if(isp.eq.5) mion=mo
               if(isp.eq.6) mion=mne
               if(isp.eq.7) mion=mna
               if(isp.eq.8) mion=mmg
               if(isp.eq.9) mion=mal
               if(isp.eq.10) mion=msi
               if(isp.eq.11) mion=ms
               if(isp.eq.12) mion=mar
               if(isp.eq.13) mion=mca
               if(isp.eq.14) mion=mfe
               if(isp.eq.15) mion=mni
               abxa=abn(isp)*xion(ik,isp,mion-1)
               if(abxa.lt.1.d-10) goto 5002
               exch=10.1985d0*(dble(mion-1))**2
               if(ej.lt.exch) then
                  exchtev=exch/tev
                  if(exchtev.gt.1.d2) goto 5003
                  y=ej/exch
                  fi=2.623d0*dsqrt(dcos(pi*(y-0.5d0)))
                  f2h2(isp)=(ctwo*abxa*gh*fh*fi*dexp(-exch/tev))/
     &                 exch**2 
                  ca=1.d0+(0.0125d0*xel*dens)/(dsqrt(te)*
     &                 (dble(mion-1))**8)
                  f2h2(isp)=f2h2(isp)/ca
                  ftwo(j)=ftwo(j)+0.5d0*(f2h2(isp)+f1h2(isp))
 5003             continue
               endif
               if(j.ne.jmin) then
                  if(exch.gt.e(j-1).and.exch.le.e(j)) then
                     f2h2(isp)=0.0d0
                     ftwo(j)=ftwo(j)+(0.5d0*(f2h2(isp)+f1h2(isp))*
     &                    (exch-e(j-1)))/(e(j)-e(j-1))
                  endif
               endif
 5002          continue
               if(isp.eq.1) goto 5001
               abxa=abn(isp)*xion(ik,isp,mion-2)
               if(abxa.lt.1.d-10) goto 5001
               exch=exche(isp-1)
               if(ej.lt.exch) then
                  exchtev=exch/tev
                  if(exchtev.gt.1.d2) goto 5004
                  y=ej/exch
                  fi=2.623d0*dsqrt(dcos(pi*(y-0.5d0)))
                  fhee=fhe*(1.d0-1.34d0/dble(mion-1))
                  f2he2(isp)=(ctwo*abxa*ghe*fhee*fi*dexp(-exch/tev))/
     &                 exch**2 
                  ca=1.d0+(1.875d-4*xel*dens)/(dsqrt(te)*
     &                 (dble(mion-2))**8)
                  f2he2(isp)=f2he2(isp)/ca
                  ftwo(j)=ftwo(j)+0.5d0*(f2he2(isp)+f1he2(isp))
 5004             continue
               endif
               if(j.eq.jmin) goto 5001
               if(exch.gt.e(j-1).and.exch.le.e(j)) then
                  f2he2(isp)=0.0d0
                  ftwo(j)=ftwo(j)+(0.5d0*(f2he2(isp)+f1he2(isp))*
     &                 (exch-e(j-1)))/(e(j)-e(j-1))
               endif
 5001       continue
         enddo
      enddo
C     

C
C- Calculate emissivity and opacity -
      xelde=xel*dens**2
      do j=jmin+1,jj
         if(free(j).lt.1.d-200) free(j)=1.d-200
         if((e1(j-1)/tev).lt.2.d2) then
            planck=(5.0403d10*e1(j-1)**3)/(dexp(e1(j-1)/tev)-1.d0)
         else
            planck=6.975d-77*e1(j-1)**3
         endif
c!! skip free-bound opac
         goto 1111
         opff=1.d-50+free(j)/(4.d0*pi*planck)
         opbf=abn(1)*dens*xion(ik,1,1)*sikh(1,j-1)
         do  ii=1,mhe-1
            opbf=opbf+abn(2)*dens*xion(ik,2,ii)*sikhe(ii,j-1)
         enddo
         do  ii=1,mc-1
            opbf=opbf+abn(3)*dens*xion(ik,3,ii)*sikc(ii,j-1)
         enddo
         do ii=1,mn-1
            opbf=opbf+abn(4)*dens*xion(ik,4,ii)*sikn(ii,j-1)
         enddo
         do ii=1,mo-1
            opbf=opbf+abn(5)*dens*xion(ik,5,ii)*siko(ii,j-1)
         enddo
         do ii=1,mne-1
            opbf=opbf+abn(6)*dens*xion(ik,6,ii)*sikne(ii,j-1)
         enddo
         do ii=1,mna-1
            opbf=opbf+abn(7)*dens*xion(ik,7,ii)*sikna(ii,j-1)
         enddo
         do ii=1,mmg-1
            opbf=opbf+abn(8)*dens*xion(ik,8,ii)*sikmg(ii,j-1)
         enddo
         do  ii=1,mal-1
            opbf=opbf+abn(9)*dens*xion(ik,9,ii)*sikal(ii,j-1)
         enddo
         do ii=1,msi-1
            opbf=opbf+abn(10)*dens*xion(ik,10,ii)*siksi(ii,j-1)
         enddo
         do ii=1,ms-1
            opbf=opbf+abn(11)*dens*xion(ik,11,ii)*siks(ii,j-1)
         enddo
         do ii=1,mar-1
            opbf=opbf+abn(12)*dens*xion(ik,12,ii)*sikar(ii,j-1)
         enddo
         do ii=1,mca-1
            opbf=opbf+abn(13)*dens*xion(ik,13,ii)*sikca(ii,j-1)
         enddo
         do ii=1,mfe-1
            opbf=opbf+abn(14)*dens*xion(ik,14,ii)*sikfe(ii,j-1)
         enddo
c!! skip ni
c         do ii=1,mni-1
c            opbf=opbf+abn(15)*dens*xion(ik,15,ii)*sikni(ii,j-1)
c         enddo
         cappa(j)=opbf+opff
 1111    continue
         if(e1(j-1).gt.7.d0) then
            jik=1
         endif
         if(fbound(j).lt.1.d-200) fbound(j)=1.d-200
         if(ftwo(j).lt.1.d-200) ftwo(j)=1.d-200
         if(cappa(j).lt.1.d-50) cappa(j)=1.d-50
         emco=free(j)+fbound(j)+ftwo(j)
         emlog=dlog10(emco/xelde)
         coolco=coolco+(e(j)-e(j-1))*(free(j)+fbound(j)+ftwo(j))
         coolff=coolff+(e(j)-e(j-1))*free(j)
         cool=cool+(e(j)-e(j-1))*emco

c add em from emiss incl. lines

c!!         em(ik,j)=emco/(4.d0*pi*xelde) + em(ik,j)

         em(ik,j)=emco/(4.d0*pi*dens**2) + em(ik,j)

         if (isnan(em(ik,j))) then
            write(6,9187)ik,j,free(j),fbound(j),ftwo(j),emco
 9187       format('ik,j,free(j),fbound(j),ftwo(j),emco ',2i5,1pe12.3,10e12.3)
         endif

c save cont. emission for rad. transfer.

c!!         emc(ik,j)=emco/(4.d0*pi*xelde) + emc(ik,j)

         emc(ik,j)=emco/(4.d0*pi*dens**2) + emc(ik,j)

c         if(j.eq.-36.or.j.eq.-28) then
c            write(6,9995)j,sqrt(e(j)*e(j-1)),free(j),fbound(j),ftwo(j),
c     &           emco/(4.d0*pi*dens**2),emc(ik,j)
 9995       format('emcon ',i5,1pe12.3,10e12.3)
c         endif
            if(em(ik,j).lt.1.d-45) then
               em(ik,j)=1.d-45
c               write(6,*)' em < 1e-45 ',ik,j,e(j),em(ik,j)
            endif
         enddo
      RETURN
      END





      double precision function ffgau(e,q,z) 
C   *****************************************
C   Karzas & Latter's free-free gaunt factors
C   (Kellogg et al. Ap.J. 199, 299 (1975)).
C   e=E(ev), q=kT(eV), z=Z(eff).
C   *****************************************
      IMPLICIT REAL*8(A-H,O-Z)                                          
C
      dimension a(6,7,3),gam2(6),gam3(6)
C
      data (gam2(j),j=1,6)/.7783d0,1.2217d0,2.6234d0,4.3766d0,
     &20.d0,70.d0/
      data (gam3(j),j=1,6)/1.d0,1.7783d0,3.d0,5.6234d0,10.d0,
     &30.d0/
      data a/1.001,1.004,1.017,1.036,1.056,1.121,1.001,1.005,
     &1.017,1.046,1.073,1.115,.9991,1.005,1.03,1.055,1.102,
     &1.176,.997,1.005,1.035,1.069,1.134,1.186,.9962,1.004,
     &1.042,1.1,1.193,1.306,.9874,.9962,1.047,1.156,1.327,
     &1.485,.9681,.9755,.8363,1.208,1.525,1.965,.3029,
     &.1616,.04757,.013,.0049,-.0032,.4905,.2155,.08357,
     &.02041,.00739,.00029,.654,.2833,.08057,.03257,.00759,
     &-.00151,1.029,.391,.1266,.05149,.01274,.00324,.9569,
     &.4891,.1764,.05914,.01407,-.00024,1.236,.7579,.326,
     &.1077,.028,.00548,1.327,1.017,1.398,.205,.0605,.00187,
     &-1.323,-.254,-.01571,-.001,-.000184,.00008,-4.762,-.3386
     &,-.03571,-.001786,-.0003,.00001,-6.349,-.4206,-.02571,
     &-.003429,-.000234,.00005,-13.231,-.59,-.04571,-.005714,
     &-.000445,-.00004,-7.672,-.6852,-.0643,-.005857,-.00042,
     &.00004,-7.143,-.9947,-.12,-.01007,-.000851,-.00004,
     &-3.175,-1.116,-.8414,-.01821,-.001729,.00023/
C
      gam=(0.01358d0*z**2)/q
c??? change compared to PL!!
      gam1=1.d3*gam
      if(gam1.gt.1.d2) goto 602
      u=e/q
      u2=u**2
C -- Born approximation gaunt factors
      u1=u/2.d0
      t=u1/3.75d0
      if(u1.gt.2.d0) goto 299
      ai=1.d0+3.5156229d0*t**2+3.089942d0*t**4+1.2067492*t**6+
     &0.23069756d0*t**8+0.0360768*t**10+0.0045813d0*t**12
      ak=-dlog(u1/2.d0)*ai-0.57721566d0+0.4227842d0*(u1/2.d0)**2+
     &0.23069756d0*(u1/2.d0)**4+0.0348859d0*(u1/2.d0)**6+0.00262698
     &d0*(u1/2.d0)**8+0.0001075*(u1/2.d0)**10+0.0000074d0*(u1/2.d0)
     &**12
      goto 297
  299 ak=1.25331414d0-0.07832358d0*(2.d0/u1)+0.02189568d0*(2.d0/u1)
     &**2-0.01062446d0*(2.d0/u1)**3+0.00587872d0*(2.d0/u1)**4-
     &0.0025154d0*(2.d0/u1)**5+0.00053208d0*(2.d0/u1)**6
      ak=ak/(dexp(u1)*dsqrt(u1))
  297 born=0.5513d0*dexp(u1)*ak
C -- Polynomial factor to multiply Born approximation
      if(gam1.lt.1.d0) g=born
      if(gam1.lt.1.d0) goto 700
      if(u.lt.0.003d0) goto 401
      if(u.le.0.03d0) N=1
      if(u.le.0.3d0.and.u.gt.0.03d0) N=2
      if(u.le.1.d0.and.u.gt.0.3d0) N=3
      if(u.le.5.d0.and.u.gt.1.d0) N=4
      if(u.le.15.d0.and.u.gt.5.d0) N=5
      if(u.gt.15.d0) N=6
      if(gam1.le.1.7783d0) M=1
      if(gam1.le.3.d0.and.gam1.gt.1.7783d0) M=2
      if(gam1.le.5.6234d0.and.gam1.gt.3.d0) M=3
      if(gam1.le.10.d0.and.gam1.gt.5.6234d0) M=4
      if(gam1.le.30.d0.and.gam1.gt.10.d0) M=5
      if(gam1.le.100.d0.and.gam1.gt.30.d0) M=6
      m1=m+1
      if(n.gt.6.or.n.le.0) then
         write(6,*)' n, u',n,e,q,z,u
      endif
      g1=born*(a(n,m,1)+a(n,m,2)*u+a(n,m,3)*u2)
      g2=born*(a(n,m1,1)+a(n,m1,2)*u+a(n,m1,3)*u2)
      p=(gam1-gam3(m))/gam2(m)
      g=(1.d0-p)*g1+p*g2
      goto 700
  602 g=1.d0
      goto 700
  401 g=born
  700 ffgau=g
C
      return 
      end
