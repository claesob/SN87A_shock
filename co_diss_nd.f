c reads photodiss. osc. strengths and calculates the optical depth in the
c sobolev approx. 

c photo data from Visser et al 2009

      subroutine photo_diss_b(imax,pop,tday,vmax,rshell,phr)
      implicit real*8(a-h,o-z)
      parameter(nl=1000)
      real*8 kb,nu0,linelum
      real*8 atot,fvlvu,nvlvu,jint
      common/co_elec/atot(50),fvlvu(50),nvlvu(50),jint(50),eta(50)
      character*8 char
      integer iul(50)
      common/intensities/wlext(300),linelum(300)
      data h/6.6256e-27/,c/2.997925e10/,kb/1.38e-16/
      data oe/2169.836/,oexe/13.295/,oeye/0.0115/,be/1.931285/,
     &     ae/0.017535/,ge/1.01e-5/,de/6.12e-6/,betae/1.0e-9/
      pi=3.14159d0
      pisqrt=dsqrt(pi)
c      open(12,file='CO_diss_vDB.dat')

      open(12,file='CO_photo_diss_visser_vd_blac09.dat')

c
c read dissociation data from van Dishoeck and Black
c
      read(12,902)char
 902  format(a)
      do i=1,37
c         read(12,*,end=11)no,wl,fvlvu(i),atot(i),eta,nvlvu(i),bp,dp,
c     &        ilv,ils,ill,iuv,ius,iul(i)

c New 2009 list

c         read(12,*,end=11)no,ilv,ils,ill,iuv,ius,iul(i),wl,nvlvu(i),
c     &        fvlvu(i),vp,atot(i),eta,bp,dp,ome,omxe

         read(12,*,end=11)no,ilv,ils,ill,iuv,ius,iul(i),wl,nvlvu(i),vp,
     &        fvlvu(i),atot(i),eta,bp,dp,ome,omxe

      enddo
 11   continue

      close (12)

      pie2_over_mc=0.0265

c check which lines fall inside the velocity range of a given vib. band

      do i=1,37
         wl=1./nvlvu(i)
c         write(6,*)i,nvlvu(i),wl
         wlaa=1.e8/nvlvu(i)
         hnu=h*c/(wl)
         wlmin=(1.-vmax/3.e5)*wlaa
         wlmax=(1.+vmax/3.e5)*wlaa
         dnuline=2.*(vmax/3.e5)*c/(wl)
         jint(i)=0.
         do j=1,imax
 922        format(2i5,1pe12.3,10e12.3)
            if(wlext(j).gt.wlmin.and.wlext(j).lt.wlmax) then
               jint(i)=jint(i)+linelum(j)/(16.*pi**2*rshell**2
     &              *hnu*dnuline)
c
c               didwlaa=1.e5
c               jint(i)=wl*wl*1.e8*didwlaa/c
c conversion factor from per Hz to per A
               conv=wl*wl*1.e8/c
               flq=jint(i)/conv
c               write(6,922)i,j,wlext(j),wlmin,wlmax,wl,linelum(j),
c     &              jint(i),hnu,dnuline,rshell,linelum(j),flq
            endif
         enddo
      enddo
c
c calculate photo diss. rates
c
      phdiss=0.
      write(6,*)'   i,   wl,      jint(j),     tau,',
     &'       phdissvlvu,  phdiss'

      do i=1,37
         wl=1./nvlvu(i)
c optical depth for band
         drdv=tday*8.64e4
         tau=pie2_over_mc*fvlvu(i)*wl*pop*drdv
         phdissvlvu=pie2_over_mc*fvlvu(i)*4.*pi*jint(i)*
     &        (1.-exp(-tau))/tau
         phdiss=phdiss+phdissvlvu
         write(6,912)i,wl*1.e8,jint(i),tau,phdissvlvu,phdiss
 912     format(' sob1 ',i5,f10.2,1pe12.3,10e12.3)
      enddo

      hc_over_k=h*c/kb
      wljujl=1.e8/nvlvu(n)
      aco=12.+16.
c
c doppler width in wavelengths (A)
c
      wl_dopp=wl_0*12.85e5*sqrt(t4/aco)/c
      x=dabs(wl_0-wl)/wl_dopp
c     damping constants
      gamma=atot(n)
      av=gamma*wl_0**2*1.e-8/(4.d0*pi*wl_dopp*c)
      return
      end


      subroutine photo_diss(te,den)
      implicit real*8(a-h,o-z)
      real*8 atot,fvlvu,nvlvu,jint
      common/co_elec/atot(50),fvlvu(50),nvlvu(50),jint(50),eta(50)
      parameter(nspec=91)
      common/colmol/coldnmol(nspec),coldmol(nspec),coltotmol(nspec)
      common/intphd/initphdiss
      COMMON/FRE/NINQ,jfmin,jfmax
      include 'param'
      COMMON/INT/FL(2,NE1:NE2),SI(14,27,NE1:NE2),E1(NE1:NE3),E(NE1:NE3)
      COMMON/contspec/fc(MD,NE1:NE2),flc(2,ne1:ne2)
      COMMON/IND/Ik
      common/cophoto/cophr
      real*8 nu0,nui,nur,nup,nu(10000),r(10000),el(10000)
      real*8 wlp(10000),wlr(10000)
      real*8 pop(0:40),fjlju(50,0:40,0:40),
     &     nujlju(50,0:40,0:40)
      real*8 kb
      character*9 q1i,q2i
      character*9 q1(10000),q2(10000)
      integer inu1(10000),inu2(10000),iul(50)
      character*8 char
      integer v,j,vu,vl,vp
      dimension evj(0:10,0:100)
      dimension jcont(50,0:40,0:40)
      data h/6.6256e-27/,c/2.997925e10/,kb/1.38e-16/
      data oe/2169.836/,oexe/13.295/,oeye/0.0115/,be/1.931285/,
     &     ae/0.017535/,ge/1.01e-5/,de/6.12e-6/,betae/1.0e-9/
      data eev/1.602e-12/,hpl/6.6260755e-27/
      save jcont,fjlju,nujlju
c      open(12,file='CO_diss_vDB.dat')

      ico = 26
      col_co = coldnmol(ico)

      write(6,*)'col_co ',col_co,den,te
c      write(6,*)'istat ',istat

      write(6,*)' flc  i diss ',flc(2,1)

      t4=te/1.e4

      jmax=38

      if(initphdiss.eq.0) then

         open(12,file='CO_photo_diss_visser_vd_blac09.dat')

c
c read dissociation data from van Dishoeck and Black
c

         read(12,902)char
 902     format(a)

c ground state
         bg=1.9225288
         dg=6.120e-6

         do i=1,37


c New 2009 list

            read(12,*,end=11)no,ilv,ils,ill,iuv,ius,iul(i),wl,nvlvu(i),
     &           vp,fvlvu(i),atot(i),eta(i),bp,dp,ome,omxe

c            write(6,9381)no,ilv,ils,ill,iuv,ius,iul(i),wl,nvlvu(i),
c     &           vp,fvlvu(i),atot(i),eta(i),bp,dp,ome,omxe

 9381       format(' read ',7i5,1pe12.3,e12.3,i5,10e12.3)
            
            if(iul(i).eq.0) then

               jmin = 0

            elseif(iul(i).eq.1) then

               jmin = 1

            endif

            do ju=jmin,jmax

               do k=1,3

c for sigma -> sigma (L=0->0) D J = 0 not allowed
                  if(k.eq.2.and.iul(i).eq.0) goto 55

c P branch

                  if(k.eq.1) then
                     jl=ju+1
                  endif

c Q branch

                  if(k.eq.2) then
                     jl=ju
                  endif

c R branch
                  if(k.eq.3.and.ju.ge.1) then
                     jl=ju-1
                  elseif(k.eq.3.and.ju.eq.0) then
                     goto 44
                  endif

                  nu0 = nvlvu(i)

                  nujlju(i,jl,ju)=nu0+bp*ju*(ju+1.)-dp*ju**2*(ju+1)**2-
     &                 (bg*jl*(jl+1.)-dg*jl**2*(jl+1)**2)

c                  write(6,9287)ju,jl,nu0,bp,dp,bg,dg,nujlju(i,jl,ju)
 9287             format(' nu ',2i5,1pe13.5,10e12.3)

c
c Hoenl-London factors
c

c lower state 1 Sigma i.e. L = 0

                  ll = 0

                  call hoenl_up(jl,ju,ll,iul(i),kk,sj)

                  fjlju(i,jl,ju)=fvlvu(i)*nujlju(i,jl,ju)*sj/
     &                 (nu0*(2.*jl+1.))

c correct for sum rule of delta L = 1, where sum is (2J+1)/2
c while it is 2J+1 for delta l = 0 (see Hertzberg p 208 who, however, does 
c  not say this explicitly)

                  if(iul(i).eq.1) then
                     
                     fjlju(i,jl,ju)=2. * fjlju(i,jl,ju)

                  endif

c                  write(6,913)no,ilv,jl,iuv,ju,kk,(ju-jl),
c     &                 1.e8/nujlju(i,jl,ju)
c     &                 ,sj,fjlju(i,jl,ju),
c     &                 log10(fjlju(i,jl,ju)*1.e8/nujlju(i,jl,ju))
 913              format(' fosc ',7i5,f10.2,f10.2,1pe12.3,0pf12.4)
 55               continue
               enddo
 44            continue
            enddo
         enddo
 11   continue

      endif

c
c calculate rotational populations of ground state, assuming LTE
c
      part=0.
      hc_over_k=h*c/kb
      do jl=0,jmax+1
         ej=bg*jl*(jl+1.)-dg*jl**2*(jl+1)
         pop(jl)=exp(-ej*hc_over_k/te)/(2.*jl+1.)
         part=part+pop(jl)
c         write(0,*)jl,pop(jl),part,ej*hc_over_k/te
      enddo

      do jl=0,jmax+1

         pop(jl)=pop(jl)/part

      enddo

      isob=0
      
      istat = 1

      if(isob.eq.1) then

         do n=1,37
            
            tau_tot = 0.

            if(iul(n).eq.0) then

               jmin = 0

            elseif(iul(n).eq.1) then

               jmin = 1

            endif

            do ju=jmin,jmax

               do k=1,3

c for sigma -> sigma (L=0->0) D J = 0 not allowed

                  if(k.eq.2.and.iul(n).eq.0) goto 561

c P branch
                  if(k.eq.1) then

                     jl=ju+1

                  endif
c Q branch
                  if(k.eq.2) then

                     jl=ju

                  endif
c R branch
                  if(k.eq.3.and.ju.ge.1) then

                     jl=ju-1

                  elseif(k.eq.3.and.ju.eq.0) then

                     goto 451

                  endif

                  wljujl=1.e8/nujlju(n,jl,ju)

                  wl_0=wljujl

c optical depth

                  wl_jl_ju = 1./nujlju(n,jl,ju)

                  drdv=tday*8.64e4  
                
                  tau = 0.0265 * fjlju(n,jl,ju) * wl_jl_ju * pop(jl) *
     &                 den * drdv

                  tau_tot = tau + tau_tot

                  if(n.le.3333) then

c                     write(6,944)n,iul(n),ju,ll,jl,k,wl_jl_ju*1.e8,
c     &                    fjlju(n,jl,ju),pop(jl)*den,tau,tau_tot
 944                 format(' sobol ',6i4,f12.3,1pe12.3,10e12.3)
                        
                  endif

 561              continue

               enddo

 451           continue

            enddo            

            flux=exp(-tau_tot)
            wl = 1.e8/nujlju(n,0,1)
c            write(18,9812)n,wl,opac,tau_tot,flux
c            write(6,9812)n,wl,opac,tau_tot,flux
 9812       format(' total tau ',i5,f12.3,1pe12.3,10e12.3)

         enddo

      elseif(istat.eq.1) then

c
c calculate the psimple spectrum
c      

         pd_rate = 0.

         do n=1,37

            if(iul(n).eq.0) then

               jmin = 0

            elseif(iul(n).eq.1) then

               jmin = 1

            endif

            do ju=jmin,jmax

               do k=1,3

c for sigma -> sigma (L=0->0) D J = 0 not allowed

                  if(k.eq.2.and.iul(n).eq.0) goto 56

c P branch
                  if(k.eq.1) then
                     jl=ju+1
                  endif
c Q branch
                  if(k.eq.2) then
                     jl=ju
                  endif
c R branch
                  if(k.eq.3.and.ju.ge.1) then
                     jl=ju-1
                  elseif(k.eq.3.and.ju.eq.0) then
                     goto 45
                  endif
                  wljujl=1.e8/nujlju(n,jl,ju)

                  wl_0=wljujl
                  pi=3.14159d0
                  pisqrt=dsqrt(pi)
                  aco=12.+16.
c
c     doppler width in wavelengths (A)
c
                  wl_dopp=wl_0*12.85e5*sqrt(t4/aco)/c

c damping width in cm/s ND eq 18

                  deltav =  1.e5 * (wl_0/1.e3) * (atot(n)/1.e10)/
     &                 (4.d0*pi)
                     
                  pie2_over_mc=0.0265

c ND eq 23

                  tau0 = pie2_over_mc * wl_0 * 1.e-8 * pop(jl) * 
     &                 col_co * fjlju(n,jl,ju)/ deltav
                  
c approx of ND89 
                  self_abs = 1./(2./3. + sqrt(tau0))

                  if(initphdiss.eq.0) then                     

                     enline = hpl * 3.e10 /(wl_0 * 1.e-8 * eev)

                     do j=jfmin,jfmax


                        if(enline.gt.e(j).and.enline.le.e(j+1)) then

                           jcont(n,jl,ju) = j

                        endif

                     enddo

                  endif
                  
                  jc = jcont(n,jl,ju)

c flc in erg/cm2 s ster eV, e1 in eV, sint in 1/cm2 s ster Hz                  

                  sint = (hpl /eev) * flc(2,jc)/(e1(jc)*eev)
                  
                  pd_rate = eta(n) * pie2_over_mc * fjlju(n,jl,ju) * 
     &                 pop(jl) * self_abs * sint + pd_rate
                     
                  if(n.le.399.or.n.ge.37) then
                     write(6,933)n,ju,jl,k,jc,e1(jc),wl_0,eta(n),
     &                    fjlju(n,jl,ju),pop(jl),tau0,sint,pd_rate
 933                 format(' pd-rate ',5i4,2f12.3,1pe12.3,10e12.3)
                  endif

 56            continue

               enddo
               
 45            continue

            enddo

         enddo

      endif

      cophr = pd_rate

      write(6,*)' Total pd-rate ',pd_rate

      do i=1,37

         wl=1./nvlvu(i)

c optical depth for band

         drdv=tday*8.64e4

         pie2_over_mc=0.0265
            
         tau=pie2_over_mc * fvlvu(i) * wl * den * drdv

c         write(6,9288)i,wl*1.e8,fvlvu(i),den,tau
 9288    format(' sob ',i5,1pe12.3,10e12.3)

      enddo

      initphdiss = 1

      end





      subroutine hoenl_low(jl,ju,ll,lu,k,s)

c This routine gives the Hoenl-London factor for a given LOWER j = jup
c See Hertzberg 1950 IV.81-83 p 208

c delta j = j_upper - j_lower

c R delta j = +1

c Q delta j = 0

c P delta j = -1

      implicit none

      real*8 sr,sp,sq,s

      integer delta_j,ju,jl,k
      integer delta_l,lu,ll

      delta_j = ju - jl

      delta_l = lu - ll


c R branch

      if(delta_j.eq.1) then

         k = 1

c Q branch

      elseif(delta_j.eq.0) then

         k = 2

c P branch

      elseif(delta_j.eq.-1) then

         k = 3

      endif

      if(delta_l.eq.0) then

         sr = real((jl + 1 + ll)*(jl + 1 - ll))/real(jl + 1)

         sq = real((2 * jl + 1)*ll**2)/real(jl *(jl + 1))

         sp = real((jl + ll)*(jl - ll))/real(jl)

      elseif(delta_l.eq.1) then

         if(ll.eq.0) then

            sr = real((jl + 2))/real(4)

            sq = real((2 * jl +1))/
     &           real(4)

            sp = real((jl -1))/real(4)

         else

            sr = real((jl + 2 + ll)*(jl + 1 + ll))/real(4 * (jl + 1))

            sq = real((jl + 1 + ll)*(jl - ll)*(2 * jl +1))/
     &           real(4 * jl *(jl + 1))
            
            sp = real((jl -1 - ll)*(jl - ll))/real(4 *jl)
            
         endif

      elseif(delta_l.eq.-1) then

         sr = real((jl + 2 - ll)*(jl + 1 - ll))/real(4 * (jl + 1))

         sq = real((jl + 1 - ll)*(jl + ll)*(2 * jl +1))/
     &        real(4 * jl *(jl + 1))

         sp = real((jl -1 + ll)*(jl + ll))/real(4 *jl)

      endif

      if(k.eq.1) then

         s = sr

      elseif(k.eq.2) then

         s = sq

      elseif(k.eq.3) then

         s = sp

      endif

c      write(6,9)k,ju,jl,delta_j,delta_l,s
c      write(0,9)k,ju,jl,delta_j,delta_l,s
 9    format('Hoenl factor ',5i5,1pe12.3)

      return

      end

      subroutine hoenl_up(jl,ju,ll,lu,k,s)

c This routine gives the Hoenl-London factor for a given UPPER j = jup
c See Hertzberg 1950 IV.81-83 p 208

c delta j = j_upper - j_lower

c R delta j = +1

c Q delta j = 0

c P delta j = -1

      real*8 sr,sp,sq,s

      integer delta_j,ju,jl,k
      integer delta_l,lu,ll

      character*1 lab

      delta_j = ju - jl

      delta_l = lu - ll


c R branch

      if(delta_j.eq.1) then

         k = 1

c Q branch

      elseif(delta_j.eq.0) then

         k = 2

c P branch

      elseif(delta_j.eq.-1) then

         k = 3

      endif

      if(delta_l.eq.0) then

         sr = real((ju + lu)*(ju - lu))/real(ju)

         sq = real((2 * ju + 1)*lu**2)/real(ju *(ju + 1))

         sp = real((ju + 1 + lu)*(ju + 1 - lu))/real(ju + 1)

      elseif(delta_l.eq.1) then

            sr = real((ju + lu))/real(4 * ju)

            sq = real((ju + lu)* (ju + 1 -lu) * (2 * ju +1))/
     &           real(4 * ju * (ju + 1))

            sp = real((ju + 1 - lu) * (ju + 2 - lu))/real(4 * (ju + 1))

      elseif(delta_l.eq.-1) then

         sr = real((ju - lu)*(ju - 1 - lu))/real(4 * ju)

         sq = real((ju - lu)*(ju + 1 + lu)*(2 * ju +1))/
     &        real(4 * ju *(ju + 1))

         sp = real((ju + 1 + lu)*(ju + 2 + lu))/real(4 *(ju + 1))

      endif

      if(k.eq.1) then

         s = sr

         lab = 'R'

      elseif(k.eq.2) then

         s = sq

         lab = 'Q'

      elseif(k.eq.3) then

         s = sp

         lab = 'P'

      endif

c      write(6,9)lab,jl,ju,delta_j,delta_l,s
c      write(0,9)lab,jl,ju,delta_j,delta_l,s
 9    format('Hoenl factor ',a1,4i3,1pe12.3)

      return

      end








