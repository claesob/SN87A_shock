      SUBROUTINE ATDATchianti(iel,ion,te)
      IMPLICIT REAL*8(A-H,O-Z)
      character*1 dum
      character*6 conf
      character*6 lang
      character*120 clab

      PARAMETER (NL=340,NLP1=NL+1)
      COMMON/A14/C(NL,NL),CI(NL),G(NL),E(NL),A(NL,NL),
     &     WL(NL,NL),DCDT(NL,NL),DCIDT(NL)
      COMMON/A27/OM(NL,NL)
      COMMON/NBA/NBACK
      COMMON/A16/SIG(NL),GAP(NL),BP(NL)
      COMMON/NLEV/e00,NION,NH,NP1H,NHMAX,NHMIN
      common/npoint/ifivepoint
      common/timecheck/time,itime
      common/colls/u
      common/debug/ideb
      integer n1,i2jp1(nl),i2sp1
      real*8 eijr(nl)
      TEV=Te/1.1609E4
      
      rewind (86)

      nc1=13
      nc2=20
c!!! mod 050807
      nc2 = 19
      nmax=3000
      if(iel.eq.3.and.ion.eq.5) then
         open(86,file='./ATDAT/c5_lev_coll_rad.dat',status='old')
         n=49
         ninepoint=0
      elseif(iel.eq.3.and.ion.eq.6) then
         open(86,file='./ATDAT/c6_lev_coll_rad.dat',status='old')
         n=25
         ninepoint=0
      elseif(iel.eq.4.and.ion.eq.6) then
         open(86,file='./ATDAT/n6_lev_coll_rad.dat',status='old')
         n=49
         ninepoint=0
      elseif(iel.eq.4.and.ion.eq.7) then
         open(86,file='./ATDAT/n7_lev_coll_rad.dat',status='old')
         n=25
         ninepoint=0
      elseif(iel.eq.5.and.ion.eq.7) then
         open(86,file='./ATDAT/o7_lev_coll_rad.dat',status='old')
         n=49
         ninepoint=0
      elseif(iel.eq.5.and.ion.eq.8) then
         open(86,file='./ATDAT/o8_lev_coll_rad.dat',status='old')
         n=25
         ninepoint=0
      elseif(iel.eq.6.and.ion.eq.6) then
         open(86,file='./ATDAT/ne6_lev_coll_rad.dat',status='old')
         n=180
         ninepoint=10000
      elseif(iel.eq.6.and.ion.eq.7) then
         open(86,file='./ATDAT/ne7_lev_coll_rad.dat',status='old')
         n=46
         ninepoint=38
      elseif(iel.eq.6.and.ion.eq.8) then
         open(86,file='./ATDAT/ne8_lev_coll_rad.dat',status='old')
         n=40
         ninepoint=0
      elseif(iel.eq.6.and.ion.eq.9) then
         open(86,file='./ATDAT/ne9_lev_coll_rad.dat',status='old')
         n=49
         ninepoint=0
      elseif(iel.eq.6.and.ion.eq.10) then
         open(86,file='./ATDAT/ne10_lev_coll_rad.dat',status='old')
         n=36
         ninepoint=0
      elseif(iel.eq.8.and.ion.eq.8) then
         open(86,file='./ATDAT/mg8_lev_coll_rad.dat',status='old')
         n=125
         ninepoint=0
      elseif(iel.eq.8.and.ion.eq.9) then
         open(86,file='./ATDAT/mg9_lev_coll_rad.dat',status='old')
         n=46
         ninepoint=0
         nc2=26
c!! 
         nc2=25
      elseif(iel.eq.8.and.ion.eq.10) then
         open(86,file='./ATDAT/mg10_lev_coll_rad.dat',status='old')
         n=40
         ninepoint=0
      elseif(iel.eq.8.and.ion.eq.11) then
         open(86,file='./ATDAT/mg11_lev_coll_rad.dat',status='old')
         n=49
         ninepoint=0
      elseif(iel.eq.8.and.ion.eq.12) then
         open(86,file='./ATDAT/mg12_lev_coll_rad.dat',status='old')
         n=25
         ninepoint=0
      elseif(iel.eq.10.and.ion.eq.8) then
         open(86,file='./ATDAT/si8_lev_coll_rad.dat',status='old')
         n=72
         ninepoint=0
      elseif(iel.eq.10.and.ion.eq.9) then
         open(86,file='./ATDAT/si9_lev_coll_rad.dat',status='old')
         n=46
         ninepoint=0
      elseif(iel.eq.10.and.ion.eq.10) then
         open(86,file='./ATDAT/si10_lev_coll_rad.dat',status='old')
         n=125
         ninepoint=60
      elseif(iel.eq.10.and.ion.eq.11) then
         open(86,file='./ATDAT/si11_lev_coll_rad.dat',status='old')
         n=46
         ninepoint=0
         nc2=25
c!!
         nc2=24
      elseif(iel.eq.10.and.ion.eq.12) then
         open(86,file='./ATDAT/si12_lev_coll_rad.dat',status='old')
         n=40
         ninepoint=0
      elseif(iel.eq.10.and.ion.eq.13) then
         open(86,file='./ATDAT/si13_lev_coll_rad.dat',status='old')
         n=49
         ninepoint=0
      elseif(iel.eq.10.and.ion.eq.14) then
         open(86,file='./ATDAT/si14_lev_coll_rad.dat',status='old')
         n=25
         ninepoint=0
      elseif(iel.eq.11.and.ion.eq.2) then
         open(86,file='./ATDAT/s2_lev_coll_rad.dat',status='old')
         n=43
         ninepoint=220
      elseif(iel.eq.11.and.ion.eq.3) then
         open(86,file='./ATDAT/s3_lev_coll_rad.dat',status='old')
         n=49
         ninepoint=0
      elseif(iel.eq.11.and.ion.eq.13) then
         open(86,file='./ATDAT/s13_lev_coll_rad.dat',status='old')
         n=46
         ninepoint=11
      elseif(iel.eq.11.and.ion.eq.14) then
         open(86,file='./ATDAT/s14_lev_coll_rad.dat',status='old')
         n=40
         ninepoint=0
      elseif(iel.eq.11.and.ion.eq.15) then
         open(86,file='./ATDAT/s15_lev_coll_rad.dat',status='old')
         n=49
         ninepoint=0
      elseif(iel.eq.11.and.ion.eq.16) then
         open(86,file='./ATDAT/s16_lev_coll_rad.dat',status='old')
         n=25
         ninepoint=0
      elseif(iel.eq.13.and.ion.eq.9) then
         open(86,file='./ATDAT/ca9_lev_coll_rad.dat',status='old')
         n=16
         ninepoint=0
      elseif(iel.eq.13.and.ion.eq.10) then
         open(86,file='./ATDAT/ca10_lev_coll_rad.dat',status='old')
         n=21
         ninepoint=0
      elseif(iel.eq.13.and.ion.eq.11) then
         open(86,file='./ATDAT/ca11_lev_coll_rad.dat',status='old')
         n=89
         ninepoint=0
      elseif(iel.eq.13.and.ion.eq.12) then
         open(86,file='./ATDAT/ca12_lev_coll_rad.dat',status='old')
         n=3
         ninepoint=0
      elseif(iel.eq.13.and.ion.eq.13) then
         open(86,file='./ATDAT/ca13_lev_coll_rad.dat',status='old')
         n=86
         ninepoint=0
      elseif(iel.eq.13.and.ion.eq.14) then
         open(86,file='./ATDAT/ca14_lev_coll_rad.dat',status='old')
         n=91
         ninepoint=0
      elseif(iel.eq.14.and.ion.eq.10) then
         open(86,file='./ATDAT/fe10_lev_coll_rad.dat',status='old')
         n=172
         ninepoint=672
         nc2=22
c!!
         nc2=21
      elseif(iel.eq.14.and.ion.eq.11) then
         open(86,file='./ATDAT/fe11_lev_coll_rad.dat',status='old')
         n=47
         ninepoint=0
         nc1=32
         nc2=39
c!!
         nc2=38
      elseif(iel.eq.14.and.ion.eq.12) then
         open(86,file='./ATDAT/fe12_lev_coll_rad.dat',status='old')
         n=143
         ninepoint=2222
      elseif(iel.eq.14.and.ion.eq.13) then
         open(86,file='./ATDAT/fe13_lev_coll_rad.dat',status='old')
         n=27
         ninepoint=0
      elseif(iel.eq.14.and.ion.eq.14) then
         open(86,file='./ATDAT/fe14_lev_coll_rad.dat',status='old')
         n=40
         ninepoint=0
      elseif(iel.eq.14.and.ion.eq.15) then
         open(86,file='./ATDAT/fe15_lev_coll_rad.dat',status='old')
         n=283
         ninepoint=0
      elseif(iel.eq.14.and.ion.eq.16) then
         open(86,file='./ATDAT/fe16_lev_coll_rad.dat',status='old')
         n=49
         ninepoint=8
      elseif(iel.eq.14.and.ion.eq.17) then
         open(86,file='./ATDAT/fe17_lev_coll_rad.dat',status='old')
         n=267
         ninepoint=36
      elseif(iel.eq.14.and.ion.eq.18) then
         open(86,file='./ATDAT/fe18_lev_coll_rad.dat',status='old')
         n=337
         ninepoint=212
      elseif(iel.eq.14.and.ion.eq.19) then
         open(86,file='./ATDAT/fe19_lev_coll_rad.dat',status='old')
         n=636
c note! truncated 
         nmax=200
         ninepoint=212
      elseif(iel.eq.14.and.ion.eq.20) then
         open(86,file='./ATDAT/fe20_lev_coll_rad.dat',status='old')
         n=741
c note! truncated 
         nmax=302
         ninepoint=1052
      elseif(iel.eq.14.and.ion.eq.21) then
         open(86,file='./ATDAT/fe21_lev_coll_rad.dat',status='old')
         n=620
c note! truncated 
         nmax=200
         ninepoint=598
         nc1=18
         nc2=20
      elseif(iel.eq.14.and.ion.eq.22) then
         open(86,file='./ATDAT/fe22_lev_coll_rad.dat',status='old')
         n=548
c note! truncated 
         nmax=200
         ninepoint=405
         nc1=19
         nc2=20
      elseif(iel.eq.14.and.ion.eq.23) then
         open(86,file='./ATDAT/fe23_lev_coll_rad.dat',status='old')
         n=58
         ninepoint=0
         nc1=23
         nc2=31
      elseif(iel.eq.14.and.ion.eq.24) then
         open(86,file='./ATDAT/fe24_lev_coll_rad.dat',status='old')
         n=40
         ninepoint=0
         nc1=13
         nc2=15
      elseif(iel.eq.14.and.ion.eq.25) then
         open(86,file='./ATDAT/fe25_lev_coll_rad.dat',status='old')
         n=49
         nc1=10
         nc2=16
      elseif(iel.eq.14.and.ion.eq.26) then
         open(86,file='./ATDAT/fe26_lev_coll_rad.dat',status='old')
         n=25
         nc1=9
         nc2=10
      endif

      do i=1,nl
         g(i)=0.
         e(i)=0.
         do j=1,nl
            a(j,i)=0.
            c(i,j)=0.
         enddo
      enddo

      eijcm=0.

c      write(6,*)' ionq ',ion,n

      do i=1,n
c         read(86,9646)iq,ch29,jlev,wn
 99      format(i3,a34,i3,f15.3,f15.6,f15.3,f15.6)
         eijcm_old=eijcm
         read(86,99)n1,conf,i2jp1q,eijcmo,eijro,eijcm,eijrq
         if(n1.le.nmax) then
            i2jp1(i)=i2jp1q
            if(eijro.gt.0..and.n1.ne.0) then
               eijr(i)=eijro
            elseif(eijro.le.0..and.n1.ne.0) then
               eijr(i)=eijrq

               if(eijcm.eq.eijcm_old.and.n1.gt.1) then
                  eijr(i)=eijr(i-1)+0.000001
               endif

               if(eijcm.eq.0..and.n1.gt.1) then
                  eijr(i)=eijr(i-1)+0.000001
               endif

            endif

c            write(0,99)n1,conf,
c     &           i2jp1(i),eijcm,eijr(i)
            e(i)=eijr(i)*13.5978
            e(i)=eijr(i)*13.605698

            if(iel.eq.14.and.ion.ge.19) then
               wlq=12398.54/E(I)
               wn=1.e8/wlq
c               write(6,938)i,ion,eijcmo,eijcm,e(i),eijr(i),wlq
 938           format('fe lev ',2i5,2f15.3,1pe15.6,10e15.6)
            endif


            g(i)=2.*jlev+1.
c     e(i)=wn/8065.46d0

         endif
      enddo

      do i=1,nc1
         read(86,987)  dum
         if(itime.le.2.and.iel.eq.14) then
c            write(6,987)dum
         endif
      enddo
987   format(a)

      if(dum.ne.'q') then
         write(6,*)' warning in cianti at trans! ',iel,ion,nc1
      endif

      do i=1,100000
         read(86,*,err=11,end=11)il,iu,wlq,gf,a21
c         write(0,*)il,iu,nmax
         if(iu.le.nmax) then
            a(iu,il)=a21
         endif
         
         if(iel.eq.14.and.ion.ge.18) then
            if(iu.le.nmax) then
               de=e(iu)-e(il)
               wll=12398./de
c               write(6,928)iu,il,wll,wlq,e(iu)*8065.46,e(il)*8065.46,a21
 928           format('Fe rates >18 ',2i5,2f10.3,1pe12.3,10e12.3)

            endif
         endif


      enddo

 11   continue

      do i=1,nc2
         read(86,987)  dum
         if(iel.eq.14) then
c            write(6,987)dum
         endif
      enddo
      if(dum.ne.'q') then
         write(6,*)' warning in cianti! ',iel,ion,nc2
      endif
c      write(6,*)' coll strengths',iel,ion,ninepoint
      do i=1,10000
         if(i.le.ninepoint) then
            read(86,92,end=12,err=12)iz,ionq,nlo,nu,k,gf,eijryd,cc,
     &           p1,p2,p3,p4,p5,p6,p7,p8,p9
            if(k.eq.6) k=1
            if(k.eq.7) k=2
            if(k.eq.8) k=3
            if(k.eq.9) k=4
            ifivepoint=0
            if(iel.eq.14.and.ion.ge.19) then
c               write(6,923)iz,ionq,nlo,nu,k,gf,eijryd,cc,
c     &              p1,p2,p3,p4,p5,p6,p7,p8,p9
 923           format(' Fe coll>19 ',5i3,1pe10.3,20e10.3)
c               write(0,*)' nmax ',nmax
            endif
c            if(iel.eq.10.and.ion.eq.10.and.itime.le.3) then
c               write(6,92)iz,ionq,nlo,nu,k,gf,eijryd,cc,
c     &              p1,p2,p3,p4,p5,p6,p7,p8,p9
c            endif
         else
            read(86,92,end=12,err=12)iz,ionq,nlo,nu,k,gf,eijryd,cc,
     &           p1,p2,p3,p4,p5
            if(iel.eq.14.and.ion.ge.19.and.itime.le.3) then
c               write(6,924)iz,ionq,nlo,nu,k,gf,eijryd,cc,
c     &              p1,p2,p3,p4,p5
 924           format(' Fe coll>19 5p ',5i5,1pe10.3,20e10.3)
            endif
 92         format(5i3,1pe10.3,20e10.3)
            ifivepoint=1
         endif
         if(k.eq.5.and.ideb.eq.1) then
            write(6,*)'diel: iz,ionq,nlo,nu ',iz,ionq,nlo,nu
         endif
         eijev=eijryd*13.597
         if(nlo.le.nmax.and.nu.le.nmax) then
            wi=1.0*i2jp1(nlo)
            c(nlo,nu)= ratech(k,eijryd,wi,cc,p1,p2,p3,p4,p5,p6,p7,p8,
     &           p9,te)
            eijryd=-eijryd
            wi=1.0*i2jp1(nu)
            c(nu,nlo)= ratech(k,eijryd,wi,cc,p1,p2,p3,p4,p5,p6,p7,p8,
     &           p9,te)

            if(iel.eq.14.and.(ion.ge.19)) 
     &           then
c               write(6,925)nlo,nu,i,te,u,c(nu,nlo),c(nlo,nu)
 925           format('Fe >19 OM ',3i5,1pe12.3,10e12.3)
            endif

            if(c(nu,nlo).lt.0..or.c(nlo,nu).lt.0.) then

               write(6,948)iel,ion,nlo,nu,te,c(nu,nlo),c(nlo,nu)
 948           format(' error in col str ',4i5,1pe12.3,10e12.3)
               c(nu,nlo)=0.
               c(nlo,nu)=0.
            endif

            if(i.eq.1.and.nlo.ne.1.and.nu.ne.2.and.ninepoint.eq.0) then
               write(6,*)' error in read ',iel,ion,nl,nu
            endif

            if(iel.eq.14.and.ion.ge.19) then
c               write(6,968)nu,nlo,te,c(nu,nlo),c(nlo,nu),a(nu,nlo),u
 968           format(' fe coll ',2i5,1pe12.3,10e12.3)
            endif

         endif
      enddo    
 12   continue
      close (86)

c      if(ideb.eq.1) write(6,*)' out of chianti data',iel,ion,n,nlo,nu

      if(n.gt.nmax) then
         n=nmax
      endif
      
c      write(0,*)' out of chianti data',iel,ion,n,nmax

      return

      end


      real*8 function spline(p1,p2,p3,p4,p5,p6,p7,p8,p9,xa)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 yp1,ypn,x(10),y(10),y2(10)
      PARAMETER (NMAX=500)
      INTEGER i,k
      common/npoint/ifivepoint


      if(ifivepoint.eq.1) then
         do i=1,9
            x(i)=0.
            y(i)=0.
         enddo
         x(1)=0.
         x(2)=0.25
         x(3)=0.5
         x(4)=0.75
         x(5)=1.0

         n=5

         y(1)=p1
         y(2)=p2
         y(3)=p3
         y(4)=p4
         y(5)=p5
      else

         x(1)=0.
         x(2)=0.125
         x(3)=0.25
         x(4)=0.375
         x(5)=0.5
         x(6)=0.625
         x(7)=0.75
         x(8)=0.875
         x(9)=1.0

         n=9

         y(1)=p1
         y(2)=p2
         y(3)=p3
         y(4)=p4
         y(5)=p5
         y(6)=p6
         y(7)=p7
         y(8)=p8
         y(9)=p9

      endif



      yp1=1.e31
      ypn=1.e31


      call splinec(x,y,n,yp1,ypn,y2)      

      call splint(x,y,y2,n,xa,ya)

      spline = ya

      if(spline<0.) then
c         write(6,9)ifivepoint,xa,(x(i),i=1,9),(y(i),i=1,9),y2,ya
 9       format('spline ',i5,1pe12.3,30e12.3)
      endif

      return
      end


      SUBROUTINE splint(xa,ya,y2a,n,x,y)
      implicit real*8(a-h,o-z)
      INTEGER n
      REAL*8 x,y,xa(10),y2a(10),ya(10)
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
      if (h.eq.0.) pause 'bad xa input in splint'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**
     *2)/6.
      return
      END

      SUBROUTINE splinec(x,y,n,yp1,ypn,y2)
      implicit real*8(a-h,o-z)
      INTEGER n,NMAX
      REAL*8 yp1,ypn,x(10),y(10),y2(10)
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


      real*8 function ratech(k,eij,wi,c,p1,p2,p3,p4,p5,p6,p7,p8,p9,t)
      IMPLICIT REAL*8(A-H,O-Z)
      common/colls/u
      common/upsi/y,y1
      r=1
      if(eij.le.0.0) goto 1
      r=0
      if(t.ge.1.e3*eij) r=exp(-eij*1.57888e5/t)
 1    s=1.d70
      if(t.ge.1.d-70) s=sqrt(1.57888e5/t)
      u=upsil(k,eij,c,p1,p2,p3,p4,p5,p6,p7,p8,p9,t)
      if(u.lt.0.) then
c         write(6,9)k,eij,c,wi,p1,p2,p3,p4,p5,p6,p7,p8,p9,t,u
 9       format(' ups < 0 ',i5,1pe12.3,20e12.3)
      endif
      ratech=2.1716e-8*s*r*u/wi
c      write(0,91)k,t,u,ratech,y,y1
 91   format('ratech ',i5,1pe12.3,20e12.3)
      return
      end


      real*8 function upsil(k,eij,c,p1,p2,p3,p4,p5,p6,p7,p8,p9,t)
      IMPLICIT REAL*8(A-H,O-Z)
      common/upsi/y,y1
c      write(0,*)'k,eij,c,p1,p2,p3,p4,p5,p6,p7,p8,p9,t',k,eij,c,p1,p2,p3,p4,p5,p6,p7,p8,p9,t
      e=dabs(t/(1.57888e5*eij))
c!!! if k=5?????
      if((k.eq.1).or.(k.eq.4)) x=log((e+c)/c)/log(e+c)
      if((k.eq.2).or.(k.eq.3)) x=e/(e+c)
c!!!! this is wrong!!!
      if(k.eq.5) x=log((e+c)/c)/log(e+c)

      y=spline(p1,p2,p3,p4,p5,p6,p7,p8,p9,x)
      y1=y
      if(k.eq.1) y=y*log(e+2.71828)
      if(k.eq.3) y=y/(e+1.)
      if(k.eq.4) y=y*log(e+c)
c!!!! this is wrong!!!
      if(k.eq.5) y=y*log(e+c)

      upsil=y
      if(upsil.lt.0.) then
c         write(6,91)k,eij,c,p1,p2,p3,p4,p5,p6,p7,p8,p9,t,x,y,
c     &        e,y1,upsil
      endif
      if(upsil.lt.0.) then
         upsil=0.
      endif
 91   format('ups less 0 ',i5,1pe12.3,20e12.3)
      return
      end

      real*8 function spline_old(p1,p2,p3,p4,p5,x)

c 5 point spline from Burgess and Tully
c not used since replaced by spline which is good also for 9 point

      IMPLICIT REAL*8(A-H,O-Z)
      s=1./30.
      s2=32.0*s*(19.*p1-43.*p2+30.*p3-7.*p4+p5)
      s3=160.*s*(-p1+7.*p2-12.*p3+7.*p4-p5)
      s4=32.*s*(p1-7.*p2+30.*p3-43.*p4+19.*p5)
      if(x.gt.0.25) goto 1
      x0=x-0.125
      t3=0.
      t2=0.5*s2
      t1=4.*(p2-p1)
      t0=0.5*(p1+p2)-0.015625*t2
      goto 4
 1    if(x.gt.0.5) goto 2
      x0=x-0.375
      t3=20.*s*(s3-s2)
      t2=0.25*(s2+s3)
      t1=4.*(p3-p2)-0.015625*t3
      t0=0.5*(p2+p3)-0.015625*t2
      goto 4
 2    if(x.gt.0.75) goto 3
      x0=x-0.625
      t3=20.*s*(s4-s3)
      t2=0.25*(s3+s4)
      t1=4.*(p4-p3)-0.015625*t3
      t0=0.5*(p3+p4)-0.015625*t2
      goto 4
 3    x0=x-0.875
      t3=0.
      t2=0.5*s4
      t1=4.*(p5-p4)
      t0=0.5*(p4+p5)-0.015625*t2
 4    spline=t0+x0*(t1+x0*(t2+x0*t3))
      return
      end

