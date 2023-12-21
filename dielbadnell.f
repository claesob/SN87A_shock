      subroutine badnell_et_al
c initiates and reads data for the Badnell et al fits      
      implicit none
      integer initdiel,initrr
      common/init_diel/initdiel,initrr
      integer i,ion,iel
      real*8 te,diel,rrec
      te=1.e4
      initdiel=1
c for initial reading iel & ion arbitrary!      
      iel=12
      ion=6
      call diel_fits_badnell(11,iel,ion,te,diel)
      call diel_fits_badnell(12,iel,ion,te,diel)
      call diel_fits_badnell(13,iel,ion,te,diel)
      call diel_fits_badnell(14,iel,ion,te,diel)      
      call diel_fits_badnell(15,iel,ion,te,diel)
      initdiel=0
      initrr=1
      
      call rr_fits_badnell(iel,ion,te,rrec)
      initrr=0
      
      do i=1,10
         call rr_diel_badnell(te)
      enddo
      return
      end
      
      subroutine rr_diel_badnell(te)
      implicit none
      real*8 te,drec,dielbadn,rrec,alrecbadn,altotbadn
      integer iel,ion
      common/rr_diel_badn/dielbadn(30,30),alrecbadn(30,30)
     &     ,altotbadn(30,30)
c S II -> S I
      call diel_fits_badnell(15,16,2,te,drec)
      dielbadn(16,2)=drec
c Ar IV -> Art III
      call diel_fits_badnell(15,18,4,te,drec)
      dielbadn(18,4)=drec
c Ca VI -> Ca V
      call diel_fits_badnell(15,20,6,te,drec)      
      dielbadn(20,6)=drec
c S III -> S II
      call diel_fits_badnell(14,16,3,te,drec)
      dielbadn(16,3)=drec
c Ar V -> Ar IV
      call diel_fits_badnell(14,18,5,te,drec)
      dielbadn(18,5)=drec
c Ca VII -> Ca VI
      call diel_fits_badnell(14,20,7,te,drec)
      dielbadn(20,7)=drec
c Si II -> Si I      
      call diel_fits_badnell(13,14,2,te,drec)
      dielbadn(14,2)=drec
c S IV -> S III      
      call diel_fits_badnell(13,16,4,te,drec)
      dielbadn(16,4)=drec
c Ar VI -> Ar V
      call diel_fits_badnell(13,18,6,te,drec)
      dielbadn(18,6)=drec
Ca VIII -> Ca VII
      call diel_fits_badnell(13,20,8,te,drec)               
      dielbadn(20,8)=drec
c Si III -> Si II      
      call diel_fits_badnell(12,14,3,te,drec)
      dielbadn(14,3)=drec
c S V -> S IV      
      call diel_fits_badnell(12,16,5,te,drec)
      dielbadn(16,5)=drec
c Ar VII -> Ar VI
      call diel_fits_badnell(12,18,7,te,drec)
      dielbadn(18,7)=drec
c Ca IX -> Ca VIII
      call diel_fits_badnell(12,20,9,te,drec)                  
      dielbadn(20,9)=drec
c Si IV -> Si III      
      call diel_fits_badnell(11,14,4,te,drec)
      dielbadn(14,4)=drec
c S VI -> S V      
      call diel_fits_badnell(11,16,6,te,drec)
      dielbadn(16,6)=drec
c Ar VIII -> Ar VII
      call diel_fits_badnell(11,18,8,te,drec)
      dielbadn(18,8)=drec
Ca X -> Ca IX
      call diel_fits_badnell(11,20,10,te,drec)
      dielbadn(20,10)=drec
      do   iel=14,20,2
         do ion=2,iel+1
            call rr_fits_badnell(iel,ion,te,rrec)
            if(rrec > 0.) then
               alrecbadn(iel,ion)=rrec
               altotbadn(iel,ion)=alrecbadn(iel,ion)+dielbadn(iel,ion)
            endif
         enddo
      enddo
      return
      end
      
      subroutine diel_fits_badnell(iso,iel,ion,te,drec)
      implicit none
      real*8 c(30,30,8),e(30,30,8),drec,te
      character*5 ionel,dum
      integer i,k,iel,ion,iso,initdiel,initrr,nel,kmax(30),npar(30,30)
      common/init_diel/initdiel,initrr
      save kmax,c,e,npar
      if(initdiel==1) then
         if(iso==11) then
            open(11,file='./ATDAT/diel_na_like_althun2006.txt',status='old')
            kmax(iso)=7
            nel=17
         elseif(iso==12) then
            open(11,file='./ATDAT/diel_mg_like_althun2007.txt',status='old')
            kmax(iso)=6
            nel=14
         elseif(iso==13) then
            open(11,file='./ATDAT/diel_al_like_abdel_naby2012.txt',status='old')
            kmax(iso)=5
            nel=13
c            open(12,file='rr_al_like_abdel_naby2012.txt',status='old')
         elseif(iso==14) then
            open(11,file='./ATDAT/diel_si_like_kaur2018.txt',status='old')
            kmax(iso)=8
            nel=14
         elseif(iso==15) then
            open(11,file='./ATDAT/diel_p_like_bleda2022.txt',status='old')
            kmax(iso)=7
            nel=13
         endif
         read(11,*)dum
         do i=1,nel
            iel=iso+i
            ion=iel-iso+1
            npar(iel,ion)=kmax(iso)
            read(11,*,end=99)ionel,(c(iel,ion,k),k=1,kmax(iso))
         enddo
         read(11,*)dum
         do i=1,nel
            iel=iso+i
            ion=iel-iso+1
            read(11,*,end=99)ionel,(e(iel,ion,k),k=1,kmax(iso))
         enddo
 99      continue
      else
         drec=0.
         do k=1,npar(iel,ion)
            drec=drec + c(iel,ion,k)*exp(-e(iel,ion,k)/te)
         enddo
         drec=drec/te**1.5
      endif
      return
      end
      
      subroutine rr_fits_badnell(iel,ion,te,rrec)
c     read radiative rec rates from Badnell 2006 fits
c ion = ionization stage of recombining ion, e.g., H II = 2      
      implicit none
      integer i,iel,ion,initdiel,initrr
      real*8 a(30,30),b(30,30),t0(30,30),t1(30,30),c(30,30),
     &     t2(30,30),te,rrec
      real*8  aa,bb,t00,t11,cc,t22
      integer z,nel,iz0,k
      character*5 ionel,dum
      common/init_diel/initdiel,initrr
      real*8 qq0,qq1,qq2
      save a,b,t0,t1,c,t2
      if(initrr==1) then
         open(11,file='./ATDAT/rr_rec_badnell2006.dat',status='old')
         do i=1,19
            read(11,*)dum
         enddo
         do i=1,100000
            read(11,*,end=99)z,nel,aa,bb,t00,t11,cc,t22
            if(z>=3.and.z<=26) then
               a(z,z-nel+1)=aa
               b(z,z-nel+1)=bb
               t0(z,z-nel+1)=t00
               t1(z,z-nel+1)=t11
               c(z,z-nel+1)=cc
               t2(z,z-nel+1)=t22
            endif
         enddo
 99      close(11)
         do k=1,4
            if(k==1) then
               open(11,file='./ATDAT/rr_rec_p_like_bleda2022.dat',status='old')
               iz0=14
            elseif(k==2) then
               open(11,file='./ATDAT/rr_rec_al_like_abdel_naby2012.dat',
     &              status='old')
               iz0=13
            elseif(k==3) then
               open(11,file='./ATDAT/rr_rec_si_like_kaur2018.dat',
     &              status='old')
               iz0=15
            elseif(k==4) then
               open(11,file='./ATDAT/rr_rec_mg_like_althun2007.dat',
     &              status='old')
               iz0=12               
            endif
            do i=1,1
               read(11,*)dum
            enddo
            do i=1,13
               read(11,*)dum,aa,bb,t00,t11,cc,t22
               iel=iz0+i
               ion=i+1
               if(iel < 30) then
                  a(iel,ion)=aa
                  b(iel,ion)=bb
                  t0(iel,ion)=t00
                  t1(iel,ion)=t11
                  c(iel,ion)=cc
                  t2(iel,ion)=t22
               endif
            enddo
            close(11)
         enddo         
      endif
      if(c(iel,ion) > 0.d0) then
         bb=b(iel,ion) + c(iel,ion)*exp(-t2(iel,ion)/te)
      else
         bb=b(iel,ion)
      endif

      qq0= a(iel,ion) * (t0(iel,ion)/te)**0.5      
      qq1= ( 1.+(te/t0(iel,ion))**0.5 )**(1.-bb) 
      qq2=    ( 1.+(te/t1(iel,ion))**0.5 )**(1.+bb)
      rrec=qq0/(qq1*qq2)

      return
      end


      
