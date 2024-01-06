      subroutine rec_rr_dr(iel,ion,te,alpha)
c interpoltes total recombination rate and direct recob. to 10 lowest levels =  al_gs      
      implicit none
      real*8 altot,algs,tetab,te,alpha,al_gs,fac
      common/gsreco/al_gs(26,26,10)
      common/adasrec/altot(26,26,20),algs(26,26,10,20),tetab(26,26,20)
      integer ind,i,j,iel,ion
      if(te <= tetab(iel,ion,1)) then
         ind=1
      elseif(te >= tetab(iel,ion,19)) then
         ind=19
      else
         do i=1,18
            if(te > tetab(iel,ion,i) .and. te <= tetab(iel,ion,i+1)) then
               ind=i
               fac=(te-tetab(iel,ion,i))/(tetab(iel,ion,i+1)-tetab(iel,ion,i))
               alpha=altot(iel,ion,i)+fac*(altot(iel,ion,i+1)-altot(iel,ion,i))
               do j=1,10
                  al_gs(iel,ion,j)=algs(iel,ion,j,i)+
     &                 fac*(algs(iel,ion,j,i+1)-algs(iel,ion,j,i))
               enddo
            endif
         enddo
      endif
      return
      end


c     subroutine recomb_adas(iel,ion,te,algs,altot)
      subroutine recomb_adas
c     calculate total recombination rate (DR and RR) from Banells compilation
c     Now only Si II-VII and S IV-IX, Ar VI-XI      .
c      
c     Note that Badnell uses the recombined ion as index and CF the recombining

      implicit none
      character*11 dum
      character*120 flrr,fldr,fl
      real*8 al(5000,19),totrec(20)
      real*8 altot,algs,tetab,te,alpha,al_gs(10),fac,enlev,e00
      common/adasrec/altot(26,26,20),algs(26,26,10,20),tetab(26,26,20)
      real*8 encont
      common/conten/encont(26,26,10)      
      integer i,j,k,m,ii,n,npart,iel,ion,npmax,irec,nlev


      do iel=1,26
         do ion=1,26
            do n=1,19
               altot(iel,ion,n)=0.
               do j=1,10
                  algs(iel,ion,j,n)=0.
               enddo
            enddo
            if(iel==8.or.iel==10.or.iel==14.or.iel==16.or.iel==18) then               
               irec=0
               do m=1,3

                  if(iel==8) then
                     npmax=1
                     if(ion==2) then
                        irec=1
                        if(m==1) then
                           fl='./ATDAT/o1ic22.dat'
                        elseif(m==2) then
                           fl='./ATDAT/o1ic23.dat'
                        elseif(m==3) then
                           fl='./ATDAT/o1ic223.dat'            
                        endif
                     elseif(ion==3) then
                        irec=1
                        if(m==1) then
                           fl='./ATDAT/o2ic22.dat'
                        elseif(m==2) then
                           fl='./ATDAT/o2ic23.dat'
                        elseif(m==3) then
                           goto 22
                        endif
                     elseif(ion==4) then
                        irec=1
                        if(m==1) then
                           fl='./ATDAT/o3ic22.dat'
                        elseif(m==2) then
                           fl='./ATDAT/o3ic23.dat'
                        elseif(m==3) then
                           goto 22
                        endif
                     elseif(ion==5) then
                        irec=1
                        npmax=0
                        if(m==1) then
                           fl='./ATDAT/o4ic12.dat'
                        elseif(m==2) then
                           fl='./ATDAT/o4ic23.dat'
                        elseif(m==3) then
                           goto 22
                        endif
                     elseif(ion==6) then
                        irec=1
                        npmax=0
                        if(m==1) then
                           fl='./ATDAT/o5ic12.dat'
                        elseif(m==2) then
                           fl='./ATDAT/o5ic23.dat'
                        elseif(m==3) then
                           goto 22
                        endif
                     elseif(ion==7) then
                        irec=1
                        npmax=0
                        if(m==1) then
                           fl='./ATDAT/o6ic12.dat'
                        elseif(m==2) then
                           fl='./ATDAT/o6ic23.dat'
                        elseif(m==3) then
                           goto 22
                        endif
                     elseif(ion==8) then
                        irec=1
                        npmax=0
                        if(m==1) then
                           fl='./ATDAT/o7ic12.dat'
                        elseif(m==2) then
                           goto 22
                        elseif(m==3) then
                           goto 22
                        endif                       
                     endif
                  elseif(iel==10) then
                     npmax=1
                     if(ion==2) then
                        irec=1
                        if(m==1) then
                           fl='./ATDAT/ne1ic22.dat'
                        elseif(m==2) then
                           fl='./ATDAT/ne1ic23.dat'
                        elseif(m==3) then
                           goto 22
                        endif
                     elseif(ion==3) then
                        irec=1
                        if(m==1) then
                           fl='./ATDAT/ne2ic22.dat'
                        elseif(m==2) then
                           fl='./ATDAT/ne2ic23.dat'
                        elseif(m==3) then
                           goto 22
                        endif
                     elseif(ion==4) then
                        irec=1
                        if(m==1) then
                           fl='./ATDAT/ne3ic22.dat'
                        elseif(m==2) then
                           fl='./ATDAT/ne3ic23.dat'
                        elseif(m==3) then
                           goto 22
                        endif
                     elseif(ion==5) then
                        irec=1
                        npmax=0
                        if(m==1) then
                           fl='./ATDAT/ne4ic22.dat'
                        elseif(m==2) then
                           fl='./ATDAT/ne4ic23.dat'
                        elseif(m==3) then
                           goto 22
                        endif
                     elseif(ion==6) then
                        irec=1
                        npmax=0
                        if(m==1) then
                           fl='./ATDAT/ne5ic22.dat'
                        elseif(m==2) then
                           fl='./ATDAT/ne5ic23.dat'
                        elseif(m==3) then
                           goto 22
                        endif
                     elseif(ion==7) then
                        irec=1
                        npmax=0
                        if(m==1) then
                           fl='./ATDAT/ne6ic12.dat'
                        elseif(m==2) then
                           fl='./ATDAT/ne6ic23.dat'
                        elseif(m==3) then
                           goto 22
                        endif
                     elseif(ion==8) then
                        irec=1
                        npmax=0
                        if(m==1) then
                           fl='./ATDAT/ne7ic12.dat'
                        elseif(m==2) then
                           fl='./ATDAT/ne7ic23.dat'
                        elseif(m==3) then
                           goto 22
                        endif                       
                     elseif(ion==9) then
                        irec=1
                        npmax=0
                        if(m==1) then
                           fl='./ATDAT/ne8ic12.dat'
                        elseif(m==2) then
                           fl='./ATDAT/ne8ic23.dat'
                        elseif(m==3) then
                           goto 22
                        endif                       
                     elseif(ion==10) then
                        irec=1
                        npmax=0
                        if(m==1) then
                           fl='./ATDAT/ne9ic12.dat'
                        elseif(m==2) then
                           goto 22
                        elseif(m==3) then
                           goto 22
                        endif                       
                     endif

                  elseif(iel==14) then
                     npmax=1
                     if(ion==2) then
                        irec=1
                        if(m==1) then
                           fl='./ATDAT/si1ic33.dat'
                        elseif(m==2) then
                           fl='./ATDAT/si1ic34.dat'
                        elseif(m==3) then
                           fl='./ATDAT/si1ic.dat'            
                        endif
                     elseif(ion==3) then
                        irec=1
                        if(m==1) then
                           fl='./ATDAT/si2ic23.dat'
                        elseif(m==2) then
                           fl='./ATDAT/si2ic334.dat'
                        elseif(m==3) then
                           fl='./ATDAT/si2ic.dat'            
                        endif
                     elseif(ion==4) then
                        irec=1
                        if(m==1) then
                           fl='./ATDAT/si3icm23.dat'
                        elseif(m==2) then
                           fl='./ATDAT/si3icm33.dat'
                        elseif(m==3) then
                           fl='./ATDAT/si3ic.dat'            
                        endif
                     elseif(ion==5) then
                        irec=1
                        npmax=0
                        if(m==1) then
                           fl='./ATDAT/si4ic23.dat'
                        elseif(m==2) then
                           fl='./ATDAT/si4ic.dat'
                        elseif(m==3) then
                           goto 22
                        endif
                     elseif(ion==6) then
                        irec=1
                        npmax=0
                        if(m==1) then
                           fl='./ATDAT/si5ic22.dat'
                        elseif(m==2) then
                           fl='./ATDAT/si5ic23.dat'
                        elseif(m==3) then
                           fl='./ATDAT/si5ic.dat'
                        endif
                     elseif(ion==7) then
                        irec=1
                        npmax=0
                        if(m==1) then
                           fl='./ATDAT/si6ic22.dat'
                        elseif(m==2) then
                           fl='./ATDAT/si6ic23.dat'
                        elseif(m==3) then
                           fl='./ATDAT/si6ic.dat'
                        endif                              

                     endif
                  elseif(iel==16) then
                     npmax=1
                     if(ion==3) then
                        irec=1
                        if(m==1) then
                           fl='./ATDAT/s2ic33.dat'
                        elseif(m==2) then
                           fl='./ATDAT/s2ic34.dat'
                        elseif(m==3) then
                           fl='./ATDAT/s2ic.dat'            
                        endif
                     elseif(ion==4) then
                        irec=1
                        if(m==1) then
                           fl='./ATDAT/s3ic33.dat'
                        elseif(m==2) then
                           fl='./ATDAT/s3ic34.dat'

                        endif
                     elseif(ion==5) then
                        irec=1
                        if(m==1) then
                           fl='./ATDAT/s4ic23.dat'
                        elseif(m==2) then
                           fl='./ATDAT/s4ic334.dat'
                        elseif(m==3) then
                           fl='./ATDAT/s4ic.dat'            
                        endif
                     elseif(ion==6) then
                        irec=1
                        if(m==1) then
                           fl='./ATDAT/s5icm23.dat'
                        elseif(m==2) then
                           fl='./ATDAT/s5icm33.dat'
                        elseif(m==3) then
                           fl='./ATDAT/s5ic.dat'            
                        endif
                     elseif(ion==7) then
                        irec=1
                        npmax=0
                        if(m==1) then
                           fl='./ATDAT/s6ic23.dat'
                        elseif(m==2) then
                           fl='./ATDAT/s6ic.dat'
                        elseif(m==3) then
                           goto 22
                        endif
                     elseif(ion==8) then
                        irec=1
                        if(m==1) then
                           fl='./ATDAT/s7ic22.dat'
                        elseif(m==2) then
                           fl='./ATDAT/s7ic23.dat'
                        elseif(m==3) then
                           fl='./ATDAT/s7ic.dat'            
                        endif
                     elseif(ion==9) then
                        irec=1
                        if(m==1) then
                           fl='./ATDAT/s8ic22.dat'
                        elseif(m==2) then
                           fl='./ATDAT/s8ic23.dat'
                        elseif(m==3) then
                           fl='./ATDAT/s8ic.dat'            
                        endif               
                     endif

                  elseif(iel==18) then
                     npmax=1
                     if(ion==5) then
                        irec=1
                        if(m==1) then                           
                           fl='./ATDAT/ar4ic33.dat'
                        elseif(m==2) then
                           fl='./ATDAT/ar4ic34.dat'
                        elseif(m==3) then
                           fl='./ATDAT/ar4ic.dat'            
                        endif
                     elseif(ion==6) then
                        irec=1
                        if(m==1) then                           
                           fl='./ATDAT/ar5ic33.dat'
                        elseif(m==2) then
                           fl='./ATDAT/ar5ic34.dat'
                        elseif(m==3) then
                           fl='./ATDAT/ar5ic.dat'            
                        endif
                     elseif(ion==7) then
                        irec=1
                        if(m==1) then
                           fl='./ATDAT/ar6ic23.dat'
                        elseif(m==2) then
                           fl='./ATDAT/ar6ic334.dat'
                        elseif(m==3) then
                           fl='./ATDAT/ar6ic.dat'            
                        endif
                     elseif(ion==8) then
                        irec=1
                        if(m==1) then
                           fl='./ATDAT/ar7icm23.dat'
                        elseif(m==2) then
                           fl='./ATDAT/ar7icm33.dat'
                        elseif(m==3) then
                           fl='./ATDAT/ar7ic.dat'            
                        endif
                     elseif(ion==9) then
                        irec=1
                        npmax=0
                        if(m==1) then
                           fl='./ATDAT/ar8ic23.dat'
                        elseif(m==2) then
                           fl='./ATDAT/ar8ic.dat'
                        elseif(m==3) then
                           goto 22
                        endif
                     elseif(ion==10) then
                        irec=1
                        if(m==1) then
                           fl='./ATDAT/ar9ic22.dat'
                        elseif(m==2) then
                           fl='./ATDAT/ar9ic23.dat'
                        elseif(m==3) then
                           fl='./ATDAT/ar9ic.dat'            
                        endif
                     elseif(ion==11) then
                        irec=1
                        if(m==1) then
                           fl='./ATDAT/ar10ic22.dat'
                        elseif(m==2) then
                           fl='./ATDAT/ar10ic23.dat'
                        elseif(m==3) then
                           fl='./ATDAT/ar10ic.dat'            
                        endif               
                     endif
                  endif
                  if(irec==1) then
c                     write(6,*)'File ',fl
                     open(11,file=fl,status='old')
                     npart=0
                     do k=1,20000
                        read(11,91,end=22)dum
 91                     format(a11)
                        if(dum(1:11) =='   INDX TE=') then
                           if(npart <= npmax) then
                              npart=npart+1
                              backspace(11)
                              read(11,92)dum,(tetab(iel,ion,i),i=1,10)
                              read(11,92)dum,(tetab(iel,ion,i),i=11,19)
 92                           format(a11,1pe10.2,9e10.2)
                              do n=1,19
                                 totrec(n)=0.
                              enddo
                              do j=1,5000
c     write(6,*)'m,npart,j ',m,npart,j
                                 read(11,*,end=11,err=11)ii,(al(j,n),n=1,10)
                                 read(11,*,end=11,err=11)(al(j,n),n=11,19)
                                 do n=1,19
                                    totrec(n)=totrec(n)+al(j,n)
                                 enddo
                              enddo
 11                           continue

c     add rad and diel recomb for each temp.                                       
                              do n=1,19
                                 altot(iel,ion,n)=totrec(n)+altot(iel,ion,n)
                                 do j=1,10
                                    algs(iel,ion,j,n)=algs(iel,ion,j,n)+al(j,n)
                                 enddo
                              enddo
c                              write(6,98)iel,ion,m,k,npart,totrec(4),totrec(7),totrec(10),altot(iel,ion,10)
 98                           format('Total rec ',5i5,1pe12.3,20e12.3)
                              goto 44
                           else
                              goto 22
                           endif
                        endif
 44                     continue
                     enddo
 22                  continue
                  endif
               enddo                                    
            endif
         enddo
      enddo
      close(11)
      open(11,file='./ATDAT/elev_rec.dat',status='old')

      do iel=1,26
         do ion=1,26
            do j=1,10
               encont(iel,ion,j)=0.
            enddo
         enddo
      enddo
            
      do i=1,1000
         read(11,*,end=88)iel,ion,nlev,e00
         do j=1,nlev
            read(11,*,end=88)enlev
            encont(iel,ion,j)=e00-enlev
c            write(6,*)'e00,enlev,encont(iel,ion,j) ',iel,ion,j,e00,enlev,encont(iel,ion,j)
         enddo
      enddo
 88   continue
      close(11)
      return
      end
      
            
      
