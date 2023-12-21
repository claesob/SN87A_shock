      subroutine auger_fr
c Note thta iel is the atomic number, ielcf the CF numbering used for fr_aug()
      implicit none
      integer nz,nion,nshell,i,k,n
      parameter(nz=30,nion=27,nshell=10)
      integer ieln(30),iz,iel,ielcf,ion,shell,epsi
      integer ns,kmax,init_augfrac
      real*8 fr_aug,eion,en_aug,en_augi,eioni,fr(10)
      data ieln/1,2,0,0,0,3,4,5,0,6,7,8,9,10,0,11,0,12,0,13,0,0,0,0,0,
     &     14,0,0,0,0/
      common/augerfrac/eion(nz,nion,nshell),en_aug(nz,nion,nshell),
     &     fr_aug(nz,nion,nshell,10),kmax(nz,nion,nshell),ns(nz,nion),
     &     init_augfrac
      write(6,*)' reading auger data',init_augfrac
      if(init_augfrac==1) then
         open(86,file='./ATDAT/auger_frac_cf.dat',status='old')
         write(6,*)' reading auger data'
         do i=1,10000
            read(86,*,end=99)ielcf,ion,shell,(fr(k),k=1,7)
            ns(ielcf,ion)=shell
            do k=1,7
               fr_aug(ielcf,ion,shell,k)=0.
               if(fr(k)>0 .or. k==1) then
                  fr_aug(ielcf,ion,shell,k)=fr(k)
c$$$c only for test without auger!                  
c$$$                  if(k==1) then
c$$$                     fr_aug(ielcf,ion,shell,k)=1.0
c$$$                  else
c$$$                     fr_aug(ielcf,ion,shell,k)=0.
c$$$                  endif
                  kmax(ielcf,ion,shell)=k
c$$$                  write(6,88)i,ielcf,ion,k,shell,fr_aug(ielcf,ion,shell,k)
c$$$ 88               format('i,ielcf,ion,k,shell,fr ',5i4,f10.6)
               endif
            enddo
         enddo
 99      continue
         close(86)
         init_augfrac=0
         
c$$$         do iel=1,27
c$$$            do ion=1,iel
c$$$               if(ns(iel,ion) >= 1) then
c$$$                  do n=1,ns(iel,ion)
c$$$c                     write(6,*)'iel,ion,shell,kmax(iel,ion,n) ',iel,ion,
c$$$c     &                    n,kmax(iel,ion,n)
c$$$                  enddo
c$$$               else
c$$$                  ns(iel,ion)=1
c$$$                  kmax(iel,ion,1)=1
c$$$                  fr_aug(iel,ion,1,1)=1.d0
c$$$c                  write(6,*)'iel,ion,shell,ns(iel,ion),kmax(iel,ion,n)',
c$$$c     &                    iel,ion,ns(iel,ion),kmax(iel,ion,ns(iel,ion))
c$$$               endif
c$$$            enddo
c$$$         enddo

      endif
c$$$      ielcf=10
c$$$      do ion=1,14
c$$$         write(6,*)'ns(ielcf,ion)',ielcf,ion,ns(ielcf,ion)
c$$$         do i=1,ns(ielcf,ion)
c$$$            write(6,*)'kmax(ielcf,ion,i) ',kmax(ielcf,ion,i)
c$$$         enddo
c$$$      enddo
      return
      end
      
               
         
