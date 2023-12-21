      subroutine num_shell
      implicit none
      integer i,j,iel,ion,ne,ns,ielcf,nscf(14,27)
      common/shells/ns(30,27)

      do i=1,30
         do j=1,i
            ns(i,1)=0
         enddo
      enddo
      
      ns(1,1)=1
      ns(2,1)=1
      ns(2,2)=1
      ns(6,1)=3
      ns(6,2)=3
      ns(6,3)=2
      ns(6,4)=2
      ns(6,5)=1
      ns(6,6)=1
c O      
      ns(8,1)=3
      ns(8,2)=3
      ns(8,3)=3
      ns(8,4)=3
      ns(8,5)=2
      ns(8,6)=2
      ns(8,7)=1
      ns(8,8)=1

c     Ne
      ns(10,1)=3
      ns(10,2)=3
      ns(10,3)=3
      ns(10,4)=3
      ns(10,5)=3
      ns(10,6)=3
      ns(10,7)=2
      ns(10,8)=2
      ns(10,9)=1
      ns(10,10)=1

c     Si
      ns(14,1)=5
      ns(14,2)=5
      ns(14,3)=4
      ns(14,4)=4
      ns(14,5)=3
      ns(14,6)=3
      ns(14,7)=3
      ns(14,8)=3
      ns(14,9)=3
      ns(14,10)=3
      ns(14,11)=2
      ns(14,12)=2
      ns(14,13)=1
      ns(14,14)=1

c     S
      ns(16,1)=5
      ns(16,2)=5
      ns(16,3)=5
      ns(16,4)=5
      ns(16,5)=4
      ns(16,6)=4
      ns(16,7)=3
      ns(16,8)=3
      ns(16,9)=3
      ns(16,10)=3
      ns(16,11)=3
      ns(16,12)=3
      ns(16,13)=2
      ns(16,14)=2
      ns(16,15)=1
      ns(16,16)=1
c     Ar
      ns(18,1)=5
      ns(18,2)=5
      ns(18,3)=5
      ns(18,4)=5
      ns(18,5)=5
      ns(18,6)=5
      ns(18,7)=4
      ns(18,8)=4
      ns(18,9)=3
      ns(18,10)=3
      ns(18,11)=3
      ns(18,12)=3
      ns(18,13)=3
      ns(18,14)=3
      ns(18,15)=2
      ns(18,16)=2
      ns(18,17)=1
      ns(18,18)=1
      do iel=1,26
         do i=1,iel
            ne=iel-i+1
            if(ne<=2) then
               ns(iel,i)=1
            elseif(ne<=4) then
               ns(iel,i)=2
            elseif(ne<=10) then
               ns(iel,i)=3
            elseif(ne<=12) then
               ns(iel,i)=4
            elseif(ne<=18) then
               ns(iel,i)=5
            elseif(ne<=24) then
               ns(iel,i)=6
               if(iel==20) then
                  ns(iel,i)=7
               endif
            elseif(ne<=26) then
               ns(iel,i)=7
            endif
            if(iel==18.or.iel==26.or.iel==8.or.iel==10..or.iel==20) then
c               write(6,*)iel,i,ne,ns(iel,i)
            endif
         enddo
      enddo
            

c     Fe
      ns(26,1)=7
      ns(26,2)=7
      ns(26,3)=6
      ns(26,4)=6
      ns(26,5)=6
      ns(26,6)=6
      ns(26,7)=6
      ns(26,8)=6
      ns(26,9)=5
      ns(26,10)=5
      ns(26,11)=5
      ns(26,12)=5
      ns(26,13)=5
      ns(26,14)=5
      ns(26,15)=4
      ns(26,16)=4
      ns(26,17)=3
      ns(26,18)=3
      ns(26,19)=3
      ns(26,20)=3
      ns(26,21)=3
      ns(26,22)=3
      ns(26,23)=3
      ns(26,24)=3
      ns(26,25)=1
      ns(26,26)=1

c     translate of CF element ordering

      do ielcf=1,14
         if(ielcf<=2) then
            iel=ielcf
         elseif(ielcf<=5) then
            iel=ielcf+3
         elseif(ielcf<=10) then
            iel=ielcf+4
         elseif(ielcf==11) then
            iel=16
         elseif(ielcf==12) then
            iel=18
         elseif(ielcf==13) then
            iel=20
         elseif(ielcf==14) then
            iel=26
         endif         
         do ion=1,27
            nscf(ielcf,ion)=ns(iel,ion)
         enddo
      enddo

      do ielcf=1,14
         do ion=1,27
            ns(ielcf,ion)=nscf(ielcf,ion)
         enddo
      enddo

      
      return
      end
