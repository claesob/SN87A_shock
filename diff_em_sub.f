      SUBROUTINE diffuse_sph_aan(imode,lambda,IMAX,ishock,te,totfl)
      IMPLICIT REAL*8(A-H,O-Z)
      integer step
      include 'param'
      PARAMETER (MPAR=1500,MPP1=1501)
      real*8 mu,mumin
      dimension f(100001)
      parameter(nang=21)
      COMMON/FRE/NINT,JMIN,JJ
      common/frek/jf
      COMMON /CTRAN/X(MPAR),S(MPAR),XJ(MPAR),XH(MPAR),XK(MPAR)
     & ,FJ(MPAR),SOURCE(MPAR),JTAU0,JTAU1
      COMMON /TAUCQ/TAUREF(MPAR),TT(MPP1,MPAR),DTAULN(MPAR),JTAU
      COMMON /CSPHER/NCORE,NATMOS,TAUM,RR(MPAR)

      common/opacity/TAop(NE1:NE2),Ss(MD,NE1:NE2),copac(md,ne1:ne2)
      COMMON/DIF/EM(MD,NE1:NE2),TAU(MD,NE1:NE2),TAUTOT(MD,NE1:NE2),
     &           EMC(MD,NE1:NE2)
       common/trq/rtr(2,md),drtr(2,md),dentr(2,md)
      COMMON/INT/FL(2,NE1:NE2),SI(14,27,NE1:NE2),E1(NE1:NE3),E(NE1:NE3)
      COMMON/SPECT/TEL,FD(MD,NE1:NE2),F0(NE1:NE2),IPAR
      COMMON/contspec/fc(MD,NE1:NE2),flc(2,ne1:ne2)
      common/taushok/taush(md)
      dimension r(md),denm(md)
      dimension flj(md),flh(md),flmu(nang),opac(md),emm(md)
      dimension specj(md,ne1:ne2),spech(md,ne1:ne2)
      data pi/3.141592/

c      write(6,*)' imax ',imax,md,lambda


C     I     1         2         3         4                            *
C           +         +         +         +             +         +   
C          R(1)      R(2)      R(3)      R(4)        R(imax-1)  R(imax)         
C         taut(k)   taut(2)  taut(3)    taut(4)    taut(imax-1) taut(imax)       
C          EM(1)    EM(2)     EM(3)      EM(4)      EM(imax-1)  EM(imax)    
C        copac(1)   copac(2)  copac(3)  copac(4)   copac(imax-1) copac(imax) 

C                    J(2)      J(3)       J(4)       j(imax-1)   j(imax)  
C                                                                      
C                +         +          +       +     +         +   
C              DX(2)     DX(3)      DX(4)       DX(imax-1)  DX(imax)     
c              emm(2)    emm(3)     emm(4)      emm(imax-1) emm(imax)
c             opac(2)    opac(3)    opac(4)     opac(imax-1 opac(imax)
C               S(2)      S(3)       S(4)       s(imax-1)   s(imax)

C             tmid(1)    tmid(3)    tmid(i)                            *
C                                                                      *

c  We want mean intensities at the grid points for the calc. of ionization and 
c     temperature at each point. To calc. this we need emissivities and opacities 
c     at the midpoints. This gives the constribution from each shell to the mean int. 

c first calc. emissivities, densities  and opacities at midpoints

c note that emm(i+1), etc is evaluated at the midpoint between i and i+1


c mu = 1

c dl = r(k+1) - r(k) 

c s(k), opac(k) denm(k) alla borde vara k+1!


c Note!

c After second and higher iteration input has i=1 for preshock boundary and i=imax at 
c     postshock boundary. After first iteratio i=1 at post-shock boundary and i=imax 
c     at preshock. 

c this program calculates the radiation transport from outer boundary (i=1) to 
c     inner (i=imax)

c we want the postshock gas to be at the outer boundary and preshock in center. 
c     After the second and higher iteration we therefore need to renumber from 
c     inside out
c   



c first renumber so that i=1 is at the outer boundary and i = imax closest to the center.

c  from first sol. (lambda = 2)


c     i=1              ishock           imax

c  cool shell          shock           preshock


c 2nd and higher


c     i=1              ishock           imax

c  pre-shock          shock           cool shell



c required for transfer solution

c        i=1     2                     imax

c   outer boundary      shock       inner boundary

c      post-shock                      pre-shock  



c radius of first shell

c!!      rmin = 1.0d16
      rmin = 1.0d17

      r(imax)=rmin

      do m=imax,2,-1

         if(lambda.eq.2) then

            r(m-1) = r(m) + dabs(drtr(2,m))

c            write(6,9238)m-1,drtr(2,m),r(m),r(m-1)
 9238       format(' rr ',i5,1pe25.15,10e25.15)

         else

            r(m-1) = r(m) + dabs(drtr(2,imax-m+2))

c            write(6,9238)m-1,drtr(2,imax-m+2),r(m),r(m-1)


         endif

         rr(m) = r(m)

      enddo 

c!! 

      rr(1) = rr(2) + drtr(2,2)


c density at mid-points. for lambda > 2 r and denm numbered from 
c        outside to inside 

c r    surface                center

c      r(1)        r(2)       r(imax)     r(imax+1)
c          denm(1)                 denm(imax)


c  
c         drtr(imax)               drtr(1)

      do i=1,imax-1

         if(lambda.eq.2) then

            denm(i+1) = (dentr(2,i) + dentr(2,i+1))/2.

         else

            denm(imax-i+1) = (dentr(2,i) + dentr(2,i+1))/2.

         endif
            
      enddo

c add a point inside of rmin

      if(lambda.ge.3) then

         r(imax+1) = r(imax) - dabs(drtr(2,2))
         
         denm(imax+1) = denm(imax)

         write(6,*)' last shell added',imax+1,drtr(2,2),r(imax+1),
     &        denm(imax+1)
         
      endif


c optical depth scale (ie column density)

c density at mid-points. for lambda > 2 r and denm numbered from 
c        outside to inside 

c r    surface                center

c      r(1)        r(2)       r(imax)        r(imax+1)

c          denm(1)                  denm(imax)


c  
c         drtr(imax)                 drtr(1)

c    tauref(1)=0  tauref(2)  tauref(imax)  tauref(imax+1)

      tauref(1) = 0.

      do i=2,imax

         if(lambda.eq.2) then

            if(i.ge.2) then
               tauref(i) = abs(drtr(2,i))*denm(i) + tauref(i-1)
            else
               tauref(i) = abs(drtr(2,i))*denm(i)
            endif

         else

            if(i.ge.2) then
               tauref(i) = abs(drtr(2,imax-i+2))*denm(i) + tauref(i-1)
            else
               tauref(i) = abs(drtr(2,imax-i+2))*denm(i)
            endif

         endif

      enddo

      jtau = imax

      jfq = 3

c      write(6,*)' parameters for spectrum imode',imode,ishock,jfq,imax

      do i=1,imax+1

         sjfq = em(i,jfq)*dentr(2,i)/copac(i,jfq)
         sjfqc = emc(i,jfq)*dentr(2,i)/copac(i,jfq)
         s10 = em(i,-10)*dentr(2,i)/copac(i,-10)
c         write(6,9288)i,rr(i),tauref(i),drtr(2,i),denm(i),
c     &        em(i,jfq),emc(i,jfq),copac(i,jfq),em(i,3),copac(i,3),
c     &        s10,sjfq,sjfqc
c         write(82,*)i,rr(i),tauref(i),drtr(2,i),denm(i),
c     &        em(i,jfq),copac(i,jfq),em(i,3),copac(i,3),
c     &        s10,s3
 9288    format(' inp ',i5,1pe22.14,1e22.14,20e12.3)

      enddo


c      write(6,*)' parameters for spectrum imode',imode,ishock,jfq

      do i=1,imax
         
         if(i.le.3.or.i.ge.imax-3) then

c            write(6,9288)i,rr(i),tauref(i),drtr(2,i),denm(i)

            sjfq = em(i,jfq)*dentr(2,i)/copac(i,jfq)
            sjfqc = emc(i,jfq)*dentr(2,i)/copac(i,jfq)
            s10 = em(i,-10)*dentr(2,i)/copac(i,-10)
c            write(6,9283)(em(i,jfx),jfx=jmin,jj)
 9283       format('em ',1pe12.3,9e12.3)
c            write(6,9284)(emc(i,jfx),jfx=jmin,jj)
 9284       format('emc ',1pe12.3,9e12.3)
c            write(6,9285)(copac(i,jfx),jfx=jmin,jj)
 9285       format('opa ',1pe12.3,9e12.3)
         endif
      enddo


      do j=jmin,jj

         jf = j


c replace the emissivity of the first bin with that of the second

         emc(1,j) = emc(2,j)

c!! mod 070807

         em(1,j) = em(2,j)


c replace the emissivity of the first bin with that of the second

         emc(imax,j) = emc(imax-1,j)

c!! mod 070807

         em(imax,j) = em(imax-1,j)


         do i=1,imax-1

            if(lambda.eq.2) then

               if(imode.eq.2) then
                  emm(i) = em(i,j) 
               elseif(imode.eq.3) then
                  emm(i) = emc(i,j) 
               endif
            
               opac(i) = copac(i,j)

               source(i) = emm(i)*dentr(2,i)/opac(i)
            
            else

               if(imode.eq.2) then
                  emm(imax-i+1) = em(i,j) 
               elseif(imode.eq.3) then
                  emm(imax-i+1) = emc(i,j)
               endif
            
               opac(imax-i+1) = copac(i,j)

               source(imax-i+1) = emm(imax-i+1)*dentr(2,i)/
     &              opac(imax-i+1)

            endif

         enddo

c now calculate the transfer with i=1 at outer bondary and imax at center

         do i=1,imax-1
               
            if(imode.eq.2) then
               FD(I,J)=0.0E0
            elseif(imode.eq.3) then
               fc(I,J)=0.0E0
            endif
 
         enddo
         
         if(opac(1).eq.0.) then
            opac(1) = opac(3)
         endif

         emm(1) = emm(2)

c!!!         emm(2) = emm(3)

         opac(1) = opac(2)

c         opac(1) = opac(3)

         emm(imax) = emm(imax-1)

         opac(imax) = opac(imax-1)

         do kk=1,2

            if(kk.eq.1) then

               i=1

            elseif(kk.eq.2) then

               i = imax

            endif

            if(lambda.eq.2) then

               source(i) = emm(i)*dentr(2,i)/opac(i)

            else

               source(imax-i+1) = emm(imax-i+1)*dentr(2,i)/
     &             opac(imax-i+1)

            endif

         enddo

         do i=1,imax

            x(i) = opac(i)

            s(i) = 0.

         enddo

         if(j.eq.jfq) then

c            write(6,*)' parameters for spectrum imax',j,imax,imode
c
            do i=1,imax

c            if(i.le.3.or.i.ge.imax-3) then
c               write(6,9297)j,i,imode,rr(i),dentr(2,i),tauref(i),x(i),
c     &              emm(i),source(i)
 9297          format('par ',3i5,1pe15.7,10e12.3)
               
c            endif

            enddo

         endif




c calculate optical depth at lyman edge from shock

         if(j.eq.3) then

            do k=1,2

               if(lambda.eq.2) then

                  if(k.eq.1) then

c  i < ishock
            
                     i1 = ishock-1

                     i2 = 1

                     step = -1

                  elseif(k.eq.2) then
            
c  i > ishock
                     i1 = ishock+1

                     i2 = imax
                  
                     step = 1
            
                  endif

                  taui = 0.

                  do i=i1,i2,step
                  
                     if(k.eq.1) then

                        dl = abs(r(i)-r(i+1))

                     else

                        dl = abs(r(i)-r(i-1))

                     endif

                     dtaui = dl*opac(i)*denm(i)

                     taui = taui + dtaui

c     taush(i) = optical depth from shock to r(i)

                     taush(i) = taui

                  enddo

               else

                  if(k.eq.1) then

c  i < ishock
            
                     i1 = imax-ishock

                     i2 = imax

                     step = -1

                  elseif(k.eq.2) then
            
c  i > ishock
                     i1 = imax-ishock

                     i2 = 1
                  
                     step = 1
            
                  endif

                  if(k.eq.1) then

c  i < ishock
            
                     i1 = ishock-1

                     i2 = 1

                     step = -1

                  elseif(k.eq.2) then
            
c  i > ishock
                     i1 = ishock+1

                     i2 = imax
                  
                     step = 1
            
                  endif

                  taui = 0.

                  do i=i1,i2,step
                  
                     if(k.eq.1) then

                        dl = abs(r(imax-i+1)-r(imax-i))

                     else

                        dl = abs(r(imax-i+1)-r(imax-i+2))

                     endif

                     dtaui = dl*opac(imax-i+1)*denm(imax-i+1)

                     taui = taui + dtaui

c taush(i) = optical depth from shock to r(i)

                     taush(i) = taui

                  enddo

               endif
                     
            enddo

            taush(ishock) = 0.

            do i=1,imax

               if(lambda.eq.2) then

                  write(6,922)i,r(i),opac(i),denm(i),taush(i)

               else

                  write(6,922)i,opac(imax-i+1),denm(imax-i+1),taush(i)
 922              format(i5,1pe12.3,10e12.3)
               endif

            enddo

         endif

         r15=1.e-15
c         call SPHTRP(IMAX,J,R15,TEINT,RM,IMAXP)

         call traneq

         do i=1,imax
            fd(i,j) = xj(i)
         enddo

         if(j.eq.3.or.j.eq.-10) then
c            write(6,*)' spec for j=',j
            do i=1,imax
c               write(6,9278)i,rr(i),xj(i),xh(i),fd(i,j)
 9278          format(' spec ',i5,1pe16.8,10e12.3)
            enddo

         endif

      enddo

c      write(57,*)' new spectrum for lambda = ',lambda,ishock
      do i=2,imax,10
         do j=jmin,jj
c            write(57,956)i,j,r(i),e1(j),specj(i,j),spech(i,j)
 956        format(2i5,1pe15.3,20e12.3)
         enddo
      enddo

      ispec = 1

      if((imode.eq.2.or.imode.eq.3).and.ispec.eq.1) then

c         write(76,*)' lambda,imode,ishock ',lambda,imode,ishock,imax
         xtot=0.
         do i=1,imax
            xtot = xtot + drtr(2,i)
            do j=jmin,jj
               if(imode.eq.2) then
c                  write(76,9240)i,j,rtr(2,i),drtr(2,i),dentr(2,i),
c     &                 e1(j),em(i,j),copac(i,j),fd(i,j)
               elseif(imode.eq.3) then
c                  write(76,9240)i,j,rtr(2,i),drtr(2,i),dentr(2,i),
c     &                 e1(j),emc(i,j),copac(i,j),fc(i,j)
               endif
 9240           format(2i5,1pe15.7,3e12.4,e24.16,10e12.4)
            enddo
         enddo
      endif

      if(imode.eq.3) then

         rewind (76)

      endif

c      stop
      return
      
      end


      SUBROUTINE diffuse_sph_imp(imode,lambda,IMAX,ishock,te,totfl)
      IMPLICIT REAL*8(A-H,O-Z)
      include 'param'
      dimension f(100001)
      COMMON/FRE/NINT,JMIN,JJ
      common/opacity/TAop(NE1:NE2),Ss(MD,NE1:NE2),copac(md,ne1:ne2)
      COMMON/DIF/EM(MD,NE1:NE2),TAU(MD,NE1:NE2),TAUTOT(MD,NE1:NE2),
     &           EMC(MD,NE1:NE2)
       common/trq/rtr(2,md),drtr(2,md),dentr(2,md)
      COMMON/INT/FL(2,NE1:NE2),SI(14,27,NE1:NE2),E1(NE1:NE3),E(NE1:NE3)
      COMMON/SPECT/TEL,FD(MD,NE1:NE2),F0(NE1:NE2),IPAR
      COMMON/contspec/fc(MD,NE1:NE2),flc(2,ne1:ne2)
      common/taushok/taush(md)
      dimension r(md),denm(md),p(md)
      dimension flj(md),flh(md),flp(md),flm(md),s(md),opac(md),emm(md)
      dimension specj(md,ne1:ne2),spech(md,ne1:ne2)
      data pi/3.141592/

c fi = max angle for spherical geometry counted from the center. 

c fi = pi/2 for half sphere

c fi = pi for whole sphere

      fi = pi/2.d0

c      imax = 101

      iqq=-99999

      write(6,*)' imax ',imax,md
      write(6,*)' imax ',lambda


C     I     1         2         3         4                            *
C           +         +         +         +             +         +   
C          R(1)      R(2)      R(3)      R(4)        R(imax-1)  R(imax)         
C         taut(k)   taut(2)  taut(3)    taut(4)    taut(imax-1) taut(imax)       
C          EM(1)    EM(2)     EM(3)      EM(4)      EM(imax-1)  EM(imax)    
C        copac(1)   copac(2)  copac(3)  copac(4)   copac(imax-1) copac(imax) 

C                    J(2)      J(3)       J(4)       j(imax-1)   j(imax)  
C                                                                      
C                +         +          +       +     +         +   
C              DX(2)     DX(3)      DX(4)       DX(imax-1)  DX(imax)     
c              emm(2)    emm(3)     emm(4)      emm(imax-1) emm(imax)
c             opac(2)    opac(3)    opac(4)     opac(imax-1 opac(imax)
C               S(2)      S(3)       S(4)       s(imax-1)   s(imax)

C             tmid(1)    tmid(3)    tmid(i)                            *
C                                                                      *

c  We want mean intensities at the grid points for the calc. of ionization and 
c     temperature at each point. To calc. this we need emissivities and opacities 
c     at the midpoints. This gives the constribution from each shell to the mean int. 

c first calc. emissivities, densities  and opacities at midpoints

c note that emm(i+1), etc is evaluated at the midpoint between i and i+1


c mu = 1

c dl = r(k+1) - r(k) 

c s(k), opac(k) denm(k) alla borde vara k+1!


c Note!

c     After second and higher iteration input has i=1 for preshock
c     boundary and i=imax at postshock boundary. After first iteration
c     i=1 at post-shock boundary and i=imax at preshock.

c     this program calculates the radiation transport from outer
c     boundary (i=1) to inner (i=imax)

c     we want the postshock gas to be at the outer boundary and preshock
c     in center for that of a reverse shock situation. For 87A it is the
c     opposite case After the second and higher iteration we therefore
c     need to renumber from inside out
c   



c     first renumber so that i=1 is at the outer boundary and i = imax
c     closest to the center.  eg 87A case

c for a reverse shock case opposite situation

c radius of first shell

c!!! here input      

      rmin = 1.0d17

      r(imax)=rmin

      do m=imax,2,-1

         if(lambda.eq.2) then

            r(m-1) = r(m) + dabs(drtr(2,m))

         else

            r(m-1) = r(m) + dabs(drtr(2,imax-m+2))

         endif

      enddo 

c density at mid-points


c density at mid-points

      do i=1,imax-1

         if(lambda.eq.2) then

            denm(i+1) = (dentr(2,i) + dentr(2,i+1))/2.

         else

            denm(imax-i+1) = (dentr(2,i) + dentr(2,i+1))/2.

         endif
            
      enddo


      write(6,*)' parameters for spectrum imode',imode,ishock

      do i=1,imax

         write(6,911)i,r(i),denm(i),copac(i,3),emc(i,3),em(i,3),
     &        copac(i,-100),emc(i,-100),em(i,-100)
 911     format(i5,1pe15.7,10e12.3)

      enddo

      write(0,*)' jmin,jj',jmin,jj

      do j=jmin,jj

c         if(j.eq.3) then

         write(6,*)' j= ',j
c first calc. emissivities, densities  and opacities at midpoints

c note that emm(i+1), etc is evaluated at the midpoint between i and i+1
            
         emc(1,j) = emc(2,j)

c renumber and calculate the emissivities and opacties at the midpoints

         do i=1,imax-1

            if(lambda.eq.2) then

               if(imode.eq.2) then
                  emm(i+1) = (em(i,j) + em(i+1,j))/2.
               elseif(imode.eq.3) then
                  emm(i+1) = (emc(i,j) + emc(i+1,j))/2.
               endif
                  
               opac(i+1) = (copac(i,j) + copac(i+1,j))/2.
                  
            else
                  
               if(imode.eq.2) then
                  emm(imax-i+1) = (em(i,j) + em(i+1,j))/2.
               elseif(imode.eq.3) then
                  emm(imax-i+1) = (emc(i,j) + emc(i+1,j))/2.
               endif
                  
               opac(imax-i+1) = (copac(i,j) + copac(i+1,j))/2.
c     write(6,9288)i,imax-i+1,opac(imax-i+1),
c     &              copac(i,j),copac(i+1,j)
 9288          format('aa ',2i5,1pe12.3,10e12.3)
            endif

         enddo
            
c     now calculate the transfer with i=1 at outer bondary and imax at center
            
         do i=1,imax-1
            
            if(imode.eq.2) then
               FD(I,J)=0.0E0
            elseif(imode.eq.3) then
               fc(I,J)=0.0E0
            endif
            
            s(i+1) = emm(i+1)*denm(i+1)/opac(i+1)
            
            if(j.eq.3) then
c     write(6,9298)i,r(i),denm(i),opac(i),emm(i)
 9298          format(' qq ',i5,1pe13.5,10e12.3)
            endif
            
         enddo
         
         if(opac(2).eq.0.) then
            opac(2) = opac(3)
         endif
         
         if(j.eq.3.or.j.eq.114) then
            
            write(6,*)' parameters for spectrum',j,imode

            do i=1,imax
               write(6,9297)i,r(i),drtr(2,i),denm(i),emm(i),opac(i),
     &              s(i)
 9297          format(i5,1pe15.7,10e12.3)
            enddo
         endif

c         write(16,*)'i,j,e1(j),r(i),taush(i),flj(i),flh(i)'

         do i=1,imax

            flj(i) = 0.

            flh(i) = 0.

            nimp = imax

            do m=i,nimp
            
               p(m) = r(m)

c max distance from r(i) until ray hits the cone of the max angle 
c     from the center, fi


               flp(m) = 0.

c I+
               tauk = 0.

c upper half
               do k=i,m-1

                  dl = sqrt(r(k)**2-p(m)**2) - sqrt(r(k+1)**2-p(m)**2)
                  
                  dtauk = dl*opac(k+1)*denm(k+1)
                  
                  t2 = tauk

                  if(tauk.gt.300.) goto 44
                  
                  tauk = tauk + dtauk
                  
                  if(dtauk.lt.1.e-5) then
                     
                     dfl = s(k+1)*exp(-abs(t2))*dtauk
                     
                  else

                     dfl = s(k+1)*(-exp(-abs(dtauk))+1.d0)*
     &                    exp(-abs(t2))
                     
                  endif
                  
                  flp(m) = flp(m) + dfl

c                  if(j.eq.3.and.i.eq.iqq) then
c                     write(6,93)i,m,k,r(i),p(m),tauk,s(k+1),flp(m)
 93                  format(' flpm ',3i5,1pe15.6,10e15.6)
c                  endif
                  
               enddo


c lower half
               do k=m-1,1,-1

                  dl = sqrt(r(k)**2-p(m)**2) - 
     &                 sqrt(r(k+1)**2-p(m)**2)

                  dtauk = dl*opac(k+1)*denm(k+1)

                  t2 = tauk

                  if(tauk.gt.300.) goto 44
                  
                  tauk = tauk + dtauk

                  if(dtauk.lt.1.e-5) then

                     dfl = s(k+1)*exp(-abs(t2))*dtauk

                  else

                     dfl = s(k+1)*(-exp(-abs(dtauk))+1.d0)*
     &                       exp(-abs(t2))
                     
                  endif

                  flp(m) = flp(m) + dfl

c                  if(j.eq.3.and.i.eq.iqq) then
c                     write(6,97)i,m,k,r(i),p(m),tauk,s(k+1),flp(m)
c 97                  format(' flpml ',3i5,1pe15.6,10e15.6)
c                  endif

               enddo

 44            continue


c I-
               tauk = 0.

               flm(m) = 0.

               do k=i-1,1,-1

                  dl = sqrt(r(k)**2-p(m)**2) - 
     &                 sqrt(r(k+1)**2-p(m)**2)
                  
                  dtauk = dl*opac(k+1)*denm(k+1)
                  
                  t2 = tauk

                  if(tauk.gt.300.) goto 55
                  
                  tauk = tauk + dtauk
                  
                  if(dtauk.lt.1.e-5) then
                     
                     dfl = s(k+1)*exp(-abs(t2))*dtauk
                     
                  else

                     dfl = s(k+1)*(-exp(-abs(dtauk))+1.d0)*
     &                    exp(-abs(t2))
                        
                  endif
                  
                  flm(m) = flm(m) + dfl

                  if(j.eq.3.and.i.eq.iqq) then
                     write(6,98)i,m,k,r(i),p(m),tauk,s(k+1),flp(m)
 98                  format(' flpmm ',3i5,1pe15.6,10e15.6)
                  endif
                  
               enddo

 55            continue

            enddo


            do m=i,nimp

               dp = p(m)-p(m+1)

               if(m.eq.nimp) then
                  dp = p(nimp-1)-p(nimp)
               endif

               if(m.gt.i.and.m.lt.nimp) then

                   flj(i) = flj(i) + 0.5*(flp(m)+flm(m))*
     &                 p(m)*dp/(r(i)*sqrt(r(i)**2-p(m)**2))

                  flh(i) = flh(i) + 0.5*(flp(m)-flm(m))
     &                 *p(m)*dp/r(i)**2

               else

c end points

                  pmean = 0.5*(p(i)+p(i+1))

                  flj(i) = flj(i) + 0.25*(flp(m) + flm(m))*
     &                 p(m)*dp/(r(i)*sqrt(r(i)**2-pmean**2))

                  flh(i) = flh(i) + 0.25*(flp(m) - flm(m))*
     &                 p(m)*dp/r(i)**2

               endif

               if(j.eq.3.and.i.eq.iqq) then
                  write(6,9)i,m,r(i),p(m),flp(m),flm(m),flj(i),
     &                 flh(i)
 9                format(' flp ',2i5,1pe15.6,10e15.6)
               endif
                           
            enddo

            specj(i,j) = flj(i)
            spech(i,j) = flh(i)

            if(j.eq.3.or.j.eq.-100) then
               write(36,96)i,r(i),s(i),flj(i),flh(i)
 96            format(i5,1pe15.6,10e15.6)
            endif

            if(imode.eq.2) then
               if(lambda.eq.2) then
                  FD(I,J)=flj(i)
               else
                  FD(imax-I+1,J)=flj(i)
               endif
               
            elseif(imode.eq.3) then
               
               if(lambda.eq.2) then
                  fc(I,J)=flj(i)
               else
                  fc(imax-I+1,J)=flj(i)
               endif
            endif
            
         enddo

c         endif

      enddo

      write(27,*)' new spectrum for lambda = ',lambda,ishock
      do i=2,imax,10
         do j=jmin,jj
            write(27,956)i,j,r(i),e1(j),specj(i,j),spech(i,j)
 956        format(2i5,1pe15.3,20e12.3)
         enddo
      enddo


      if(imode.eq.2.or.imode.eq.3) then

         write(26,*)' lambda,imode,ishock ',lambda,imode,ishock
         xtot=0.
         do i=1,imax
            xtot = xtot + drtr(2,i)
            do j=jmin,jj
               if(imode.eq.2) then
                  write(26,9240)i,j,rtr(2,i),drtr(2,i),dentr(2,i),
     &                 e1(j),em(i,j),copac(i,j),fd(i,j)
               elseif(imode.eq.3) then
                  write(26,9240)i,j,rtr(2,i),drtr(2,i),dentr(2,i),
     &                 e1(j),emc(i,j),copac(i,j),fc(i,j)
               endif
 9240           format(2i5,1pe15.7,3e12.4,e24.16,10e12.4)
            enddo
         enddo
      endif


      return
      
      end


      SUBROUTINE diffuse_sph_imp_eff(imode,lambda,IMAX,ishock,
     &     te,totfl)
      IMPLICIT REAL*8(A-H,O-Z)
      include 'param'
      integer step
      real*8 mu,mup1
      real*8 intp(md,md),intm(md,md),dcol(md,md)
      COMMON/FRE/NINT,JMIN,JJ
      common/taushok/taush(md)
      common/opacity/TAop(NE1:NE2),Ss(MD,NE1:NE2),copac(md,ne1:ne2)
      COMMON/DIF/EM(MD,NE1:NE2),TAU(MD,NE1:NE2),TAUTOT(MD,NE1:NE2),
     &           EMC(MD,NE1:NE2)
       common/trq/rtr(2,md),drtr(2,md),dentr(2,md)
      COMMON/INT/FL(2,NE1:NE2),SI(14,27,NE1:NE2),E1(NE1:NE3),E(NE1:NE3)
      COMMON/SPECT/TEL,FD(MD,NE1:NE2),F0(NE1:NE2),IPAR
      COMMON/EDDFL/FH(MD,NE1:NE2),FHC(MD,NE1:NE2)
      COMMON/contspec/fc(MD,NE1:NE2),flc(2,ne1:ne2)
      common/rinner/rmin,rmax
      dimension r(md),denm(md),p(md)
      dimension flj(md),flh(md),s(md),opac(md),emm(md)
      dimension specj(md,ne1:ne2),spech(md,ne1:ne2)
      data pi/3.141592/

      write(6,*)' imax ',imax,md
      write(6,*)' imax ',lambda

      iqq = -999

C     I     1         2         3         4                            *
C           +         +         +         +             +         +   
C          R(1)      R(2)      R(3)      R(4)        R(imax-1)  R(imax)         
C         taut(k)   taut(2)  taut(3)    taut(4)    taut(imax-1) taut(imax)       
C          EM(1)    EM(2)     EM(3)      EM(4)      EM(imax-1)  EM(imax)    
C        copac(1)   copac(2)  copac(3)  copac(4)   copac(imax-1) copac(imax) 

C                    J(2)      J(3)       J(4)       j(imax-1)   j(imax)  
C                                                                      
C                +         +          +       +     +         +   
C              DX(2)     DX(3)      DX(4)       DX(imax-1)  DX(imax)     
c              emm(2)    emm(3)     emm(4)      emm(imax-1) emm(imax)
c             opac(2)    opac(3)    opac(4)     opac(imax-1 opac(imax)
C               S(2)      S(3)       S(4)       s(imax-1)   s(imax)

C             tmid(1)    tmid(3)    tmid(i)                            *
C                                                                      *

c  We want mean intensities at the grid points for the calc. of ionization and 
c     temperature at each point. To calc. this we need emissivities and opacities 
c     at the midpoints. This gives the constribution from each shell to the mean int. 

c first calc. emissivities, densities  and opacities at midpoints

c note that emm(i+1), etc is evaluated at the midpoint between i and i+1


c mu = 1

c dl = r(k+1) - r(k) 

c s(k), opac(k) denm(k) alla borde vara k+1!


c Note!

c     After second and higher iteration input has i=1 for preshock
c     boundary and i=imax at postshock boundary. After first iteration
c     i=1 at post-shock boundary and i=imax at preshock.

c     this program calculates the radiation transport from outer
c     boundary (i=1) to inner (i=imax)

c     we want the postshock gas to be at the outer boundary and preshock
c     in center for that of a reverse shock situation. For 87A it is the
c     opposite case After the second and higher iteration we therefore
c     need to renumber from inside out
c   



c     first renumber so that i=1 is at the outer boundary and i = imax
c     closest to the center.  eg 87A case

c for a reverse shock case opposite situation

c radius of first shell



      r(imax)=rmin

      do m=imax,2,-1

         if(lambda.eq.2) then

            r(m-1) = r(m) + dabs(drtr(2,m))

         else

            r(m-1) = r(m) + dabs(drtr(2,imax-m+2))

         endif

      enddo 

      nimp = imax

      do m=1,nimp
            
         p(m) = r(m)

      enddo

c add a central ray

      p(nimp+1) = 0.

c density at mid-points

      do i=1,imax-1

         if(lambda.eq.2) then

            denm(i+1) = (dentr(2,i) + dentr(2,i+1))/2.

         else

            denm(imax-i+1) = (dentr(2,i) + dentr(2,i+1))/2.

         endif
            
      enddo


      write(6,*)' parameters for spectrum imode',imode,ishock

      rmax = r(1)

      do i=1,imax

         write(6,911)i,r(i),denm(i),copac(i,3),emc(i,3)
 911     format(i5,1pe15.7,10e12.3)

      enddo

      write(0,*)' jmin,jj',jmin,jj

      do j=jmin,jj

c         if(j.eq.3) then

            write(6,*)' j= ',j
c first calc. emissivities, densities  and opacities at midpoints

c note that emm(i+1), etc is evaluated at the midpoint between i and i+1
            
         emc(1,j) = emc(2,j)

c renumber and calculate the emissivities and opacties at the midpoints

         do i=1,imax-1

            if(lambda.eq.2) then

               if(imode.eq.2) then
                  emm(i+1) = (em(i,j) + em(i+1,j))/2.
               elseif(imode.eq.3) then
                  emm(i+1) = (emc(i,j) + emc(i+1,j))/2.
               endif
            
               opac(i+1) = (copac(i,j) + copac(i+1,j))/2.
            
            else

               if(imode.eq.2) then
                  emm(imax-i+1) = (em(i,j) + em(i+1,j))/2.
               elseif(imode.eq.3) then
                  emm(imax-i+1) = (emc(i,j) + emc(i+1,j))/2.
               endif
            
               opac(imax-i+1) = (copac(i,j) + copac(i+1,j))/2.
c               write(6,9288)i,imax-i+1,opac(imax-i+1),
c     &              copac(i,j),copac(i+1,j)
 9288          format('aa ',2i5,1pe12.3,10e12.3)
            endif

         enddo

c now calculate the transfer with i=1 at outer bondary and imax at center

         do i=1,imax-1
               
            if(imode.eq.2) then
               FD(I,J)=0.0E0
            elseif(imode.eq.3) then
               fc(I,J)=0.0E0
            endif
 
            s(i+1) = emm(i+1)*denm(i+1)/opac(i+1)

            if(j.eq.3) then
c               write(6,9298)i,r(i),denm(i),opac(i),emm(i)
 9298          format(' qq ',i5,1pe13.5,10e12.3)
            endif
            
         enddo
         
         if(opac(2).eq.0.) then
            opac(2) = opac(3)
         endif

         if(j.eq.3) then

            write(6,*)' parameters for spectrum j= ',j

            do i=1,imax
               write(6,9297)i,r(i),drtr(2,i),denm(i),emm(i),opac(i),
     &              s(i)
 9297          format(i5,1pe15.7,10e12.3)
            enddo
         endif

c first all distances between shells

         if(j.eq.jmin) then

            do m=1,nimp
            
               do i=1,m

c distance between shell i and i+1 along ray m

                  dl = sqrt(r(i)**2-p(m)**2) - sqrt(r(i+1)**2-p(m)**2)

c column density of d:o. note that first shell is dcol(m,2) and last 
c   dcol(m,m)
                  
                  dcol(m,i+1) = dl*denm(i+1)
                  
               enddo

            enddo

         endif

c now the intensities along all points along each p-ray

c first the int in the inward direction, intm

         do m=1,nimp
            
            intm(m,1) = 0.

            do i=1,m-1

               dtau = dcol(m,i+1)*opac(i+1)
                  
               if(dtau.lt.1.e-5) then
                     
                  dep = dtau
                     
               else

                  dep = 1.d0 - exp(-abs(dtau))
                  
               endif

               intm(m,i+1) = s(i+1)*dep + intm(m,i)*exp(-dtau)
                  
               if(j.eq.3.and.i.eq.iqq) then
                  write(6,977)m,s(i+1),dtau,dep,intm(m,i+1)
 977              format('intm ',i5,1pe12.3,10e12.3)
               endif
            enddo

c now the int in the outward direction, intp

            intp(m,m) = intm(m,m)

            do i=m,2,-1

               dtau = dcol(m,i)*opac(i)
                  
               if(dtau.lt.1.e-5) then
                     
                  dep = dtau
                     
               else

                  dep = 1.d0 - exp(-abs(dtau))
                  
               endif

               intp(m,i-1) = s(i)*dep + intp(m,i)*exp(-dtau)

               if(j.eq.3.and.i.eq.iqq) then
                  write(6,978)m,s(i),dtau,dep,intp(m,i+1)
 978              format('intp ',i5,1pe12.3,10e12.3)
               endif
                  
            enddo
               
         enddo

c and now the mean intensity and flux

         do i=1,imax

            flj(i) = 0.

            flh(i) = 0.

            intm(nimp+1,i) = intm(nimp,i) 

            intp(nimp+1,i) = intp(nimp,i) 

            do m=i,nimp

               mu = sqrt(1.-(p(m)/r(i))**2)

               mup1 = sqrt(1.-(p(m+1)/r(i))**2)

               dmu = mup1 - mu

               flj(i) = flj(i) + 0.25*(intp(m,i) + intm(m,i) +
     &              intp(m+1,i) + intm(m+1,i))*dmu

               flh(i) = flh(i) + 0.25*(intp(m,i) - intm(m,i) +
     &              intp(m+1,i) - intm(m+1,i))*(mu+mup1)*dmu/2.

               if(j.eq.jpq.and.i.eq.iqq) then

                  write(6,9)i,m,r(i),p(m),mu,mup1,dmu,intp(m,i),
     &                 intm(m,i),intp(m+1,i),intm(m+1,i),flj(i),flh(i)
 9                format(' flp ',2i5,1pe15.6,10e15.6)

               endif
                           
            enddo

            specj(i,j) = flj(i)

            spech(i,j) = flh(i)

            write(36,96)j,i,r(i),s(i),flj(i),flh(i)
            if(j.eq.3) then
               write(6,96)j,i,r(i),s(i),flj(i),flh(i)
 96            format(2i5,1pe15.6,10e15.6)
            endif

            if(imode.eq.2) then
               if(lambda.eq.2) then
                  FD(I,J)=flj(i)
                  FH(I,J)=flh(i)
               else
                  FD(imax-I+1,J)=flj(i)
                  FH(imax-I+1,J)=flh(i)
               endif
               
            elseif(imode.eq.3) then
               
               if(lambda.eq.2) then
                  fc(I,J)=flj(i)
                  fhc(I,J)=flh(i)
               else
                  fc(imax-I+1,J)=flj(i)
                  fhc(imax-I+1,J)=flh(i)
               endif
            endif
            
         enddo

c      endif

c calculate optical depth at lyman edge from shock

         if(j.eq.3) then

            do k=1,2

               if(lambda.eq.2) then

                  if(k.eq.1) then

c  i < ishock
            
                     i1 = ishock-1

                     i2 = 1

                     step = -1

                  elseif(k.eq.2) then
            
c  i > ishock
                     i1 = ishock+1

                     i2 = imax
                  
                     step = 1
            
                  endif

                  taui = 0.

                  do i=i1,i2,step
                  
                     if(k.eq.1) then

                        dl = abs(r(i)-r(i+1))

                     else

                        dl = abs(r(i)-r(i-1))

                     endif

                     dtaui = dl*opac(i)*denm(i)

                     taui = taui + dtaui

c     taush(i) = optical depth from shock to r(i)

                     taush(i) = taui

                  enddo

               else

                  if(k.eq.1) then

c  i < ishock
            
                     i1 = imax-ishock

                     i2 = imax

                     step = -1

                  elseif(k.eq.2) then
            
c  i > ishock
                     i1 = imax-ishock

                     i2 = 1
                  
                     step = 1
            
                  endif

                  if(k.eq.1) then

c  i < ishock
            
                     i1 = ishock-1

                     i2 = 1

                     step = -1

                  elseif(k.eq.2) then
            
c  i > ishock
                     i1 = ishock+1

                     i2 = imax
                  
                     step = 1
            
                  endif

                  taui = 0.

                  do i=i1,i2,step
                  
                     if(k.eq.1) then

                        dl = abs(r(imax-i+1)-r(imax-i))

                     else

                        dl = abs(r(imax-i+1)-r(imax-i+2))

                     endif

                     dtaui = dl*opac(imax-i+1)*denm(imax-i+1)

                     taui = taui + dtaui

c taush(i) = optical depth from shock to r(i)

                     taush(i) = taui

                  enddo

               endif
                     
            enddo

            taush(ishock) = 0.

            do i=1,imax

               if(lambda.eq.2) then

                  write(6,922)i,r(i),opac(i),denm(i),taush(i)

               else

                  write(6,922)i,r(imax-i+1),opac(imax-i+1),
     &                 denm(imax-i+1),taush(i)
 922              format(' lym ',i5,1pe16.8,10e12.3)
               endif

            enddo

         endif


      enddo

      write(27,*)' new spectrum for lambda = ',lambda,ishock
      do i=2,imax,10
         do j=jmin,jj
c            write(27,956)i,j,r(i),e1(j),specj(i,j),spech(i,j)
 956        format(2i5,1pe15.3,20e12.3)
         enddo
      enddo

      if(imode.eq.2.or.imode.eq.3) then

         write(26,*)' lambda,imode,ishock ',lambda,imode,ishock
         xtot=0.
         do i=1,imax
            xtot = xtot + drtr(2,i)
            do j=jmin,jj
               if(imode.eq.2) then
                  write(26,9240)i,j,rtr(2,i),drtr(2,i),dentr(2,i),
     &                 e1(j),em(i,j),copac(i,j),fd(i,j),fh(i,j)
               elseif(imode.eq.3) then
                  write(26,9240)i,j,rtr(2,i),drtr(2,i),dentr(2,i),
     &                 e1(j),emc(i,j),copac(i,j),fc(i,j),fhc(i,j)
               endif
 9240           format(2i5,1pe15.7,3e12.4,e24.16,10e12.4)
            enddo
         enddo
      endif


      return
      
      end
