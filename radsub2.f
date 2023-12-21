      subroutine feiiiatom(te,den,z,rtot)
C
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(NMAX=30)
      SAVE
      character*80 lab
      character*11 lev
      common/wfeiv/wl(45,45)
      COMMON/FE3/WLFE3(30,30),AFE3(30,30),EMFE3(30,30)
      DIMENSION RL(45,45),X(30)
      dimension e(45),g(45),a(45,45),c(45,45)
      if(initfe3.eq.0) then
      initfe3=1
      do i=1,45
         do j=1,45
            a(i,j)=0.
            c(i,j)=0.
         enddo
      enddo
      REWIND 28
      read(28,1921)lab
      read(28,1921)lab
1921   format(a)
      n=nmax
      do 34 i=1,n
34    read(28,912)ix,e(i),g(i),lev
912   format(i3,f12.2,f7.1,a11)
      fe4=0.
      read(28,1921)lab
      read(28,1921)lab
      do 19 i=2,n
      im1=i-1
      do 19 j=1,im1
19    read(28,*)iq,jq,dlair,a(i,j),c(i,j)
      endif
      t4=te/1.e4
      CALL FORBN(NMAX,E,G,A,C,TE,Z,DEN,RL,X)
      rtot=0.
      do 22 i=2,nmax
      do 22 j=1,i-1
      wlfe3(i,j)=wl(i,j)
      emfe3(i,j)=z*rl(i,j)
 22   rtot=rtot+rl(i,j)
      return
      END

      subroutine feivatom(te,den,z,rtot)
c!!cfq hela rutinen och naesta!
C
C       THIS PROGRAM CALCULATES THE EMISSIVITY OF FORBIDDEN LINES
C       USING A 10-LEVEL SCHEME FOR CIII
C
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(NMAX=22)
      SAVE
      character*80 lab
      character*11 lev
      common/wfeiv/wl(45,45)
      COMMON/FE4/WLFE4(30,30),AFE4(30,30),EMFE4(30,30)
      DIMENSION RL(45,45),X(30)
      dimension e(45),g(45),a(45,45),c(45,45)
      if(initfe.eq.0) then
         do i=1,45
            do j=1,45
               a(i,j)=0.
               c(i,j)=0.
            enddo
         enddo
      initfe=1
      REWIND 27
      read(27,1921)lab
      read(27,1921)lab
1921   format(a)
      n=nmax
      do 34 i=1,n
34    read(27,912)ix,e(i),g(i),lev
912   format(i3,f12.2,f7.1,a11)
      fe4=0.
      read(27,1921)lab
      read(27,1921)lab
      do 19 i=2,n
      im1=i-1
      do 19 j=1,im1
19    read(27,*)iq,jq,dlair,a(i,j),c(i,j)
      endif
      t4=te/1.e4
c!!
      CALL FORBN(NMAX,E,G,A,C,TE,Z,DEN,RL,X)
      rtot=0.

      do i=2,nmax
         do j=1,i-1
            wlfe4(i,j)=wl(i,j)
            emfe4(i,j)=z*rl(i,j)
            rtot=rtot+rl(i,j)
         enddo
      enddo
      return
      END

      SUBROUTINE FORBN(N,E,G,A,CQ,TE,Z,DEN,RL,X)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER N
c      PARAMETER (MD=350,MDP1=MD+1)
      include 'PARAM'
c      PARAMETER (NL=30,NLP1=NL+1)
      COMMON/PHY/DENS(MD)
      COMMON/ELEC/DEL(MD)
      COMMON/TPAR/RIN,DRQ,R1Q,TDAY
      COMMON/IND/IK
      common/wfeiv/wl(45,45)
      DIMENSION  C(45,45),CQ(45,45),G(45),E(45),EEV(45),A(45,45)
      DIMENSION AA(31,31),X(30),RL(45,45),BE(45,45)
c      DIMENSION AA(nl+1,nl+1),X(nl),RL(45,45),BE(45,45)

      DATA PI/3.1415926d0/
      DEN=DENS(IK)
      DENEL=DEN*DEL(IK)
      NM1=N-1
      NP1=N+1
      DO  I=1,NP1
         DO  J=1,NP1
            AA(I,J)=0.
            c(i,j)=0.
         enddo
      enddo
      DO 220 I=1,N
 220  EEV(I)=E(I)/8065.46
      DO 5395 I=1,N
      DO 5394 J=1,N
      WL(I,J)=0.0
      IF(I.EQ.J) GOTO 5394
      IF(E(I).EQ.E(J)) GOTO 5394
      WL(I,J)=DABS(12398.54/(EEV(I)-EEV(J)))
5394  CONTINUE
5395  CONTINUE
      TEV=TE/1.1609D4
      TIME=8.64E4*TDAY
      DO 5 I=1,NM1
      IP1=I+1
      DO 5 J=IP1,N
      C(J,I)=8.63D-6*DEN*CQ(J,I)/(DSQRT(TE)*G(J))
5     C(I,J)=C(J,I)*DEXP(-(EEV(J)-EEV(I))/TEV)*G(J)/G(I)
      DO L=1,4 
      DO I=1,N
         DO J=1,N
            IF(i.ne.j) then
               DEI=X(I)*Z*DEN
C
C     CALCULATE THE ESCAPE PROBABILITY.
C
               T0=1.e-24*WL(J,I)**3*A(J,I)*G(J)*DEI*TIME/(8.*PI*G(I))
               IF(L.EQ.1) T0=.0
               IF(T0.LT.0.01) BE(J,I)=1.-0.5*T0+T0**2/6.
               IF(T0.GE.0.01) BE(J,I)=(1.-EXP(-T0))/T0
               IF(I.EQ.N) then
                  AA(I,J)=1.D0
               else
                  AA(I,J)=BE(J,I)*A(J,I)+C(J,I)
               endif
            elseif(i.eq.j) then
               be(i,j)=0.   
            endif
         enddo
      enddo
 100  CONTINUE
      do i=1,n
         S=0.
         DO  j=1,N
            S=S+C(I,j)+BE(I,j)*A(I,j)
         enddo
         AA(I,I)=-S
      enddo
      NRC=N+1
      AA(N,NRC)=1.
      EPS=1.D-30
      NA=N
      DI=SIMULdi(NA,AA,X,EPS,1,NRC)
      ENDDO
      RTOT=0.
      DO 10 I=1,NM1
      IP1=I+1
      DO 10 J=IP1,N
      IF(A(J,I).GT.0.) THEN
      RL(J,I)=1.602E-12*(EEV(J)-EEV(I))*X(J)*BE(J,I)*A(J,I)/DEN
      RTOT=RTOT+RL(J,I)
      ENDIF
 10   CONTINUE
25    CONTINUE
      RETURN
      END

      DOUBLE PRECISION FUNCTION SIMULdi(N,A,X,EPS,II,NRC)
      IMPLICIT REAL*8(A-H,O-Z)
c      PARAMETER (NL=340,NLP1=NL+1)
C     SOLVE A*X = B
C     B CONTAINS AS INPUT THE R.H. SIDE OF THE SYSTEM, AND AS OUTPUT THE
C     SOLUTION
C     EPS, II, NRC NOT USED
      DIMENSION A(31,31),X(30),B(31)
      DIMENSION LASTN(31),ASAVE(31)                                     
C                                                                       
C -- Find the column index of the last non-zero element in each row --  
C -- Speeds up the solution of loose matrices considerably --           
C                                                                       
      DO I=1,N
      B(I)=A(I,N+1)
      ENDDO
      DO 100 I=1,N                                                      
      DO 101 J=N,I,-1                                                   
      JSAVE=J                                                           
      IF(A(I,J).NE.0.0D0) GOTO 102                                      
  101 CONTINUE                                                          
  102 LASTN(I)=JSAVE                                                    
  100 CONTINUE                                                          
C                                                                       
C -- ************************** --                                      
C -- Forward elimination scheme --                                      
C -- ************************** --                                      
C                                                                       
      DO 200 I=1,N-1                                                    
C                                                                       
C -- Partial pivoting routine --                                        
C                                                                       
      AI=0.0D0                                                          
      DO 201 L=I,N                                                      
      ASAVE(L)=0.0D0                                                    
      IF(DABS(A(L,I)).GT.AI) LSAVE=L                                    
      IF(DABS(A(L,I)).GT.AI) AI=DABS(A(L,I))                            
  201 CONTINUE                                                          
      IF(LSAVE.EQ.I) GOTO 202                                           
C                                                                       
C -- Interchange rows A(I,?) and A(LSAVE,?) if LSAVE.NE.I --            
C                                                                       
      LAST=MAX0(LASTN(I),LASTN(LSAVE))                                  
      DO 203 L=I,LAST                                                   
      ASAVE(L)=A(I,L)                                                   
      A(I,L)=A(LSAVE,L)                                                 
  203 A(LSAVE,L)=ASAVE(L)                                               
      BSAVE=B(I)                                                        
      B(I)=B(LSAVE)                                                     
      B(LSAVE)=BSAVE                                                    
      LASTNI=LASTN(I)                                                   
      LASTN(I)=LASTN(LSAVE)                                             
      LASTN(LSAVE)=LASTNI                                               
  202 CONTINUE                                                          
C                                                                       
C -- Elimination routine --                                             
C                                                                       
      BI=B(I)                                                           
      AII=A(I,I)                                                        
      DO 301 J=I+1,N                                                    
      IF(A(J,I).EQ.0.0D0.OR.(I+1).GT.LASTN(I)) GOTO 301                 
      AJIAII=A(J,I)/A(I,I)                                              
      DO 302 K=I+1,LASTN(I)                                             
  302 A(J,K)=A(J,K)-AJIAII*A(I,K)                                       
      B(J)=B(J)-AJIAII*BI                                               
      LASTN(J)=MAX0(LASTN(J),LASTN(I))                                  
  301 CONTINUE                                                          
C                                                                       
C -- ************************** --                                      
C -- End of forward elimination --                                      
C -- ************************** --                                      
C                                                                       
  200 CONTINUE                                                          
C                                                                       
C -- *************** --                                                 
C -- Back-substitute --                                                 
C -- *************** --                                                 
C                                                                       
      DO 400 K=N,1,-1                                                   
      BK=B(K)                                                           
      IF((K+1).GT.N) GOTO 401                                           
      DO 402 L=K+1,N                                                    
      BK=BK-A(K,L)*B(L)                                                 
  402 CONTINUE                                                          
  401 B(K)=BK/A(K,K)                                                    
  400 CONTINUE                                                          
C                                                                       
      DO I=1,N
      X(I)=B(I)
      ENDDO
      SIMUL=0.
      RETURN                                                            
      END



      subroutine hrec(t,aln,ala,alb)
c
c     Hydrogen recombination rates taken from Osterbrock at 1E4 K
c     powerlaw fit between 5000 and 10000 K
c
      implicit real*8(a-h,o-z)
      dimension all(4,4),aln(4)
      t4=t/1.e4
      do k=1,4
         do l=1,4
            all(k,l)=0.d0     
         enddo
      enddo
      all(1,1)=1.58e-13/t4**.529
      all(2,1)=2.34e-14/t4**.365
      all(2,2)=5.35e-14/t4**.639
      all(3,1)=7.81e-15/t4**.532
      all(3,2)=2.04e-14/t4**.643
      all(3,3)=1.73e-14/t4**.809
      all(4,1)=3.59e-15/t4**.543
      all(4,2)=9.66e-15/t4**.644
      all(4,3)=1.08e-14/t4**.815
      all(4,4)=5.54e-15/t4**.976
      do k=1,4
         aln(k)=0.
         do l=1,4
            aln(k)=all(k,l)+aln(k)
         enddo
      enddo
      ala=4.18e-13/t4**.706
      alb=2.59e-13/t4**.810
      return
      end


      subroutine hrecem(t,jem)
c
c     Hydrogen recombination emission rates taken from Osterbrock at 1E4 K
c     Table 4.4 T=1.e4 n=1.e6
c     case B
c     powerlaw fit between 5000 and 10000 K.
c
      implicit real*8(a-h,o-z)
      real*8 jem,jhb
      dimension jem(7,4)
      t4=t/1.e4
      jem(3,2)=2.81
      jem(4,2)=1.
      jem(5,2)=0.471
      jem(6,2)=0.262
      jem(7,2)=0.163
      jem(4,3)=0.317
      jem(5,3)=0.335
      jem(6,3)=0.339
      jem(7,3)=0.333
      jem(5,4)=0.154
      jem(6,4)=0.163
      jem(7,4)=0.163
      jhb=1.25e-25/t4**0.828
      do k=2,4
      do l=k+1,5
      jem(l,k)=jhb*jem(l,k)
      enddo
      enddo
      return
      end

