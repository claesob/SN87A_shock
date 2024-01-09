C 900228
      DOUBLE PRECISION FUNCTION EXP1(X)
      IMPLICIT REAL*8(A-H,O-Z)
      IF(X.LT.1.) THEN
            E1=-LOG(X)-.577216+X-.24991*X*X+.0552*X*X*X-.00976
     &            *X*X*X*X
      ELSE
            E1=EXP(-X)*(X*X+2.334733*X+.250621)/((X*X+3.330657*X
     &            +1.681534)*X)
      ENDIF
      EXP1=E1
      RETURN
      END

      DOUBLE PRECISION FUNCTION E2(X)
      IMPLICIT REAL*8(A-H,O-Z)
      IF(X.GT.0.) GOTO 500
      E2=1.
      GOTO 400
500   IF(X.GT.300.) GOTO 100
      IF(X.GT.1.) GOTO 200
      E1=-LOG(X)-.577216+X-.24991*X*X+.0552*X*X*X-.00976
     &*X*X*X*X
      GOTO 300
200   E1=EXP(-X)*(X*X+2.334733*X+.250621)/((X*X+3.330657*X+1.681534
     &)*X)
300   E2=EXP(-X)-X*E1
      GOTO 400
100   E2=0.
400   CONTINUE
      RETURN
      END

      SUBROUTINE SOLVE(M1,M2,N1,XR,XC)
      IMPLICIT REAL*8(A-H,O-Z)
C     **************************************************************
C     *****
C     THIS ROUTINE SOLVES THE COUPLED IONIZATION EQUIL. EQUATIONS.
C      M1 = LOWEST IONIZATION STAGE
C      M2 = HIGHEST IONIZATION STAGE = NUMBER OF EQUATIONS-M1+1
C      N1 IS NOT USED
C     XR(I)=N(+I)/N(+I+1)
C     XO(I+1)=N(+I)/N(JK)
C     XC(I+1)=N(+I)/N
C     *****
C     **************************************************************
      DIMENSION XO(30),Q(30),XR(30),XC(30)
      DO L=1,30
         XC(L)=1.E-20
      enddo
      DO 725 L=M1,M2
      IF(XR(L).GE.0.0) GOTO 726
  726 IF(XR(L).GT.0.0E0) GOTO 725
      XR(L)=1.E-25
  725 Q(L)=LOG10(XR(L))
      M3=M2-M1
      M4=M2-1
      M8=M2+1
      JK=0
      IF(Q(M2).GT.0.0) GOTO 730
      NN=1
      XSUM=0.0E0
      DO 740 L=1,M3
      MA=M2-L
      IF(Q(MA).LT.0.0) GOTO 740
      IF(NN.GT.1) GOTO 740
      JK=M2-L
      NN=2
  740 CONTINUE
      IF(JK.GT.0) GOTO 720
      XO(M1)=1.0E0
      XSUM=0.0E0
      M9=M1+1
      DO 100 J=M9,M8
  100 XO(J)=XO(J-1)*XR(J-1)
      XSUM=0.E0
      DO 150 J=M1,M8
  150 XSUM=XSUM+XO(J)
      XC(M1)=1./XSUM
      JK=M1
      GOTO 770
  720 XO(JK+1)=1.0E0
      M5=JK-M1+1
      M6=JK+2
      DO 600 I=1,M5
      MA=JK-I+1
       MB=JK-I+2
  600 XO(MA)=XO(MB)/XR(MA)
      DO 700 I=M6,M8
  700 XO(I)=XO(I-1)*XR(I-1)
      M7=M1
      JK=JK+1
      GOTO 750
  730 XO(M2+1)=1.0E0
      M7=M3+1
      DO 735 I=1,M7
      MA=M2+1-I
      MB=M2-I+2
  400 IF(ABS(XR(MA)).GT.1.D-60) XO(MA)=XO(MB)/XR(MA)
      IF(ABS(XR(MA)).LT.1.D-60) XO(MA)=0.
  735 CONTINUE
      XSUM=0.E0
      DO 760 K=M1,M8
  760 XSUM=XSUM+XO(K)
      JK=M8
      XC(M8)=1./XSUM
      M7=M1
      GOTO 770
 750  XSUM=0.E0
      DO 775 J=M1,M8
      XSUM=XSUM+XO(J)
  775 CONTINUE
      XC(JK)=1./(XSUM)
  770 DO 780 J=M1,M8
      IF(J.EQ.JK) GOTO 780
      XC(J)=XC(JK)*XO(J)
  780 CONTINUE
      RETURN
      END

      DOUBLE PRECISION FUNCTION SIMUL_ion(N,A,X,EPS,II,NRC)
      IMPLICIT REAL*8(A-H,O-Z)
c      PARAMETER (NL=130,NLP1=NL+1)
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

      DO I=1,N                                                      
         DO J=N,I,-1                                                   
            JSAVE=J                                                           
            IF(A(I,J).NE.0.0D0) GOTO 102                                      
         enddo
 102     LASTN(I)=JSAVE
      enddo

C     
C -- ************************** --                                      
C -- Forward elimination scheme --                                      
C -- ************************** --                                      
C                 
      DO I=1,N-1    
C     
C     -- Partial pivoting routine --                                        
C     
         AI=0.0D0                                                          
         DO L=I,N                                                      
            ASAVE(L)=0.0D0                                                    
            IF(DABS(A(L,I)).GT.AI) then
               LSAVE=L                                    
            endif
            IF(DABS(A(L,I)).GT.AI) AI=DABS(A(L,I))                            
         enddo
         IF(LSAVE.ne.I) then
C     
C     -- Interchange rows A(I,?) and A(LSAVE,?) if LSAVE.NE.I --            
C     
            LAST=MAX0(LASTN(I),LASTN(LSAVE))                                  
            DO L=I,LAST                                                   
               ASAVE(L)=A(I,L)                                                   
               A(I,L)=A(LSAVE,L)                                                 
               A(LSAVE,L)=ASAVE(L)
            enddo
            BSAVE=B(I)                                                        
            B(I)=B(LSAVE)                                                     
            B(LSAVE)=BSAVE                                                    
            LASTNI=LASTN(I)                                                   
            LASTN(I)=LASTN(LSAVE)                                             
            LASTN(LSAVE)=LASTNI
         endif                                             
C     
C     -- Elimination routine --                                             
C     
            BI=B(I)                                                           
            AII=A(I,I)                                                        
            DO  J=I+1,N                                                    
               IF(A(J,I).EQ.0.0D0.OR.(I+1).GT.LASTN(I)) GOTO 301                
               AJIAII=A(J,I)/A(I,I)                                             
               DO K=I+1,LASTN(I)                                             
                  A(J,K)=A(J,K)-AJIAII*A(I,K)
               enddo
               B(J)=B(J)-AJIAII*BI                                              
               LASTN(J)=MAX0(LASTN(J),LASTN(I))                                
 301           CONTINUE
            enddo
C     
C     -- ************************** --                                      
C     -- End of forward elimination --                                      
C     -- ************************** --                                      
C     
         enddo
C     
C -- *************** --                                                 
C -- Back-substitute --                                                 
C -- *************** --                                                 
C                                                                       
      DO K=N,1,-1                                                   
         BK=B(K)                                                           
         IF((K+1).GT.N) GOTO 401                                           
         DO L=K+1,N                                                    
            BK=BK-A(K,L)*B(L)                                                 
         enddo
 401     B(K)=BK/A(K,K)                                                    
      enddo
C                                                                       
      DO I=1,N
      X(I)=B(I)
      ENDDO
      SIMUL=0.
      RETURN                                                            
      END

      DOUBLE PRECISION FUNCTION SIMUL(N,A,X,EPS,II,NRC)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NL=340,NLP1=NL+1)
C     SOLVE A*X = B
C     B CONTAINS AS INPUT THE R.H. SIDE OF THE SYSTEM, AND AS OUTPUT THE
C     SOLUTION
C     EPS, II, NRC NOT USED
      DIMENSION A(NLP1,NLP1),X(Nl),B(NLP1)
      DIMENSION LASTN(NLP1),ASAVE(NLP1)                                     
C     
C     -- Find the column index of the last non-zero element in each row --  
C     -- Speeds up the solution of loose matrices considerably --           
C     
      DO I=1,N
         B(I)=A(I,N+1)
      ENDDO
      DO  I=1,N                                                      
         DO J=N,I,-1                                                   
            JSAVE=J                                                           
            IF(A(I,J).NE.0.0D0) GOTO 102                                      
         enddo
 102     LASTN(I)=JSAVE
      enddo

C     
C     -- ************************** --                                      
C     -- Forward elimination scheme --                                      
C     -- ************************** --                                      
C     
      DO I=1,N-1                                                    
C     
C     -- Partial pivoting routine --                                        
C     
         AI=0.0D0                                                          
         DO L=I,N                                                      
            ASAVE(L)=0.0D0                                                    
            IF(DABS(A(L,I)).GT.AI) LSAVE=L                                    
            IF(DABS(A(L,I)).GT.AI) AI=DABS(A(L,I))                            
         enddo
         IF(LSAVE.EQ.I) GOTO 202                                           
C     
C     -- Interchange rows A(I,?) and A(LSAVE,?) if LSAVE.NE.I --            
C     
         LAST=MAX0(LASTN(I),LASTN(LSAVE))                                  
         DO  L=I,LAST                                                   
            ASAVE(L)=A(I,L)                                                   
            A(I,L)=A(LSAVE,L)                                                 
            A(LSAVE,L)=ASAVE(L)
         enddo
         BSAVE=B(I)                                                        
         B(I)=B(LSAVE)                                                     
         B(LSAVE)=BSAVE                                                    
         LASTNI=LASTN(I)                                                   
         LASTN(I)=LASTN(LSAVE)                                             
         LASTN(LSAVE)=LASTNI                                               
 202     CONTINUE                                                          
C     
C     -- Elimination routine --                                             
C     
         BI=B(I)                                                           
         AII=A(I,I)                                                        
         DO 301 J=I+1,N                                                    
            IF(A(J,I).EQ.0.0D0.OR.(I+1).GT.LASTN(I)) GOTO 301                 
            AJIAII=A(J,I)/A(I,I)                                              
            DO  K=I+1,LASTN(I)                                             
               A(J,K)=A(J,K)-AJIAII*A(I,K)
            enddo
            B(J)=B(J)-AJIAII*BI                                               
            LASTN(J)=MAX0(LASTN(J),LASTN(I))                                  
 301     CONTINUE                                                          
C     
C     -- ************************** --                                      
C     -- End of forward elimination --                                      
C     -- ************************** --                                      
C     
      enddo
C     
C     -- *************** --                                                 
C     -- Back-substitute --                                                 
C     -- *************** --                                                 
C     
      DO K=N,1,-1                                                   
         BK=B(K)                                                           
         IF((K+1).GT.N) GOTO 401                                           
         DO L=K+1,N                                                    
            BK=BK-A(K,L)*B(L)
         enddo
 401     B(K)=BK/A(K,K)
      enddo
C     
      DO I=1,N
         X(I)=B(I)
      ENDDO
      SIMUL=0.
      RETURN                                                            
      END


      DOUBLE PRECISION FUNCTION SIMULT(N,A,X,EPS,II,NRC)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NL=340,NLP1=NL+1)
C     SOLVE A*X = B
C     B CONTAINS AS INPUT THE R.H. SIDE OF THE SYSTEM, AND AS OUTPUT THE
C     SOLUTION
C     EPS, II, NRC NOT USED
      DIMENSION A(NLP1,NLP1),X(Nlp1),B(NLP1)
      DIMENSION LASTN(NLP1),ASAVE(NLP1)                                     
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
      SIMULt=0.
      RETURN                                                            
      END

      DOUBLE PRECISION FUNCTION SIMULB(N,A,X,EPS,INDIC,NRC)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION IROW(200),A(21,21),JCOL(200),JORD(200),Y(200),X(N)
      MAX=N
      IF(INDIC.GE.0) MAX=N+1
      IF(N.LE.50) GOTO 5
      SIMULB=0.
      RETURN
5     DETER=1.
      DO 18 K=1,N
      KM1=K-1
      PIVOT=0.
      DO 11 I=1,N
      DO 11 J=1,N
      IF(K.EQ.1) GOTO 9
      DO 8 ISCAN=1,KM1
      DO 8 JSCAN=1,KM1
      IF(I.EQ.IROW(ISCAN)) GOTO 11
      IF(J.EQ.JCOL(JSCAN)) GOTO 11
8     CONTINUE
9     IF(ABS(A(I,J)).LE.ABS(PIVOT)) GOTO 11
      PIVOT=A(I,J)
      IROW(K)=I
      JCOL(K)=J
11    CONTINUE
777   FORMAT(1X,'STOPP I POP',3I4,2E12.5)
      IF(ABS(PIVOT).GT.EPS) GOTO 13
      SIMULB=0.
      WRITE(6,777)K,IROW(K),JCOL(K),PIVOT,A(IROW(K),JCOL(K))
      STOP
13    IROWK=IROW(K)
      JCOLK=JCOL(K)
      DETER=DETER*PIVOT
      DO 14 J=1,MAX
14    A(IROWK,J)=A(IROWK,J)/PIVOT
      A(IROWK,JCOLK)=1./PIVOT
      DO 19 I=1,N
      AIJCK=A(I,JCOLK)
      IF(I.EQ.IROWK) GOTO 19
      A(I,JCOLK)=-AIJCK/PIVOT
      DO 17 J=1,MAX
17    IF(J.NE.JCOLK) A(I,J)=A(I,J)-AIJCK*A(IROWK,J)
19    CONTINUE
18    CONTINUE
      DO 20 I=1,N
      IROWI=IROW(I)
      JCOLI=JCOL(I)
      JORD(IROWI)=JCOLI
20    IF(INDIC.GE.0) X(JCOLI)=A(IROWI,MAX)
      INTCH=0
      NM1=N-1
      DO 22 I=1,NM1
      IP1=I+1
      DO 22 J=IP1,N
      IF(JORD(J).GE.JORD(I)) GOTO 22
      JTEMP=JORD(J)
      JORD(J)=JORD(I)
      JORD(I)=JTEMP
      INTCH=INTCH+1
22    CONTINUE
      IF(INTCH/2*2.NE.INTCH) DETER=-DETER
      IF(INDIC.LE.0) GOTO 26
      SIMULB=DETER
      RETURN
26    CONTINUE
      DO 28J=1,N
      DO 27 I=1,N
      IROWI=IROW(I)
      JCOLI=JCOL(I)
27    Y(JCOLI)=A(IROWI,J)
      DO 28 I=1,N
28    A(I,J)=Y(I)
      DO 30 I=1,N
         DO J=1,N
            IROWJ=IROW(J)
            JCOLJ=JCOL(J)
            Y(IROWJ)=A(I,JCOLJ)
         enddo
         DO 30 J=1,N
 30         A(I,J)=Y(J)
      SIMULB=DETER
      RETURN
  200  FORMAT(1X,'TOO BIG!')
      END

      SUBROUTINE SIMPA(N,F,H,S)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION F(N)
C     OBS! N MUST BE AN ODD NUMBER
      NC=N-1
      ND=N-2
      NA=N+1
      NB=N-3
      SUM4=0.
      SUM2=0.
      DO I=2,NC,2
         SUM4=SUM4+F(I)
      enddo
      DO I=3,ND,2
         SUM2=SUM2+F(I)
      enddo
      S=H*(4.*SUM4+2.*SUM2+F(1)+F(N))/3.
      RETURN
      END



      SUBROUTINE SNPAR(TDAY,RPH,TEFF,RLUM)
      IMPLICIT REAL*8(A-H,O-Z)
C     PARAMETERS FOR SN 1980K D=8MPC EB-V=0.35 T SINCE MAX
      IF(TDAY.LT.90.) BM0=0.0582*TDAY-19.49
      IF(TDAY.GT.90.) BM0=0.01005*TDAY-14.87
      IF(TDAY.LT.60.) BV0=-.05+1.3333E-2*TDAY
      IF(TDAY.GT.60.) BV0=0.75
      IF(TDAY.LT.60.) TEFF=7300./(.55+1.3333E-2*TDAY)
      IF(TDAY.GT.60.) TEFF=5407.
      TEFF=6400.
      BC=-42.54+10.*LOG10(TEFF)+2.9E4/TEFF
      BOL=BM0-BC
      RPH=6.96E10*10.**(8.472-2.*LOG10(TEFF)-0.2*BOL)
      RPH=2.E14
      RLUM=4.*3.1415*RPH**2*5.67E-5*TEFF**4
      RETURN
      END

      SUBROUTINE ABUND(INI,RM)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 MTOT,MNI
      common/abl/abn(15)
      COMMON/HYDROS/HSCALE,DEN1,XEL,AMEAN
      COMMON/MPAR/MTOT,MNI,VEXP
      COMMON/THER/A1,B1,TIN,E10,E20
      DATA SOLMA/1.989E33/
      RSUN=MTOT/SOLMA
      GOTO 110
      IF(RM.GT.0.3*RSUN/1.4) GOTO 120
      ABn(1)=0.1E-10
      ABn(2)=0.00000001
      ABn(5)=0.990
      ABN(3)=0.005
      ABN(4)=0.00000001
      ABN(10)=0.00001
      ABN(8)=0.001
      ABN(14)=0.00001
      ABN(9)=1.E-11
      ABN(7)=0.0001
      ABN(13)=0.0001
      ABN(11)=0.0010
      GOTO 110
 120  CONTINUE
      ABN(1)=0.1E-10
      ABN(2)=0.00000001
      ABN(5)=0.900
      ABN(3)=0.10
      ABN(4)=0.00000001
      ABN(10)=0.00001
      ABN(8)=0.0010
      ABN(14)=4.E-9
      ABN(9)=1.E-11
      ABN(7)=0.0001
      ABN(13)=0.0001
      ABN(11)=0.0010
 110  CONTINUE
      DO IP=1,11
         ABn(IP)=1.E-6
      enddo
      ABN(14)=1.
C     NOMOTO ET AL MIXED C6 MODEL
      GOTO 130
      IF(RM.LT.0.65) GOTO 130
       IF(INI.LE.0) A1=3300.
       IF(INI.LE.0) B1=7000.
       INI=1
      ABn(1)=0.1E-10
      ABn(2)=0.00000001
      ABn(5)=0.410
      ABn(3)=0.350
      ABn(4)=0.00000001
      ABn(10)=0.12
      ABn(8)=0.018
      ABn(14)=0.00001
      ABn(9)=1.E-11
      ABn(7)=0.0001
      ABn(13)=0.0022
      ABn(11)=0.0550
 130  CONTINUE
C     GOTO 150
C     NOMOTO ET AL UNMIXED C6 MODEL
      IF(RM.LT.0.65) GOTO 140
       IF(INI.LE.0) A1=5000.
       IF(INI.LE.0) B1=10000.
       INI=2
      ABn(1)=0.1E-10
      ABn(2)=0.00000001
      ABn(5)=0.372
      ABn(3)=0.000
      ABn(4)=0.00000001
      ABn(10)=0.311
      ABn(8)=0.047
      ABn(14)=0.0001
      ABn(9)=1.E-11
      ABn(7)=0.0001
      ABn(13)=0.059
      ABn(11)=0.146
 140  CONTINUE
      IF(RM.LT.1.034) GOTO 150
       IF(INI.LE.2) A1=5000.
       IF(INI.LE.2) B1=10000.
       INI=3
      ABn(1)=0.1E-10
      ABn(2)=0.00000001
      ABn(5)=0.429
      ABn(3)=0.569
      ABn(4)=0.00000001
      ABn(10)=4.56E-4
      ABn(8)=3.69E-4
      ABn(14)=5.48E-4
      ABn(9)=1.E-11
      ABn(7)=0.0001
      ABn(13)=2.72E-5
      ABn(11)=2.20E-4
 150  CONTINUE
C     NORMALIZE ABUNDANCES
      SAB=0.
      DO  I=1,14
         SAB=SAB+ABn(I)
      enddo
      DO I=1,14
         ABn(I)=ABn(I)/SAB
      enddo
C     MEAN ATOMIC WHEIGHT
      AMEAN=ABn(1)+4.*ABn(2)+12.*ABN(3)+14.*ABN(4)+16.*ABN(5)+
     &28.*ABN(10)+56.*ABN(14)+26.*ABN(9)+24.*ABN(8)+23.*ABN(7)+
     &40.*ABN(13)+32.*ABN(11)
      RETURN
      END


      SUBROUTINE STRUC(IREV,ISTRU,RM,DENI,RI,RINNER,RM1CGS)
C
C     READS THE STRUCTURE OF THE SUPERNOVA FROM FILE 12
C     DESIGNED FOR THE ENSMAN-WOOSLEY MODELS, BUT MAY BE USED IN GENERAL
C     RM= MASS FROM CENTER
C     IF ISTRU=1 USE DENSITY FROM MODEL, IF NOT USE A UNIFORM
C     DENSITY SHELL
C

c abi = usual order H, He, C, N, O....

c ab = 'cf order' H, He, O, C, N...

      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 MR,MSUN,MTOT,MTOTS,MNI,MHE,MOX
c      PARAMETER (MD=350,MDP1=MD+1)
      include 'param'
      SAVE
      COMMON/IND/II
      COMMON/RADIE/RA(0:MD)
      COMMON/PHY/DEN(MD)
      COMMON/ELEC/DEL(MD)
      COMMON/TEM/TE(MD)
      common/abl/abn(15)
      COMMON/MASSES/RMI(MD)
      COMMON/QSOM/qso,IAGN,ISTAT,ISOBOL
      COMMON/AGNPARA/DECON,GAMMION,IPRESS,IDENS,ICDEN
      COMMON/MPAR/MTOTS,MNI,VEXP
      COMMON/HYDROS/HSCALE,DEN1,XEL,AMEAN
      COMMON/TPAR/RIN,DRQ,R1Q,TDAY
      common/initstr/init
      common/revs/terev,tecs
      COMMON/CHEMCOMP/ABIN(12),ISOLAR,ISNCOMP,INPUTCOMP
      common/abun_ord/abseq(20)
      DIMENSION ABU(15,0:1200),MR(0:1200),R(0:1200),DENS(0:1200),
     &            ABI(15)
      DIMENSION ABIS(15),ABSUN(12),ABSUNH(12)
      DATA MSUN/1.989d33/
      DATA PI/3.1415926d0/,ELCH/1.60219d-12/,AMU/1.660531d-24/
c      DATA ABSUN/7.500d-01,2.500d-01,3.981d-03,1.288d-03,1.047d-02,
c     &      1.660d-03,6.457d-04,9.333d-04,5.128614d-04,2.512d-04,
c     &      7.943d-05,2.239d-03/
c     Ferland and Netzer
c      DATA ABSUN/0.909,.0909,3.36d-4,1.09d-4,6.18d-4,1.0d-4,0.30d-4,
c     &0.29d-4,0.145d-4,1.d-30,1.d-30,0.236d-4/
c     NRL test
      DATA ABSUN/0.909d0,.0909d0,2.73d-4,0.91d-4,5.45d-4,1.4d-4,0.27d-4,
     &0.27d-4,0.135d-4,1.d-30,1.d-30,0.91d-7/
C     SOLAR COMPOSITION FROM GREVESSE AND ANDERS 1989, COSMIC ABUNDANCES OF
C     MATTER, AIP CONFERENCE PROC 183.
C
      DATA ABSUNH/1.000d0,.098d0,3.63d-4,1.12d-4,8.51d-4,1.23d-4,
     &     3.80d-5,3.55d-5,1.62d-5,3.63d-6,2.29d-6,4.68d-5/
      dimension ab87a(14)
      data ab87a/1.,0.25,4.2e-5,1.9e-4,1.9e-4,6.2e-5,1.0e-6,1.5e-5,
     &     1.2e-6,9.7e-6,5.6e-6,3.5e-6,1.1e-6,3.4e-5/
      if(iagn.eq.1) t0=1.
c!!
      ipwn=1
      if(ipwn==1) then
         write(6,*)'Core composition '
         open(21,file='average_abundances_by_num_19M.dat',status='old')
         
c     izone: O-Ne-Mg = 1; O-Si-Mg = 2; O(0.61)-Si-S = 3;  O(0.56)-Si-S = 4;
c                 O(0.46)-Si-S = 5; O(0.16)-Si-S = 6; Si-S-Ar = 7
         izone=3
         do iz=1,7
            read(21,*)zone_mass,(abn(i),i=1,14)
            if(iz==izone) then
               goto 333
            endif
         enddo
 333     write(6,9275)zone_mass,(abn(i),i=1,14)
 9275    format(' Zone mass and input abundances by number ',1pe11.3,15e11.3)
      endif

C     MEAN ATOMIC WEIGHT PER ION
      AMEAN=ABN(1)+4.*ABN(2)+12.*ABN(3)+14.*ABN(4)+16.*ABN(5)+
     &     28.*ABN(10)+56.*ABN(14)+26.*ABN(9)+24.*ABN(8)+23.*ABN(7)+
     &     40.*ABN(13)+32.*ABN(11)+20.*ABN(6)+36.*ABN(12)
      write(6,*)' TDAY, T0  ',TDAY,T0

      RETURN
      END
      
      DOUBLE PRECISION FUNCTION FMEAN(M,E1)
c
c calculate mean intensity J
c
c  fsh = mean inte. from file with the input en. bins.
c  fshm = mean intensity with the energy bins used in the code
c
      IMPLICIT REAL*8(A-H,O-Z)
      include 'param'
c      PARAMETER (NE1=-100,NE2=250,NE3=NE2+1)
      SAVE
C
C     IONIZING SPECTRUM
C
      COMMON/TPAR/RIN,DR,R1,TDAYS
      COMMON/QSPEC/GAMMA,ALQSO
      COMMON/FNORM/FNORM
      COMMON/HYDROS/HSCALE,DEN1,XEQ,AMEAN
      COMMON/RQW/TEFF,RQ
      COMMON/PULSAR/RLTOT,ALFA,EMIN,EMAX
      common/obstar/rltotob,teffob
      COMMON/CINOUT/INOUT,IPULS
      COMMON/CICS/RSHOCK,ICS
      COMMON/QSOM/qso,IAGN,ISTAT,ISOBOL
      COMMON/INSPEC/IBLACKB,IPULSSP,IFERNET,INLR,IFERMAT,ICRAB,
     &     IFREEFREE,imewe,iobstar
      COMMON/FRE/NINT,JMIN,JJ
      COMMON/INT/FL(2,NE1:NE2),SI(14,27,NE1:NE2),EM(NE1:NE3),E(NE1:NE3)
      common/specyn/ispec
      common/initmw/initmw,initr2
      common/amw/emw(1000),fj(1000),jmaxmw
      common/contpoint/otscon(100),ncont,jcont(100),iioncont(100)
      dimension flnorm(1000)
      DIMENSION EI(10000),FSH(10000),FSHM(NE1:NE3)
      DATA PI/3.1415926d0/,ELCH/1.60219d-12/,AMU/1.660531d-24/
      DIMENSION RYCRB(10),FNUCRB(10)
*     Crab Nebula continuum from Davidson and Fesen, 1985 Ann Rev
*     second vector is F_nu as seen from Earth
      DATA RYCRB/-5, -3.284d0, -2.823d0, -0.959d0, -0.137d0, 0.863d0, 
     &     1.863d0, 3.863d0, 6.176d0, 6.869d0 /
     &  FNUCRB/-20.423d0, -20.860d0, -20.677d0, -22.308d0, -22.721d0,
     &     -23.650d0, -25.192d0, -27.395d0, -30.682d0, -31.780d0/
C
C     FMEAN = MEAN INTENSITY IN ERGS / CM**2 S EV
C
      FMEAN=0.
      IF(IPULSSP.EQ.1.AND.ICRAB.EQ.0) THEN
c       PULSAR SPECTRUM
        CON=RLTOT*(1.-ALFA)/(EMIN*((EMAX/EMIN)**(1.-ALFA)-1.))
        RL0=CON*(EMIN/E1)**ALFA
        if(iobstar.eq.1) then
c OB-star instead of pulsar in the ceter
c          rltotob,teffob
c use pure black-body
           hnukt=1.602e-12*e1/(1.38e-16*teffob)
           if(hnukt.ge.300.) then
              hnukt=300.
           endif
           rl0=2.80e15*rltotob*e1**3/(teffob**4*(exp(hnukt)-1.))
        endif
        FMEAN=RL0/(4.*PI**2*(R1*1.e15)**2)
      ENDIF
      IF(ICS.EQ.1.AND.IPULSSP.NE.1 .and.ifreefree.ne.1) THEN
c     in which bin is e1?
        DO J=JMIN,JJ
          IF(E1.GT.E(J).AND.E1.LT.E(J+1)) THEN
            FMEAN=FSHM(J)
            GOTO 66
          ENDIF
        ENDDO
 66     CONTINUE
      ENDIF
      IF(INOUT.EQ.1) THEN
        W=0.5*(1.-SQRT(1.-1./RQ**2))
      ELSE
        W=1.
      ENDIF
      IF(ICS.EQ.1.AND.IPULSSP.NE.1.and.ifreefree.eq.0) W=1.
      FMEAN=W*FMEAN
C
C     D:O IN ERGS / CM**2 S HZ
C
      IF(M.EQ.2) FMEAN=FMEAN*4.138E-15
C
C     OCCUPATION NUMBER
C
      IF(M.EQ.3) FMEAN=2.00E-11*FMEAN/E1**3
      RETURN
      END


      SUBROUTINE SPEC(TEMP)
C     ***********************************************************
C     *****
C     FOR THE FIRST ITERATION (NITSPH=1) SOLVE THE EQUATION OF TRANSFER
C         IN THE OUTWARD ONLY APPROXIMATION. FOR THE NEXT ITERATIONS
C         (NITSPH > 1) CALCULATE THE CONTINOUS OPACITIES FOR THE SHELL.
C     *****
C     ***********************************************************
      IMPLICIT REAL*8(A-H,O-Z)
c      PARAMETER (MD=350,MDP1=MD+1)
      PARAMETER (NL=340,NLP1=NL+1)
      PARAMETER (NFEL=3000)
      include 'param'
c      PARAMETER (NE1=-100,NE2=250,NE3=NE2+1)
      SAVE
      common/ionx/xion(md,14,27)
      COMMON/INSPEC/IBLACKB,IPULSSP,IFERNET,INLR,IFERMAT,ICRAB,
     &     IFREEFREE,imewe,iobstar
      COMMON/UVBL/IUVBLANK
      parameter (nlp=30000)
      COMMON/EQUIV/RADF,FR(2,nlp),W(nlp),CIN(nlp),WLI(nlp+NFEL),FB(300),tauline(nlp+nfel)
      COMMON/IND/I
      COMMON/SPH/ISPH
      COMMON/CINOUT/INOUT,IPULS
      COMMON/FRE/NINT,JMIN,JJ
      COMMON/PHY/DEN(MD)
      COMMON/DXA/DX(MD)
      COMMON/INT/FL(2,NE1:NE2),SI(14,27,NE1:NE2),E1(NE1:NE3),E(NE1:NE3)
      COMMON/HPOP/XN(6,NL),XN1,XN2,XN3
      COMMON/HYPOP/XNH(NL)
      COMMON/LITER/N
      COMMON/A26/SIGO(20,NE1:NE2),SIGCA(10,NE1:NE2),SIGHE(16,NE1:NE2),
     &SIGFE(NL,NE1:NE2),sigfei(nl,ne1:ne2),SIGH(50,NE1:NE2),SIGHE2(50,NE1:NE2)
      COMMON/DIF/EM(MD,NE1:NE2),TAU(MD,NE1:NE2),TAUTOT(MD,NE1:NE2),
     &           EMC(MD,NE1:NE2)
      COMMON/TEM/TE(MD)
      COMMON/ELEC/DEL(MD)
      COMMON/SIK/SK(14,27,NE1:NE2)
      common/abl/abn(15)
      COMMON/SPECT/TEL,FD(MD,NE1:NE2),F0(NE1:NE2),IPAR
      common/sorc/so(NE1:NE2),KM(NE1:NE2)
      COMMON/A16/SIG(NL),GA(NL),BP(NL)
      COMMON/NLEV/e000,NION,NH,NP1H,NHMAX,NHMIN
      COMMON/A14/C(NL,NL),CI(NL),G(NL),EN(NL),A(NL,NL),
     &WL(NL,NL),DCDT(NL,NL),DCIDT(NL)
      COMMON/CICS/RSHOCK,ICS
      common/linecont/opaclhe(100)
      common/fill/filling
      common/linepoint/nlines,jline(nlp),iion(nlp)
      common/lineopac/totopl(nlp)
      common/contpoint/otscon(100),ncont,jcont(100),iioncont(100)
      COMMON/RADIE/R(0:MD)
      common/ize/izero
      DIMENSION TA(NE1:NE2),XP(100),S(MD,NE1:NE2),opac(105,ne1:ne2)
      DIMENSION SIGM(5,10,NE1:NE2)
      DATA PI/3.1415926E0/,ELCH/1.60219E-12/,AMU/1.660531E-24/
      dimension nionel(14)
      data nionel/1,2,6,7,8,10,11,12,13,14,16,18,20,26/
      K=91
      IF(XN1.LE.0.) XN1=1.E-37
      K=92
      IF(XN2.LE.0.) XN2=1.E-37
      K=93
      IF(XN3.LE.0.) XN3=1.E-37

      RADF=0.

      izero=0

      CALL ATDATO
      DO J=JMIN,JJ
        DO K=4,9
            SIGM(2,K,J)=1.E-18*SIGOX(K,E1(J))
        ENDDO
      ENDDO
      ftotq=0.

      DO J=JMIN,JJ
C     ***********************************************************
C     *****
C     CALCULATE OPTICAL DEPTHS AT EACH ENERGY
C     *****
C     ***********************************************************
         TA(J)=0.0E0
         DTAM=0.

         do iel=3,14
            do ion=1,nionel(iel)
               dta = abn(iel)*xion(i,iel,ion)*
     &              (si(iel,ion,j)+sk(iel,ion,j))
               IF(DTA.GT.DTAM) THEN
                  DTAM=DTA
               ENDIF
               TA(J)=TA(J)+DTA
            enddo
         enddo

C     EXCITED O I CONTINUA
         DO K=4,9
            DTAB=ABn(5)*XN(2,K)*SIGM(2,K,J)
            TA(J)=TA(J)+DTAB

            IF(DTAB.GT.DTAM) THEN
               DTAM=DTAB
            ENDIF
         ENDDO
C     

C     TAX=TAU(I)=OPTICAL DEPTH OF SLAB BETWEEN I AND I+1
C
      
         TAX=DX(I)*TA(J)*DEN(I)*filling
         TAU(I,J)=TAX
         TAUTOT(I,J)=TAUTOT(I-1,J)+TAU(I,J)
         IF(TA(J).GT.0.) S(I,J)=EM(I,J)*DEN(I)/TA(J)
         IF(TA(J).LE.0.) S(I,J)=0.
         SO(J)=S(I,J)
         CALL RADTRANPULS(0,J,TAX,S)      
         ftotq=ftotq+fl(2,j)*(e(j+1)-e(j))
      ENDDO
c     jline(l) = energy bin of line l 
c     k = ion number
      do l=1,nlines
         totopl(l)=0.
c     if energy of line e(jmin) = jline(l) = -999 -> put opacity = 0 
         if(jline(l).gt.-999) then
            do k=1,105
               totopl(l)=totopl(l)+opac(k,jline(l))
c               opacl(l,k)=opac(k,jline(l))
            enddo
         endif
      enddo
      
      RETURN
      END


      SUBROUTINE opac_source(iqq,imode,dx,den,TEMPq,xel)
C     ***********************************************************
C     *****
C     FOR THE FIRST ITERATION (NITSPH=1) SOLVE THE EQUATION OF TRANSFER
C         IN THE OUTWARD ONLY APPROXIMATION. FOR THE NEXT ITERATIONS
C         (NITSPH > 1) CALCULATE THE CONTINOUS OPACITIES FOR THE SHELL.
c
c  imode = 1 : update tautot etc
c          0   only calculate opacity
C     *****
C     ***********************************************************
      IMPLICIT REAL*8(A-H,O-Z)
c      PARAMETER (MD=350,MDP1=MD+1)
      PARAMETER (NL=340,NLP1=NL+1)
      PARAMETER (NFEL=3000)
      include 'param'
c      PARAMETER (NE1=-100,NE2=250,NE3=NE2+1)
      SAVE
      COMMON/INSPEC/IBLACKB,IPULSSP,IFERNET,INLR,IFERMAT,ICRAB,
     &     IFREEFREE,imewe,iobstar
      COMMON/UVBL/IUVBLANK
      parameter (nlp=30000)
      COMMON/EQUIV/RADF,FR(2,nlp),W(nlp),CIN(nlp),WLI(nlp+NFEL),FB(300),tauline(nlp+nfel)
      COMMON/IND/I
      COMMON/mdIND/kradius
      COMMON/SPH/ISPH
      COMMON/CINOUT/INOUT,IPULS
      COMMON/FRE/NINT,JMIN,JJ
      COMMON/INT/FL(2,NE1:NE2),SI(14,27,NE1:NE2),E1(NE1:NE3),E(NE1:NE3)
      COMMON/HPOP/XN(6,NL),XN1,XN2,XN3
      COMMON/HYPOP/XNH(NL)
      COMMON/LITER/N
      COMMON/A26/SIGO(20,NE1:NE2),SIGCA(10,NE1:NE2),SIGHE(16,NE1:NE2),
     &SIGFE(NL,NE1:NE2),sigfei(nl,ne1:ne2),SIGH(50,NE1:NE2),SIGHE2(50,NE1:NE2)
      COMMON/DIF/EM(MD,NE1:NE2),TAU(MD,NE1:NE2),TAUTOT(MD,NE1:NE2),
     &           EMC(MD,NE1:NE2)
      COMMON/TEM/TE(MD)
      COMMON/ELEC/DEL(MD)
      COMMON/SIK/SK(14,27,NE1:NE2)
      common/abl/abn(15)
      common/ionx/xion(md,14,27)
      COMMON/SPECT/TEL,FD(MD,NE1:NE2),F0(NE1:NE2),IPAR
      common/sorc/so(NE1:NE2),KM(NE1:NE2)
      COMMON/A16/SIG(NL),GA(NL),BP(NL)
      COMMON/NLEV/e000,NION,NH,NP1H,NHMAX,NHMIN
      COMMON/A14/C(NL,NL),CI(NL),G(NL),EN(NL),A(NL,NL),
     &     WL(NL,NL),DCDT(NL,NL),DCIDT(NL)
      COMMON/CICS/RSHOCK,ICS
      common/linecont/opaclhe(100)
      common/fill/filling
      common/timecheck/time,itime
      common/linepoint/nlines,jline(nlp),iion(nlp)
      common/lineopac/totopl(nlp)
      common/contpoint/otscon(100),ncont,jcont(100),iioncont(100)
      COMMON/RADIE/R(0:MD)
      common/ize/izero
      common/opacity/TA(NE1:NE2),S(MD,NE1:NE2),copac(md,ne1:ne2)
      common/debug/ideb
      dimension nionel(14)
      data nionel/1,2,6,7,8,10,11,12,13,14,16,18,20,26/
      DIMENSION XP(100),opac(ne1:ne2)
      DIMENSION SIGM(5,10,NE1:NE2)
      DATA PI/3.1415926E0/,ELCH/1.60219E-12/,AMU/1.660531E-24/
      K=91
      IF(XN1.LE.0.) XN1=1.E-37
      K=92
      IF(XN2.LE.0.) XN2=1.E-37
      K=93
      IF(XN3.LE.0.) XN3=1.E-37
      RADF=0.
      CALL ATDATO
      DO J=JMIN,JJ
        DO K=4,9
            SIGM(2,K,J)=1.E-18*SIGOX(K,E1(J))
        ENDDO
      ENDDO
      ftotq=0.
      DO J=JMIN,JJ
C     ***********************************************************
C     *****
C     CALCULATE OPTICAL DEPTHS AT EACH ENERGY
C     *****
C     ***********************************************************
         TA(J)=0.0E0
         DTAM=0.
         opac(j)=0.

         do iel=3,14
            do ion=1,nionel(iel)
               if(xion(i,iel,ion).gt.1.e-15) then
                  dtfb = abn(iel)*xion(i,iel,ion)*
     &                 (si(iel,ion,j)+sk(iel,ion,j))
               else
                  dtfb = 0.
               endif
               opac(j)=dtfb+opac(j)
              IF(DTfb.GT.DTAM) THEN
                  DTAM=DTfb
                  dtfbm=dtfb
               ENDIF
               TA(J)=TA(J)+DTfb
            enddo
         enddo
C     EXCITED O I CONTINUA
         DO K=4,9
            DToi=ABn(5)*XN(2,K)*SIGM(2,K,J)
            TA(J)=TA(J)+DToi
            IF(DToi.GT.DTAM) THEN
               DTAM=DToi
            ENDIF
         ENDDO
C     
         if(imode.eq.1) then

C     TAX=TAU(I)=OPTICAL DEPTH OF SLAB BETWEEN I AND I+1
C
      
            TAX=abs(dx)*TA(J)*DEN*filling
            TAU(I,J)=TAX
            if(i.ge.2) then
               TAUTOT(I,J)=TAUTOT(I-1,J)+TAU(I,J)
            else
               TAUTOT(I,J)=TAU(I,J)
            endif

            IF(TA(J).GT.0.) then
               S(I,J)=EM(I,J)*DEN/TA(J)
            else
               S(I,J)=0.
            endif
            SO(J)=S(I,J)
         endif
      ENDDO

c
c     calculate the total and partial (due to ion k) opacity at line l
c     jline(l) = energy bin of line l 
c     k = ion number

      do l=1,nlines
         totopl(l)=0.
         if(jline(l).gt.ne2) then
         endif
c     if energy of line e(jmin) = jline(l) = -999 -> put opacity = 0 
c!!         if(l.eq.88) then
c         if(jline(l).ge.2) then
         if(jline(l).gt.-999) then
            totopl(l)=opac(jline(l))
         else
            totopl(l)= 0.
         endif
      enddo
      RETURN
      END

      SUBROUTINE diffuse_pl(imode,icont,IMAX,ishock,te,totfl)
      IMPLICIT REAL*8(A-H,O-Z)
      include 'param'
c      PARAMETER (NE1=-100,NE2=250,NE3=NE2+1)
C      *************************************************************
C     THIS ROUTINE CALCULATES THE DIFFUSE EMISSION FOR GIVEN EMISSION
C     RATE (PER STERADIAN) AND OPTICAL DEPTHS OF THE DIFFERENT SLABS.
C     THE METHOD IS THE SAME AS THAT USED BY WILLIAMS (1967).
C     FL0(MA,J), MA=1,2 IS THE FLUX FROM THE SURFACE AND THE INNER
C     BOUNDARY. THE OPTICAL DEPTH POINTS ARE CALCULATED FROM THE SUR-  
C     FACE WITH TAU(1,J)=0.
C     CENTERING OF VARIABLES                                           *
c     FD = mean intensity J
C     *
C     I     1         2         3         4                            *
C           +         +         +         +                            *
C          R(1)      R(2)      R(3)      R(4)                          *
C         taut(1)   taut(2)  taut(3)    taut(4)                        *
C                                                                      *
C                +         +          +       +                        *
C              DX(2)     DX(3)      DX(4)                              *
C              EM(2)     EM(3)      EM(4)                              *
C               S(2)      S(3)       S(4)                              *
C               J(2)      J(3)       J(4)                              *
C             tmid(2)    tmid(3)    tmid(4)                            *
C                                                                      *
C      *****
C      *************************************************************
      COMMON/FRE/NINT,JMIN,JJ
      COMMON/DIF/EM(MD,NE1:NE2),TAU(MD,NE1:NE2),TAUTOT(MD,NE1:NE2),
     &           EMC(MD,NE1:NE2)
      COMMON/SPECT/TEL,FD(MD,NE1:NE2),F0(NE1:NE2),IPAR
      COMMON/contspec/fc(MD,NE1:NE2),flc(2,ne1:ne2)
      COMMON/PHY/DEN(MD)
      COMMON/DXA/DX(MD)
      COMMON/ELEC/DEL(MD)
      common/fill/filling
      common/trq/rtr(2,md),drtr(2,md),dentr(2,md)
      common/diffsh/fdshock(ne1:ne2)
      COMMON/INT/FL(2,NE1:NE2),SI(14,27,NE1:NE2),E1(NE1:NE3),E(NE1:NE3)
      COMMON/SIK/SK(14,27,NE1:NE2)
      common/timecheck/time,itime
      common/preion/ipre
      common/tauelec/tau_tot_es,tau_elec_sh,tau_elec(md)
      common/taushok/taush(md)
      common/opacity/TAop(NE1:NE2),Ss(MD,NE1:NE2),copac(md,ne1:ne2)
      COMMON/mdIND/kradius
      data pi/3.14159/
      DIMENSION FL0(2,NE1:NE2)
      dimension sf(ne1:ne2),taj(ne1:ne2),tmid(md)
      dimension emm(md),copm(md),denm(md)

      boost=1.
      
      DO J=JMIN,JJ
         DO  MA=1,2
            FL0(MA,J)=0.
         enddo
      enddo

      totfl=0.
      
      DO J=JMIN,JJ

         do k=1,imax

            tmid(k+1) = (TAUTOT(k,J)+TAUTOT(K+1,J))/2.

         enddo

         if(imode.eq.1) then
            
            i = imax

            if(icont.eq.0) then

               FD(I,J)=0.0E0

            elseif(icont.eq.1) then

               fc(I,J)=0.0E0

            endif
            
            IMM1=I-1
c!!!           
c            DO K=2,IMM1
           DO K=2,i
C
C     CALCULATE THE OPTICAL DEPTH BETWEEN SHELL I AND SHELL K
C
               T2=abs(TAUTOT(I,J)-TAUTOT(K,J))
               if(t2.eq.0.) then
                  t2=1.e-30
               endif
               T3=abs(TAUTOT(I,J)-TAUTOT(K-1,J))

                  
               IF(DRTR(2,K).LE.0..and.ipre.eq.0) then
                  do ii=1,imax+1
                     write(6,*)ii,drtr(2,ii),dentr(2,ii)
                  enddo
                  stop
               endif
                  
               IF(DENTR(2,K).LE.0.) WRITE(6,*)' DIFN',K,DENTR(2,K)
c TA = opacity (mean cross section cm^2 )               
               TA=abs(TAUTOT(K,J)-TAUTOT(K-1,J))/(filling*
     &              abs(DRTR(2,K))*DENTR(2,K))

               TA=ABS(TA)
C     
C     
C     SOURCE FUNCTION EVALUATED IN THE MIDDLE BETWEEN K-1 AND K
C     
               IF(TA.LE.0.) then
                  EMTA=0.
               else
                  if(icont.eq.0) then
                     EMTA=EM(K,J)/TA
                  elseif(icont.eq.1) then
                     EMTA=EMc(K,J)/TA
                  endif
               endif
               
c     one more density in j (=em*den)
               
               SOURCe = DENTR(2,K)*EMTA

               sf(j)=source
               taj(j)=ta

               icase = 0

c               if(t2.gt.1.e-3.or.(t3-t2).gt.1.e-12) then
               if((t3-t2).gt.0.01) then

                  FDIF=(E2(T2)-E2(T3))*SOURCe/2.

                  icase = 1

               elseif((t3-t2).ge.0..and.(t3-t2).le.0.01) then

                  tm=(t2+t3)/2.

                  if(icont.eq.0) then
                     
c use correction for spherical effects (eg Cox 1972, Shull & McKee 1979)

                     f1 = exp1(tm)/2.

                     if(f1.gt.3.) then

                        f1 = 3.

                     endif
                     
                     fdif = em(k,j)*dentr(2,k)**2*abs(drtr(2,k))* f1

                     icase = 2

                  elseif(icont.eq.1) then

c use correction for spherical effects (eg Cox 1972, Shull & McKee 1979)

                     f1 = exp1(tm)/2.

                     if(f1.gt.3.) then

                        f1 = 3.

                     endif

                     fdif = emc(k,j)*dentr(2,k)**2*abs(drtr(2,k))* f1

                     icase = 3

                  endif
               else

                  fdif=0.
                  
                  icase = 4

               endif

               if(k.eq.i) then

c for the emission from the shell itself include the emission 
c     from both directions

                  fdif = 2.*fdif

                  if(t3.gt.100.) then

                     fdif = source
                     
                     icase  = 5

                  endif
                  
               endif

               tm=(t2+t3)/2.
               fdif2 = source*exp1(tm)*(t3-t2)/2.

c use correction for spherical effects (eg Cox 1972, Shull & McKee 1979)

               f1 = exp1(tm)/2.

               if(f1.gt.3.) then
                  
                  f1 = 3.
                  
               endif


               if(icont.eq.0) then

                  fdif3 = em(k,j)*dentr(2,k)**2*abs(drtr(2,k))* f1

               elseif(icont.eq.1) then

                  fdif3 = emc(k,j)*dentr(2,k)**2*abs(drtr(2,k))* f1

               endif


               fapp=(1.d0-.577216d0)*(t3-t2)-t3*log(t3)+t2*log(t2)-
     &              0.5*(t2**2-t3**2)+
     &              (.24991d0-1.d0/6.d0)*(t2**3-t3**3)-
     &              (.0552d0-1.d0/24.d0)*(t2**4-t3**4)
     

               FDIF1=fapp*SOURCe/2.

               if(icont.eq.0) then
                  FD(I,J)=FD(I,J)+FDIF

               elseif(icont.eq.1) then
                  Fc(I,J)=fc(I,J)+FDIF
               endif
               

               ma=1

               E32=0.
               E33=0.
               IF(T2.LT.300.) then
                  E32=(EXP(-T2)-T2*E2(T2))/2.
               endif
               IF(T3.LT.300.) then
                  E33=(EXP(-T3)-T3*E2(T3))/2.
               endif
c     if(t2.gt.1.e-3.or.(t3-t2).gt.1.e-12) then
               if(t2.gt.1..or.t3.gt.1.) then
                  flux=(E32-E33)*source/2.
               elseif(t3-t2.gt.0.d0) then
                  fapp=(t3-t2)-t3**2*(1.5-.577216-log(t3))/2.+
     &                 t2**2*(1.5-.577216-log(t2))/2.
                  flux=fapp*source/2.
               else
                  flux=0.
               endif

C*****MULTPLY BY AN ARBITARY FACTOR 2
c               FLUX=2.*FLUX
               if(icont.eq.0) then
                  FL0(MA,J)=FL0(MA,J)+FLUX
               endif             
            enddo

c add primary flux from shock (only for 1st iteration of the pre shock gas)

            if(ipre.eq.1) then
               if(icont.eq.0) then
                  fd(i,j) = fd(i,j) + fdshock(j)*e2(tautot(i,j))
               endif
            endif

            if(icont.eq.0) then
               fl(2,j) = fd(i,j)

            elseif(icont.eq.1) then
               flc(2,J)=fc(I,J)
            endif

            if(fdif.lt.0.) then
               write(6,9247)i,j,e1(j),drtr(2,i),dentr(2,i),tautot(i,j),em(i,j),
     &              ta,source,fdif,fd(i,j),fl0(1,j),fdshock(j)
 9247          format(' fd lt 0 ',2i5,1pe12.3,20e12.3)
            endif

c     correction of mean intensity for electron scattering
            flc(2,j)=boost*flc(2,j)
            fl(2,j)=boost*fl(2,j)

         elseif(imode.eq.2.or.imode.eq.3) then

C     I     1         2         3         4                            *
C           +         +         +         +                            *
C          R(1)      R(2)      R(3)      R(4)                          *
C         taut(k)   taut(2)  taut(3)    taut(4)                        *
C          EM(1)    EM(2)     EM(3)      EM(4)                        *
C        copac(1)   copac(2)  copac(3)  copac(4)                        *

C                    J(2)      J(3)       J(4)                         *
C                                                                      *
C                +         +          +       +                        *
C              DX(2)     DX(3)      DX(4)                              *
c              emm(2)    emm(3)     emm(4)
C               S(2)      S(3)       S(4)                              *

C             tmid(1)    tmid(3)    tmid(i)                            *
C                                                                      *

c  We want mean intensities at the grid points for the calc. of ionization and 
c     temperature at each point. To calc. this we need emissivities and opacities 
c     at the midpoints. This gives the contribution from each shell to the mean int. 

c first calc. emissivities, densities  and opacities at midpoints

c note that emm(i+1), etc is evaluated at the midpoint between i and i+1


            emc(1,j) = emc(2,j)

            do i=1,imax
               
               if(imode.eq.2) then
                  emm(i+1) = (em(i,j) + em(i+1,j))/2.
               elseif(imode.eq.3) then
                  emm(i+1) = (emc(i,j) + emc(i+1,j))/2.
               endif

               copm(i+1) = (copac(i,j) + copac(i+1,j))/2.
               denm(i+1) = (dentr(2,i) + dentr(2,i+1))/2.
               
            enddo

            if(copm(2).eq.0.) then
               copm(2) = copm(3)
            endif

            DO I=2,IMAX

               if(imode.eq.2) then
                  FD(I,J)=0.0E0
               elseif(imode.eq.3) then
                  fc(I,J)=0.0E0
               endif
               IMM1=IMAX-1
               DO K=2,imax
C
C     CALCULATE THE OPTICAL DEPTH BETWEEN SHELL I AND SHELL K
C

                  T2=abs(tmid(i)-TAUTOT(K,J))
                  T3=abs(tmid(i)-TAUTOT(K-1,J))



                  t2 = 0.


                  if(i.le.k) then

c t2 = optical depth from i to k-1

c t3 = optical depth from i to k

                     do kk=i+1,k-1
                        t2 = t2 + abs(drtr(2,kk))*denm(kk)*copm(kk)
                     enddo

                     t3 = t2 + abs(drtr(2,k))*denm(k)*copm(k)

                     t3mt2 = abs(drtr(2,k))*denm(k)*copm(k)

                     if(i.eq.ishock.and.j.eq.3) then

                        taush(k) = t3

                     endif

                  elseif(k.lt.i) then

c t2 = optical depth from i to k

c t3 = optical depth from i to k-1

                     do kk=k+1,i
                        t2 = t2 + abs(drtr(2,kk))*denm(kk)*copm(kk)
                     enddo

                     t3 = t2 + abs(drtr(2,k))*denm(k)*copm(k)

                     t3mt2 = abs(drtr(2,k))*denm(k)*copm(k)

                     if(i.eq.ishock.and.j.eq.3) then

                        taush(k) = t2

                     endif

                  endif

                  IF(DRTR(2,K).LE.0..and.ipre.eq.0) then
                     WRITE(6,*)' DIFX',K,DRTR(2,K)
                     do ii=1,imax+1
                        write(6,*)ii,drtr(2,ii),dentr(2,ii)
                     enddo
                     stop
                  endif
                  
                  IF(DENTR(2,K).LE.0.) WRITE(6,*)' DIFN',K,DENTR(2,K)
                  TA=abs(TAUTOT(K,J)-TAUTOT(K-1,J))/(filling*
     &                 abs(DRTR(2,K))*DENTR(2,K))
                  TA=ABS(TA)

C     
C     
C     SOURCE FUNCTION EVALUATED IN THE MIDDLE BETWEEN K-1 AND K
C     
                  IF(TA.LE.0.) then
                     EMTA=0.
                  else
                     EMTA=EMm(K)
                  endif
                  SOURCe= DENm(K)*EMTA/copm(k)

                  icase=0

                  if(t3mt2.gt.0.01) then

                     FDIF=abs(E2(T2)-E2(T3))*SOURCe/2.
                     fdif1 =fdif
                     icase=1

                  elseif(t3mt2.le.0.01.and.t3mt2.gt.0.) then

                     tm=(t2+t3)/2.
                     fdif = emm(k)*denm(k)**2*abs(drtr(2,k))*exp1(tm)/2.

                     thickness=1.e16
                     radius = 1.e18
                     corr = -0.5*log(0.5*thickness/radius)
                     fdif = emm(k)*denm(k)**2*abs(drtr(2,k))*corr/2.
                     icase=2

                  else

                     fdif=0.

                     icase=3
                  endif

                  if(k.eq.i) then

c for the emission from the shell itself include the emission 
c     from both directions

                     fdif = 2.*fdif

                     t2 = 1.e-20
                     
                     t3 = abs(tautot(i,j)-tmid(i))

                     FDIF=(E2(T2)-E2(T3))*SOURCe

                     icase = 4

                     if(abs(t3).gt.0.01) then

                        FDIF=abs(1.-E2(T3))*SOURCe
                        fdif1 =fdif
                        icase=41

                     elseif(abs(t3).le.0.01.and.abs(t3).gt.0.) then

c Note that t3 = half the total optical depth over the shell

                        fac = 1.-0.577216-log(t3)
                        fdif = emm(k)*denm(k)**2*abs(drtr(2,k))*fac/2.
                        icase=42

                     else

                        fdif=0.

                        icase=43
                     endif


                     if(t3.gt.100.) then

                        fdif = source

                     endif
                     
                     fdif4 = fdif
c!!!

                     fdif=0.
                  
                  endif

                  if(imode.eq.2) then
                     FD(I,J)=FD(I,J)+FDIF
                     fl(2,j) = fd(i,j)
                  elseif(imode.eq.3) then
                     fc(I,J)=fc(I,J)+FDIF
                     flc(2,J)=fc(I,J)
                  endif

                  if(i.eq.-100.and.(j.eq.-48.or.j.eq.-49)) then

                     write(6,9237)i,k,icase,drtr(2,k),dentr(2,k),denm(k),
     &                    tautot(k,j),t2,
     &                    t3,t3-t2,tm,em(k,j),emc(k,j),emm(k),
     &                    ta,copac(k,j),fdif,source,fc(i,j),fd(i,j)
 9237          format(' fd ',3i5,1pe12.3,20e12.3)
                  endif
                  if(icase.eq.0) then
                     write(6,9247)i,k,icase,drtr(2,k),dentr(2,k),tautot(k,j),t2,
     &                    t3,t3-t2,tmid(i),em(k,j),ta,fdif,source,
     &                    fd(i,j)
                     stop
                  endif
               enddo
            enddo
         endif
      enddo
      
         xtot=0.
         do i=1,imax
            xtot = xtot + drtr(2,i)
            do j=jmin,jj
               if(i.eq.imax-1) then
c include these for ionizing spectrum
c     write(76,9240)i,j,xtot,drtr(2,i),dentr(2,i),
c     &                 e1(j),tautot(i,j),emc(i,j),em(i,j),copac(i,j),
c     &                 fc(i,j),fd(i,j)
 9240             format(2i5,1pe15.7,3e12.4,12e12.4)
               endif
            enddo
         enddo

         if(imode.eq.2.or.imode.eq.3) then


         do i=1,imax
            totfl=0.
            pht=0.
            pht2=0.
            pht3=0.
            do j=3,jj
               totfl=totfl+fd(i,j)*(e(j)-e(j-1))
               pht=pht+4.*pi*fd(imax,j)*sk(1,1,j)*(e(j+1)-e(j))/
     &              (1.602e-12*e1(j))
               pht2=pht2+4.*pi*fd(imax,j)*sk(2,1,j)*(e(j+1)-e(j))/
     &              (1.602e-12*e1(j))
               pht3=pht3+4.*pi*fd(imax,j)*sk(2,2,j)*(e(j+1)-e(j))/
     &              (1.602e-12*e1(j))
            enddo
         enddo
      elseif(imode.eq.1) then
         totfl=0.
         pht=0.
         pht2=0.
         pht3=0.
         fltot=0.


         if(icont.eq.0) then
            do j=jmin,jj
               totfl=totfl+fd(imax,j)*(e(j)-e(j-1))
               if(j.ge.2) then 
                  fltot=fltot+fl0(1,j)*(e(j)-e(j-1))
               endif
c            write(76,992)j,e1(j),1.24e4/e1(j),em(imax,j),
c     &           tautot(imax,j),fl(2,j),fl0(1,j),fltot,sf(j),taj(j)
 992           format(i5,1pe12.3,0pf13.4,1pe12.3,10e12.3)
               if(j.ge.3) then
                  pht=pht+4.*pi*fd(imax,j)*sk(1,1,j)*(e(j+1)-e(j))/
     &                 (1.602e-12*e1(j))
                  pht2=pht2+4.*pi*fd(imax,j)*sk(2,1,j)*(e(j+1)-e(j))/
     &                 (1.602e-12*e1(j))
                  pht3=pht3+4.*pi*fd(imax,j)*sk(2,2,j)*(e(j+1)-e(j))/
     &                 (1.602e-12*e1(j))
               endif
            enddo
         endif
      endif

c convert to Eddington flux

      totfl = 4.*pi*fltot

      RETURN



      
  997 WRITE(6,996)
  996 FORMAT(1X,'STOPP')
      STOP
      END


      SUBROUTINE BIS(TE,XEL,TMIN,IFAIL,IFPOP,ICONV,ISTOP)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/THER/A1,B1,TIN,TOL,E20
      COMMON/HRAT/RHYD,ZEL,XELQ,HYR,HEAT,COOL
C     IF ICONV = 1 AT EXIT THE FUNCTION AT THE ZERO IS LARGER THAN TOLERATED
C     IF IFAIL = 1 AT EXIT THE FUNCTION HAS THE SAME SIGN AT BOTH LIMITS
C     E20 = TOL. IN FUNCTION VALUE
      a2=0.
      a3=0.
      a4=0.
      b2=0.
      b3=0.
      b4=0.
      tme=0.
      iter=0
      init=0
      eps=1.
      eps4=1.e-2
      ICONV=1
      IFAIL=0
      XE=XEL
      MM=1
      A=A1
      B=B1
      IF(B1.LT.TMIN) GOTO 44
      FB=RAD(B1,XE,IFPOP)
      IF(IFPOP.EQ.1) GOTO 4
c first round redo the calc. to get a more stable result
      FB=RAD(B1,XE,IFPOP)
      IF(IFPOP.EQ.1) GOTO 4
      iit=0
      mw=0
C
C     LEFT BOUNDARY 
C
 9    CONTINUE
      IF(A1.LT.TMIN) GOTO 44
      FA=RAD(A1,XE,IFPOP)
      if(iit.eq.0) then
c first round redo the calc. to get a more stable result
         FA=RAD(A1,XE,IFPOP)
         iit=1
      endif
      if(ifpop.eq.1) a1=(a1+b1)/2.       
      if(ifpop.eq.1) mw=mw+1
      if(mw.ge.10) goto 4
      IF(IFPOP.EQ.1) GOTO 9
c 6    IF(FA*FB.GT.0.) WRITE(6,*)' SAME SIGN',MM,FA,FB
c     SAME SIGN: IFAIL=1
6     IF(FA*FB.GT.0.) IFAIL=1
      IF(MM.GT.40) GOTO 4
C     IF OK CONTINUE
      IF(IFAIL.NE.1) GOTO 1
      a4=a3
      b4=b3
      a3=a2
      b3=b2
      a2=a1
      b2=b1
      if(abs(a2-a3)/a2.lt.eps4.and.abs((fa2-fa)/fa).lt.eps4) then
      a1=0.95*a1
      A=A1
      fa2=fa
      IF(A1.LT.TMIN) GOTO 44
      FA=RAD(A1,XE,IFPOP)
      IF(IFPOP.EQ.1) GOTO 4
      IFAIL=0
      MM=MM+1
      GOTO 6
      endif
      if(abs(b2-b3)/b2.lt.eps4.and.abs((fb2-fb)/fb).lt.eps4) then
C     INCREASE UPPER LIMIT
      b1=1.075*b1
      b=b1
      fb2=fb
      IF(B1.LT.TMIN) GOTO 44
      Fb=RAD(b1,XE,IFPOP)
      IF(IFPOP.EQ.1) GOTO 4
      IFAIL=0
      MM=MM+1
      GOTO 6
      endif
      if(abs(a2-a4).lt.eps.and.abs(b2-b4).lt.eps) then
      delt=log10(a1/1.e2)/21.
      a10=a1
      do i=1,21
      b1=a1
      a1=a10*10.**(-(i-1)*delt)
      fb=fa
      IF(A1.LT.TMIN) GOTO 44
      fa=rad(a1,XE,IFPOP)
      if(fa*fb.lt.0.) goto 55
      enddo
55    a=a1
      b=b1
      goto 1
      else
      endif
      IF(ABS(FB).LT.ABS(FA)) GOTO 7
      B1=A1
      B=B1
      FB=FA
      IF(A1.LE.3500.) A1=0.95*A1
      IF(A1.GT.3500.) A1=0.95*A1
      A=A1
      fa2=fa
      IF(A1.LT.TMIN) GOTO 44
      FA=RAD(A1,XE,IFPOP)
      IF(IFPOP.EQ.1) GOTO 4
      IFAIL=0
      MM=MM+1
      GOTO 6
7     A1=B1
      A=A1
      FA=FB
      IF(A1.LE.3500.) B1=1.05*B1
      IF(A1.GT.3500.) B1=1.051*B1
      B=B1
      fb2=fb
      IF(B1.LT.TMIN) GOTO 44
      FB=RAD(B1,XE,IFPOP)
      IF(IFPOP.EQ.1) GOTO 4
      IFAIL=0
      MM=MM+1
      GOTO 6
C     NEW TEMP. BETWEEN A AND B
1     TME=(A+B)/2.
      mw=0
      a4=a3
      b4=b3
      a3=a2
      b3=b2
      a2=a
      b2=b
10    CONTINUE
      IF(TME.LT.TMIN) GOTO 44
      FM=RAD(TME,XE,IFPOP)
      if(ifpop.eq.1) tme=(a+tme)/2.       
      if(ifpop.eq.1) mw=mw+1
      if(mw.ge.10) goto 4
      IF(IFPOP.EQ.1) GOTO 10
      IF(FA*FM.GE.0.) THEN
C     IF FM HAS THE SAME SIGN AS FA THEN A=TME
c!            IF(ABS(FM).GT.ABS(FA)) THEN
c!                  A=A-ABS(TME-A)
c!                  FA=RAD(A,XE,IFPOP)
c!            ELSE
                  A=TME
                  FA=FM
c!            ENDIF
            GOTO 3
C     IF FM HAS A DIFFERENT SIGN FROM FA THEN B=TME
      ELSE
c!            IF(ABS(FM).GT.ABS(FB)) THEN
c!                  B=B+ABS(TME-B)
c!                  FB=RAD(B,XE,IFPOP)
c!            ELSE
                  B=TME
                  FB=FM
c!            ENDIF
      ENDIF
      if(abs(a2-a3)/a2.lt.eps4.and.abs((fa2-fa)/fa).lt.eps4) then
            a1=0.95*a1
            A=A1
            fa2=fa
            IF(A1.LT.TMIN) GOTO 44
            FA=RAD(A1,XE,IFPOP)
            IF(IFPOP.EQ.1) GOTO 4
            IFAIL=0
            MM=MM+1
            GOTO 6
      endif
      if(abs(b2-b3)/b2.lt.eps4.and.abs((fb2-fb)/fb).lt.eps4) then
C     INCREASE UPPER LIMIT
            b1=1.075*b1
            b=b1
            fb2=fb
            IF(B1.LT.TMIN) GOTO 44
            Fb=RAD(b1,XE,IFPOP)
            IF(IFPOP.EQ.1) GOTO 4
            IFAIL=0
            MM=MM+1
            GOTO 6
      endif
      if(abs(a2-a4).lt.1.e-10.and.abs(b2-b4).lt.eps) then
            delt=log10(a1/1.e2)/21.
            a10=a1
            do i=1,21
                  a1=a10*10.**(-(i-1)*delt)
                  IF(A1.LT.TMIN) GOTO 44
                  fa=rad(a1,XE,IFPOP)
                  if(fa*fb.lt.0.) goto 57
            enddo
57          a=a1
            goto 1
      endif
3     continue
      deltta=abs(a-b)/tme
c      IF(ABS(FM).LE.E20*ABS(HEAT).and.deltta.lt.tol) then
      IF(ABS(FM).LE.E20*ABS(HEAT)) then
         ICONV=0
         IFAIL=0
         GOTO 4
      endif
      MM=MM+1
      ITERM=2
      IF(MM.GE.20) then 
            if(abs(a-tme).gt.abs(b-tme)) then
                  a1=0.99*a1
                  IF(A1.LT.TMIN) GOTO 44
                  fa=rad(a1,XE,IFPOP)
                  a=a1
            else
                  b1=1.01*b1
                  IF(B1.LT.TMIN) GOTO 44
                  fb=rad(b1,xe,ifpop)
                  b=b1
            endif
            mm=1
            goto 1
      endif
      ITERM=1
C     GO ON WITH A NEW TEMPERATURE
      GOTO 1
4     TE=TME
      IF(ITER.EQ.2) WRITE(6,*)'ITERM'
      ISTOP=0
      RETURN
44    ISTOP=1
      END

      SUBROUTINE RADTRANPULS(IH,J,TAX,S)
C     ***********************************************************
C     *****
C     FOR THE FIRST ITERATION (NITSPH=1) SOLVE THE EQUATION OF TRANSFER
C         IN THE OUTWARD ONLY APPROXIMATION. FOR THE NEXT ITERATIONS
C         (NITSPH > 1) CALCULATE THE CONTINOUS OPACITIES FOR THE SHELL.
C     WORKS BOTH IN OUT AND OUT IN
C     HERE OUT TO IN
C     USED FOR ALL CS-INT. CALC.
C     IF IH = 1 CALCULATE FLUX & MEAN INTENSITY
C     IF IH = 1 CALCULATE MEAN INTENSITY
C     *****
C     ***********************************************************
      IMPLICIT REAL*8(A-H,O-Z)
c      PARAMETER (MD=350,MDP1=MD+1)
      include 'param'
      PARAMETER (NL=340,NLP1=NL+1)
      PARAMETER (NFEL=3000)
c      PARAMETER (NE1=-100,NE2=250,NE3=NE2+1)
      SAVE
      COMMON/IND/I
      COMMON/CINOUT/INOUT,IPULS
      COMMON/FRE/NINT,JMIN,JJ
      COMMON/INT/FL(2,NE1:NE2),SI(14,27,NE1:NE2),E1(NE1:NE3),E(NE1:NE3)
      COMMON/ITSPH/NITSPH
      COMMON/LITER/N
      COMMON/DIF/EM(MD,NE1:NE2),TAU(MD,NE1:NE2),TAUTOT(MD,NE1:NE2),
     &           EMC(MD,NE1:NE2)
      COMMON/RADIE/R(0:MD)
      COMMON/SPECT/TEL,FD(MD,NE1:NE2),F0(NE1:NE2),IPAR
      COMMON/DTAU/FLUX(NE1:NE2)
      COMMON/DIFFP/FH0(NE1:NE2),FHD(NE1:NE2)
      DIMENSION S(MD,NE1:NE2)
      DATA PI/3.1415926E0/,ELCH/1.60219E-12/,AMU/1.660531E-24/
C
C     PRIMARY FLUX
C
      IF(IPULS.EQ.1) THEN
            FPRIM=EXP(-TAUTOT(I,J))*FMEAN(1,E1(J))
      ELSE
            FPRIM=0.
      ENDIF
      IF(IH.EQ.1) GOTO 888
C     ***********************************************************
C     *****
C
C     FOR SECOND AND HIGHER ITERATIONS THE DIFFUSE FLUX IS 
C     ALREADY CALCULATED
C     MEAN INTENSITY = PRIMARY FLUX + DIFFUSE FLUX
c!    IF TAU > 1 USE J = SOURCE FCN.
C     *****
C     ***********************************************************
C     ***********************************************************
C     *****
C     USE THE HATCHETT, BUFF & MCCRAY TRANSPORT EQ. FOR SPHERICAL
C     GEOMETRY.
C     *****
C     ***********************************************************
      IF(N.EQ.1) THEN
        IF(TAX.LT.1.E-5) THEN
          EXTAU=1.-TAX+TAX**2/2.
        ELSEIF(TAX.LT.700.) THEN
          EXTAU=EXP(-TAX)
        ELSEIF(TAX.GE.700.) THEN
          EXTAU=0.
        ENDIF
        FL(2,J)=(FL(1,J)*EXTAU+S(I,J)*(1.-EXTAU))
     &                  /(1.+(R(I)-R(I-1))/R(I-1))**2.
        FD(I,J)=FL(2,J)-FPRIM
      ELSE
C
C       FOR SECOND AND HIGHER ITERATIONS THE DIFFUSE FLUX IS 
C       ALREADY CALCULATED
C       MEAN INTENSITY = PRIMARY FLUX + DIFFUSE FLUX
C
        FL(2,J)=FPRIM+FD(I,J)
      ENDIF
888   FHD(J)=FD(I,J)
      FH0(J)=FPRIM
      FLUX(J)=FL(2,J)
      RETURN
      END




