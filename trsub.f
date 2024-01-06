



      SUBROUTINE IONLABEL(LAB1)
      CHARACTER*8 LAB(200),LAB1(200)
      DATA LAB/'H I     ','HE I    ','HE II   ','O VI    ','O VII   '
     &        ,'O VIII  ','O V     ','C III   ','C IV    ','C V     '
     &        ,'C VI    ','C I     ','C II    ','O I     ','O II    '
     &        ,'O III   ','O IV    ','N I     ','N II    ','N III   '
     &        ,'N IV    ','N V     ','N VI    ','N VII   ','SI I    '
     &        ,'SI II   ','SI III  ','SI IV   ','SI V    ','SI VI   '
     &        ,'SI VII  ','SI VIII ','SI IX   ','SI X    ','SI XI   '
     &        ,'SI XII  ','SI XIII ','SI XIV  ','MG I    ','MG II   '
     &        ,'MG III  ','FE I    ','FE II   ','FE III  ','FE IV   '
     &        ,'AL I    ','AL II   ','AL III  ','AL IV   ','CA I    '
     &        ,'CA II   ','CA III  ','NA I    ','NA II   ','NA III  '
     &        ,'S I     ','S II    ','S III   ','S IV    ','S V     '
     &        ,'S VI    ','S VII   ','S VII   ','S IX    ','S X     '
     &        ,'S XI    ','S XII   ','S XIII  ','S XIV   ','S XV    '
     &        ,'S XVI   ','NE I    ','NE II   ','NE III  ','NE IV   '
     &        ,'NE V    ','NE VI   ','NE VII  ','NE VIII ','NE IX   '
     &        ,'NE X    ','AR I    ','AR II   ','FE V    ','FE VI   '
     &        ,'FE VII  ','FE VIII ','FE IX   ','FE X    ','FE XI   '
     &        ,'FE XII  ','FE XIII ','FE XIV  ','FE XV   ','O III 1D'
     &      ,5*'        ','BALMER  '
     &        ,'PASCHEN ','BRACKETT','PFUND   ','O I 4   ','O I 5   '
     &        ,'O I 6   ','O I 7   ','O I 8   ','O I 9   ','FF      '
     &        ,89*' '/
      DO K=1,111
            LAB1(K)=LAB(K)
c            write(0,9)k,lab(k)
 9          format('lab ',i5,a8)
      ENDDO
      RETURN
      END

      SUBROUTINE WLFORB(WLF)
      IMPLICIT REAL*8(A-H,O-Z)
      common/forbwl/wlft(300)
      common/wlforbk/wlfm(300)
      DIMENSION WLF(300),WLF1(71)
C     WAVELENGTHS OF ALL FORBIDDEN LINES
      DATA WLF1/5007.,2321.,4363.,3726.,2470.,7320.,6548.,3063.,
     &         5755.,6300.,2964.,5581.,5201.,3468.,10406.,0.,631800.,
     &         440560.,1455600.,9818.,4619.,8729.,6718.,4071.,10320.,
     &         252460.,563060.,1297016.,684932.,16360.,10995.,
     &         3704000.,6091000.,9069.,3722.,6312.,3869.,1815.,3342.,
     &         2423.,1602.,4714.,1575.,2975.,3346.,155500.,108640.,
     &         360200.,883300.,518000.,11287.,7002.,9264.,6157.,
     &         89910.,2.184e5,7135.8,7751.1,3006.1,3109.1,5191.8,
     &         4711.3,4740.2,2853.7,2868.2,7.741e5,564721.,7237.3,
     &         7170.6,7333.4,7262.8/
      nforb=110
      do k=1,nforb
         wlf(k)=wlfm(k)
      enddo
      DO K=1,71
            WLF(K)=WLF1(K)
      ENDDO
c!!!
      goto 1111
c  All lines which should not be included in summation
      do k=1,71
         if(k.ge.1.and.k.le.6) then
            wlf(k)=-1.
         elseif(k.ge.13.and.k.le.16) then
            wlf(k)=-1.
         elseif(k.ge.23.and.k.le.25) then
            wlf(k)=-1.
         elseif(k.ge.40.and.k.le.45) then
            wlf(k)=-1.
         elseif(k.ge.49.and.k.le.50) then
            wlf(k)=-1.
         endif
      enddo    
 1111 continue
      RETURN
      END
