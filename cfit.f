      subroutine cfit(iz,in,t,c)
* Version 2, March 24, 1997
******************************************************************************
*** This subroutine calculates rates of direct collisional ionization 
*** for all ionization stages of all elements from H to Ni (Z=28)
*** by use of the fits from G. S. Voronov, 1997, ADNDT, 65, 1
*** Input parameters:  iz - atomic number 
***                    in - number of electrons from 1 to iz 
***                    t  - temperature, K
*** Output parameter:  c  - rate coefficient, cm^3 s^(-1)
******************************************************************************
      implicit real*8(a-h,o-z)
      common/CF/CF(5,30,30)
      c=0.0
      if(iz.lt.1.or.iz.gt.28)return
      if(in.lt.1.or.in.gt.iz)return
      te=t*8.617385e-05
      u=cf(1,iz,in)/te
      if(u.gt.80.0)return
      c=cf(3,iz,in)*(1.0+cf(2,iz,in)*sqrt(u))/(cf(4,iz,in)+u)*
     &     u**cf(5,iz,in)*exp(-u)
      return
      end
*********************************************
