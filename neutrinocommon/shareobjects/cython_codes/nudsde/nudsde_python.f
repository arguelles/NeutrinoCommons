      function dsde(E1,E2,neutype,targettype,intertype)

***********************************************************************
*** Program nusigma calculates the neutrino-nucleon cross sections
*** with a given parton distribution function (currently CTEQ6).
*** Both charged and neutral current cross sections on neutrons
*** and protons respectively are calculated. The cross sections on neutrons
*** and protons are calculated with routines written by Joakim Edsjo.
*** The isoscalar cross sections are calculated written by Joakim Edsjo,
*** but based on java-routines from Dima Chirkin and Wolfgang Rodhe (MMC).
*** The isoscalar routines are only implemented for testing purposes.
***
*** This routine returns the differntial cross section dsigma/dEmuon
*** where Emuon is the energy of the final state lepton.
***
*** Date: 2005-10-21
*** Joakim Edsjo, edsjo@physto.se
*** Modified by C. Arguelles to calculate 1 point.
*** ModDate : 2011-06-11
***********************************************************************
      implicit none

      include 'nupar.h'
      
      real*8 dsde
      real*8 NuCrossDiffl
      real*8 E1,E2
      integer neutype
      character targettype
      character*2 intertype

      call nusetup

      dsde = Nucrossdiffl(E1,E2,neutype,targettype,intertype)

      end


