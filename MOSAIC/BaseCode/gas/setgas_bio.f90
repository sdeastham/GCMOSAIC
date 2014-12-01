
!***********************************************************************
! subroutine SetGas_Bio: sets up gas-phase species indices for
! the selected mechanism.
!
! author: Rahul A. Zaveri
! date  : february 1996
!-------------------------------------------------------------------------
      subroutine SetGas_Bio(ilast)
      use module_data_mosaic_main
      use module_data_mosaic_gas

      implicit none

      integer ilast

      iisop	= ilast + 1
      iisoprd	= ilast + 2
      iisopp	= ilast + 3
      iisopn	= ilast + 4
      iisopo2	= ilast + 5
      iapi	= ilast + 6
      ilim	= ilast + 7

      ilast	= ilim

      return
      end subroutine SetGas_Bio

