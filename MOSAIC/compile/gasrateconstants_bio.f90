!**************************************************************************
! subroutine GasRateConstants_Bio: generates thermal rate coefficients
!                   for the selected mechanism
! nomenclature:
! rk_bio    = reaction rate constants for hc2 mechanism    (molec-cc-s)
! te        = ambient atmospheric temperature (K)
!
! author: Rahul A. Zaveri
! date  : february 1996
!
!-------------------------------------------------------------------------
      subroutine GasRateConstants_Bio
      use module_data_mosaic_main
      use module_data_mosaic_gas

      implicit none

      real(r8) :: ARR

      rk_bio(1)  = ARR(2.6d-11, 409.d0)
      rk_bio(2)  = ARR(1.2d-14, -2013.d0)
      rk_bio(3)  = ARR(3.0d-12, -446.d0)
      rk_bio(4)  = rk_photo(jphoto_isoprd)
      rk_bio(5)  = 3.3e-11
      rk_bio(6)  = 7.0e-18
      rk_bio(7)  = 1.0e-15
      rk_bio(8)  = 4.0e-12
      rk_bio(9)  = 4.0e-12
      rk_bio(10) = 4.0e-12
      rk_bio(11) = ARR(1.7d-13, 1300.d0)
      rk_bio(12) = ARR(1.7d-13, 1300.d0)
      rk_bio(13) = ARR(1.7d-13, 1300.d0)
      rk_bio(14) = rk_param(jisopp)
      rk_bio(15) = rk_param(jisopn)
      rk_bio(16) = rk_param(jisopo2)

! a-pinene
      rk_bio(17) = ARR(12.1d-12, 444.d0)
      rk_bio(18) = ARR(1.01d-15, -732.d0)
      rk_bio(19) = ARR(1.19d-12, 490.d0)

! limonene
      rk_bio(20) = 1.71e-10
      rk_bio(21) = 2.00e-16
      rk_bio(22) = 1.22e-11

      return
      end subroutine GasRateConstants_Bio



