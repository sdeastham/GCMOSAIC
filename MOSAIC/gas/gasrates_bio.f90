      subroutine GasRates_Bio(s)
      use module_data_mosaic_main
      use module_data_mosaic_gas

      implicit none

      real(r8) :: s(ntot_max)

! isoprene
      r_bio(1)  = rk_bio(1) *s(iisop)*s(ioh)
      r_bio(2)  = rk_bio(2) *s(iisop)*s(io3)
      r_bio(3)  = rk_bio(3) *s(iisop)*s(ino3)
      r_bio(4)  = rk_bio(4) *s(iisoprd)
      r_bio(5)  = rk_bio(5) *s(iisoprd)*s(ioh)
      r_bio(6)  = rk_bio(6) *s(iisoprd)*s(io3)
      r_bio(7)  = rk_bio(7) *s(iisoprd)*s(ino3)
      r_bio(8)  = rk_bio(8) *s(iisopp)*s(ino)
      r_bio(9)  = rk_bio(9) *s(iisopn)*s(ino)
      r_bio(10) = rk_bio(10)*s(iisopo2)*s(ino)
      r_bio(11) = rk_bio(11)*s(iisopp)*s(iho2)
      r_bio(12) = rk_bio(12)*s(iisopn)*s(iho2)
      r_bio(13) = rk_bio(13)*s(iisopo2)*s(iho2)
      r_bio(14) = rk_bio(14)*s(iisopp)
      r_bio(15) = rk_bio(15)*s(iisopn)
      r_bio(16) = rk_bio(16)*s(iisopo2)

! a-pinene
      r_bio(17) = rk_bio(17)*s(iapi)*s(ioh)
      r_bio(18) = rk_bio(18)*s(iapi)*s(io3)
      r_bio(19) = rk_bio(19)*s(iapi)*s(ino3)
! limonene
      r_bio(20) = rk_bio(20)*s(ilim)*s(ioh)
      r_bio(21) = rk_bio(21)*s(ilim)*s(io3)
      r_bio(22) = rk_bio(22)*s(ilim)*s(ino3)

      return
      end subroutine GasRates_Bio



