!
!***********************************************************************
      subroutine ode_bio
      use module_data_mosaic_main
      use module_data_mosaic_gas

      implicit none

      real(r8) :: yield_factor
!
!
      p_bio(ino)= 0.0
      d_bio(ino)= r_bio(8)+r_bio(9)+r_bio(10)
!
      p_bio(ino2)= .91*r_bio(8)+1.2*r_bio(9)+r_bio(10)
      d_bio(ino2)= 0.0
!
      p_bio(ino3)= 0.0
      d_bio(ino3)= r_bio(3)+r_bio(7)
!
      p_bio(ihno3)= .07*r_bio(7)
      d_bio(ihno3)= 0.0
!
      p_bio(io3)= 0.0
      d_bio(io3)= r_bio(2)+r_bio(6)+r_bio(18)+r_bio(21)
!
      p_bio(ioh)= .27*r_bio(2)+.27*r_bio(6)
      d_bio(ioh)= r_bio(1)+r_bio(5)
!
      p_bio(iho2)= .07*r_bio(2)+.33*r_bio(4)+.1*r_bio(6)   &
              +.93*r_bio(7)+.91*r_bio(8)+.8*r_bio(9)+r_bio(10)
      d_bio(iho2)= r_bio(11)+r_bio(12)+r_bio(13)
!
      p_bio(ih2o2)= 0.0
      d_bio(ih2o2)= 0.0
!
      p_bio(ico)= .07*r_bio(2)+.33*r_bio(4)+.16*r_bio(6)   &
             +.64*r_bio(7)+.59*r_bio(10)
      d_bio(ico)= 0.0
!
      p_bio(ihcho)= .6*r_bio(2)+.2*r_bio(4)+.15*r_bio(6)   &
               +.28*r_bio(7)+.63*r_bio(8)+.25*r_bio(10)
      d_bio(ihcho)= 0.0
!
      p_bio(iald2)= .15*r_bio(2)+.07*r_bio(4)+.02*r_bio(6)   &
               +.28*r_bio(7)+.8*r_bio(9)+.55*r_bio(10)+r_bio(15)   &
               +.5*r_bio(16)
      d_bio(iald2)= 0.0
!
      p_bio(ipar)= 1.86*r_bio(7)+.18*r_bio(8)+1.6*r_bio(9)+   &
                   2*r_bio(12)+2*r_bio(15)
      d_bio(ipar)= 0.0
!
      p_bio(iaone)= .03*r_bio(4)+.09*r_bio(6)+.63*r_bio(10)   &
               +.5*r_bio(16)
      d_bio(iaone)= 0.0
!
      p_bio(imgly)= .85*r_bio(6)+.34*r_bio(10)
      d_bio(imgly)= 0.0
!
      p_bio(ionit)= .93*r_bio(7)+.09*r_bio(8)+.8*r_bio(9)+r_bio(12)   &
               +r_bio(15)
      d_bio(ionit)= 0.0
!
      p_bio(ircooh)= .39*r_bio(2)+.46*r_bio(6)
      d_bio(ircooh)= 0.0
!
      p_bio(irooh)= r_bio(11)+r_bio(13)
      d_bio(irooh)= 0.0
!
      p_bio(ich3o2)= .7*r_bio(4)+.05*r_bio(6)
      d_bio(ich3o2)= 0.0
!
      p_bio(ic2o3)= .2*r_bio(2)+.97*r_bio(4)+.5*r_bio(5)   &
               +.11*r_bio(6)+.07*r_bio(7)
      d_bio(ic2o3)= 0.0
!
      p_bio(ixo2)= .08*r_bio(1)+.2*r_bio(2)+.2*r_bio(5)+.07*r_bio(6)   &
              +.93*r_bio(7)
      d_bio(ixo2)= 0.0
!
      p_bio(iisop)= 0.0
      d_bio(iisop)= r_bio(1)+r_bio(2)+r_bio(3)
!
      p_bio(iisoprd)= .65*r_bio(2)+.91*r_bio(8)+.2*r_bio(9)+r_bio(14)
      d_bio(iisoprd)= r_bio(4)+r_bio(5)+r_bio(6)+r_bio(7)
!
      p_bio(iisopp)= r_bio(1)
      d_bio(iisopp)= r_bio(8)+r_bio(11)+r_bio(14)
!
      p_bio(iisopn)= r_bio(3)
      d_bio(iisopn)= r_bio(9)+r_bio(12)+r_bio(15)
!
      p_bio(iisopo2)= .5*r_bio(5)
      d_bio(iisopo2)= r_bio(10)+r_bio(13)+r_bio(16)
!
!
! SORGAM
      yield_factor = 2.0 ! 1.0 = original SORGAM

      p_bio(iapi)= 0.0
      d_bio(iapi)= r_bio(17)+r_bio(18)+r_bio(19)
!
      p_bio(ilim)= 0.0
      d_bio(ilim)= r_bio(20)+r_bio(21)+r_bio(22)
!
      p_bio(iapi1)= yield_factor*(foh *0.028*r_bio(17) +   &
                    fo3 *0.028*r_bio(18)     +   &
                    fno3*0.028*r_bio(19))

      d_bio(iapi1)= 0.0
!
      p_bio(iapi2)= yield_factor*(foh *0.241*r_bio(17) +   &
                    fo3 *0.241*r_bio(18)     +   &
                    fno3*0.241*r_bio(19))
      d_bio(iapi2)= 0.0
!
      p_bio(ilim1)= yield_factor*(foh *0.163*r_bio(20) +   &
                    fo3 *0.163*r_bio(21)     +   &
                    fno3*0.163*r_bio(22))
      d_bio(ilim1)= 0.0
!
      p_bio(ilim2)= yield_factor*(foh *0.247*r_bio(20) +   &
                    fo3 *0.247*r_bio(21)     +   &
                    fno3*0.247*r_bio(22))
      d_bio(ilim2)= 0.0

      return
      end subroutine ode_bio






