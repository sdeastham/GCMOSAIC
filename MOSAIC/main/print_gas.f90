! Gas-phase species
! subroutine for printing output at iprint time steps

      subroutine print_gas
      use module_data_mosaic_main
      use module_data_mosaic_gas

      implicit none

      integer l, it !BALLI, 'it' shld be subr arg
      real(r8) :: cppb(ntot_max)


      if (iwrite_gas <= 0) return


      do l=1,ngas_max
        cppb(l) = (cnn(l)/cair_mlc)*ppb	! converting (molecules/cc) to (ppb)
      enddo


      if(it.eq.0)write(lun_gas,201)(species(l), l=1, ngas_max)
201   format(' UTC(hr)    t(hr)   Temp(K)    Pr(atm)   RH(%)',   &
             '  AIR(molec/cc) H2O(molec/cc)', 100(2x,a11))



      write(lun_gas,202)time_UTC,time_hrs,te,pr_atm,RH,cair_mlc,h2o,   &
                   (cppb(l),l=1,ngas_max)
202   format(f8.4,4(2x,f8.4),100(2x,e11.5))

!      close(20)


      return
      end subroutine print_gas

