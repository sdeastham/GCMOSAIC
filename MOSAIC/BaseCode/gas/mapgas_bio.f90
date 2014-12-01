
!***********************************************************************
! subroutine MapGas_Bio: maps cnn to and fro stot for the biogenic
!                        gas-phase mechanism.
!
! nomenclature:
! cnn       = full species concentration array.
! stot      = subset of cnn. species concentration array to be supplied to
!             lsodes. length of stot depends on the selected mechanism
! iregime   = selected chemical regime (1-6)
! imap      = 0 : map cnn to stot
!           = 1 : map stot to cnn
!
! author: Rahul A. Zaveri
! date  : february 1996
!
!-------------------------------------------------------------------------
      subroutine MapGas_Bio(stot,imap)
      use module_data_mosaic_main
      use module_data_mosaic_gas

      implicit none

      integer imap
      real(r8) :: stot(ntot_max)

      emit(iisop)	=emission(kisop)
      emit(iisoprd)	=emission(kisoprd)
      emit(iisopp)	=emission(kisopp)
      emit(iisopn)	=emission(kisopn)
      emit(iisopo2)	=emission(kisopo2)
      emit(iapi)	=emission(kapi)
      emit(ilim)	=emission(klim)

      if(imap.eq.0)then    ! map cnn into stot
      stot(iisop)	=cnn(kisop)
      stot(iisoprd)	=cnn(kisoprd)
      stot(iisopp)	=cnn(kisopp)
      stot(iisopn)	=cnn(kisopn)
      stot(iisopo2)	=cnn(kisopo2)
      stot(iapi)	=cnn(kapi)
      stot(ilim)	=cnn(klim)
!
      else                 ! map stot back into cnn
      cnn(kisop)	=stot(iisop)
      cnn(kisoprd)	=stot(iisoprd)
      cnn(kisopp)	=stot(iisopp)
      cnn(kisopn)	=stot(iisopn)
      cnn(kisopo2)	=stot(iisopo2)
      cnn(kapi)		=stot(iapi)
      cnn(klim)		=stot(ilim)
      endif

      return
      end subroutine MapGas_Bio
