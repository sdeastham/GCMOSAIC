
!***********************************************************************
! subroutine MapGas_Com: maps cnn to and fro stot for the common
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
      subroutine MapGas_Com(stot,imap)
      use module_data_mosaic_main
      use module_data_mosaic_gas

      implicit none

      integer imap
      real(r8) :: stot(ntot_max)

      emit(ih2so4)	=emission(kh2so4)
      emit(ihno3)	=emission(khno3)
      emit(ihcl)	=emission(khcl)
      emit(inh3)	=emission(knh3)
      emit(ino)		=emission(kno)
      emit(ino2)	=emission(kno2)
      emit(ino3)	=emission(kno3)
      emit(in2o5)	=emission(kn2o5)
      emit(ihono)	=emission(khono)
      emit(ihno4)	=emission(khno4)
      emit(io3)		=emission(ko3)
      emit(io1d)	=emission(ko1d)
      emit(io3p)	=emission(ko3p)
      emit(ioh)		=emission(koh)
      emit(iho2)	=emission(kho2)
      emit(ih2o2)	=emission(kh2o2)
      emit(ico)		=emission(kco)
      emit(iso2)	=emission(kso2)
      emit(ich4)	=emission(kch4)
      emit(ic2h6)	=emission(kc2h6)
      emit(ich3o2)	=emission(kch3o2)
      emit(iethp)	=emission(kethp)
      emit(ihcho)	=emission(khcho)
      emit(ich3oh)	=emission(kch3oh)
      emit(ianol)	=emission(kanol)
      emit(ich3ooh)	=emission(kch3ooh)
      emit(iethooh)	=emission(kethooh)
      emit(iald2)	=emission(kald2)
      emit(ihcooh)	=emission(khcooh)
      emit(ircooh)	=emission(krcooh)
      emit(ic2o3)	=emission(kc2o3)
      emit(ipan)	=emission(kpan)

      emit(iaro1)	=emission(karo1)
      emit(iaro2)	=emission(karo2)
      emit(ialk1)	=emission(kalk1)
      emit(iole1)	=emission(kole1)
      emit(iapi1)	=emission(kapi1)
      emit(iapi2)	=emission(kapi2)
      emit(ilim1)	=emission(klim1)
      emit(ilim2)	=emission(klim2)


      if(imap.eq.0)then    ! map cnn into stot
      stot(ih2so4)	=cnn(kh2so4)
      stot(ihno3)	=cnn(khno3)
      stot(ihcl)	=cnn(khcl)
      stot(inh3)	=cnn(knh3)
      stot(ino)		=cnn(kno)
      stot(ino2)	=cnn(kno2)
      stot(ino3)	=cnn(kno3)
      stot(in2o5)	=cnn(kn2o5)
      stot(ihono)	=cnn(khono)
      stot(ihno4)	=cnn(khno4)
      stot(io3)		=cnn(ko3)
      stot(io1d)	=cnn(ko1d)
      stot(io3p)	=cnn(ko3p)
      stot(ioh)		=cnn(koh)
      stot(iho2)	=cnn(kho2)
      stot(ih2o2)	=cnn(kh2o2)
      stot(ico)		=cnn(kco)
      stot(iso2)	=cnn(kso2)
      stot(ich4)	=cnn(kch4)
      stot(ic2h6)	=cnn(kc2h6)
      stot(ich3o2)	=cnn(kch3o2)
      stot(iethp)	=cnn(kethp)
      stot(ihcho)	=cnn(khcho)
      stot(ich3oh)	=cnn(kch3oh)
      stot(ianol)	=cnn(kanol)
      stot(ich3ooh)	=cnn(kch3ooh)
      stot(iethooh)	=cnn(kethooh)
      stot(iald2)	=cnn(kald2)
      stot(ihcooh)	=cnn(khcooh)
      stot(ircooh)	=cnn(krcooh)
      stot(ic2o3)	=cnn(kc2o3)
      stot(ipan)	=cnn(kpan)
      stot(iaro1)	=cnn(karo1)
      stot(iaro2)	=cnn(karo2)
      stot(ialk1)	=cnn(kalk1)
      stot(iole1)	=cnn(kole1)
      stot(iapi1)	=cnn(kapi1)
      stot(iapi2)	=cnn(kapi2)
      stot(ilim1)	=cnn(klim1)
      stot(ilim2)	=cnn(klim2)
!
      else                 ! map stot back into cnn
      cnn(ih2so4)	=stot(kh2so4)
      cnn(ihno3)	=stot(khno3)
      cnn(ihcl)		=stot(khcl)
      cnn(inh3)		=stot(knh3)
      cnn(ino)		=stot(kno)
      cnn(ino2)		=stot(kno2)
      cnn(ino3)		=stot(kno3)
      cnn(in2o5)	=stot(kn2o5)
      cnn(ihono)	=stot(khono)
      cnn(ihno4)	=stot(khno4)
      cnn(io3)		=stot(ko3)
      cnn(io1d)		=stot(ko1d)
      cnn(io3p)		=stot(ko3p)
      cnn(ioh)		=stot(koh)
      cnn(iho2)		=stot(kho2)
      cnn(ih2o2)	=stot(kh2o2)
      cnn(ico)		=stot(kco)
      cnn(iso2)		=stot(kso2)
      cnn(ich4)		=stot(kch4)
      cnn(ic2h6)	=stot(kc2h6)
      cnn(ich3o2)	=stot(kch3o2)
      cnn(iethp)	=stot(kethp)
      cnn(ihcho)	=stot(khcho)
      cnn(ich3oh)	=stot(kch3oh)
      cnn(ianol)	=stot(kanol)
      cnn(ich3ooh)	=stot(kch3ooh)
      cnn(iethooh)	=stot(kethooh)
      cnn(iald2)	=stot(kald2)
      cnn(ihcooh)	=stot(khcooh)
      cnn(ircooh)	=stot(krcooh)
      cnn(ic2o3)	=stot(kc2o3)
      cnn(ipan)		=stot(kpan)
      cnn(karo1)	=stot(iaro1)
      cnn(karo2)	=stot(iaro2)
      cnn(kalk1)	=stot(ialk1)
      cnn(kole1)	=stot(iole1)
      cnn(kapi1)	=stot(iapi1)
      cnn(kapi2)	=stot(iapi2)
      cnn(klim1)	=stot(ilim1)
      cnn(klim2)	=stot(ilim2)

      endif

      return
      end subroutine MapGas_Com
