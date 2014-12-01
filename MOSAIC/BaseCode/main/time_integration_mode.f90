subroutine time_integration_mode( gam_ratio, iter_MESA, aH2O_a)
  use module_data_mosaic_kind,  only: r8
  use module_data_mosaic_main,  only: mgas, dt_min, dt_aeroptic_min,&
       maer, maeroptic, nstep, msolar, ipmcmos, iprint, pblh, rh,   &
       pr_atm, te, time_sec, lun_sect_190, cair_mlc, ppb, ngas_max, &
       cnn
  use module_data_mosaic_aero,  only: nbin_a_max, nsalt, mYES

  use module_pmcmos_subaa,      only:  pmcmos_update_met, pmcmos_emit_dilu
  use module_pmcmos_aer,        only:  pmcmos_print

  implicit none

  ! Subroutine Arguments
  integer,  intent(out), dimension(nbin_a_max) :: iter_MESA
  real(r8), intent(out), dimension(nbin_a_max) :: aH2O_a, gam_ratio

  !Local variables
  integer :: mfreq_aeroptic, it
  integer, dimension(nbin_a_max) :: jaerosolstate

  real(r8) :: cair_mol_m3,cair_mol_cc  
  real(r8), dimension(nbin_a_max) :: dp_wet_a


  jaerosolstate(:) = 0      !BSINGH - Temporarily initialized, please modify if required

  it = 0                        ! initialize time step counter
  if(mgas == mYES)then
     call print_gas             ! print inital gas concentrations
  endif

  !------------------------------------------------------------------
  !
  ! main "time" loop begins...

  mfreq_aeroptic = max( 1, nint(dt_aeroptic_min/dt_min) )
  if ((maeroptic > 0) .and. (maer == myes)) then
     !BSINGH - The following is not implemented correctly in the code
     !         In this scenario, the code should print out an error message and stop!
     print*, 'Stopping in time_integration_mode.f90 before aerosol_optical call...'
     stop     
     !call aerosol_optical(                    &
     !     cair_mol_m3, dens_comp_a, mw_comp_a,  & ! all in
     !     dens_aer_mac, mw_aer_mac, ref_index_a ) ! all in
  endif

  call pmcmos_print(0,jaerosolstate)

  do  it = 1, nstep
     
     call UpdateTime
     
     if(msolar.eq.1)then
        call SolarZenithAngle
     endif
     
     if (ipmcmos <=  0) then
        call UpdateMetFields
        call UpdateEmissions
     else
        call pmcmos_update_met( 1 )
        call pmcmos_emit_dilu
        
        !        write(*,'(1p,a,3e14.6)') 'cos_sza,tcur,tmid', cos_sza, tcur_sec, tmid_sec
        !        write(*,'(a,1p,5e13.5)') 'co,so2,ch4,c2h6', cnn(17:20)*ppb/cair_mlc
     endif
     
     call IntegrateChemistry(it,                                & !intent-ins
          jaerosolstate, dp_wet_a,                              & !intent-inouts
          cair_mol_m3, cair_mol_cc, aH2O_a, gam_ratio,iter_mesa ) !intent-outs
     
     if ((ipmcmos >  0) .and. &
          (mod(it,iprint) == 0)) then
        write(lun_sect_190,'(/a,5(1x,f14.6))') 'time(h)', time_sec/3600.0, &
             te, pr_atm*1013.25, rh, pblh
        write(lun_sect_190,'(1p,5e14.6)') cnn(1:ngas_max)*ppb/cair_mlc
     endif
     
     !      call DoMassBalance
     
     if ((maeroptic > 0) .and. (maer == myes)) then
        !BSINGH - The following is not implemented correctly in the code
        !         In this scenario, the code should print out an error message and stop!
        print*, 'Stopping in time_integration_mode.f90 before aerosol_optical call (2nd call)...'
        stop  
        
        !if (mod(it,mfreq_aeroptic) .eq. 0) &
        !     call aerosol_optical(                    &
        !     cair_mol_m3, dens_comp_a, mw_comp_a,  & ! all in 
        !     dens_aer_mac, mw_aer_mac, ref_index_a ) ! all in
     endif
     
     
     if(mgas == mYES)then
        call print_gas
     endif
     
     !!      if(maer == mYES)then               ! UNCOMMENT THIS LINE
     !!        call print_aer(1)                ! UNCOMMENT THIS LINE
     !!      endif                              ! UNCOMMENT THIS LINE
     
     call pmcmos_print(it,jaerosolstate)
     
  enddo   ! time loop
  
  !------------------------------------------------------------------
  
  
  return
end subroutine time_integration_mode




