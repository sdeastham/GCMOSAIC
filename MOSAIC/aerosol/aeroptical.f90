! aeroptical.f90
!-----------------------------------------------------------------------
subroutine aerosol_optical(                 &
      cair_mol_m3, dens_comp_a, mw_comp_a,  & ! all in
      dens_aer_mac, mw_aer_mac, ref_index_a,& ! all in
      jaerosolstate, jphase, tot_cl_in,     &
      tot_nh4_in, tot_no3_in, XT, gas, aer, &
      water_a, num_a, Dp_dry_a, Dp_wet_a,   &
      dp_core_a, area_dry_a, area_wet_a,    &
      dens_dry_a, dens_wet_a, mass_dry_a,   &
      mass_wet_a, vol_dry_a, vol_wet_a,     &
      total_species, electrolyte,           &
      ri_shell_a, ri_avg_a, ri_core_a,      &
      sigmag_a, jhyst_leg                   )
      
!
! rce 2013-07-30
! *** this is note fully implemented yet ***
! *** the call from time_integration_mode is incorrect ***
!
! ??? question ???
! what is the purpose of the "aer" argument
!     will the aerosol information be passed in through aer ??
!     or will aerosol info be passed in through cnn (as was done previously)
!         and need to be mapped into aer ??
!

  use module_data_mosaic_kind, only: r8
  use module_data_mosaic_main, only: cnn, lun_aeroptic, time_hrs, time_utc
  use module_data_mosaic_aero, only: naer, naercomp, nbin_a, nbin_a_max, ngas_volatile, nelectrolyte
  use module_mosaic_box, only: aerosol_optical_properties

  implicit none
  !   subr arguments
  integer, intent(inout), dimension(nbin_a_max) :: jaerosolstate, jphase 
  integer, intent(inout), dimension(nbin_a_max) :: jhyst_leg

  real(r8), intent(inout) :: cair_mol_m3
  real(r8), intent(inout), dimension(nbin_a_max) :: sigmag_a

  
  real(r8), intent(in), dimension(naer)       :: dens_aer_mac, mw_aer_mac
  real(r8), intent(in), dimension(naercomp)   :: dens_comp_a,mw_comp_a
  
  real(r8), intent(inout) :: tot_cl_in, tot_nh4_in, tot_no3_in, XT
  real(r8), intent(inout), dimension(ngas_volatile) :: gas
  real(r8), intent(inout), dimension(naer,3,nbin_a_max) :: aer
  real(r8), intent(inout), dimension(nbin_a_max) :: water_a
  real(r8), intent(inout), dimension(nbin_a_max) :: num_a
  real(r8), intent(inout), dimension(nbin_a_max) :: Dp_dry_a
  real(r8), intent(inout), dimension(nbin_a_max) :: Dp_wet_a, dp_core_a
  real(r8), intent(inout), dimension(nbin_a_max) :: area_dry_a, area_wet_a
  real(r8), intent(inout), dimension(nbin_a_max) :: dens_dry_a, dens_wet_a
  real(r8), intent(inout), dimension(nbin_a_max) :: mass_dry_a, mass_wet_a
  real(r8), intent(inout), dimension(nbin_a_max) :: vol_dry_a, vol_wet_a
  real(r8), intent(inout), dimension(ngas_volatile) :: total_species    
  real(r8), intent(inout), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte
  
  
  complex,  intent(in), dimension(naercomp)   :: ref_index_a
  complex,  intent(inout), dimension(nbin_a_max) :: ri_shell_a, ri_avg_a, ri_core_a
  
  !   local variables
  integer ibin, it!BALLI its wring to mention 'it' here..it should be subr arg
  
  real(r8), save ::      &
       vlambc,   &  ! wavelength in microns
       qextc,            &  ! extinction efficiency, dimensionless
       qscatc,           &  ! scattering efficiency, dimensionless
       gscac,            &  ! asymmetry parameter for given size, dimensionless
       extc,             &  ! extinction cross section, cm^2
       scatc                ! scattering cross section, cm^2

  real(r8), dimension(nbin_a_max)     :: ext_cross, scat_cross, asym_particle

  !
  !   NOTES
  !
  !   mshellcore = 0 - do non-shell/core mie calc
  !              = 1 - do shell/core mie calc, and shell = { bc }
  !             ??? what about insoluble oin and dust species ???
  !
  !   core volume == 0 if either
  !      ri_core_a(ibin) == (0.0, 0.0)
  !      dp_core_a(ibin) == 0.0
  !   shell volume == 0 if either
  !      ri_shell_a(ibin) == (0.0, 0.0)
  !

  write( lun_aeroptic, '(a,2f12.3)' )   &
       'in aerosol_optical - run hrs, utc hrs = ', time_hrs, time_utc

  if (it == 0) then
     write( lun_aeroptic, '(a)' )   &
          '   *** currently cannot do optical calcs when it=0'
     write( lun_aeroptic, '(a)' )
     return
  endif

  ! map variables from cnn arrays to mosaic aerchem working arrays
  call map_mosaic_species_BOX( 0, cnn, aer, gas, jhyst_leg, num_a, Dp_dry_a,      &
            sigmag_a, water_a, cair_mol_m3 )

  !
  ! calc optical properties of each bin
  !
  write( lun_aeroptic, '(a)' )   &
       '   calling aerosol_optical_properties'
  call aerosol_optical_properties(                                       &
          gas, aer, num_a, water_a,                                      &  ! inout, inout, inout, in
          dens_comp_a, mw_comp_a, dens_aer_mac, mw_aer_mac, ref_index_a, &  ! all in
          Dp_dry_a, Dp_wet_a, dp_core_a,                                 &  ! inout, out, out
          ri_shell_a, ri_core_a, ri_avg_a, jaerosolstate, jphase,        &
          tot_cl_in, tot_nh4_in, tot_no3_in, XT, area_dry_a, area_wet_a, &
          dens_dry_a, dens_wet_a, mass_dry_a, mass_wet_a, vol_dry_a,     &
          vol_wet_a, total_species, electrolyte                          )  ! all out

  !
  ! do mie calculations
  !
  write( lun_aeroptic, '(a)' )   &
       '   doing mie calculations'
  do ibin = 1, nbin_a

     ! set vlambc to convention "550 nm" wavelength
     vlambc= 0.550  ! wavelength in microns

     call miedriver( dp_wet_a(ibin), dp_core_a(ibin), &
          ri_shell_a(ibin), ri_core_a(ibin), &
          vlambc, qextc, qscatc, gscac, extc, scatc )

     ! save cross sections & asymmetry parameter for each particle
     ext_cross(ibin) = extc
     scat_cross(ibin) = scatc
     asym_particle(ibin) = gscac

  end do

  !
  ! output
  !
  write( lun_aeroptic, '(3a)' )   &
       '   bin, num (#/cm3), dpwet, dpcore (cm)',   &
       ', ri_tot, ri_shell, ri_core',   &
       ', ext_cross, scat_cross, asymmetry_particle'
  do ibin = 1, nbin_a
     write( lun_aeroptic, '(i7,1p,e10.2,4(2x,2e10.2)2x,3e10.2)' )   &
          ibin, num_a(ibin), dp_wet_a(ibin), dp_core_a(ibin),   &
          ri_avg_a(ibin), ri_shell_a(ibin), ri_core_a(ibin),   &
          ext_cross(ibin), scat_cross(ibin), asym_particle(ibin)
  end do
  write( lun_aeroptic, '(a)' )

  return
end subroutine aerosol_optical

