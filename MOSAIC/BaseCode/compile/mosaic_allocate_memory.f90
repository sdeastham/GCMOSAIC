
      subroutine mosaic_allocate_memory( initialize_flag )
!
! allocate arrays whose dimensions are set at run time
!
! initialize them to zero when initialize_flag > 0 (or is not present)
!

      use module_data_mosaic_main
      use module_data_mosaic_aero
      use module_data_mosaic_asect

      implicit none

      integer, intent(in), optional :: initialize_flag


      allocate( lun_aer(naerbin) )
      allocate( lun_aer_status(naerbin) )

      allocate( aer_output(naerbin) )

      allocate( species(ntot_max) )

      allocate( cnn(ntot_max) )

      allocate( emission(ntot_max), emit(ntot_max) )

      !allocate( jaerosolstate(nbin_a_max) )
      !allocate( jaerosolstate_bgn(nbin_a_max) )
      !allocate( jphase(nbin_a_max) )
      !allocate( jhyst_leg(nbin_a_max) )

      !allocate( num_a(nbin_a_max) )
      !allocate( Dpgn_a(nbin_a_max) )
      !allocate( Dp_dry_a(nbin_a_max) )
      !allocate( Dp_wet_a(nbin_a_max) )
      !allocate( Dp_core_a(nbin_a_max) )
      !allocate( area_dry_a(nbin_a_max) )
      !allocate( area_wet_a(nbin_a_max) )
      !allocate( mass_dry_salt(nbin_a_max) )
      !allocate( mass_dry_a_bgn(nbin_a_max) )
      !allocate( mass_dry_a(nbin_a_max) )
      !allocate( mass_wet_a(nbin_a_max) )
      !allocate( mass_soluble_a(nbin_a_max) )
      !allocate( vol_dry_a(nbin_a_max) )
      !allocate( vol_wet_a(nbin_a_max) )
      !allocate( dens_dry_a_bgn(nbin_a_max) )
      !allocate( dens_dry_a(nbin_a_max) )
      !allocate( dens_wet_a(nbin_a_max) )
      !allocate( sigmag_a(nbin_a_max) )
      !allocate( water_a(nbin_a_max) )
      !allocate( water_a_hyst(nbin_a_max) )
      !allocate( water_a_up(nbin_a_max) )
      !allocate( pH(nbin_a_max) )
      !allocate( aer(naer,3,nbin_a_max) )
      !allocate( aer_sum(3,nbin_a_max) )
      !allocate( aer_percent(naer,3,nbin_a_max) )
      !allocate( electrolyte(nelectrolyte,3,nbin_a_max) )
      !allocate( electrolyte_sum(3,nbin_a_max) )
      !allocate( epercent(nelectrolyte,3,nbin_a_max) )
      !allocate( aH2O_a(nbin_a_max) )
      !allocate( DpmV(nbin_a_max) )
      !allocate( volume_a(nbin_a_max) )
      !allocate( volume_bin(nbin_a_max) )
      !allocate( kelvin(nbin_a_max) )
      !allocate( kel(ngas_volatile,nbin_a_max) )
      !allocate( ext_cross(nbin_a_max) )
      !allocate( scat_cross(nbin_a_max) )
      !allocate( asym_particle(nbin_a_max) )

      !allocate( idry_case3a(nbin_a_max) )
      !allocate( ieqblm_bin(nbin_a_max) )
      !allocate( integrate(ngas_volatile,3,nbin_a_max) )

      !allocate( Heff(ngas_volatile,nbin_a_max) )
      !allocate( kg(ngas_volatile,nbin_a_max) )
      !allocate( df_gas_s(ngas_volatile,nbin_a_max) )
      !allocate( df_gas_l(ngas_volatile,nbin_a_max) )
      !allocate( df_gas_o(ngas_volatile,nbin_a_max) )
      !allocate( df_gas(ngas_volatile,nbin_a_max) )
      !allocate( flux_s(ngas_volatile,nbin_a_max) )
      !allocate( flux_l(ngas_volatile,nbin_a_max) )
      !allocate( flux_o(ngas_volatile,nbin_a_max) )
      !allocate( flux(ngas_volatile,nbin_a_max) )
      !allocate( delta_nh3_max(nbin_a_max) )
      !allocate( delta_hno3_max(nbin_a_max) )
      !allocate( delta_hcl_max(nbin_a_max) )
      !allocate( volatile_s(ngas_volatile,nbin_a_max) )
      !allocate( phi_volatile_s(ngas_volatile,nbin_a_max) )
      !allocate( phi_volatile_l(ngas_volatile,nbin_a_max) )
      !allocate( phi_volatile_o(ngas_volatile,nbin_a_max) )
      !allocate( h_s_i_m(ngas_volatile,nbin_a_max) )

      !allocate( iter_MESA(nbin_a_max) )

      !allocate( growth_factor(nbin_a_max) )
      !allocate( MDRH(nbin_a_max) )
      !allocate( sigma_soln(nbin_a_max) )

      !allocate( ri_avg_a(nbin_a_max) )
      !allocate( ri_shell_a(nbin_a_max) )
      !allocate( ri_core_a(nbin_a_max) )

      !allocate( mc(Ncation,nbin_a_max) )
      !allocate( ma(Nanion,nbin_a_max) )
      !allocate( gam(nelectrolyte,nbin_a_max) )
      !allocate( gam_ratio(nbin_a_max) )
      !allocate( activity(nelectrolyte,nbin_a_max) )


      if (m_partmc_mosaic <= 0) then
! following only used for mosaic box-model (not for partmc-mosaic)

         allocate( nsize_aer( maxd_atype ) )
         allocate( ncomp_aer( maxd_atype ) )
         allocate( ncomp_plustracer_aer( maxd_atype ) )
         allocate( mastercompptr_aer(maxd_acomp, maxd_atype) )
         allocate( massptr_aer( maxd_acomp, maxd_asize, maxd_atype, maxd_aphase ) )
         
         allocate( waterptr_aer( maxd_asize, maxd_atype ) )
         allocate( hyswptr_aer( maxd_asize, maxd_atype ) )
         allocate( numptr_aer( maxd_asize, maxd_atype, maxd_aphase ) )
         
         allocate( mprognum_aer(maxd_asize,maxd_atype,maxd_aphase) )
         
         allocate( ibin_of_isize_itype( maxd_asize, maxd_atype ) )
         allocate( isize_of_ibin( nbin_a_max ) )
         allocate( itype_of_ibin( nbin_a_max ) )
         allocate( itype_of_itype_md1md2(maxd_atype_md1,maxd_atype_md2) )
         allocate( itype_md1_of_itype(maxd_atype) )
         allocate( itype_md2_of_itype(maxd_atype) )
         
         allocate( dens_aer( maxd_acomp, maxd_atype ) )
         allocate( mw_aer( maxd_acomp, maxd_atype ) )
         allocate( hygro_aer( maxd_acomp, maxd_atype ) )
         
         allocate( volumcut_sect( 0:maxd_asize, maxd_atype ) )
         allocate( volumcen_sect(   maxd_asize, maxd_atype ) )
         allocate( volumlo_sect(    maxd_asize, maxd_atype ) )
         allocate( volumhi_sect(    maxd_asize, maxd_atype ) )
         allocate( dcut_sect( 0:maxd_asize, maxd_atype ) )
         allocate( dcen_sect(   maxd_asize, maxd_atype ) )
         allocate( dlo_sect(    maxd_asize, maxd_atype ) )
         allocate( dhi_sect(    maxd_asize, maxd_atype ) )
         allocate( sigmag_aer(maxd_asize, maxd_atype) )
         
         allocate( xcut_atype_md1( 0:maxd_atype_md1 ) )
         allocate( xcut_atype_md2( 0:maxd_atype_md2 ) )
         
         allocate( name_aer( maxd_acomp, maxd_atype ) )
         
         allocate( lptr_so4_aer(maxd_asize, maxd_atype, maxd_aphase) )
         allocate( lptr_msa_aer(maxd_asize, maxd_atype, maxd_aphase) )
         allocate( lptr_no3_aer(maxd_asize, maxd_atype, maxd_aphase) )
         allocate( lptr_cl_aer( maxd_asize, maxd_atype, maxd_aphase) )
         allocate( lptr_co3_aer(maxd_asize, maxd_atype, maxd_aphase) )
         allocate( lptr_nh4_aer(maxd_asize, maxd_atype, maxd_aphase) )
         allocate( lptr_na_aer( maxd_asize, maxd_atype, maxd_aphase) )
         allocate( lptr_ca_aer( maxd_asize, maxd_atype, maxd_aphase) )
         allocate( lptr_oin_aer(maxd_asize, maxd_atype, maxd_aphase) )
         allocate( lptr_oc_aer( maxd_asize, maxd_atype, maxd_aphase) )
         allocate( lptr_bc_aer( maxd_asize, maxd_atype, maxd_aphase) )
         
      end if
      
      
      if ( present( initialize_flag ) ) then
         if (initialize_flag <= 0) return
      end if
      ! initialize the allocated arrays
      
      lun_aer = 0
      lun_aer_status = 0
      
      aer_output = '???'

      species = '???'

      cnn = 0.0

      emission = 0.0

      !jaerosolstate = 0
      !jaerosolstate_bgn = 0
      !jphase = 0
      !jhyst_leg = 0

      !num_a = 0.0
      !Dpgn_a = 0.0
      !Dp_dry_a = 0.0
      !Dp_wet_a = 0.0
      !Dp_core_a = 0.0
      !area_dry_a = 0.0
      !area_wet_a = 0.0
      !mass_dry_salt = 0.0
      !mass_dry_a_bgn = 0.0
      !mass_dry_a = 0.0
      !mass_wet_a = 0.0
      !mass_soluble_a = 0.0
      !vol_dry_a = 0.0
      !vol_wet_a = 0.0
      !dens_dry_a_bgn = 0.0
      !dens_dry_a = 0.0
      !dens_wet_a = 0.0
      !sigmag_a = 0.0
      !water_a = 0.0
      !water_a_hyst = 0.0
      !water_a_up = 0.0
      !pH = 0.0
      !aer = 0.0
      !aer_sum = 0.0
      !aer_percent = 0.0
      !electrolyte = 0.0
      !electrolyte_sum = 0.0
      !epercent = 0.0
      !aH2O_a = 0.0
      !DpmV = 0.0
      !volume_a = 0.0
      !volume_bin = 0.0
      !kelvin = 0.0
      !kel = 0.0
      !ext_cross = 0.0
      !scat_cross = 0.0
      !asym_particle = 0.0

      !idry_case3a = 0
      !ieqblm_bin = 0
      !integrate = 0

      !Heff = 0.0
      !kg = 0.0
      !df_gas_s = 0.0
      !df_gas_l = 0.0
      !df_gas_o = 0.0
      !df_gas = 0.0
      !flux_s = 0.0
      !flux_l = 0.0
      !flux_o = 0.0
      !flux = 0.0
      !delta_nh3_max = 0.0
      !delta_hno3_max = 0.0
      !delta_hcl_max = 0.0
      !volatile_s = 0.0
      !phi_volatile_s = 0.0
      !phi_volatile_l = 0.0
      !phi_volatile_o = 0.0
      !h_s_i_m = 0.0

      !iter_MESA = 0

      !growth_factor = 0.0
      !MDRH = 0.0
      !sigma_soln = 0.0

      !ri_avg_a = ( 0.0, 0.0 )
      !ri_shell_a = ( 0.0, 0.0 )
      !ri_core_a = ( 0.0, 0.0 )

      !mc = 0.0
      !ma = 0.0
      !gam = 0.0
      !gam_ratio = 0.0
      !activity = 0.0


      if (m_partmc_mosaic <= 0) then
! following only used for mosaic box-model (not for partmc-mosaic)

      nsize_aer = 0
      ncomp_aer = 0
      ncomp_plustracer_aer = 0
      mastercompptr_aer = 0
      massptr_aer = 0

      waterptr_aer = 0
      hyswptr_aer = 0
      numptr_aer = 0

      mprognum_aer = 0

      ibin_of_isize_itype = 0
      isize_of_ibin = 0
      itype_of_ibin = 0
      itype_of_itype_md1md2 = 0
      itype_md1_of_itype = 0
      itype_md2_of_itype = 0

      dens_aer = 0.0
      mw_aer = 0.0
      hygro_aer = 0.0

      volumcut_sect = 0.0
      volumcen_sect = 0.0
      volumlo_sect = 0.0
      volumhi_sect = 0.0
      dcut_sect = 0.0
      dcen_sect = 0.0
      dlo_sect = 0.0
      dhi_sect = 0.0
      sigmag_aer = 0.0

      xcut_atype_md1 = 0.0
      xcut_atype_md2 = 0.0

      name_aer = '???'

      lptr_so4_aer = 0
      lptr_msa_aer = 0
      lptr_no3_aer = 0
      lptr_cl_aer = 0
      lptr_co3_aer = 0
      lptr_nh4_aer = 0
      lptr_na_aer = 0
      lptr_ca_aer = 0
      lptr_oin_aer = 0
      lptr_oc_aer = 0
      lptr_bc_aer = 0

      end if


      return
      end subroutine mosaic_allocate_memory
