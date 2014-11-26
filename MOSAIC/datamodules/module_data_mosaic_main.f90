      module module_data_mosaic_main

      use module_data_mosaic_kind, only:  r8
      use module_data_mosaic_constants, only:  &
          avogad, deg2rad, pi, piover4, piover6, third


      implicit none

      integer, parameter ::   &
                ngas_com = 40,   &
                ngas_urb = 19,   &
                ngas_bio =  7,   &
                ngas_mar = 11
      !BSINGH - 05/28/2013(RCE updates)
      integer, parameter ::   &
                  naer_tot = 24 		      ! total num of 3-D variables per bin

      integer, save ::   &
                naerbin  = -999888777         ! number of bins (set at run time)
      !BSINGH - 05/28/2013(RCE updates ENDS)
!               naerbin  = 41760  	      ! ( 48 size)*(29 wbc)*(30 kappa)
!               naerbin  = 3240  	      ! ( 24 size)*(15 wbc)*( 9 kappa)
!               naerbin  = 90000 	      ! (100 size)*(30 wbc)*(30 kappa)

      integer, parameter ::   &
                ncld_tot = 13,		   &  ! + 8 = total num of 3-D variables/bin
                ncldbin  =  4,		   &  ! num of cloud bins
                ncld     = 22		! num of dynamic cloud species/bin

      integer, parameter :: ngas_max = ngas_com + ngas_urb + ngas_bio + ngas_mar
      
      integer, parameter :: ncld_max = ncld_tot*ncldbin
      !BSINGH - 05/28/2013(RCE updates)
      integer, save :: naer_max = -999888777  ! set at run time to naer_tot*naerbin

      integer, save :: ntot_max = -999888777  ! set at run time to (ngas_max + naer_max + ncld_max)
      !BSINGH - 05/28/2013(RCE updates ENDS)

      integer, parameter :: kmaxd = 100		! num of vertical levels in a column mode

      integer, save ::   &
      		naerbin_used=0,   &   ! num of aerosol bins being used
      		ncldbin_used=0,   &   ! num of  cloud  bins being used
      		ntot_used=ngas_max    ! portion of cnn array being used

      integer, save ::   &
      		iwrite_gas,   &
      		iwrite_aer_bin,   &
      		iwrite_aer_dist,   &
      		iwrite_aer_species

      integer, save, allocatable ::        &!BSINGH - 05/28/2013(RCE updates)
      		lun_aer(:),   &
      		lun_aer_status(:)

      integer, save ::   &
      		lun_inp,   &
      		lun_gas,   &      		
      		lun_aeroptic,   &
      		lun_drydist,   &
      		lun_wetdist,   &
      		lun_species,   &
      		lun_fullout

      integer, save ::   &
                lun_sect_170,   & 
                lun_sect_171,   & 
                lun_sect_172,   & 
                lun_sect_180,   & 
                lun_sect_183,   & 
                lun_sect_184,   & 
                lun_sect_185,   & 
                lun_sect_186,   & 
                lun_sect_188,   & 
                lun_sect_190

      integer, save ::   &
      		ipmcmos = 0,   &   ! if > 0, do emissions, dilution, air density,
      		                   ! and relative humidity as in partmc_mosaic 
      		istate_pblh = 0    ! used for dilution from pblh changes

      real(r8), parameter :: press0_pa = 1.01325d5  ! pressure of 1 atm [Pa]
      real(r8), parameter :: mw_air = 28.966d0      ! dry-air mean molecular weight [g/mol]

      character(len=64), save, allocatable ::   &!BSINGH - 05/28/2013(RCE updates)
      		aer_output(:)

      character(len=64), save ::   &
                inputfile,   &
      		gas_output,   &
      		aeroptic_output,   &
      		drydist_output,   &
      		wetdist_output,   &
      		species_output,   &
      		fullout_fname

!------------------------------------------------------------------------
! Global Species Indices
!
      integer, save ::   &
       kh2so4,      khno3,       khcl,        knh3,        kno,   &
       kno2,        kno3,        kn2o5,       khono,       khno4,   &
       ko3,         ko1d,        ko3p,        koh,         kho2,   &
       kh2o2,       kco,         kso2,        kch4,        kc2h6,   &
       kch3o2,      kethp,       khcho,       kch3oh,      kanol,   &
       kch3ooh,     kethooh,     kald2,       khcooh,      krcooh,   &
       kc2o3,       kpan,   &
       karo1,       karo2,       kalk1,       kole1,       kapi1,   &
       kapi2,       klim1,       klim2,   &
       kpar,        kaone,       kmgly,       keth,        kolet,   &
       kolei,       ktol,        kxyl,        kcres,       kto2,   &
       kcro,        kopen,       konit,       krooh,       kro2,   &
       kano2,       knap,        kxo2,        kxpar,   &
       kisop,       kisoprd,     kisopp,      kisopn,      kisopo2,   &
       kapi,        klim,   &
       kdms,        kmsa,        kdmso,       kdmso2,      kch3so2h,   &
       kch3sch2oo,  kch3so2,     kch3so3,     kch3so2ch2oo,kch3so2oo,   &
       ksulfhox

      integer, save ::   &
       knum_a,      kdpdry_a,    ksigmag_a,  kjhyst_a,   &
       kwater_a,    kso4_a,      kno3_a,     kcl_a,       knh4_a,   &
       koc_a,       kmsa_a,      kco3_a,     kna_a,       kca_a,   &
       kbc_a,       koin_a,      karo1_a,    karo2_a,     kalk1_a,   &
       kole1_a,     kapi1_a,     kapi2_a,    klim1_a,     klim2_a

      integer, save ::   &
       knum_c,      kwater_c,    kso4_c,     kno3_c,     kcl_c,   &
       kmsa_c,      kco3_c,      knh4_c,     kna_c,      kca_c,   &
       koc_c,       kbc_c,       koin_c,   &
       karo1_c,     karo2_c,     kalk1_c,    kole1_c,    kapi1_c,   &
       kapi2_c,     klim1_c,     klim2_c



!-------------------------------------------------------------------------

      character(len=40), save, allocatable :: species(:)!BSINGH - 05/28/2013(RCE updates)

      integer, save ::	m_partmc_mosaic  ! >0 for partmc_mosaic, <=0 for mosaic box model!BSINGH - 05/28/2013(RCE updates)

      integer, save ::	mmode, mgas, maer, mcld

      integer, save ::	maeroptic, mshellcore

      real(r8), save, allocatable :: cnn(:)!BSINGH - 05/28/2013(RCE updates)

      real(r8), save :: tNi, tSi, tCli, tNH4i, DN, DS, DCl, DNH4

      real(r8), save, allocatable :: emission(:), emit(:)!BSINGH - 05/28/2013(RCE updates)

      integer, save ::   &
       tbeg_dd,   tbeg_mo,   &
       tbeg_hh,   tbeg_mm,   tbeg_ss,   tmar21_sec,   &
       trun_dd,   trun_hh,   trun_mm,   trun_ss,   &
       tbeg_sec,  trun_sec,   &
       nstep,     iprint

      real(r8), save ::   &
       tcur_sec,  tcur_min,  tcur_hrs,		   &  ! time since beginning of year (UTC)
       time_sec,  time_min,  time_hrs,		   &  ! time since beginning of simulation
       time_UTC,  time_UTC_beg,			   &  ! time of day in hrs (UTC)
       tmid_sec,  tsav_sec,  told_sec, time_sec_old,   &
       dt_sec,    dt_min,    dt_aeroptic_min,      &
       rlon,      rlat,                            &
       zalt_m,    cos_sza

      real(r8), save ::   &
       cair_mlc,  cair_molm3,     h2o,       o2,   &
       h2,        ppb,            speed_molec,     &
       te,        pr_atm,         RH,        pblh, &
       cair_mlc_old, cair_molm3_old,               &
       te_old,    pr_atm_old,     RH_old,    pblh_old

      integer, save ::   &
       idaytime,  msolar, mphoto


!      real, save ::   &
!       gas(ngas_max), flux(ngas_max,naerbin)

      ! SDE: OpenMP
!$OMP THREADPRIVATE( cnn                                  )
!$OMP THREADPRIVATE( cair_mlc, cair_molm3, h2o, o2, o2    )
!$OMP THREADPRIVATE( h2, speed_molec, te, pr_atm, RH      )
!$OMP THREADPRIVATE( pblh, cair_mlc_old, cair_molm3_old   )
!$OMP THREADPRIVATE( te_old, pr_atm_old, RH_old, pblh_old )

!------------------------------------------------------------------------

      end module module_data_mosaic_main
