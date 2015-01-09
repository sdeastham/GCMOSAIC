!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!     
! !MODULE: mosaic
!     
! !DESCRIPTION: Module MOSAIC\_MOD is the interface between GEOS-Chem 
!  and the MOSAIC aerosol thermodynamics and microphysics module.
!\\   
!\\   
! !INTERFACE: 
!
      MODULE MOSAIC_MOD
!
! !USES:
!
      USE HCO_ERROR_MOD  ! For real precisions (hp)

      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: CLEANUP_MOSAIC
      PUBLIC  :: DO_MOSAIC
!
! !PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: INIT_MOSAIC
      PRIVATE :: INIT_DATA_MODULES_INBUILT
      PRIVATE :: INIT_AEROSOL_INBUILT
      PRIVATE :: WATER_VAPOR
      PRIVATE :: SO4_PARTITION
!
! !REMARKS:
!  Original Author:
!**********************************************************************************  
! This computer software was prepared by Battelle Memorial Institute, hereiNa_after
! the Contractor, under Contract No. DE-AC05-76RL0 1830 with the Department of 
! Energy (DOE). NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY WARRANTY,
! EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
! Copyright (c) 2007 Battelle Memorial Institute
! Written by Rahul A. Zaveri
!
!********************************************************************************************
! Model         : MOSAIC (Model for Simulating Aerosol Interactions & Chemistry)
!
! Last Update   : October 2007
!
! Sub-Modules   : CBM-Z  (gas-phase photochemistry)
!                 ASTEM  (Adaptive Step Time-split Euler Method)
!                 MESA   (Multicomponent Equilibrium Solver for Aerosols)
!                 MTEM   (Multicomponent Taylor Expansion Method for activity coefficients)
!
! Purpose       : CBM-Z simulates lower tropospheric trace gas photochemistry for
!                 background, urban, biogenic, and marine (DMS) sub-regimes.
!
!                 ASTEM solves the dynamic partitioning of semi-volatile
!                 species between gas and partiCl_ae phases.
!
!                 MESA solves the multicomponent solid-liquid equilibria within
!                 the aerosol phase.
!
!                 MTEM computes the multicomponent activity coefficients of
!                 electrolytes in aqueous atmospheric aerosols.
!
! Author        : Rahul A. Zaveri, PhD
!                 Senior Research Scientist
!                 Pacific Northwest National Laboratory
!                 Atmospheric Sciences TechniCa_al Group
!                 P.O. Box 999, MSIN K9-30
!                 Richland, WA 99352
!                 Phone: (509) 372-6159, Fax: (509) 372-6168
!                 Email: Rahul.Zaveri@pnl.gov
!
! Bugs/Problems : Please report any bugs or problems to Rahul.Zaveri@pnl.gov
!
! Terms of Use  : (1) MOSAIC and its submodules CBM-Z, ASTEM, MESA, and MTEM may not be
!                     included in any commercial package, or used for any commercial
!                     appliCa_ations without prior authorization from the author.
!                 (2) The MOSAIC code may be used for educational or non-profit purposes
!                     only. Any other usage must be first approved by the author.
!                 (3) The MOSAIC code cannot be modified in any way or form or distributed
!                     without the author's prior consent.
!                 (4) No portion of the MOSAIC source code can be used in other codes
!                     without the author's prior consent.
!                 (5) The MOSAIC code is provided on an as-is basis, and the author
!                     bears no liability from its usage.
!                 (6) PubliCa_ations resulting from the usage of MOSAIC must cite
!                     the references below for proper acknowledgment.
!
! References    : Zaveri R.A., R.C. Easter, J.D. Fast, and L.K. Peters, Model
!                   for simulating aerosol interactions and chemistry (MOSAIC),
!                   J. Geophys. Res., in review.
!
!                 Zaveri R.A., R.C. Easter, and L.K. Peters (2005a) A computationally
!                   efficient multicomponent equilibrium solver for aerosols (MESA),
!                   J. Geophys. Res, 110, D24203, doi:10.1029/2004JD005618.
!
!                 Zaveri R.A., R.C. Easter, and A.S. Wexler (2005b) A new method for
!                   multicomponent activity coefficients of electrolytes in aqueous
!                   atmospheric aerosols, J. Geophys. Res., 110, D02201,
!                   doi:10.1029/2004JD004681.
!
!                 Zaveri R.A. and L.K. Peters (1999) A new lumped structure photoChemical
!                   mechanism for large-scale appliCa_ations
!
! Support       : Funding for the development and evaluation of MOSAIC and
!                 its sub-modules was provided by:
!                 (a) the U.S. Department of Energy (DOE) under the auspices of the
!                     Atmospheric Science Program (ASP) of the Office of BiologiCa_al and
!                     Environmental Research
!                 (b) the NASA Aerosol Program and NASA Earth Science Enterprise
!                 (c) the U.S. Environmental Protection Agency (EPA) Aerosol Program
!                 (d) PNNL Laboratory Directed Research and Development (LDRD) Program
!--------------------------------------------------------------------------------------------
!  Original v4.1a MOSAIC implementation into GEOS-Chem by
!  Sebastian D. Eastham and Jessica Kunke (SDE, JK, 2014-11-13)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
!
! PRIVATE TYPES:
!
      !=================================================================
      ! MODULE VARIABLES:
      !
      ! Arrays
      !
      ! X                  : ???
      !=================================================================

      ! Arrays

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS
!
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_mosaic
!
! !DESCRIPTION: Subroutine DO\_MOSAIC is the interface between the 
!  GEOS-Chem model and the aerosol model MOSAIC.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DO_MOSAIC( am_I_Root, Input_Opt, 
     &                      State_Met, State_Chm, RC )
!
! !USES:
!
      USE CHEMGRID_MOD,       ONLY : ITS_IN_THE_STRATMESO
      USE CMN_SIZE_MOD
      USE CMN_GCTM_MOD,       ONLY : AIRMW, PI
      USE ERROR_MOD,          ONLY : DEBUG_MSG
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE ERROR_MOD,          ONLY : SAFE_DIV
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE PRESSURE_MOD,       ONLY : GET_PCENTER
      USE HCO_EMISLIST_MOD,   ONLY : HCO_GetPtr
      USE TIME_MOD,           ONLY : GET_MONTH
      USE TIME_MOD,           ONLY : ITS_A_NEW_MONTH
      ! Gas-phase species
      USE TRACERID_MOD,       ONLY : IDTSO4
      USE TRACERID_MOD,       ONLY : IDTHNO3
      USE TRACERID_MOD,       ONLY : IDTHCl
      USE TRACERID_MOD,       ONLY : IDTNH3
      USE TRACERID_MOD,       ONLY : IDTMSA

      ! Bulk only
      USE TRACERID_MOD,       ONLY : IDTNH4, IDTNIT

      ! Aerosol species
      USE TRACERID_MOD,       ONLY : IDTaNUM,  IDTaHyst, IDTaH2O
      USE TRACERID_MOD,       ONLY : IDTaSO4,  IDTaNO3,  IDTaNH4
      USE TRACERID_MOD,       ONLY : IDTaCl,   IDTaCO3,  IDTaMSA
      USE TRACERID_MOD,       ONLY : IDTaNa,   IDTaCa,   IDTaOC
      USE TRACERID_MOD,       ONLY : IDTaBC,   IDTaOIN,  IDTSALA
      ! SOA not yet implemented
      USE TRACERID_MOD,       ONLY : IDTaAro1, IDTaAro2, IDTaAlk1
      USE TRACERID_MOD,       ONLY : IDTaOle1, IDTaApi1, IDTaApi2
      USE TRACERID_MOD,       ONLY : IDTaLim1, IDTaLim2
      ! MOSAIC indexing
      use module_data_mosaic_main, only: r8,naer_tot,ngas_max, !Parameters
     &   ntot_max,avogad,kh2so4,khno3,khcl,knh3,kmsa,knum_a,piover6,
     &   ksigmag_a,kjhyst_a,kwater_a,ntot_used,cnn,kdpdry_a,dt_sec,
     &   species

      use module_data_mosaic_aero, only: nbin_a_max,naer,ngas_volatile,
     &   jtotal,nbin_a
    
      ! CNN aerosol mass indices 
      use module_data_mosaic_aero, only: 
     &   iSO4_a,iNO3_a,iCl_a,iNH4_a,iMSA_a,iCO3_a,iNa_a,iCa_a,iOC_a,
     &   iBC_a,iOIN_a 

      ! Switches
      use module_data_mosaic_aero, only:
     &   jhyst_lo, jhyst_up, jhyst_undefined,mhyst_method,
     &   mhyst_uporlo_waterhyst,all_solid,all_liquid,mixed,no_aerosol

      use module_data_mosaic_main, only: pr_atm, RH, te, o2, h2, h2o,
     &   pr_atm_old, RH_old, te_old, speed_molec, pblh, pblh_old,
     &   cair_mlc, cair_molm3, cair_mlc_old, cair_molm3_old, ppb
     
      use module_data_mosaic_asect, only:
     &    dcen_sect, sigmag_aer 

!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS: 
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!  Original v4.1a MOSAIC implementation into GEOS-Chem by
!  Sebastian D. Eastham and Jessica Kunke (SDE, JK, 2014-11-13)
! 
! !REVISION HISTORY: 
!  13 Nov 2014 - S. D. Eastham - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! Array dimensions
      INTEGER, PARAMETER       :: NTEMP = 1

      ! Constants
      REAL*8, PARAMETER        :: R_UNIV = 8.31448d0
      REAL*8, PARAMETER        :: AVOGADRO = 6.022d23
!
! !LOCAL VARIABLES:
!
      INTEGER,DIMENSION(nbin_a_max)  :: iter_MESA
      REAL(r8),DIMENSION(nbin_a_max) :: gam_ratio, aH2O_a
      REAL(r8),DIMENSION(nbin_a_max) :: dp_wet_a

      LOGICAL, SAVE            :: FIRST = .TRUE.
      INTEGER                  :: I,    J,    L,    N, NOFFSET
      CHARACTER(LEN=255)       :: X
      real(r8)                 :: molecPerCm3 ! Molecules of air per cm3
      real(r8)                 :: boxPerCm3 ! 1/(Box volume in cm3)
      real(r8)                 :: cmpConv   ! kg -> umol/m3 (needs /(g/mol))
      real(r8)                 :: VOL       ! Box volume in m3

      ! Unused outputs from aerchemistry
      real(r8)                 :: cair_mol_m3, cair_mol_cc

      LOGICAL                  :: LPRT
      INTEGER                  :: N_TRACERS
      INTEGER,DIMENSION(nbin_a_max) :: jaerosolstate
      REAL*8                   :: XNUMOL(Input_Opt%N_TRACERS)
      REAL*8                   :: TRACER_MW_KG(Input_Opt%N_TRACERS)
      REAL*8                   :: TSO4, ASO4, GSO4
      REAL*8                   :: DRYDIAM

      ! Bulk calculation
      REAL(r8)                 :: MOLSO4, MOLNH4, MOLNO3, MOLMSA
      REAL(r8)                 :: MOLH, AVOGDIV, VOLDRY, P_hPa
                  
      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL*8, POINTER          :: STT(:,:,:,:)

      ! Debug
      INTEGER,PARAMETER        :: ITARG = 0 
      INTEGER,PARAMETER        :: JTARG = 0 
      INTEGER,PARAMETER        :: LTARG = 0

      ! For getting fields from HEMCO 
      LOGICAL                  :: aIR

      ! If TRUE, sea salt is also handled in the bulk model
      LOGICAL, PARAMETER       :: BULK_SS = .FALSE.

      !=================================================================
      ! DO_MOSAIC begins here!
      !=================================================================

      ! Assume success
      RC                   = GIGC_SUCCESS

      ! am I root? 
      aIR                  = am_I_Root

      ! Copy fields from INPUT_OPT to local variables for use below
      LPRT        = Input_Opt%LPRT  
      N_TRACERS   = Input_Opt%N_TRACERS
      XNUMOL      = Input_Opt%XNUMOL
      TRACER_MW_KG= Input_Opt%TRACER_MW_KG(1:N_TRACERS)
      DRYDIAM     = Input_Opt%DRYDIAM

      ! Location string
      X = 'DO_MOSAIC (mosaic_mod.F)'

      ! First-time initialization
      IF ( FIRST ) THEN

         ! Make sure certain tracers are defined
         ! Gas-phase first
         IF ( IDTSO4 == 0 ) 
     &       CALL ERROR_STOP( 'Gas-phase SO4 is undefined!', X)
         IF ( IDTHNO3 == 0 ) 
     &       CALL ERROR_STOP( 'Gas-phase HNO3 is undefined!', X)
         IF ( IDTNH3 == 0 ) 
     &       CALL ERROR_STOP( 'Gas-phase NH3 is undefined!', X)
         IF ( IDTHCl == 0 ) 
     &       CALL ERROR_STOP( 'Gas-phase HCl is undefined!', X)
         IF ( IDTMSA == 0 ) 
     &       CALL ERROR_STOP( 'Gas-phase MSA is undefined!', X)

#if defined( MOSAIC) 
#if defined( MBULK ) 
         IF ( IDTSALA == 0 )
     &       CALL ERROR_STOP( 'Seasalt is undefined!', X)
   
#else
         ! Check aerosol phase - only look at final bin
         IF ( IDTaSO4(NBINMOS)  == 0 ) 
     &       CALL ERROR_STOP( 'Aerosol-phase SO4 is undefined!', X)
         IF ( IDTaNO3(NBINMOS)  == 0 ) 
     &       CALL ERROR_STOP( 'Aerosol-phase NO3 is undefined!', X)
         IF ( IDTaNH4(NBINMOS)  == 0 ) 
     &       CALL ERROR_STOP( 'Aerosol-phase NH4 is undefined!', X)
         IF ( IDTaCl(NBINMOS)  == 0 ) 
     &       CALL ERROR_STOP( 'Aerosol-phase Cl is undefined!', X)
         IF ( IDTaCO3(NBINMOS)  == 0 ) 
     &       CALL ERROR_STOP( 'Aerosol-phase CO3 is undefined!', X)
         IF ( IDTaNa(NBINMOS)  == 0 ) 
     &       CALL ERROR_STOP( 'Aerosol-phase Na is undefined!', X)
         IF ( IDTaCa(NBINMOS)  == 0 ) 
     &       CALL ERROR_STOP( 'Aerosol-phase Ca is undefined!', X)
         IF ( IDTaOC(NBINMOS)  == 0 ) 
     &       CALL ERROR_STOP( 'Aerosol-phase OC is undefined!', X)
         IF ( IDTaBC(NBINMOS)  == 0 ) 
     &       CALL ERROR_STOP( 'Aerosol-phase BC is undefined!', X)
         IF ( IDTaOIN(NBINMOS)  == 0 ) 
     &       CALL ERROR_STOP( 'Aerosol-phase OIN is undefined!', X)
         IF ( IDTaNUM(NBINMOS)  == 0 ) 
     &       CALL ERROR_STOP( 'Number concentration is undefined!', X)
         IF ( IDTaMSA(NBINMOS)  == 0 ) 
     &       CALL ERROR_STOP( 'Aerosol-phase MSA is undefined!', X)
#endif
#endif

         ! Aerosol data necessary to both bulk and sectional
         IF ( IDTaHyst(NBINMOS)  == 0 ) 
     &       CALL ERROR_STOP( 'Hysteresis tracer is undefined!', X)
         IF ( IDTaH2O(NBINMOS)  == 0 ) 
     &       CALL ERROR_STOP( 'Aerosol-phase H2O is undefined!', X)

         ! Initialize arrays
         CALL INIT_MOSAIC

         ppb = 1.d9

         ! Reset first-time flag
         FIRST = .FALSE. 
      ENDIF
 
      ! Factor needed for bulk calculation
!      AVOGDIV = 1.d0/AVOGADRO
     
      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      ! (mpayer, 12/6/12)
      STT => State_Chm%Tracers

      !=================================================================
      ! Loop over grid boxes and call MOSAIC
      !=================================================================

!$OMP PARALLEL
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,    J,      L,     N,                  VOL )
!$OMP+PRIVATE(      NOFFSET,boxPerCm3,molecPerCm3,cmpConv  )
!$OMP+PRIVATE( iter_MESA, aH2O_a, gam_ratio, jaerosolstate  )
!$OMP+PRIVATE( dp_wet_a, cair_mol_m3, cair_mol_cc, TSO4     )
!$OMP+PRIVATE( GSO4, ASO4, MOLSO4, MOLNH4, MOLNO3        )
!$OMP+PRIVATE( MOLMSA, MOLH, VOLDRY, P_hPa                )

      ! Memory allocation and initialization
      allocate( cnn(ntot_max) )

!$OMP DO
!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Only applying MOSAIC in troposphere
         IF ( ITS_IN_THE_STRATMESO( I, J, L, State_Met ) ) THEN 
            ! TEMPORARY - Move all SO4 into gas-phase for UCX
#if !defined( MBULK )
            DO N = 1,NBINMOS
               IF (STT(I,J,L,IDTaSO4(N)) > 0) THEN
                   STT(I,J,L,IDTSO4) = STT(I,J,L,IDTSO4) + 
     &                                 STT(I,J,L,IDTaSO4(N))
                   STT(I,J,L,IDTaSO4(N)) = 0.d0
               ENDIF
            ENDDO
#endif
            CYCLE
         ENDIF

         ! Reinitialize
         jaerosolstate(:) = 0
         dp_wet_a(:) = 0

         ! Set MOSAIC variables
         RH       = State_Met%RH(I,J,L)   ! %
         te       = State_Met%T(I,J,L)    ! K
         te_old   = te                    ! K
!         P_hPa    = State_Met%PMID(I,J,L) ! hPa
         P_hPa    = GET_PCENTER(I,J,L)
         pr_atm   = P_hPa * 9.86923267d-4! atm
         pr_atm_old= pr_atm

!         ! DEBUG - Cycle if pressure < 200 hPa
!         IF (P_hPa .lt. 200.d0) THEN
!            CYCLE
!         ENDIF

         ! Volume of grid box [m3] 
         VOL      = State_Met%AIRVOL(I,J,L)

         ! Conversion factors
         ! Gridboxes per cm3
         boxPerCm3   = (1.d-6/VOL)
         cmpConv = 1.d6/VOL

         ! Set the MOSAIC air composition
         ! mol/m3
         ! kg/m3 * (g/kg) * (mol/g)
         cair_molm3 = State_Met%AIRDEN(L,I,J) * 1.d3
     &                            / AIRMW
         ! molec/cm3
         cair_mlc = 1.d-6*avogadro*cair_molm3
         cair_mlc_old = cair_mlc
         cair_molm3_old = cair_molm3

         h2o = WATER_VAPOR(RH, cair_mlc, te, pr_atm)
         o2  = 0.21*cair_mlc
         h2  = 0.58e-6*cair_mlc

         ! Molecules of air per cm3
         ! kg/m3 * (m3/cm3 * g/kg) * molec/mol * (mol/g)
         molecPerCm3 = cair_mlc

         ! Build the MOSAIC interface array
         CNN = 0d0

         ! STT is in kg
         ! MOSAIC wants molec/cm3
         CNN(khno3)   = XNUMOL(IDTHNO3)*boxPerCm3*STT(I,J,L,IDTHNO3)
         CNN(knh3)    = XNUMOL(IDTNH3)*boxPerCm3*STT(I,J,L,IDTNH3)

#if defined( MBULK )
         ! Assume that H2SO4 partitions according to Kulmala
         TSO4 = STT(I,J,L,IDTSO4)
         CALL SO4_PARTITION(TSO4,State_Met%AD(I,J,L),P_hPa,te,
     &                     ASO4,GSO4)
         CNN(kh2so4) = XNUMOL(IDTSO4) * boxPerCm3 * GSO4

         ! Assume no gas-phase MSA
         CNN(kmsa) = 0d0

         ! Take seasalt as 30.61% by weight Na+ and 55.04% by weight Cl,
         ! 1.16% by weight Ca2+ and ignore the rest (4.79%). At some point
         ! should put this (K+ and Mg+) into "other inorganic"...

         ! CNN aerosols should be umol/m3
         noffset = ngas_max
         CNN(noffset + kwater_a+iso4_a) = ASO4 * 
     &              cmpConv / TRACER_MW_KG(IDTSO4)
         CNN(noffset + kwater_a+ino3_a) = STT(I,J,L,IDTNIT) * 
     &              cmpConv / TRACER_MW_KG(IDTNIT)
         CNN(noffset + kwater_a+iNH4_a) = STT(I,J,L,IDTNH4) * 
     &              cmpConv / TRACER_MW_KG(IDTNH4)
         CNN(noffset + kwater_a+iMSA_a) = STT(I,J,L,IDTMSA) * 
     &              cmpConv / TRACER_MW_KG(IDTMSA)

         IF (BULK_SS) THEN
            CNN(noffset + kwater_a+iNa_a)  = STT(I,J,L,IDTSALA) * 
     &                 0.5504d0 * 1.d9 / ( 35.45d0  * VOL  )
            CNN(noffset + kwater_a+iCa_a)  = STT(I,J,L,IDTSALA) * 
     &                 0.0116d0 * 1.d9 / ( 40.08d0  * VOL  )
            CNN(noffset + kwater_a+iCl_a)  = STT(I,J,L,IDTSALA) * 
     &                 0.3061d0 * 1.d9 / ( 22.99d0  * VOL  )
         ENDIF

         ! Calculate number directly for a 1-bin model
         N = 1
         
         CNN(noffset + kdpdry_a) = DRYDIAM*1.d-4 ! Dry diameter in um - convert to cm

         ! Estimate local number density of aerosol
         ! This is probably inefficient but should be safe
         ! MOLX = mol(x) per m3(air)
         MOLSO4 = CNN(noffset + kwater_a+iso4_a)*1.d-6
         MOLNO3 = CNN(noffset + kwater_a+iNO3_a)*1.d-6
         MOLNH4 = CNN(noffset + kwater_a+iNH4_a)*1.d-6
         MOLMSA = CNN(noffset + kwater_a+iMSA_a)*1.d-6
         MOLH   = (2*MOLSO4) + MOLNO3 + MOLMSA - MOLNH4
         MOLH   = max(MOLH,0d0)

         ! Assume a density of 1000 kg/m3 for H+ (ie MW/rho = 1)
         ! Assuming 1800 kg/m3 density for SO4, NO3, NH4 and MSA
         ! (same as used in MOSAIC)
         ! VOLDRY is (cm3 aerosol)/(cm3 air)
         VOLDRY = (MOLSO4*TRACER_MW_KG(IDTSO4)/1.8d3) +
     &            (MOLNO3*TRACER_MW_KG(IDTNIT)/1.8d3) +
     &            (MOLNH4*TRACER_MW_KG(IDTNH4)/1.8d3) +
     &            (MOLMSA*TRACER_MW_KG(IDTMSA)/1.8d3) + 
     &            (MOLH/1.0d3)

         ! Number per cm3
         CNN(noffset + knum_a) = 6.d0*VOLDRY/(PI*((DRYDIAM*1.d-4)**3))
         !cnn(noffset + ksigmag_a) = sigmag_aer(1,1)
         CNN(noffset + ksigmag_a) = 1.6d0     ! Geometric SD - not used
         CNN(noffset + kwater_a) = STT(I,J,L,IDTaH2O(N))/VOL
         ! Hysteresis
         IF (mhyst_method == mhyst_uporlo_waterhyst) THEN
            CNN(noffset + kjhyst_a) = STT(I,J,L,IDTaHyst(N))/VOL
         ELSE
            CNN(noffset + kjhyst_a) = 0d0
         ENDIF

         !SDE DEBUG==================================================================
         IF ((I.eq.ITARG).and.(J.eq.JTARG).and.(L.eq.LTARG)) THEN
            WRITE(6,*) 'CNN SO4g pre MOSAIC:  ', CNN(kh2so4)
            WRITE(6,*) 'CNN SO4a pre MOSAIC:  ', CNN(noffset + kwater_a+
     &                                          iso4_a)
            WRITE(6,*) 'RH (%):               ', RH
            WRITE(6,*) 'CNN pre MOSAIC'
            WRITE(6,*) '====================================='
            DO N=1,NTOT_MAX
               IF (trim(species(N)).ne."???") THEN
                  WRITE(6,'(a,I03,x,a20,x,E10.4)') 
     &                      'I/SPC/TOT:        ', 
     &                      N,trim(species(N)),CNN(N)
               ENDIF
            ENDDO
            WRITE(6,*) '====================================='
            N = 1
            WRITE(6,*) 'TOTALS pre MOSAIC'
            WRITE(6,*) '====================================='
            WRITE(6,*) 'STT aH2O:         ', STT(I,J,L,IDTaH2O(1))
            WRITE(6,*) 'STT SO4:          ', STT(I,J,L,IDTSO4)
            WRITE(6,*) 'STT TNO3 as NO3:  ', 
     &           (STT(I,J,L,IDTHNO3)*62.d0/63.d0)+STT(I,J,L,IDTNIT)
            WRITE(6,*) ' => HNO3 as NO3:  ', STT(I,J,L,IDTHNO3)*62.d0/63.d0
            WRITE(6,*) ' => NIT:          ', STT(I,J,L,IDTNIT)
            WRITE(6,*) 'STT TNH3 as NH3:  ', STT(I,J,L,IDTNH3) + 
     &           (STT(I,J,L,IDTNH4)*17.d0/18.d0)
            WRITE(6,*) ' => NH3:          ', STT(I,J,L,IDTNH3)
            WRITE(6,*) ' => NH4 as NH3:   ', STT(I,J,L,IDTNH4)*17.d0/18.d0
            WRITE(6,*) '====================================='
         ENDIF
         !SDE DEBUG==================================================================

#else
         ! STT now stores only gas phase SO4 and MSA
         CNN(kh2so4)  = XNUMOL(IDTSO4)*boxPerCm3*STT(I,J,L,IDTSO4)
         CNN(kmsa)    = XNUMOL(IDTMSA)*boxPerCm3*STT(I,J,L,IDTMSA)
         CNN(khcl)    = XNUMOL(IDTHCl)*boxPerCm3*STT(I,J,L,IDTHCl)
         DO N = 1, NBINMOS
            NOFFSET = ngas_max + naer_tot*(n - 1)
            ! #/cm3 from "kg" - convert to moles, then molecules
            ! Divide by volume
            CNN(NOFFSET+knum_a) = (STT(I,J,L,IDTaNUM(N)) / 
     &         TRACER_MW_KG(IDTaNUM(N))) * AVOGADRO * boxPerCm3

            ! These would be used in a modal representation
            cnn(noffset + kdpdry_a) = dcen_sect(N,1)*1.d4
            cnn(noffset + ksigmag_a) = sigmag_aer(N,1)
!            CNN(NOFFSET+kdpdry_a) = AERDRYDIAM(N) ! Dry diameter in um
!            CNN(NOFFSET+ksigmag_a) = 0d0       ! Geometric SD

            ! Aerosol water content (kg/m3 air)
            CNN(NOFFSET+kwater_a) = STT(I,J,L,IDTaH2O(N))/VOL

            ! Hysteresis
            IF (mhyst_method == mhyst_uporlo_waterhyst) THEN
               CNN(NOFFSET+kjhyst_a) = STT(I,J,L,IDTaHyst(N))/VOL
            ELSE
               CNN(NOFFSET+kjhyst_a) = 0d0
            ENDIF

            ! Aerosol components (umol/m3)
            CNN(NOFFSET+kwater_a+iso4_a) = STT(I,J,L,IDTaSO4(N)) * 
     &                 cmpConv / TRACER_MW_KG(IDTaSO4(N))
            CNN(NOFFSET+kwater_a+ino3_a) = STT(I,J,L,IDTaNO3(N)) * 
     &                 cmpConv / TRACER_MW_KG(IDTaNO3(N))
            CNN(NOFFSET+kwater_a+iCl_a)  = STT(I,J,L,IDTaCl(N)) * 
     &                 cmpConv / TRACER_MW_KG(IDTaCl(N))
            CNN(NOFFSET+kwater_a+iNH4_a) = STT(I,J,L,IDTaNH4(N)) * 
     &                 cmpConv / TRACER_MW_KG(IDTaNH4(N))
            CNN(NOFFSET+kwater_a+iCO3_a) = STT(I,J,L,IDTaCO3(N)) * 
     &                 cmpConv / TRACER_MW_KG(IDTaCO3(N))
            CNN(NOFFSET+kwater_a+iMSA_a) = STT(I,J,L,IDTaMSA(N)) * 
     &                 cmpConv / TRACER_MW_KG(IDTaMSA(N))
            CNN(NOFFSET+kwater_a+iNa_a)  = STT(I,J,L,IDTaNa(N)) * 
     &                 cmpConv / TRACER_MW_KG(IDTaNa(N))
            CNN(NOFFSET+kwater_a+iCa_a)  = STT(I,J,L,IDTaCa(N)) * 
     &                 cmpConv / TRACER_MW_KG(IDTaCa(N))
            CNN(NOFFSET+kwater_a+iOC_a)  = STT(I,J,L,IDTaOC(N)) * 
     &                 cmpConv / TRACER_MW_KG(IDTaOC(N))
            CNN(NOFFSET+kwater_a+iBC_a)  = STT(I,J,L,IDTaBC(N)) * 
     &                 cmpConv / TRACER_MW_KG(IDTaBC(N))
            CNN(NOFFSET+kwater_a+iOIN_a) = STT(I,J,L,IDTaOIN(N)) * 
     &                 cmpConv / TRACER_MW_KG(IDTaOIN(N))
            ! SOA not yet implemented
         ENDDO
#endif

         ! Run MOSAIC
            CALL aerchemistry(1, dt_sec, 0.d0,                        !intent-ins
     &       jaerosolstate, dp_wet_a,                                 !intent-inouts
     &       cair_mol_m3, cair_mol_cc, aH2O_a, gam_ratio, iter_MESA ) !intent-outs

         ! Save data back into tracer array
         ! STT wants kg
         ! MOSAIC is in molec/cm3
         STT(I,J,L,IDTHNO3) = CNN(kHNO3) /(boxPerCm3*XNUMOL(IDTHNO3))
         STT(I,J,L,IDTNH3)  = CNN(kNH3)  /(boxPerCm3*XNUMOL(IDTNH3))

#if defined( MBULK )
         N = 1
         ! Add together aerosol and gas-phase SO4
         STT(I,J,L,IDTSO4)  = CNN(kh2so4)/(boxPerCm3*XNUMOL(IDTSO4)) + 
     &     (CNN(noffset+kwater_a+iso4_a)*TRACER_MW_KG(IDTSO4)/cmpConv)

         ! Add together aerosol and gas-phase MSA
         STT(I,J,L,IDTMSA)  = CNN(kMSA)/(boxPerCm3*XNUMOL(IDTMSA)) + 
     &     (CNN(noffset+kwater_a+iMSA_a)*TRACER_MW_KG(IDTMSA)/cmpConv)

         ! Return nitrate and ammonium to GEOS-Chem
         STT(I,J,L,IDTNIT)  =  
     &     (CNN(noffset+kwater_a+iNO3_a)*TRACER_MW_KG(IDTNIT)/cmpConv)

         STT(I,J,L,IDTNH4)  =  
     &     (CNN(noffset+kwater_a+iNH4_a)*TRACER_MW_KG(IDTNH4)/cmpConv)

         ! Just 1 aerosol with a specified size
         STT(I,J,L,IDTaH2O(N)) = CNN(noffset+kwater_a)*VOL

         ! Hysteresis
         IF (mhyst_method == mhyst_uporlo_waterhyst) THEN
            STT(I,J,L,IDTaHyst(N)) = CNN(noffset+kjhyst_a)*VOL
         ENDIF

         !SDE DEBUG==================================================================
         IF ((I.eq.ITARG).and.(J.eq.JTARG).and.(L.eq.LTARG)) THEN
            WRITE(6,*) 'CNN SO4g post MOSAIC: ', CNN(kh2so4)
            WRITE(6,*) 'CNN SO4a post MOSAIC: ', CNN(noffset + kwater_a+
     &                                          iso4_a)
            WRITE(6,*) 'JAEROSOLSTATE:        ', jaerosolstate(1)
            WRITE(6,*) 'RH (%):               ', RH
            WRITE(6,*) 'CNN post MOSAIC'
            WRITE(6,*) '====================================='
            DO N=1,NTOT_MAX
               IF (trim(species(N)).ne."???") THEN
                  WRITE(6,'(a,I03,x,a20,x,E10.4)') 
     &                      'I/SPC/TOT:        ', 
     &                      N,trim(species(N)),CNN(N)
               ENDIF
            ENDDO
            WRITE(6,*) '====================================='
            N = 1
            WRITE(6,*) 'TOTALS post MOSAIC'
            WRITE(6,*) '====================================='
            WRITE(6,*) 'STT aH2O:         ', STT(I,J,L,IDTaH2O(1))
            WRITE(6,*) 'STT SO4:          ', STT(I,J,L,IDTSO4)
            WRITE(6,*) 'STT TNO3 as NO3:  ', 
     &           (STT(I,J,L,IDTHNO3)*62.d0/63.d0)+STT(I,J,L,IDTNIT)
            WRITE(6,*) ' => HNO3 as NO3:  ', STT(I,J,L,IDTHNO3)*62.d0/63.d0
            WRITE(6,*) ' => NIT:          ', STT(I,J,L,IDTNIT)
            WRITE(6,*) 'STT TNH3 as NH3:  ', STT(I,J,L,IDTNH3) + 
     &           (STT(I,J,L,IDTNH4)*17.d0/18.d0)
            WRITE(6,*) ' => NH3:          ', STT(I,J,L,IDTNH3)
            WRITE(6,*) ' => NH4 as NH3:   ', STT(I,J,L,IDTNH4)*17.d0/18.d0
            WRITE(6,*) '====================================='
         ENDIF
         !SDE DEBUG==================================================================

#else
         STT(I,J,L,IDTMSA)  = CNN(kMSA)  /(boxPerCm3*XNUMOL(IDTMSA))
         STT(I,J,L,IDTSO4)  = CNN(kh2so4)/(boxPerCm3*XNUMOL(IDTSO4))
         STT(I,J,L,IDTHCl)  = CNN(kHCl)  /(boxPerCm3*XNUMOL(IDTHCl))
         ! Reverse all previous conversions
         DO N = 1, NBINMOS
            NOFFSET = ngas_max + naer_tot*(n - 1)
            ! #/cm3 from ppbv (ppbv ~ # per bn molecules of air)
            STT(I,J,L,IDTaNUM(N)) = CNN(NOFFSET+knum_a) * 
     &         TRACER_MW_KG(IDTaNUM(N)) / (AVOGADRO*boxPerCm3)
            STT(I,J,L,IDTaH2O(N)) = CNN(NOFFSET+kwater_a) * VOL
            IF (mhyst_method == mhyst_uporlo_waterhyst) THEN
               STT(I,J,L,IDTaHyst) = CNN(NOFFSET+kjhyst_a) * VOL
            ENDIF

            ! Aerosol components
            STT(I,J,L,IDTaSO4(N)) = CNN(NOFFSET+kwater_a+iso4_a) *
     &                 TRACER_MW_KG(IDTaSO4(N)) / cmpConv
            STT(I,J,L,IDTaNO3(N)) = CNN(NOFFSET+kwater_a+ino3_a) *
     &                 TRACER_MW_KG(IDTaNO3(N)) / cmpConv
            STT(I,J,L,IDTaCl(N)) = CNN(NOFFSET+kwater_a+iCl_a) *
     &                 TRACER_MW_KG(IDTaCl(N)) / cmpConv
            STT(I,J,L,IDTaNH4(N)) = CNN(NOFFSET+kwater_a+iNH4_a) *
     &                 TRACER_MW_KG(IDTaNH4(N)) / cmpConv
            STT(I,J,L,IDTaCO3(N)) = CNN(NOFFSET+kwater_a+iCO3_a) *
     &                 TRACER_MW_KG(IDTaCO3(N)) / cmpConv
            STT(I,J,L,IDTaMSA(N)) = CNN(NOFFSET+kwater_a+iMSA_a) *
     &                 TRACER_MW_KG(IDTaMSA(N)) / cmpConv
            STT(I,J,L,IDTaNa(N)) = CNN(NOFFSET+kwater_a+iNa_a) *
     &                 TRACER_MW_KG(IDTaNa(N)) / cmpConv
            STT(I,J,L,IDTaCa(N)) = CNN(NOFFSET+kwater_a+iCa_a) *
     &                 TRACER_MW_KG(IDTaCa(N)) / cmpConv
            STT(I,J,L,IDTaOC(N)) = CNN(NOFFSET+kwater_a+iOC_a) *
     &                 TRACER_MW_KG(IDTaOC(N)) / cmpConv
            STT(I,J,L,IDTaBC(N)) = CNN(NOFFSET+kwater_a+iBC_a) *
     &                 TRACER_MW_KG(IDTaBC(N)) / cmpConv
            STT(I,J,L,IDTaOIN(N)) = CNN(NOFFSET+kwater_a+iOIN_a) *
     &                 TRACER_MW_KG(IDTaOIN(N)) / cmpConv
            ! SOA not yet implemented
         ENDDO
#endif
   
      ENDDO
      ENDDO
      ENDDO 
!$OMP END DO
      if (allocated(cnn)) deallocate(cnn)
!$OMP END PARALLEL

      ! Free pointer
      NULLIFY( STT )

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### MOSAIC: a DO_MOSAIC' )

      ! Return to calling program
      END SUBROUTINE DO_MOSAIC
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: water_vapor
!
! !DESCRIPTION: Function WATER\_VAPOR stands in for the MOSAIC   
!  function WaterVapor.
!\\
!\\
! !INTERFACE:
!
      FUNCTION WATER_VAPOR(rh_arg, cair_mlc_arg, te_arg, pr_atm_arg)
     &          RESULT(WVAP)
!
! !USES:
!
      USE MODULE_DATA_MOSAIC_KIND, ONLY : r8
!
! !INPUT PARAMETERS:
!
!
! !INPUT/OUTPUT PARAMETERS: 
!
!
! !OUTPUT PARAMETERS:
!
!
! !REMARKS:
!  Original v4.1a MOSAIC implementation into GEOS-Chem by
!  Sebastian D. Eastham and Jessica Kunke (SDE, JK, 2014-11-13)
! 
! !REVISION HISTORY: 
!  17 Nov 2014 - S. D. Eastham - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(r8),INTENT(IN)  :: RH_ARG, CAIR_MLC_ARG, TE_ARG, PR_ATM_ARG
      REAL(r8)             :: T_STEAM, PR_STD, A, ARG, PR_H2O
      REAL(r8)             :: WVAP

      !=================================================================
      ! WATER_VAPOR begins here!
      !=================================================================

      t_steam = 373.15 ! steam temperature  [K]
      pr_std   = 1.0   ! standard pressure  [atm]

      a      = 1.0 - t_steam/te_arg
      arg    = (((-.1299*a -.6445)*a -1.976)*a +13.3185)*a
      pr_h2o = pr_std*exp(arg)  ! [atm]
      WVAP  = RH_arg*(pr_h2o/pr_atm_arg)*cair_mlc_arg/100. ! [molec/cc]

      return
      END FUNCTION WATER_VAPOR
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CLEANUP_MOSAIC
!
! !DESCRIPTION: Subroutine CLEANUP\_MOSAIC deallocates MOSAIC
!  arrays and variables.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_MOSAIC
!
! !USES:
!
      use module_data_mosaic_main
      use module_data_mosaic_aero
      use module_data_mosaic_asect
!
! !INPUT PARAMETERS:
!
!
! !INPUT/OUTPUT PARAMETERS: 
!
!
! !OUTPUT PARAMETERS:
!
!
! !REMARKS:
!  Original v4.1a MOSAIC implementation into GEOS-Chem by
!  Sebastian D. Eastham and Jessica Kunke (SDE, JK, 2014-11-13)
! 
! !REVISION HISTORY: 
!  17 Nov 2014 - S. D. Eastham - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

      !=================================================================
      ! CLEANUP_MOSAIC  begins here!
      !=================================================================

      IF (ALLOCATED( lun_aer )) DEALLOCATE( lun_aer )
      IF (ALLOCATED( lun_aer_status )) DEALLOCATE(lun_aer_status)

      IF (ALLOCATED( aer_output )) DEALLOCATE(aer_output)

      IF (ALLOCATED( species )) DEALLOCATE( species )

      IF (ALLOCATED( cnn )) DEALLOCATE( cnn )

      IF (ALLOCATED( emission )) DEALLOCATE( emission )

      IF (ALLOCATED( nsize_aer )) DEALLOCATE( nsize_aer )
      IF (ALLOCATED( ncomp_aer )) DEALLOCATE( ncomp_aer )
      IF (ALLOCATED( ncomp_plustracer_aer )) 
     &                     DEALLOCATE( ncomp_plustracer_aer)
      IF (ALLOCATED( mastercompptr_aer )) DEALLOCATE( mastercompptr_aer)
      IF (ALLOCATED( massptr_aer )) DEALLOCATE( massptr_aer )
      
      IF (ALLOCATED( waterptr_aer )) DEALLOCATE( waterptr_aer )
      IF (ALLOCATED( hyswptr_aer )) DEALLOCATE( hyswptr_aer )
      IF (ALLOCATED( numptr_aer )) DEALLOCATE( numptr_aer )
      
      IF (ALLOCATED( mprognum_aer )) DEALLOCATE( mprognum_aer )
      
      IF (ALLOCATED( ibin_of_isize_itype )) 
     &                    DEALLOCATE(ibin_of_isize_itype)
      IF (ALLOCATED( isize_of_ibin )) DEALLOCATE(isize_of_ibin)
      IF (ALLOCATED( itype_of_ibin )) DEALLOCATE(itype_of_ibin)
      IF (ALLOCATED( itype_of_itype_md1md2 ))
     &                    DEALLOCATE( itype_of_itype_md1md2 )
      IF (ALLOCATED( itype_md1_of_itype ))
     &                    DEALLOCATE( itype_md1_of_itype )
      IF (ALLOCATED( itype_md2_of_itype )) 
     &                    DEALLOCATE( itype_md2_of_itype )
      
      IF (ALLOCATED( dens_aer )) DEALLOCATE( dens_aer )
      IF (ALLOCATED( mw_aer )) DEALLOCATE( mw_aer )
      IF (ALLOCATED( hygro_aer )) DEALLOCATE( hygro_aer )
      
      IF (ALLOCATED( volumcut_sect )) DEALLOCATE( volumcut_sect )
      IF (ALLOCATED( volumcen_sect )) DEALLOCATE( volumcen_sect )
      IF (ALLOCATED( volumlo_sect )) DEALLOCATE( volumlo_sect )
      IF (ALLOCATED( volumhi_sect )) DEALLOCATE( volumhi_sect )
      IF (ALLOCATED( dcut_sect )) DEALLOCATE( dcut_sect )
      IF (ALLOCATED( dcen_sect )) DEALLOCATE( dcen_sect )
      IF (ALLOCATED( dlo_sect )) DEALLOCATE( dlo_sect )
      IF (ALLOCATED( dhi_sect )) DEALLOCATE( dhi_sect )
      IF (ALLOCATED( sigmag_aer )) DEALLOCATE( sigmag_aer )
      
      IF (ALLOCATED( xcut_atype_md1 )) DEALLOCATE( xcut_atype_md1 )
      IF (ALLOCATED( xcut_atype_md2 )) DEALLOCATE( xcut_atype_md2 )
      
      IF (ALLOCATED( name_aer )) DEALLOCATE( name_aer )
      
      IF (ALLOCATED( lptr_so4_aer )) DEALLOCATE( lptr_so4_aer )
      IF (ALLOCATED( lptr_msa_aer )) DEALLOCATE( lptr_msa_aer )
      IF (ALLOCATED( lptr_no3_aer )) DEALLOCATE( lptr_no3_aer )
      IF (ALLOCATED( lptr_cl_aer )) DEALLOCATE( lptr_cl_aer )
      IF (ALLOCATED( lptr_co3_aer )) DEALLOCATE( lptr_co3_aer )
      IF (ALLOCATED( lptr_nh4_aer )) DEALLOCATE( lptr_nh4_aer )
      IF (ALLOCATED( lptr_na_aer )) DEALLOCATE( lptr_na_aer )
      IF (ALLOCATED( lptr_ca_aer )) DEALLOCATE( lptr_ca_aer )
      IF (ALLOCATED( lptr_oin_aer )) DEALLOCATE( lptr_oin_aer )
      IF (ALLOCATED( lptr_oc_aer )) DEALLOCATE( lptr_oc_aer )
      IF (ALLOCATED( lptr_bc_aer )) DEALLOCATE( lptr_bc_aer )

      END SUBROUTINE CLEANUP_MOSAIC
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: INIT_MOSAIC
!
! !DESCRIPTION: Subroutine INIT\_MOSAIC initializes MOSAIC
!  arrays and variables.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_MOSAIC
!
! !USES:
!
      use module_data_mosaic_main
      use module_data_mosaic_gas
      use module_data_mosaic_aero
      use module_data_mosaic_asect
      use module_data_mosaic_cloud

      use module_data_mosaic_pmcmos, only:  
     &     temp_profile_fname, pblh_profile_fname,  
     &     aer_init_fname, aer_back_fname, aer_emit_fname,  
     &     gas_init_fname, gas_back_fname, gas_emit_fname,  
           !BSINGH - 05/28/2013(RCE updates)
     &     aer_emit_adjust_fac,  
     &     msolar_pmcmos_dtshift,  
     &     pmcmos_inputfile_version, pmcmos_inputfile_directory
      !BSINGH - 05/28/2013(RCE updates ENDS)
      use module_pmcmos_init, only:  pmcmos_init
      
      use module_sect_iface, only:  !BSINGH - 05/28/2013(RCE updates)
     &     sect_iface_allocate_memory

      USE TIME_MOD, ONLY : GET_TS_CHEM
!
! !INPUT PARAMETERS:
!
!
! !INPUT/OUTPUT PARAMETERS: 
!
!
! !OUTPUT PARAMETERS:
!
!
! !REMARKS:
!  Original v4.1a MOSAIC implementation into GEOS-Chem by
!  Sebastian D. Eastham and Jessica Kunke (SDE, JK, 2014-11-13)
! 
! !REVISION HISTORY: 
!  17 Nov 2014 - S. D. Eastham - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      integer i, ibin, idum, input_version
      integer k
      integer noffset, nsize_aer_tmp
      !BSINGH - 05/28/2013(RCE updates ENDS)
      real(r8) :: tmpa, tmpb

      !=================================================================
      ! INIT_MOSAIC  begins here!
      !=================================================================

      ! Initialize MOSAIC data modules
      CALL init_data_modules_inbuilt

      ! Set the MOSAIC time step to match the chemistry time step
      dt_min = 1.d0*GET_TS_CHEM()
      dt_sec = 60.d0*dt_min 

      ! Adapted from MOSAIC's "ReadInput"
      ! Not relevant to GEOS-Chem
      m_partmc_mosaic = 0!BSINGH - 05/28/2013(RCE updates)

      !===================================================================!
      ! SDE - Do not use multi-dimensional sectional approach             !
      ! It has not been tested and there is definitely some setup code    !
      ! not yet implemented!                                              !
      !===================================================================!
      ! Number of aerosol types (?)
      ! Related to multi-dimensional sectional approach
      ntype_md1_aer = 1
      ntype_md2_aer = 1
      ntype_aer = ntype_md1_aer*ntype_md2_aer

      ! 1 indicates a multi-dimensional sectional framework
      msectional_flag2 = 0
      if (ntype_aer > 1) msectional_flag2 = 1 

      ! Number of aerosol size bins
#if defined( MOSAIC4 )
      nsize_aer_tmp = 4
#elif defined( MOSAIC8 )
      nsize_aer_tmp = 8
#elif defined( MBULK )
      nsize_aer_tmp = 1
#endif

      ! Total number of bins
      nbin_a = ntype_aer*nsize_aer_tmp

      ! Size framework:
      !   1: Modal
      !   2: Unstructured
      !   3: Sectional
#if defined( MOSAIC4 ) || defined( MOSAIC8 )
      msize_framework = msectional
#elif defined( MBULK )
      msize_Framework = munstructured
#endif
 
      ! This has to be 1
      method_bcfrac = 1

      ! This has to be either 11 or 12
      method_kappa = 11

      ! About line 125 of readinputfile.f90
      ! About line 16 of input file itself
      ! Set cloud stuff to be ignored
      ncldbin_used = 0

      ! This is basically legacy code
      nbin_a_max = nbin_a
      naerbin = nbin_a
      naerbin_used = nbin_a
  
      ! The total number of variables
      ! naer_tot is the total number of transported variables per bin
      naer_max = naer_tot*naerbin

      ! Not sure
      maxd_atype     = ntype_aer
      maxd_atype_md1 = ntype_md1_aer
      maxd_atype_md2 = ntype_md2_aer
      maxd_asize = nsize_aer_tmp
 
      ! The total length of the "cnn" array 
      ntot_used = ngas_max + naer_max + ncld_tot*ncldbin_used
      ntot_max = ntot_used
 
      ! Allocate memory
      !call mosaic_allocate_memory( 0 )
      !call mosaic_allocate_memory_inbuilt( 0 ) 
      call mosaic_allocate_memory_inbuilt(  ) 
 
      ! Allocate memory for sectional interface
      call sect_iface_allocate_memory
 
      ! Number of size bins for each aerosol type
      ! Currently assumes uniform
      nsize_aer(:) = nsize_aer_tmp

      ! If set to 2, bin sizes are calculated and the initial dry diameter is stored
      ! for each bin into CNN. The bin EDGES are always calculated for a sectional
      ! framework, but the CNN array entry is only set if this is equal to 2.
      if (msize_framework == munstructured) then
          maersize_init_flag1 = 1
      else
          maersize_init_flag1 = 2
      endif

      ! Bin limits
      if (msize_framework == msectional) then
          ! Bin upper and lower limits in um
          dlo_aersize_init = 2.d-2
          dhi_aersize_init = 3.2d-1
          !=========================================================================!
          ! Assume for the moment that we aren't using multidimensional             !
          ! sectional approach                                                      !
          !=========================================================================!
          !if (msectional_flag2 > 0) then
          !   read(lun_inp,*) method_atype_md1_init
          !   if (method_atype_md1_init <= 1) then
          !      read(lun_inp,*) xcut_atype_md1(0:ntype_md1_aer)
          !   else
          !      read(lun_inp,*) xcutlo_atype_md1_init, xcuthi_atype_md1_init
          !      tmpa = xcutlo_atype_md1_init + &
          !           (xcuthi_atype_md1_init-xcutlo_atype_md1_init)/ntype_md1_aer
          !      tmpb = max( 0.1_r8, xcutlo_atype_md1_init+0.1 )
          !      if ( (tmpa < 0.0) .or. (xcuthi_atype_md1_init < tmpb) ) then
          !         write(*,'(2a,1p,2e14.6)') &
          !              '*** readinputfile fatal error - ', &
          !              'bad xcutlo/hi_atype_md1_init', &
          !              xcutlo_atype_md1_init, xcuthi_atype_md1_init
          !         stop
          !      end if
          !   end if

          !   read(lun_inp,*) method_atype_md2_init
          !   if (method_atype_md2_init <= 1) then
          !      read(lun_inp,*) xcut_atype_md2(0:ntype_md2_aer)
          !   else
          !      read(lun_inp,*) xcutlo_atype_md2_init, xcuthi_atype_md2_init
          !      if ( (xcutlo_atype_md2_init < 1.0e-7) .or. &
          !           (xcuthi_atype_md2_init < xcutlo_atype_md2_init*1.1) ) then
          !         write(*,'(2a,1p,2e14.6)') &
          !              '*** readinputfile fatal error - ', &
          !              'bad xcutlo/hi_atype_md2_init', &
          !              xcutlo_atype_md2_init, xcuthi_atype_md2_init
          !         stop
          !      end if
          !   end if
          !end if ! (msectional_flag2 > 0)
          !=========================================================================!
      else
         dlo_aersize_init = 0.0
         dhi_aersize_init = 0.0
      end if ! (msize_framework == msectional)

      ! Set hysteresis method - use either 2 or 3 under normal circumstances
      !    1: Hysteresis with a dedicated marker (not intended for CTMs)
      !    2: Full hysteresis with HYST tracer (Expensive, untested but acc.)
      !    3: Upper hysteresis leg only (Cheaper, mostly accurate)
      !    4: Lower hysteresis leg onl (Expensive, probably unrealistic)
      mhyst_method = mhyst_uporlo_waterhyst ! 2

      ! NOW AT LINE 18 OF CASE1.INP, LINE 201 OF READINPUTFILE.F90
      ! These flags are all relevant only to sectional representation
      ! Allow aerosol coagulation? (0 = no, 1 = yes)
      mcoag_flag1 = 1
      ! Frequency of coagulation
      ! Leave this as 1 for now (looks like code doesn't follow through for != 1)
      ifreq_coag = 1
      ifreq_coag = max( 1, ifreq_coag )
      ! Allow mass transfer between bins? (0 = no, 1 = yes)
      mmovesect_flag1 = 1

      ! Nucleation parameterization
      !     0: No nucleation allowed (!)
      !     1: h2so4-nh3-h2o ternary nuc. of merikanto et al. (2007)
      !     2: h2so4-h2o binary nuc. of vehkamaki et al. (2002)
      !     3: h2so4-h2o binary nuc. of wexler et al. (1994)
      !    11: empirical first  order boundary layer nucleation
      !    12: empirical second order boundary layer nucleation
      ! This is not necessary for the bulk mode - gas/aerosol
      ! partitioning is controlled separately
      mnewnuc_flag1 = 1

      ! Not currently used
      msectional_flag1 = 0

      ! If running in bulk mode...
#if defined( MBULK )
      ifreq_coag = 0
      mmovesect_flag1 = 0
      mnewnuc_flag1 = 0
      msectional_flag1 = 0
#endif

      ! Output print control - we want no output, so...
      iprint = 0
      iwrite_gas = 0
      iwrite_aer_bin = 0
      iwrite_aer_dist = 0
      iwrite_aer_species = 0

      ! Mode - 1 is time integration, 2 is parametric
      ! We are butchering the time integration mode
      mmode = 1
      
      ! Disable the MOSAIC box-model gas-phase chemistry
      mgas = 0

      ! Enable MOSAIC aerosol chemistry
      maer = 1
 
      ! MOSAIC cloud model not yet ready
      mcld = 0

      ! MOSAIC aerosol optics - not yet ready (!)
      maeroptic = 0
      
      ! Shell/core model of aerosols
      !    0: No shell/core
      !    1: Allow BC cores
      !    2: Allow BC and DUST cores
#if defined( MBULK )
      mshellcore = 0
#else
      mshellcore = 1
#endif

      ! Solar radiation pattern - not relevant to embedded aero model
      !    1: Diurnal pattern
      !    2: Fixed
      msolar = 2

      ! Photolysis parameterization (irrelevant as we are using FAST-JX)
      !    1: Rick's parameterization
      !    2: Yang's parameterization
      mphoto = 2

      ! Allow gas-aerosol partitioning (0 is off, 1 is on)
      mGAS_AER_XFER = 1

      mDYNAMIC_SOLVER = mASTEM

      ! ASTEM tolerances
      ! Tolerance for tau (solver tolerance)
      ! Range: 0.01 - 1.0
      alpha_ASTEM = 0.5

      ! Relative equilibrium (eqb) tolerance
      ! Range: 0.01 - 0.03
      rtol_eqb_ASTEM = 0.01

      ! Percentage mol tolerange
      ! Range: 0.01 - 1.0
      ptol_mol_ASTEM = 0.01

      ! PMCMOS - this is not relevant here
      ipmcmos = 0

      ! PMCMOS parameters (ignore)
      istate_pblh = 0
      msolar_pmcmos_dtshift = 0
      
      if (msolar >= 1000) then
         if (ipmcmos > 0) msolar_pmcmos_dtshift = 1
         msolar = mod( msolar, 1000 )
      end if

      ! Set up the CNN array
      ! Only need gas-phase SO4, HNO3, NH3, HCl and MSA
      ! SOAs will complicate this
!      cnn(:) = 0.0_r8
      emission(:) = 0.0_r8

      species(kh2so4) = 'H2SO4'
      species(khno3) = 'HNO3'
      species(khcl) = 'HCl'
      species(knh3) = 'NH3'
      species(kmsa) = 'MSA'
      
      ! Aerosol data - ignore (set online)

      ! Just in case
      te_old = te
      pr_atm_old = pr_atm
      rh_old = rh
      pblh = 1.0e3   ! default pblh (m)
      pblh_old = pblh

      ! These variables are duplicated in module_data_mosaic_main and ..._aero
      ipmcmos_aero = ipmcmos
      maeroptic_aero = maeroptic

      ! Initialize the aerosol module
      ! This sets the following values in CNN:
      !    Aerosol dry diameters
      !    Aerosol geometric standard deviations
      !    Aerosol number densities
      CALL INIT_AEROSOL_INBUILT

      END SUBROUTINE INIT_MOSAIC
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: INIT_DATA_MODULES_INBUILT
!
! !DESCRIPTION: Subroutine INIT\_DATA\_MODULES\_INBUILT initializes MOSAIC
!  arrays and variables that would usually be initialized internally.
!\\
!\\
! !INTERFACE:
!
      subroutine init_data_modules_inbuilt
!
!   place various constant or initial values into common
      use module_data_mosaic_main
      use module_data_mosaic_gas
      use module_data_mosaic_aero
      use module_data_mosaic_cloud

      implicit none

!
! !INPUT PARAMETERS:
!
!
! !INPUT/OUTPUT PARAMETERS: 
!
!
! !OUTPUT PARAMETERS:
!
!
! !REMARKS:
!  Original v4.1a MOSAIC implementation into GEOS-Chem by
!  Sebastian D. Eastham and Jessica Kunke (SDE, JK, 2014-11-13)
! 
! !REVISION HISTORY: 
!  29 Nov 2014 - S. D. Eastham - Initial version. Code taken directly
!                                from init_data_modules.f90 in the
!                                MOSAIC main function directory
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

      !=================================================================
      ! INIT_MOSAIC  begins here!
      !=================================================================

!--------------------------------------------------
!     define fundamental constants...
       avogad	= 6.02217e+23

!      pi	= 3.141592654
!      deg2rad	= 0.017453293
       pi	= 4.0_r8 * atan( 1.0_r8 )
       piover4	= pi/4.0_r8
       piover6	= pi/6.0_r8
       deg2rad	= pi/180.0_r8
       third	= 1.0_r8/3.0_r8

!---------------------------------
! define species indices
!
! species in inorganic chemistry
       kh2so4		=  1
       khno3		=  2
       khcl		=  3
       knh3		=  4
       kno		=  5
       kno2		=  6
       kno3		=  7
       kn2o5		=  8
       khono		=  9
       khno4		= 10
       ko3		= 11
       ko1d		= 12
       ko3p		= 13
       koh		= 14
       kho2		= 15
       kh2o2		= 16
       kco		= 17
       kso2		= 18
!
! species in methane, ethane, formaldehyde chemistry
       kch4		= 19
       kc2h6		= 20
       kch3o2		= 21
       kethp		= 22
       khcho		= 23
       kch3oh		= 24
       kanol		= 25
       kch3ooh		= 26
       kethooh		= 27
       kald2		= 28
       khcooh		= 29
       krcooh		= 30
       kc2o3		= 31
       kpan		= 32
       karo1		= 33	! soa prec 1
       karo2		= 34	! soa prec 2
       kalk1		= 35	! soa prec 3
       kole1		= 36	! soa prec 4
       kapi1		= 37	! soa prec 5
       kapi2		= 38	! soa prec 6
       klim1		= 39	! soa prec 7
       klim2		= 40	! soa prec 8
!
! species in hc1 mechanism. initialize indices to zero
       kpar		= 41
       kaone		= 42
       kmgly		= 43
       keth		= 44
       kolet		= 45
       kolei		= 46
       ktol		= 47
       kxyl		= 48
       kcres		= 49
       kto2		= 50
       kcro		= 51
       kopen 		= 52
       konit		= 53
       krooh		= 54
       kro2		= 55
       kano2		= 56
       knap		= 57
       kxo2		= 58
       kxpar		= 59
!
! species in hc2 mechanism. initialize indices to zero
       kisop		= 60
       kisoprd		= 61
       kisopp		= 62
       kisopn		= 63
       kisopo2		= 64
       kapi		= 65
       klim		= 66

! species in dms mechanism. initialize indices to zero
       kdms		= 67
       kmsa		= 68
       kdmso		= 69
       kdmso2		= 70
       kch3so2h		= 71
       kch3sch2oo	= 72
       kch3so2		= 73
       kch3so3		= 74
       kch3so2oo	= 75
       kch3so2ch2oo	= 76
       ksulfhox		= 77

! aerosol bin: (24 species=bin)
       knum_a		=  1
       kdpdry_a		=  2
       ksigmag_a	=  3
       kjhyst_a		=  4
       kwater_a		=  5
       kso4_a		=  6
       kno3_a		=  7
       kcl_a		=  8
       knh4_a		=  9
       kmsa_a		= 10
       karo1_a		= 11
       karo2_a		= 12
       kalk1_a		= 13
       kole1_a		= 14
       kapi1_a		= 15
       kapi2_a		= 16
       klim1_a		= 17
       klim2_a		= 18
       kco3_a		= 19
       kna_a		= 20
       kca_a		= 21
       koin_a		= 22
       koc_a		= 23
       kbc_a		= 24


! cloud mode: 1-ncldmode (13=mode)
       knum_c		=  1
       kwater_c		=  2
       kso4_c		=  3
       kno3_c		=  4
       kcl_c		=  5
       knh4_c		=  6
       koc_c		=  7
       kmsa_c		=  8
       kco3_c		=  9
       kna_c		= 10
       kca_c		= 11
       kbc_c		= 12
       koin_c		= 13
!
!
! regime-dependent chemistry definitions
!
       iregime		=  1
!
!     GAS
!
       ih2so4		=  1
       ihno3		=  2
       ihcl		=  3
       inh3		=  4
       ino		=  5
       ino2		=  6
       ino3		=  7
       in2o5		=  8
       ihono		=  9
       ihno4		= 10
       io3		= 11
       io1d		= 12
       io3p		= 13
       ioh		= 14
       iho2		= 15
       ih2o2		= 16
       ico		= 17
       iso2		= 18
!
! species in methane, ethane, formaldehyde chemistry
       ich4		= 19
       ic2h6		= 20
       ich3o2		= 21
       iethp		= 22
       ihcho		= 23
       ich3oh		= 24
       ianol		= 25
       ich3ooh		= 26
       iethooh		= 27
       iald2		= 28
       ihcooh		= 29
       ircooh		= 30
       ic2o3		= 31
       ipan		= 32
       iaro1		= 33	! soa prec 1
       iaro2		= 34	! soa prec 2
       ialk1		= 35	! soa prec 3
       iole1		= 36	! soa prec 4
       iapi1		= 37	! soa prec 5
       iapi2		= 38	! soa prec 6
       ilim1		= 39	! soa prec 7
       ilim2		= 40	! soa prec 8
!
! species in hc1 mechanism. initialize indices to zero
       ipar		= 41
       iaone		= 42
       imgly		= 43
       ieth		= 44
       iolet		= 45
       iolei		= 46
       itol		= 47
       ixyl		= 48
       icres		= 49
       ito2		= 50
       icro		= 51
       iopen		= 52
       ionit		= 53
       irooh		= 54
       iro2		= 55
       iano2		= 56
       inap		= 57
       ixo2		= 58
       ixpar		= 59
!
! species in hc2 mechanism. initialize indices to zero
       iisop		= 60
       iisoprd		= 61
       iisopp		= 62
       iisopn		= 63
       iisopo2		= 64
       iapi		= 65
       ilim		= 66

! species in dms mechanism. initialize indices to zero
       idms		= 67
       imsa		= 68
       idmso		= 69
       idmso2		= 70
       ich3so2h		= 71
       ich3sch2oo	= 72
       ich3so2		= 73
       ich3so3		= 74
       ich3so2oo	= 75
       ich3so2ch2oo	= 76
       isulfhox		= 77

! alkylperoxy radical indices for parameterized permutation reactions
       jch3o2		=  1
       jethp		=  2
       jro2		=  3
       jc2o3		=  4
       jano2		=  5
       jnap		=  6
       jisopp		=  7
       jisopn		=  8
       jisopo2		=  9
       jxo2		= 10

! photolyzing species indices
       jphoto_no2	=  1
       jphoto_no3	=  2
       jphoto_hono	=  3
       jphoto_hno3	=  4
       jphoto_hno4	=  5
       jphoto_n2o5	=  6
       jphoto_o3a	=  7
       jphoto_o3b	=  8
       jphoto_h2o2	=  9
       jphoto_hchoa	= 10
       jphoto_hchob	= 11
       jphoto_ch3ooh	= 12
       jphoto_ethooh	= 13
       jphoto_ald2	= 14
       jphoto_aone	= 15
       jphoto_mgly	= 16
       jphoto_open	= 17
       jphoto_rooh	= 18
       jphoto_onit	= 19
       jphoto_isoprd	= 20


!
!     CLOUD
!
! cloud (local): used for total and undissociated species
       iso4_c		=  1
       ino3_c		=  2
       icl_c		=  3
       inh4_c		=  4
       ioc_c		=  5
       imsa_c		=  6
       ico2_c		=  7
       ina_c		=  8
       ica_c		=  9
       ibc_c		= 10
       ioin_c		= 11
       iso2_c		= 12
       ihono_c		= 13
       ih2o2_c		= 14
       ich3ooh_c	= 15
       ihcooh_c		= 16
       ircooh_c		= 17
       ihcho_c		= 18
       io3_c		= 19
       iho2_c		= 20	! don't actually exist in cloud (surf rxn only)
       ino2_c		= 21	! don't actually exist in cloud (surf rxn only)
       ino3r_c		= 22	! don't actually exist in cloud (surf rxn only)
       in2o5_c		= 23	! don't actually exist in cloud (surf rxn only)

! cloud ionic species
       jh_c		=  1
       jnh4_c		=  2
       jna_c		=  3
       jhso4_c		=  4
       jso4_c		=  5
       jno3_c		=  6
       jcl_c		=  7
       jno2_c		=  8
       jhso3_c		=  9
       jso3_c		= 10
       jhco3_c		= 11
       jco3_c		= 12
       jho2_c		= 13
       jhcoo_c		= 14
       jrcoo_c		= 15
       jmsa_c		= 16
       joh_c		= 17
       jch2oh2_c	= 18


      return
      end subroutine init_data_modules_inbuilt
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: INIT_AEROSOL_INBUILT
!
! !DESCRIPTION: Subroutine INIT\_AEROSOL\_INBUILT initializes aerosol
!  arrays and variables that would usually be initialized internally.
!\\
!\\
! !INTERFACE:
!
      subroutine init_aerosol_inbuilt
        !
        !   does initialization of variables in module_data_mosaic_asect
        !
      use module_data_mosaic_kind, only: r8
      use module_data_mosaic_main, only: 
     &    cnn, ipmcmos, 
     &    kdpdry_a, kjhyst_a, knum_a, ksigmag_a, kwater_a, 
     &    lun_sect_171, 
     &    m_partmc_mosaic, naer_tot, ngas_max, ntot_max, 
     &    pi, species
      use module_data_mosaic_aero, only: 
     &    aer_name, 
     &    dhi_aersize_init, dlo_aersize_init, 
     &    d_mdrh_dim2, jsulf_poor_num, jsulf_rich_num, mdrh_t_num, 
     &    iso4_a,     ino3_a,     icl_a,     inh4_a,     ico3_a,  
     &    imsa_a,     ina_a,      ica_a,     ioc_a,      ibc_a,   
     &    ioin_a,     iaro1_a,    iaro2_a,   ialk1_a,    iole1_a,   
     &    iapi1_a,    iapi2_a,    ilim1_a,   ilim2_a, 
     &    maersize_init_flag1, mhyst_method, mmodal, 
     &    method_atype_md1_init, method_atype_md2_init, 
     &    msectional, msectional_flag2, msize_framework, munstructured, 
     &    nbin_a, 
     &    naer, dens_aer_mac, mw_aer_mac, kappa_aer_mac, 
     &    xcuthi_atype_md1_init, xcutlo_atype_md1_init, 
     &    xcuthi_atype_md2_init, xcutlo_atype_md2_init
      use module_data_mosaic_asect, only: 
     &    ai_phase, 
     &    dcen_sect, dcut_sect, dlo_sect, dhi_sect, 
     &    dens_aer, hygro_aer, mw_aer, name_aer, sigmag_aer, 
     &    dens_mastercomp_aer, hygro_mastercomp_aer, 
     &    mw_mastercomp_aer, name_mastercomp_aer, 
     &    ibin_of_isize_itype, isize_of_ibin, itype_of_ibin, 
     &    itype_md1_of_itype, itype_md2_of_itype, itype_of_itype_md1md2,
     &    lptr_so4_aer, lptr_msa_aer, lptr_no3_aer, lptr_cl_aer, 
     &    lptr_co3_aer, lptr_nh4_aer, lptr_na_aer, lptr_ca_aer, 
     &    lptr_oin_aer, lptr_oc_aer, lptr_bc_aer, 
     &    hyswptr_aer, massptr_aer, mastercompptr_aer, numptr_aer,
     &    waterptr_aer, 
     &    mastercompindx_so4_aer, mastercompindx_no3_aer, 
     &    mastercompindx_cl_aer, 
     &    mastercompindx_msa_aer, mastercompindx_co3_aer, 
     &    mastercompindx_nh4_aer, 
     &    mastercompindx_na_aer, mastercompindx_ca_aer, 
     &    mastercompindx_oin_aer, 
     &    mastercompindx_oc_aer, mastercompindx_bc_aer, 
     &    maxd_acomp, maxd_asize, maxd_atype, 
     &    ncomp_aer, ncomp_plustracer_aer, nphase_aer, nsize_aer, 
     &    ntot_mastercomp_aer, ntype_aer, ntype_md1_aer, ntype_md2_aer, 
     &    volumcen_sect, volumcut_sect, volumhi_sect, volumlo_sect, 
     &    xcut_atype_md1, xcut_atype_md2
      !qak
      !qak
      use module_pmcmos_aer, only: pmcmos_init_aerosol
      use module_mosaic_init, only: mosaic_init
      
      implicit none
      
      !Local Variables
      integer :: ibin, iphase, isize, itype, it1, it2
      integer ::  l, ll, lunaa, noffset
      real(r8) :: tmpa, tmpb, tmpc, tmpn
      character(len=5) :: tmpch5
      
      !BSINGH - load_mosaic_parameters is now called through mosaic_init
      !         I commented out the following call
      
      !call load_mosaic_parameters(nmax_ASTEM,b_mtem,zc,za,b_zsr,a_zsr,aw_min,        &
      !     mw_electrolyte,dens_electrolyte,partial_molar_vol,MW_c,MW_a,mw_aer_mac,   &
      !     dens_aer_mac,kappa_aer_mac,dens_comp_a,mw_comp_a,ref_index_a,rtol_mesa,   &
      !     jsalt_index,jsulf_poor,jsulf_rich,Nmax_mesa,d_mdrh)
      
      call mosaic_init 
      
      ! check for valid msize_framework
      if ((msize_framework /= mmodal        ) .and. 
     &     (msize_framework /= munstructured ) .and. 
     &     (msize_framework /= msectional    )) then
         write(*,'(2a,2(1x,i10))') 
     &        '*** init_aerosol fatal error - ', 
     &        'bad msize_framework = ', msize_framework
         stop
      end if
      
      ! check misc flags
      if ((mhyst_method < 1) .or. (mhyst_method > 4)) then !BSINGH - 05/28/2013(RCE updates)
         write(*,'(2a,2(1x,i10))') 
     &        '*** init_aerosol fatal error - ', 
     &        'bad mhyst_method = ', mhyst_method
         stop
      end if
      
      
      !BSINGH - 05/28/2013(RCE updates)
      ! check naer and naer_tot
      if ((naer /= 19) .or. (naer_tot /= 24)) then
         write(*,'(2a,2(1x,i10))') 
     &     '*** init_aerosol fatal error - ', 
     &     'expecting naer=19 & naer_tot=24 but they = ', naer, naer_tot
         stop
      end if
      
      
      ! set names for cnn array
      do l = ngas_max+1, ntot_max
         write( species(l), '(a,i9.9)' ) 's', l
      end do
      do ibin = 1, nbin_a
         noffset = ngas_max + naer_tot*(ibin - 1)
         write(tmpch5,'(a,i3.3)') '_a', mod(ibin,1000)
         species(noffset+knum_a   ) = 'num'    // tmpch5
         species(noffset+kdpdry_a ) = 'dpdry'  // tmpch5
         species(noffset+ksigmag_a) = 'sigmag' // tmpch5
         species(noffset+kjhyst_a ) = 'jhyst'  // tmpch5
         species(noffset+kwater_a ) = 'water'  // tmpch5
         do l = 1, naer
            species(noffset+5+l) = trim(aer_name(l)) // tmpch5
         end do
      end do
      
      
      !
      ! following is primarily for sectional framework,
      !    but some of it may be useful for mosaic box-model unstructured
      ! it uses the module_data_mosaic_asect arrays,
      !    so it does not work for partmc_moscic
      !
      if (m_partmc_mosaic > 0) return
      !BSINGH - 05/28/2013(RCE updates ENDS)
      
      
      
      
      
      
      
      ! check for valid ntype_aer and nsize_aer(:)
      if ((ntype_aer < 1) .or. (ntype_aer > maxd_atype)) then
         write(*,'(2a,2(1x,i10))') 
     &        '*** init_aerosol fatal error - ', 
     &        'too big/small ntype_aer = ', ntype_aer
         stop
      end if
      
      do itype = 1, ntype_aer
         if ((nsize_aer(itype) < 1) .or. 
     &        (nsize_aer(itype) > maxd_asize)) then
            write(*,'(2a,2(1x,i10))') 
     &           '*** init_aerosol fatal error - ', 
     &           'too big/small nsize_aer = ', nsize_aer(itype), itype
            stop
         end if
         if (nsize_aer(itype) /= nsize_aer(1)) then
            write(*,'(2a,2(1x,i10))') 
     &           '*** init_aerosol fatal error - ', 
     &           'non-constant nsize_aer = ', nsize_aer(itype), itype
            stop
         end if
      end do
      
      ! phases -- only 1 for now
      nphase_aer = 1
      ai_phase = 1
      
      
      ! species naer vs. maxd_acomp !BSINGH - 05/28/2013(RCE updates -Deleted an 'if' block here)
      if (naer > maxd_acomp) then
         write(*,'(2a,2(1x,i10))') 
     &        '*** init_aerosol fatal error - ', 
     &        'naer > maxd_acomp = ', naer, maxd_acomp
         stop
      end if
      
      ! set species pointers and properties !BSINGH - 05/28/2013(RCE updates - Added this comment)
      dens_mastercomp_aer(:) = 0.0
      mw_mastercomp_aer(:) = 0.0
      hygro_mastercomp_aer(:) = 0.0
      
      ntot_mastercomp_aer = naer
      do ll = 1, naer
         dens_mastercomp_aer(ll) = dens_aer_mac(ll)
         mw_mastercomp_aer(  ll) = mw_aer_mac(ll)
         name_mastercomp_aer(ll) = aer_name(ll)
      end do
      mastercompindx_so4_aer = iso4_a
      mastercompindx_no3_aer = ino3_a
      mastercompindx_cl_aer  = icl_a
      mastercompindx_msa_aer = imsa_a
      mastercompindx_co3_aer = ico3_a
      mastercompindx_nh4_aer = inh4_a
      mastercompindx_na_aer  = ina_a
      mastercompindx_ca_aer  = ica_a
      mastercompindx_oin_aer = ioin_a
      mastercompindx_oc_aer  = ioc_a
      mastercompindx_bc_aer  = ibc_a
      
      ! set hygroscopicities that are (temporarily) used
      !    by sectional routines for hysteresis water
      ! for consistency with the mosaic water uptake,
      !    only the inorganic salts/acids have non-negligible values
      hygro_mastercomp_aer(:) = 1.0e-10
      hygro_mastercomp_aer( mastercompindx_so4_aer ) = 0.5
      hygro_mastercomp_aer( mastercompindx_nh4_aer ) = 0.5
      hygro_mastercomp_aer( mastercompindx_no3_aer ) = 0.5
      hygro_mastercomp_aer( mastercompindx_msa_aer ) = 0.58
      hygro_mastercomp_aer( mastercompindx_na_aer  ) = 1.16
      hygro_mastercomp_aer( mastercompindx_cl_aer  ) = 1.16
      hygro_mastercomp_aer( mastercompindx_ca_aer  ) = 0.1
      hygro_mastercomp_aer( mastercompindx_co3_aer ) = 0.1
      
      ! 2010-02-26 - get hygro from mosaic kappa, always
      ! 2010-03-28 - only do this when ipmcmos > 0,
      !    to allow backwards compatibility with mos27a
      if (ipmcmos > 0) then
         do ll = 1, naer
            hygro_mastercomp_aer(ll) = kappa_aer_mac(ll)
         end do
      end if
      
      
      ibin_of_isize_itype(:,:) = -999888777
      isize_of_ibin(:) = -999888777
      itype_of_ibin(:) = -999888777
      itype_of_itype_md1md2(:,:) = -999888777
      itype_md1_of_itype(:) = -999888777
      itype_md2_of_itype(:) = -999888777
      
      ncomp_aer(:) = 0
      ncomp_plustracer_aer(:) = 0
      numptr_aer(:,:,:) = -999888777
      waterptr_aer(:,:) = -999888777
      hyswptr_aer(:,:) = -999888777
      massptr_aer(:,:,:,:) = -999888777
      mastercompptr_aer(:,:) = -999888777
      
      dens_aer(:,:) = 0.0
      mw_aer(:,:) = 0.0
      hygro_aer(:,:) = 0.0
      
      iphase = ai_phase
      ibin = 0
      do itype = 1, ntype_aer
         ncomp_aer(itype) = naer
         ncomp_plustracer_aer(itype) = naer
         do isize = 1, nsize_aer(itype)
            ibin = ibin + 1
            noffset = ngas_max + naer_tot*(ibin - 1)
            numptr_aer(isize,itype,iphase) = noffset + knum_a
            waterptr_aer(isize,itype) = noffset + kwater_a
            hyswptr_aer(isize,itype) = noffset + kjhyst_a
            do ll = 1, ncomp_aer(itype)
               massptr_aer(ll,isize,itype,iphase) = 
     &              noffset + (naer_tot-naer) + ll
               mastercompptr_aer(ll,itype) = ll
               name_aer( ll,itype) = name_mastercomp_aer(ll)
               dens_aer( ll,itype) = dens_mastercomp_aer(ll)
               mw_aer(   ll,itype) = mw_mastercomp_aer(ll)
               hygro_aer(ll,itype) = hygro_mastercomp_aer(ll)
            end do
            !BSINGH - 05/28/2013(RCE updates)
            lptr_so4_aer(isize,itype,iphase) = 
     &                        massptr_aer(iso4_a,isize,itype,iphase)
            lptr_msa_aer(isize,itype,iphase) = 
     &                        massptr_aer(imsa_a,isize,itype,iphase)
            lptr_no3_aer(isize,itype,iphase) = 
     &                        massptr_aer(ino3_a,isize,itype,iphase)
            lptr_cl_aer( isize,itype,iphase) = 
     &                        massptr_aer(icl_a ,isize,itype,iphase)
            lptr_co3_aer(isize,itype,iphase) = 
     &                        massptr_aer(ico3_a,isize,itype,iphase)
            lptr_nh4_aer(isize,itype,iphase) = 
     &                        massptr_aer(inh4_a,isize,itype,iphase)
            lptr_na_aer( isize,itype,iphase) = 
     &                        massptr_aer(ina_a ,isize,itype,iphase)
            lptr_ca_aer( isize,itype,iphase) = 
     &                        massptr_aer(ica_a ,isize,itype,iphase)
            lptr_oin_aer(isize,itype,iphase) = 
     &                        massptr_aer(ioin_a,isize,itype,iphase)
            lptr_oc_aer( isize,itype,iphase) = 
     &                        massptr_aer(ioc_a ,isize,itype,iphase)
            lptr_bc_aer( isize,itype,iphase) = 
     &                        massptr_aer(ibc_a ,isize,itype,iphase)
            !BSINGH - 05/28/2013(RCE updates ENDS)
      
      
            ibin_of_isize_itype(isize,itype) = ibin
            isize_of_ibin(ibin) = isize
            itype_of_ibin(ibin) = itype
         end do ! isize
      end do ! itype
      
      
      ! generate initial dry diameters automatically
      dlo_sect( :,:) = 0.0
      dhi_sect( :,:) = 0.0
      dcen_sect(:,:) = 0.0
      dcut_sect(:,:) = 0.0
      volumlo_sect( :,:) = 0.0
      volumhi_sect( :,:) = 0.0
      volumcen_sect(:,:) = 0.0
      volumcut_sect(:,:) = 0.0
      
      if ((maersize_init_flag1 == 2) .or. 
     &      (msize_framework == msectional)) then
         if ((dlo_aersize_init <= 0.0) .or. 
     &         (dlo_aersize_init >= dhi_aersize_init)) then
            write(*,'(2a,1p,2e11.3)') 
     &            '*** init_aerosol fatal error - ', 
     &            'bad dlo/hi_aersize_init =', 
     &            dlo_aersize_init, dhi_aersize_init
            stop
         end if
      
         tmpa = log(dlo_aersize_init)
         tmpb = (log(dhi_aersize_init) - tmpa)/nsize_aer(1)
         ibin = 0
         do itype = 1, ntype_aer
            do isize = 1, nsize_aer(itype)
               ! dlo_sect & dhi_sect are bin lower/upper bound dry diam (cm)
               dlo_sect( isize,itype) = exp(tmpa + (isize-1  )*tmpb)*1.0e-4_r8
               dhi_sect( isize,itype) = exp(tmpa + (isize    )*tmpb)*1.0e-4_r8
               dcen_sect(isize,itype) = exp(tmpa + (isize-0.5)*tmpb)*1.0e-4_r8
               dcut_sect(isize,itype) = dhi_sect(isize,itype)
      
               volumlo_sect( isize,itype) = (dlo_sect( isize,itype)**3)
     &                                               *pi/6.0_r8
               volumhi_sect( isize,itype) = (dhi_sect( isize,itype)**3)
     &                                               *pi/6.0_r8
               volumcen_sect(isize,itype) = (dcen_sect(isize,itype)**3)
     &                                               *pi/6.0_r8
               volumcut_sect(isize,itype) = volumhi_sect(isize,itype)
      
               ibin = ibin + 1
               ! SDe 2015-01-05: Needed to move this into OMP clause
!               if (maersize_init_flag1 == 2) then
!                  noffset = ngas_max + naer_tot*(ibin - 1)
!                  ! cnn(noffset + kdpdry_a) is bin initial dry diameter (micron)
!                  cnn(noffset + kdpdry_a) = exp(tmpa + (isize-0.5)*tmpb)
!               end if
            end do ! isize
      
            dcut_sect(0,itype) = dlo_sect(1,itype)
            volumcut_sect(0,itype) = volumlo_sect(1,itype)
         end do ! itype
      
      end if
      
      
      ! sigmag
      do itype = 1, ntype_aer
         do isize = 1, nsize_aer(itype)
            ibin = ibin_of_isize_itype(isize,itype)
            noffset = ngas_max + naer_tot*(ibin - 1)
            !sigmag_aer(isize,itype) = cnn(noffset + ksigmag_a)
            if (msize_framework == munstructured ) then
               !           assume monodisperse for each "bin"
               sigmag_aer(isize,itype) = 1.0
            else if (msize_framework == msectional) then
               !           sigmag based on width of bin - following is same as
               !           log(sigmag) = [log(dhi) - log(dlo)]/sqrt(12)
               sigmag_aer(isize,itype) = 
     &              (dhi_sect(isize,itype)/dlo_sect(isize,itype))**0.289
            end if
            ! SDE 2015-01-05: OMP prevents use of CNN outside loop
            !cnn(noffset + ksigmag_a) = sigmag_aer(isize,itype)
         end do ! isize
      end do ! itype
      
     
      !===================================================================================== 
      ! SDE 2015-01-05: This can (?) be dispensed with in the GEOS-Chem framework
      !===================================================================================== 
!      ! number
!      !    for modal or unstructured framework, the input file value should be used !BSINGH - 05/28/2013(RCE updates - Modified comment)
!      !    for sectional, initialize it correctly on the first time step
!      !*** this seems to work, but code should really be changed so that
!      !       aerosol number is initialized BEFORE the first time step
!      if (msize_framework == msectional) then !BSINGH - 05/28/2013(RCE updates)
!         do itype = 1, ntype_aer
!            do isize = 1, nsize_aer(itype)
!      !        cnn(numptr_aer(isize,itype,iphase)) = 0.0
!      ! 29-sep-2014 - only set number to zero here if
!      !     its value is negative, or
!      !     its value is below that obtained from dry-volume and volumhi_sect, or
!      !     its value is above that obtained from dry-volume and volumlo_sect
!               if (cnn(numptr_aer(isize,itype,iphase)) <= 0.0) then
!                  cnn(numptr_aer(isize,itype,iphase)) = 0.0
!               else
!                  tmpa = 0.0
!                  do ll = 1, ncomp_aer(itype)
!                     tmpa = tmpa 
!     &           + max(0.0_r8,cnn(massptr_aer(ll,isize,itype,iphase))) 
!     &           * 1.0e-12*mw_aer(ll,itype)/dens_aer(ll,itype)
!                     ! the 1.0e-12 factor converts cnn from umol/m3 to mol/cm3
!                     ! then tmpa is cm3-dry-aerosol/cm3-air
!                  end do
!                  tmpn = cnn(numptr_aer(isize,itype,iphase))
!                  tmpb = tmpa/volumlo_sect(isize,itype)
!                  tmpc = tmpa/volumhi_sect(isize,itype)
!                  if ( tmpn < 0.999*tmpc .or. tmpn > 1.001*tmpb ) then
!                     cnn(numptr_aer(isize,itype,iphase)) = 0.0
!                  end if
!               end if
!            end do ! isize
!         end do ! itype
!      endif!BSINGH - 05/28/2013(RCE updates)
!      
!      !BSINGH - 05/28/2013(RCE updates- Got rid of two do constructs here)
      !===================================================================================== 
      ! SDE 2015-01-05
      !===================================================================================== 
      
      
      
      ! 3d sectional stuff
      if (msectional_flag2 > 0) then
      
         ! mappings between itype and (itype_md1,itype_md2)
         itype = 0
         do it2 = 1, ntype_md2_aer
            do it1 = 1, ntype_md1_aer
               itype = itype + 1
               itype_of_itype_md1md2(it1,it2) = itype
               itype_md1_of_itype(itype) = it1
               itype_md2_of_itype(itype) = it2
            end do
         end do
      
         ! cut values for bc mass fraction
         ! when method_atype_md1_init <= 1, they are read from input file
         ! otherwise, they have a uniform linear spacing between the
         !    specified xcutlo & xcuthi values
         if (method_atype_md1_init >= 2) then
            tmpa = xcutlo_atype_md1_init
            tmpb = (xcuthi_atype_md1_init - tmpa) / ntype_md1_aer
            do it1 = 0, ntype_md1_aer
               xcut_atype_md1(it1) = tmpa + tmpb*it1
            end do
         else
            tmpb = 1.0
         end if
         if ( (tmpb <= 0.0) .or. 
     &         (xcut_atype_md1(1) < 0.0) .or. 
     &         (xcut_atype_md1(ntype_md1_aer) < 0.1) ) then
            write(*,'(2a)') 
     &            '*** init_aerosol fatal error - ', 
     &            'bad xcut_atype_md1 ='
            write(*,'(1p,7e11.3)') 
     &           xcut_atype_md1(0:ntype_md1_aer)
            stop
         end if
      
         ! cut values for hygroscopicity
         ! when method_atype_md2_init <= 1, they are read from input file
         ! otherwise, they have a uniform logarithmic spacing between the
         !    specified xcutlo & xcuthi values
         if (method_atype_md2_init >= 2) then
            tmpa = log(xcutlo_atype_md2_init)
            tmpb = (log(xcuthi_atype_md2_init) - tmpa) / ntype_md2_aer
            do it2 = 0, ntype_md2_aer
               xcut_atype_md2(it2) = exp( tmpa + tmpb*it2 )
            end do
         else
            tmpb = 1.0
         end if
         if ( (tmpb <= 0.0) .or. 
     &         (xcut_atype_md2(0) <= 0.0) .or. 
     &         (xcut_atype_md2(1) < 0.0001) .or. 
     &         (xcut_atype_md2(ntype_md2_aer) < 0.1) ) then
            write(*,'(2a)') 
     &            '*** init_aerosol fatal error - ', 
     &            'bad xcut_atype_md2 ='
            write(*,'(1p,7e11.3)') 
     &            xcut_atype_md2(0:ntype_md2_aer)
            stop
         end if
      
      else
         xcut_atype_md1(0:1) = (/ 0.0_r8, 1.0_r8 /)
         xcut_atype_md2(0:1) = (/ 0.1_r8, 1.0_r8 /)
      
      end if ! (msectional_flag2 > 0)
      
      
      ! diagnostic output
      lunaa = lun_sect_171
      if (lunaa > 0) then
      
         write(lunaa,'(//a,i5)') 'ntot_mastercomp_aer', 
     &                            ntot_mastercomp_aer
         do ll = 1, ntot_mastercomp_aer
            write(lunaa,'(a,i4,2x,a,1p,4e12.4)') 'll,name,dens,mw,hy', 
     &            ll, name_mastercomp_aer(ll), dens_mastercomp_aer(ll), 
     &            mw_mastercomp_aer(ll), hygro_mastercomp_aer(ll)
         end do
      
         iphase = ai_phase
         do itype = 1, ntype_aer
            write(lunaa,'(//a,2i5)') 'itype, ncomp', itype, 
     &                                     ncomp_aer(itype)
            do ll = 1, ncomp_aer(itype)
               write(lunaa,'(a,i4,2x,a,1p,4e12.4)') 'll,name,dens,mw', 
     &               ll, name_aer(ll,itype), dens_aer(ll,itype), 
     &               mw_aer(ll,itype)
            end do
            do isize = 1, nsize_aer(itype)
               write(lunaa,'(/a,5i5)')    'itype,isize,numptr,waterptr',
     &               itype, isize, numptr_aer(isize,itype,iphase), 
     &               waterptr_aer(isize,itype)
               do ll = 1, ncomp_aer(itype)
                  write(lunaa,'(a,3i5,i7,2x,a)') 
     &             'itype,isize,ll,massptr,name', 
     &            itype, isize, ll, massptr_aer(ll,isize,itype,iphase), 
     &                  name_aer(ll,itype)
               end do
            end do
         end do
      
      end if   ! (lunaa > 0)
      
      
      ! pmcmos stuff
      if (ipmcmos > 0) call pmcmos_init_aerosol(mw_aer_mac,
     &                              dens_aer_mac,kappa_aer_mac)
      
      
      return
      end subroutine init_aerosol_inbuilt
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: mosaic_allocate_memory_inbuilt
!
! !DESCRIPTION: Subroutine MOSAIC\_ALLOCATE\_MEMORY\_INBUILT is copied
!  directly from the original MOSAIC code of the same name.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE mosaic_allocate_memory_inbuilt( initialize_flag )
!
! !USES:
!
      use module_data_mosaic_main
      use module_data_mosaic_aero
      use module_data_mosaic_asect
!
! !INPUT PARAMETERS:
!
      integer, intent(in), optional :: initialize_flag
!
! !REVISION HISTORY: 
!  28 Dec 2014 - S. D. Eastham - Initial version
!  05 Jan 2015 - S. D. Eastham - Modified for OpenMP
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

      !=================================================================
      ! MOSAIC_ALLOCATE_MEMORY_INBUILT begins here!
      !=================================================================

! allocate arrays whose dimensions are set at run time
!
! initialize them to zero when initialize_flag > 0 (or is not present)

      allocate( lun_aer(naerbin) )
      allocate( lun_aer_status(naerbin) )

      allocate( aer_output(naerbin) )

      allocate( species(ntot_max) )
 
      ! SDE 2015-01-05: Needs to be declared on each machine separately
      !allocate( cnn(ntot_max) )

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
         allocate( massptr_aer( maxd_acomp, maxd_asize, maxd_atype, 
     &                                                 maxd_aphase ) )
         
         allocate( waterptr_aer( maxd_asize, maxd_atype ) )
         allocate( hyswptr_aer( maxd_asize, maxd_atype ) )
         allocate( numptr_aer( maxd_asize, maxd_atype, maxd_aphase ) )
         
         allocate( mprognum_aer(maxd_asize,maxd_atype,maxd_aphase) )
         
         allocate( ibin_of_isize_itype( maxd_asize, maxd_atype ) )
         allocate( isize_of_ibin( nbin_a_max ) )
         allocate( itype_of_ibin( nbin_a_max ) )
         allocate( itype_of_itype_md1md2(maxd_atype_md1,maxd_atype_md2))
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
      end subroutine mosaic_allocate_memory_inbuilt
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: calc_h2so4_gas
!
! !DESCRIPTION: Subroutine SO4\_PARTITION calculates the fraction of
!  SO4 mass which can be considered to be gaseous H2SO4.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SO4_PARTITION( SO4KG,AIRMASS,PCENTER,TCENTER,
     &                                         SO4AER,SO4GAS )
!
! !USES:
!
      USE CMN_GCTM_MOD,       ONLY : AIRMW
      USE ERROR_MOD,          ONLY : IS_SAFE_DIV
!
! !INPUT PARAMETERS:
!
      REAL*8,INTENT(IN)  :: SO4KG   ! Total SO4 (kg)
      REAL*8,INTENT(IN)  :: AIRMASS ! Total air (kg)
      REAL*8,INTENT(IN)  :: PCENTER ! Pressure (hPa)?
      REAL*8,INTENT(IN)  :: TCENTER ! Temperature (K)

      REAL*8,INTENT(OUT) :: SO4AER  ! Aerosol SO4 (kg)
      REAL*8,INTENT(OUT) :: SO4GAS  ! Gaseous SO4 (kg)
!
! !REVISION HISTORY: 
!  11 Apr 2013 - S. D. Eastham - Initial version
!  07 Feb 2014 - R. Yantosca   - Cosmetic changes
!  14 Feb 2014 - R. Yantosca   - Changed parallel DO loop order to L-J-I,
!                                which is the most efficient configuration
!  21 Feb 2014 - M. Sulprizio  - Now pass Input_Opt, State_Met, and State_Chm
!                                objects via the arg list
!  08 Dec 2014 - S. D. Eastham - Converted to single-box model for MOSAIC
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE      :: FIRST=.TRUE.
      REAL*8,PARAMETER   :: GF_THRESHOLD = 0.0d0
      REAL*8,PARAMETER   :: GF_RANGE     = 1.0d-8
      REAL*8,PARAMETER   :: GF_DELTAHBYR = 10156.d0
      REAL*8,PARAMETER   :: GF_T0        = 360.d0
      REAL*8,PARAMETER   :: GF_TC        = 905.d0
      REAL*8,SAVE        :: GF_LOGP0
      REAL*8,SAVE        :: GF_BFACTOR
      REAL*8,SAVE        :: GF_ATMCONV
      REAL*8,SAVE        :: GF_INVT0
      REAL*8             :: GF_INVT,GF_LOGPSULFATE,GF_CFACTOR
      REAL*8             :: GF_AFACTOR
      REAL*8             :: GF_PP,GF_PVAP,GF_DIFF

      REAL*8             :: H2SO4SUM, INVAIR, AERFRAC 

      !=================================================================
      ! SO4_PARTITION begins here!
      !=================================================================

      IF (FIRST) THEN
         FIRST = .FALSE.
         ! Calculate H2SO4 gas phase prefactors
         GF_INVT0 = 1.d0/GF_T0
         GF_LOGP0 = (-1.d0*GF_DELTAHBYR*GF_INVT0) + 16.259d0
         GF_BFACTOR = 0.38d0/(GF_TC - GF_T0)
         GF_ATMCONV = LOG(101325.d0)
      ENDIF

      INVAIR  = AIRMW / AIRMASS
      ! mol/mol
      H2SO4SUM = SO4KG*INVAIR/96.d0
      ! Use approximation from Kumala (1990)
      GF_PP = H2SO4SUM*PCENTER
      GF_INVT = 1./TCENTER
      GF_CFACTOR = 1.d0+(LOG(GF_T0*GF_INVT))-(GF_T0*GF_INVT)
      GF_LOGPSULFATE = GF_LOGP0 + (GF_DELTAHBYR*(GF_INVT0 -
     &         GF_INVT + (GF_BFACTOR*GF_CFACTOR)))
      GF_LOGPSULFATE = GF_LOGPSULFATE + GF_ATMCONV
      GF_PVAP = 1.d-2 * EXP(GF_LOGPSULFATE)
 
      GF_DIFF = GF_PP - GF_PVAP
      IF (GF_DIFF .lt. 0d0) THEN
         AERFRAC = 0.d0
      ELSEIF (.not.IS_SAFE_DIV(GF_DIFF,GF_PP)) THEN
         ! Hmm
         AERFRAC = 0.d0
      ELSE
         AERFRAC = GF_DIFF/GF_PP
      ENDIF

!      GF_DIFF = (GF_PVAP+GF_THRESHOLD) - GF_PP
!      IF (GF_DIFF .lt. 0d0) THEN
!         AERFRAC = 1.d0
!      ELSEIF (GF_DIFF .lt. GF_RANGE) THEN
!         AERFRAC = 1.d0-(GF_DIFF/GF_RANGE)
!      ELSE
!         AERFRAC = 0d0
!      ENDIF
      SO4AER = AERFRAC * SO4KG
      SO4GAS = SO4KG - SO4AER

      END SUBROUTINE SO4_PARTITION
!
!EOC

      END MODULE MOSAIC_MOD
