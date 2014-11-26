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
      PRIVATE :: WATER_VAPOR
      PRIVATE :: SET_MOSAIC_AIR
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
      USE CMN_GCTM_MOD,       ONLY : AIRMW
      USE ERROR_MOD,          ONLY : DEBUG_MSG
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE ERROR_MOD,          ONLY : SAFE_DIV
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE HCO_EMISLIST_MOD,   ONLY : HCO_GetPtr
      USE TIME_MOD,           ONLY : GET_MONTH
      USE TIME_MOD,           ONLY : ITS_A_NEW_MONTH
      ! Gas-phase species
      USE TRACERID_MOD,       ONLY : IDTSO4
      USE TRACERID_MOD,       ONLY : IDTHNO3
      USE TRACERID_MOD,       ONLY : IDTHCl
      USE TRACERID_MOD,       ONLY : IDTNH3
      USE TRACERID_MOD,       ONLY : IDTMSA
      ! Aerosol species
      USE TRACERID_MOD,       ONLY : IDTaNUM,  IDTaHyst, IDTaH2O
      USE TRACERID_MOD,       ONLY : IDTaSO4,  IDTaNO3,  IDTaNH4
      USE TRACERID_MOD,       ONLY : IDTaCl,   IDTaCO3,  IDTaMSA
      USE TRACERID_MOD,       ONLY : IDTaNa,   IDTaCa,   IDTaOC
      USE TRACERID_MOD,       ONLY : IDTaBC,   IDTaOIN
      ! SOA not yet implemented

      ! MOSAIC indexing
      use module_data_mosaic_main, only: r8,naer_tot,ngas_max, !Parameters
     &   ntot_max,avogad,kh2so4,khno3,khcl,knh3,kmsa,knum_a,piover6,
     &   ksigmag_a,kjhyst_a,kwater_a,ntot_used,cnn

      use module_data_mosaic_aero, only: nbin_a_max,naer,ngas_volatile,
     &   jtotal,nbin_a,ih2so4_g,ihno3_g,ihcl_g,inh3_g,iMSA_a_g,
     &   jhyst_lo, jhyst_up, jhyst_undefined,mhyst_method,
     &   mhyst_uporlo_waterhyst,all_solid,all_liquid,mixed,no_aerosol

      use module_data_mosaic_main, only: pr_atm, RH, te, o2, h2, h2o,
     &   pr_atm_old, RH_old, te_old, speed_molec, pblh, pblh_old,
     &   cair_mlc, cair_molm3
     

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
      INTEGER                  :: I,    J,    L,    N
      CHARACTER(LEN=255)       :: X
      real(r8)                 :: molecPerCm3 ! Molecules of air per cm3
      real(r8)                 :: boxPerCm3 ! 1/(Box volume in cm3)
      real(r8)                 :: cmpConv   ! kg -> umol/m3 (needs /(g/mol))

      ! Unused outputs from aerchemistry
      real(r8)                 :: cair_mol_m3, cair_mol_cc

      LOGICAL                  :: LPRT
      INTEGER                  :: N_TRACERS
      INTEGER,DIMENSION(nbin_a_max) :: jaerosolstate
      REAL*8                   :: XNUMOL(Input_Opt%N_TRACERS)
      REAL*8                   :: TRACER_MW_KG(Input_Opt%N_TRACERS)
                  
      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL*8, POINTER          :: STT(:,:,:,:)

      ! For getting fields from HEMCO 
      LOGICAL                  :: aIR

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
 
         ! Check aerosol phase - only look at fiNa_al bin
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
         IF ( IDTaMSA(NBINMOS)  == 0 ) 
     &       CALL ERROR_STOP( 'Aerosol-phase MSA is undefined!', X)
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
         IF ( IDTaH2O(NBINMOS)  == 0 ) 
     &       CALL ERROR_STOP( 'Aerosol-phase H2O is undefined!', X)
         IF ( IDTaNUM(NBINMOS)  == 0 ) 

         ! Aerosol data
     &       CALL ERROR_STOP( 'Number concentration is undefined!', X)
         IF ( IDTaHyst(NBINMOS)  == 0 ) 
     &       CALL ERROR_STOP( 'Hysteresis tracer is undefined!', X)

         ! Initialize arrays
         CALL INIT_MOSAIC

         ! Reset first-time flag
         FIRST = .FALSE. 
      ENDIF
     
      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      ! (mpayer, 12/6/12)
      STT => State_Chm%Tracers

      !=================================================================
      ! Loop over grid boxes and call MOSAIC
      !=================================================================

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,    J,      L,     N,                  VOL )
!$OMP+PRIVATE( CNN, NOFFSET,boxPerCm3,molecPerCm3,cmpConv  )
!$OMP+PRIVATE( iter_MESA, aH2O_a, gam_ratio, jaerosolstate  )
!$OMP+PRIVATE( dp_wet_a, cair_mol_m3, cair_mol_cc           )
!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Only applying MOSAIC in troposphere
         IF ( ITS_IN_THE_STRATMESO( I, J, L, State_Met ) ) THEN 
            ! TEMPORARY - Move all SO4 into gas-phase for UCX
            DO N = 1,NBINMOS
                STT(I,J,L,IDTSO4) = STT(I,J,L,IDTSO4) + 
     &                              STT(I,J,L,IDTaSO4(N))
                STT(I,J,L,IDTaSO4(N)) = 0.d0
            ENDDO
            CYCLE
         ENDIF

         ! Reinitialize
         jaerosolstate(:) = 0
         dp_wet_a(:) = 0

         ! Set MOSAIC variables
         RH       = State_Met%RH(I,J,L)   ! %
         te       = State_Met%T(I,J,L)    ! K
         te_old   = te                    ! K
         pr_atm   = State_Met%PMID(I,J,L) ! hPa
         pr_atm   = pr_atm * 9.86923267d-4! atm
         pr_atm_old= pr_atm

         ! Volume of grid box [m3] 
         VOL      = State_Met%AIRVOL(I,J,L)

         ! Conversion factors
         ! Gridboxes per cm3
         boxPerCm3   = (1.d-6/VOL)
         cmpConv = 1.d9/VOL

         ! Set the MOSAIC air composition
         ! mol/m3
         ! kg/m3 * (m3/cm3 * g/kg) * (mol/g)
         cair_molm3 = State_Met%AIRDEN(L,I,J) * 1.d-3 
     &                            / AIRMW
         ! molec/cm3
         cair_mlc = 1.d-6*avogadro*cair_molm3
         cair_mlc_old = cair_mlc
         cair_molm3_old = cair_molm3

         h2o = WATER_VAPOR(RH, cair_mlc, te, pr_atm)
         o2  = 0.21*cair_mlc
         h2  = 0.58e-6*cair_mlc

         ppb = 1.e+9

         ! Molecules of air per cm3
         ! kg/m3 * (m3/cm3 * g/kg) * molec/mol * (mol/g)
         molecPerCm3 = cair_mlc

         ! Build the MOSAIC interface array
         CNN = 0d0

         ! STT is in kg
         ! MOSAIC wants molec/cm3
         CNN(kh2so4)  = XNUMOL(IDTSO4)*boxPerCm3*STT(I,J,L,IDTSO4)
         CNN(khno3)   = XNUMOL(IDTHNO3)*boxPerCm3*STT(I,J,L,IDTHNO3)
         CNN(khcl)    = XNUMOL(IDTHCl)*boxPerCm3*STT(I,J,L,IDTHCl)
         CNN(knh3)    = XNUMOL(IDTNH3)*boxPerCm3*STT(I,J,L,IDTNH3)
         CNN(kmsa)    = XNUMOL(IDTMSA)*boxPerCm3*STT(I,J,L,IDTMSA)

         DO N = 1, NBINMOS
            NOFFSET = ngas_max + naer_tot*(n - 1)
            ! #/cm3 from "kg" - convert to moles, then molecules
            ! Divide by volume
            CNN(NOFFSET+knum_a) = (STT(I,J,L,IDTaNUM(N)) / 
     &         TRACER_MW_KG(IDTaNUM(N))) * AVOGADRO * boxPerCm3

            ! These would be used in a modal representation
            CNN(NOFFSET+kpdry_a) = 0d0       ! Dry diameter in um
            CNN(NOFFSET+ksigmag_a) = 0d0     ! Geometric SD

            ! Aerosol water content (kg/m3 air)
            CNN(NOFFSET+kwater_a) = STT(I,J,L,IDTaH2O(N))/VOL

            ! Hysteresis
            IF (mhyst_method == mhyst_uporlo_waterhyst) THEN
               CNN(NOFFSET+kjhyst_a) = STT(I,J,L,IDTaHyst)/VOL
            ELSE
               CNN(NOFFSET+kjhyst_a) = 0d0
            ENDIF

            ! Aerosol components
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

         ! Run MOSAIC
            CALL aerchemistry(0, 0, dt_sec,                             !intent-ins
     &       jaerosolstate, dp_wet_a,                                 !intent-inouts
     &       cair_mol_m3, cair_mol_cc, aH2O_a, gam_ratio, iter_MESA ) !intent-outs

         ! Save data back into tracer array
         ! STT wants kg
         ! MOSAIC is in molec/cm3
         STT(I,J,L,IDTSO4)  = CNN(kh2so4)/(boxPerCm3*XNUMOL(IDTSO4))
         STT(I,J,L,IDTHNO3) = CNN(kHNO3) /(boxPerCm3*XNUMOL(IDTHNO3))
         STT(I,J,L,IDTHCl)  = CNN(kHCl)  /(boxPerCm3*XNUMOL(IDTHCl))
         STT(I,J,L,IDTNH3)  = CNN(kNH3)  /(boxPerCm3*XNUMOL(IDTNH3))
         STT(I,J,L,IDTMSA)  = CNN(kMSA)  /(boxPerCm3*XNUMOL(IDTMSA))

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
            STT(I,J,L,IDTaSO4(N)) = CNN(NOFFSETT+kwater_a+iso4_a) *
     &                 TRACER_MW_KG(IDTaSO4(N)) / cmpConv
            STT(I,J,L,IDTaNO3(N)) = CNN(NOFFSETT+kwater_a+ino3_a) *
     &                 TRACER_MW_KG(IDTaNO3(N)) / cmpConv
            STT(I,J,L,IDTaCl(N)) = CNN(NOFFSETT+kwater_a+iCl_a) *
     &                 TRACER_MW_KG(IDTaCl(N)) / cmpConv
            STT(I,J,L,IDTaNH4(N)) = CNN(NOFFSETT+kwater_a+iNH4_a) *
     &                 TRACER_MW_KG(IDTaNH4(N)) / cmpConv
            STT(I,J,L,IDTaCO3(N)) = CNN(NOFFSETT+kwater_a+iCO3_a) *
     &                 TRACER_MW_KG(IDTaCO3(N)) / cmpConv
            STT(I,J,L,IDTaMSA(N)) = CNN(NOFFSETT+kwater_a+iMSA_a) *
     &                 TRACER_MW_KG(IDTaMSA(N)) / cmpConv
            STT(I,J,L,IDTaNa(N)) = CNN(NOFFSETT+kwater_a+iNa_a) *
     &                 TRACER_MW_KG(IDTaNa(N)) / cmpConv
            STT(I,J,L,IDTaCa(N)) = CNN(NOFFSETT+kwater_a+iCa_a) *
     &                 TRACER_MW_KG(IDTaCa(N)) / cmpConv
            STT(I,J,L,IDTaOC(N)) = CNN(NOFFSETT+kwater_a+iOC_a) *
     &                 TRACER_MW_KG(IDTaOC(N)) / cmpConv
            STT(I,J,L,IDTaBC(N)) = CNN(NOFFSETT+kwater_a+iBC_a) *
     &                 TRACER_MW_KG(IDTaBC(N)) / cmpConv
            STT(I,J,L,IDTaOIN(N)) = CNN(NOFFSETT+kwater_a+iOIN_a) *
     &                 TRACER_MW_KG(IDTaOIN(N)) / cmpConv
            ! SOA not yet implemented
         ENDDO
   
      ENDDO
      ENDDO
      ENDDO 
!$OMP END PARALLEL DO

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
      FUNCTION WATER_VAPOR(rh, cair_mld, te, pr_atm) 
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
      REAL(r8),INTENT(OUT) :: WATER_VAPOR
      REAL(r8),INTENT(IN)  :: RH, CAIR_MLC, TE, PR_ATM
      REAL(r8)             :: T_STEAM, PR_STD, A, ARG, PR_H2O

      !=================================================================
      ! WATER_VAPOR begins here!
      !=================================================================

      t_steam = 373.15 ! steam temperature  [K]
      pr_std   = 1.0   ! standard pressure  [atm]

      a      = 1.0 - t_steam/te
      arg    = (((-.1299*a -.6445)*a -1.976)*a +13.3185)*a
      pr_h2o = pr_std*exp(arg)  ! [atm]
      WATER_VAPOR = RH*(pr_h2o/pr_atm)*cair_mlc/100. ! [molec/cc]

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

      use module_data_mosaic_pmcmos, only:  &
           temp_profile_fname, pblh_profile_fname,  &
           aer_init_fname, aer_back_fname, aer_emit_fname,  &
           gas_init_fname, gas_back_fname, gas_emit_fname,  &
           !BSINGH - 05/28/2013(RCE updates)
           aer_emit_adjust_fac,  &
           msolar_pmcmos_dtshift,  &
           pmcmos_inputfile_version, pmcmos_inputfile_directory
      !BSINGH - 05/28/2013(RCE updates ENDS)
      use module_pmcmos_init, only:  pmcmos_init
      
      use module_sect_iface, only:  &!BSINGH - 05/28/2013(RCE updates)
           sect_iface_allocate_memory

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
      real(r8) :: tmpa, tmpb, dt_min

      !=================================================================
      ! INIT_MOSAIC  begins here!
      !=================================================================

      ! Initialize MOSAIC data modules
      CALL init_data_modules

      ! Set the MOSAIC time step to match the chemistry time step
      dt_min = 1.d0*GET_TS_CHEM()
      dt_sec = 60.*dt_min 

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
#else
      nsize_aer_tmp = 8
#endif

      ! Total number of bins
      nbin_a = ntype_aer*nsize_aer_tmp

      ! Size framework:
      !   1: Modal
      !   2: Unstructured
      !   3: Sectional
#if defined( MOSAIC4 ) || defined( MOSAIC8 )
      msize_framework = msectional
#else
      msize_Framework = munstructured
#endif
 
      ! This has to be 1
      method_bcfrac = 1

      ! This has to be either 11 or 12
      method_kappa = 12

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
      ntot_used = ngas_max + naer_max + ncldtot*ncldbin_used

      ! Allocate memory? This is a dummy routine in MOSAIC
      !call mosaic_allocate_memory( 0 )

      ! Allocate memory for sectional interface
      call sect_iface_allocate_memory
 
      ! Number of size bins for each aerosol type
      ! Currently assumes uniform
      nsize_aer(:) = nsize_aer_tmp

      ! If set to 2, bin sizes are calculated and the initial dry diameter is stored
      ! for each bin into CNN. The bin EDGES are always calculated for a sectional
      ! framework, but the CNN array entry is only set if this is equal to 2.
      maersize_init_flag1 = 2
      if (msize_framework == msectional) then
          ! Bin upper and lower limits in um
          dlo_aersize_init = 1.d-2
          dhi_aersize_init = 1.d1
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
         dlo_aersize_init = 0.0 ; dhi_aersize_init = 0.0
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
      mnewnuc_flag1 = 1

      ! Not currently used
      msectional_flag1 = 0

      ! If running in bulk mode...
#if !(defined( MOSAIC4) || defined( MOSAIC8 ))
      ifreq_coag = 0
      mmovesect_flag1 = 0
      mnewnuc_flag1 = 1
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
      mshellcore = 1

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
      cnn(:) = 0.0_r8
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
      CALL INIT_AEROSOL

      END SUBROUTINE INIT_MOSAIC

!EOC
      END MODULE MOSAIC_MOD
