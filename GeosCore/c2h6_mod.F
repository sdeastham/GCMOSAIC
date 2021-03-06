!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: c2h6_mod
!
! !DESCRIPTION: Module C2H6\_MOD contains variables and routines used for the
!  tagged C2H6 (ethane) simulation. (xyp, qli, bmy, 7/28/01, 4/5/06)
!\\
!\\
! !INTERFACE:
!
      MODULE C2H6_MOD
!
! !USES:
!
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: GET_C2H6_ANTHRO
!      PUBLIC  :: EMISSC2H6
      PUBLIC  :: CHEMC2H6
      PUBLIC  :: CLEANUP_C2H6
      PUBLIC  :: INIT_C2H6
!
! PRIVATE DATA MEMBERS:
!
      PRIVATE :: NGASC2H6
      PRIVATE :: FMOL_C2H6
      PRIVATE :: XNUMOL_C2H6
!
! !REMARKS:
!  Setting LSPLIT = T in "input.geos" will run with the following tracers:
!     (1) Total C2H6
!     (2) C2H6 from biomass burning
!     (3) C2H6 from biofuel burning
!     (4) C2H6 from natural gas leaking/venting (e.g. "anthro" C2H6)
!                                                                             .
!  Setting LSPLIT = F in "input.geos" will run w/ the following tracers:
!     (1) Total C2H6
!
! !REVISION HISTORY:
!  (1 ) Eliminated obsolete code from 1/02 (bmy, 2/27/02)
!  (2 ) Now divide module header into MODULE PRIVATE, MODULE VARIABLES, and
!        MODULE ROUTINES sections.  Updated comments (bmy, 5/28/02)
!  (3 ) Now reference BXHEIGHT and T from "dao_mod.f".  Also references
!        "error_mod.f".  Removed obsolete code.  Now references F90 module
!         tracerid_mod.f". (bmy, 11/15/02)
!  (4 ) Now references "grid_mod.f" and the new "time_mod.f" (bmy, 2/11/03)
!  (5 ) Now references "directory_mod.f", "logical_mod.f", and "tracer_mod.f".
!        (bmy, 7/20/04)
!  (6 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (7 ) Now modified
!  22 Mar 2012 - M. Payer    - Add function GET_C2H6_ANTHRO
!  22 Mar 2012 - M. Payer    - Added ProTeX headers
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  11 Apr 2014 - R. Yantosca - Now declare INIT_C2H6 public, so that it can
!                              be called from input_mod.F
!  20 Jun 2014 - R. Yantosca - Remove obsolete code now replaced by HEMCO
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! Molecular weight of C2H6 [kg/mole]
      REAL*8, PARAMETER   :: FMOL_C2H6 = 30d-3

      ! Ratio of [molec C2H6/kg C2H6]
      REAL*8, PARAMETER   :: XNUMOL_C2H6 = 6.022d+23/FMOL_C2H6
!
! PRIVATE TYPES:
!
      ! Global monthly mean natural gas C2H6 emissions
      REAL*8, ALLOCATABLE :: NGASC2H6(:,:)

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: chemc2h6
!
! !DESCRIPTION: Subroutine CHEMC2H6 performs C2H6 chemistry. Loss of C2H6 is
!  via reaction with OH. (xyp, qli, bmy, 10/19/99, 7/20/04)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CHEMC2H6( am_I_Root, Input_Opt,
     &                     State_Met, State_Chm, RC )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GLOBAL_OH_MOD,      ONLY : OH,        GET_GLOBAL_OH
      USE TIME_MOD,           ONLY : GET_MONTH, GET_TS_CHEM

! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
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
! !REVISION HISTORY:
!  (1 ) Now do chemistry all the way to the model top. 
!  (2 ) Use monthly mean OH fields for oxidation -- reference the monthly 
!        mean OH array and the routine which reads it from disk in 
!       "global_oh_mod.f" (bmy, 1/25/02)
!  (3 ) Now reference T from "dao_mod.f".  Also make FIRSTCHEM a local SAVEd
!        variable. (bmy, 11/15/02)
!  (4 ) Now use functions GET_MONTH and GET_TS_CHEM from "time_mod.f".
!  (5 ) Now reference STT & N_TRACERS from "tracer_mod.f".  Now reference 
!        LSPLIT from "logical_mod.f" (bmy, 7/20/04)
!  22 Mar 2012 - M. Payer    - Added ProTeX headers
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, State_Chm, RC
!  23 Oct 2013 - R. Yantosca - Now pass objects to GET_GLOBAL_OH routine
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! SAVEd scalars
      LOGICAL, SAVE          :: FIRSTCHEM = .TRUE.
      INTEGER, SAVE          :: LASTMONTH = -99

      ! Scalars
      LOGICAL                :: LSPLIT
      INTEGER                :: N_TRACERS
      INTEGER                :: I, J, L, N
      REAL*8                 :: DTCHEM, KRATE

      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL*8, POINTER :: STT(:,:,:,:)

      !=================================================================
      ! CHEMC2H6 begins here! 
      !=================================================================

      ! Assume success
      RC        = GIGC_SUCCESS

      ! Get fields from Input_Opt
      LSPLIT    = Input_Opt%LSPLIT
      N_TRACERS = Input_Opt%N_TRACERS


      IF ( FIRSTCHEM ) THEN
         FIRSTCHEM = .FALSE.  ! save for future use?
      ENDIF

      ! DTCHEM is the chemistry timestep in seconds
      DTCHEM = GET_TS_CHEM() * 60d0

      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      ! (mpayer, 12/6/12)
      STT => State_Chm%Tracers

      !=================================================================
      ! Read in the tropospheric OH fields of the (LMN)th month
      ! OH data will be saved into the OH array of "global_oh_mod.f" 
      !=================================================================
      IF ( GET_MONTH() /= LASTMONTH ) THEN
         CALL GET_GLOBAL_OH( am_I_Root, Input_Opt, 
     &                       State_Met, GET_MONTH(), RC )
         LASTMONTH = GET_MONTH()
      ENDIF

      !=================================================================
      ! Do C2H6 chemistry -- C2H6 Loss due to chemical reaction with OH
      !
      ! DECAY RATE: The decay rate (KRATE) is calculated by:
      !
      !    OH + C2H6 -> H2O + C2H5 (JPL '97)
      !    k = 8.7D-12 * exp(-1070/T)
      !
      ! KRATE has units of [ molec^2 C2H6 / cm6 / s ]^-1.
      !=================================================================
!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED ) 
!$OMP+PRIVATE( I, J, L, N, KRATE )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         
         ! Loss rate [molec2 C2H6/cm6/s]^-1
         KRATE = 8.7d-12 * EXP( -1070.0d0 / State_Met%T(I,J,L) )

         ! Apply loss to total C2H6 (tracer #1)
         STT(I,J,L,1) = STT(I,J,L,1) *
     &                  ( 1d0 - KRATE * OH(I,J,L) * DTCHEM )

         ! If we are running w/ tagged tracers,
         ! then also apply the loss to each of these
         IF ( LSPLIT ) THEN 
            DO N = 2, N_TRACERS
            
               ! Subtract loss of C2H6 by OH and store in STT [kg C2H6]
               ! Loss = k * [C2H6] * [OH] * dt
               STT(I,J,L,N) = STT(I,J,L,N) *
     &                        ( 1d0 - KRATE * OH(I,J,L) * DTCHEM )
            ENDDO
         ENDIF

      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Free pointer
      NULLIFY( STT )

      END SUBROUTINE CHEMC2H6
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_c2h6_anthro
!
! !DESCRIPTION: Function GET\_C2H6\_ANTHRO returns the monthly average
!  anthropogenic C2H6 emissions at GEOS-Chem grid box (I,J).  Data will
!  be returned in units of [atoms C/cm2/s].
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_C2H6_ANTHRO( I, J, N ) RESULT( C2H6_ANTHRO )
!
! !USES:
!
      USE TRACERID_MOD, ONLY : IDTC2H6
      USE CMN_SIZE_MOD       ! Size parameters
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)   :: I   ! GEOS-Chem longitude index
      INTEGER, INTENT(IN)   :: J   ! GEOS-Chem latitude index
      INTEGER, INTENT(IN)   :: N   ! GEOS-Chem tracer index
!
! !RETURN VALUE:
!
      REAL*8                :: C2H6_ANTHRO
! 
! !REVISION HISTORY: 
!  22 Mar 2012 - M. Payer    - Initial version adapted from GET_RETRO_ANTHRO
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      !  GET_C2H6_ANTHRO begins here
      !=================================================================

      IF ( N == IDTC2H6 ) THEN
         C2H6_ANTHRO = NGASC2H6(I,J)
      ENDIF

      END FUNCTION GET_C2H6_ANTHRO
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_c2h6
!
! !DESCRIPTION: Subroutine INIT\_C2H6 allocates and zeroes the NGASC2H6 array,
!  which holds global monthly mean natural gas C2H6 emissions.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_C2H6( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE ERROR_MOD,          ONLY : ALLOC_ERR
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?! 
!
! !REMARKS:
!  This routine is now called from GIGC_INIT_EXTRA (in GeosCore/input_mod.F).
!
! !REVISION HISTORY:
!  (1 ) Now references ALLOC_ERR from "error_mod.f" (bmy, 10/15/02)
!  22 Mar 2012 - M. Payer    - Added ProTeX headers
!  11 Apr 2014 - R. Yantosca - Now accept am_I_Root, Input_Opt, RC arguments
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      !=================================================================
      !  INIT_C2H6 begins here!
      !=================================================================

      ! Assume success
      RC = GIGC_SUCCESS

      ! Allocate NGASC2H6 array
      ALLOCATE( NGASC2H6( IIPAR, JJPAR ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'NGASC2H6' )

      ! Zero NGASC2H6 array
      NGASC2H6 = 0d0

      END SUBROUTINE INIT_C2H6
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_c2h6
!
! !DESCRIPTION: Subroutine CLEANUP\_C2H6 deallocates the natural gas C2H6
!  emission array.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_C2H6
! 
! !REVISION HISTORY:
!  22 Mar 2012 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! CLEANUP_C2H6 begins here!
      !=================================================================

      IF ( ALLOCATED( NGASC2H6 ) ) DEALLOCATE( NGASC2H6 )

      END SUBROUTINE CLEANUP_C2H6
!EOC
      END MODULE C2H6_MOD
