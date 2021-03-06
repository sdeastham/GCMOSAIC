!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: emissions_mod.F90
!
! !DESCRIPTION: Module emissions\_mod.F90 is a wrapper module to interface
! GEOS-Chem and HEMCO. It basically just calls the GEOS-Chem - HEMCO interface
! routines. For some specialty sims, a few additional steps are required that
! are also executed here.
!\\
!\\
! !INTERFACE:
!
MODULE EMISSIONS_MOD
!
! !USES:
!
  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: EMISSIONS_INIT
  PUBLIC :: EMISSIONS_RUN
  PUBLIC :: EMISSIONS_FINAL
!
! !REVISION HISTORY:
!  27 Aug 2014 - C. Keller   - Initial version. 
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: EMISSIONS_INIT
!
! !DESCRIPTION: Subroutine EMISSIONS\_INIT calls the HEMCO - GEOS-Chem
! interface initialization routines.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE EMISSIONS_INIT( am_I_Root, Input_Opt, State_Met, State_Chm, RC ) 
!
! !USES:
!
    USE GIGC_ErrCode_Mod
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE GIGC_State_Chm_Mod, ONLY : ChmState
    USE ERROR_MOD,          ONLY : ERROR_STOP
    USE HCOI_GC_MAIN_MOD,   ONLY : HCOI_GC_INIT
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! root CPU?
    TYPE(MetState),   INTENT(IN   )  :: State_Met  ! Met state
    TYPE(ChmState),   INTENT(IN   )  :: State_Chm  ! Chemistry state 
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(INOUT)  :: Input_Opt  ! Input opts
    INTEGER,          INTENT(INOUT)  :: RC         ! Failure or success
!
! !REVISION HISTORY: 
!  27 Aug 2014 - C. Keller    - Initial version 
!EOP
!------------------------------------------------------------------------------
!BOC

    !=================================================================
    ! EMISSIONS_INIT begins here!
    !=================================================================

    ! Assume success
    RC = GIGC_SUCCESS

    ! Initialize the HEMCO environment for this GEOS-Chem run.
    CALL HCOI_GC_Init( am_I_Root, Input_Opt, State_Met, State_Chm, RC ) 
    IF ( RC/=GIGC_SUCCESS ) RETURN 

  END SUBROUTINE EMISSIONS_INIT
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: EMISSIONS_RUN
!
! !DESCRIPTION: Subroutine EMISSIONS\_RUN calls the HEMCO - GEOS-Chem
! interface run routines.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE EMISSIONS_RUN( am_I_Root, Input_Opt, State_Met, State_Chm, RC ) 
!
! !USES:
!
    USE GIGC_ErrCode_Mod
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE GIGC_State_Chm_Mod, ONLY : ChmState
    USE ERROR_MOD,          ONLY : ERROR_STOP
    USE HCOI_GC_MAIN_MOD,   ONLY : HCOI_GC_RUN
    USE DUST_MOD,           ONLY : DUSTMIX
    USE CO2_MOD,            ONLY : EMISSCO2
    USE GLOBAL_CH4_MOD,     ONLY : EMISSCH4
    USE TRACERID_MOD,       ONLY : IDTCH4

    ! Use old mercury code for now (ckeller, 09/23/2014)
    USE MERCURY_MOD,        ONLY : EMISSMERCURY

    ! For UCX, use Seb's routines for now
#if defined( UCX )
    USE UCX_MOD,            ONLY : EMISS_BASIC
#endif
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! root CPU?
    TYPE(MetState),   INTENT(IN   )  :: State_Met  ! Met state
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState),   INTENT(INOUT)  :: State_Chm  ! Chemistry state 
    TYPE(OptInput),   INTENT(INOUT)  :: Input_Opt  ! Input opts
    INTEGER,          INTENT(INOUT)  :: RC         ! Failure or success
!
! !REVISION HISTORY: 
!  27 Aug 2014 - C. Keller    - Initial version 
!EOP
!------------------------------------------------------------------------------
!BOC
 
    !=================================================================
    ! EMISSIONS_RUN begins here!
    !=================================================================

    ! Assume success
    RC = GIGC_SUCCESS

    ! Run HEMCO
    CALL HCOI_GC_RUN( am_I_Root, Input_Opt, State_Met, State_Chm, RC ) 
    IF ( RC /= GIGC_SUCCESS ) RETURN 
  
    ! PBL mixing is not applied to dust emissions. Instead, they become 
    ! directly added to the tracer arrays.
    CALL DUSTMIX( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
 
    ! For CO2 simulation, emissions are not added to Trac_Tend and hence
    ! not passed to the Tracers array during PBL mixing. Thus, need to add 
    ! emissions explicitly to the tracers array here.
    IF ( Input_Opt%ITS_A_CO2_SIM ) THEN
       CALL EMISSCO2( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
       IF ( RC /= GIGC_SUCCESS ) RETURN 
    ENDIF

    ! For CH4 simulation or if CH4 is defined, call EMISSCH4. 
    ! This will get the individual CH4 emission terms (gas, coal, wetlands, 
    ! ...) and write them into the individual emissions arrays defined in
    ! global_ch4_mod (CH4_EMIS), from where the final emission array is
    ! assembled and passed to STT or Trac_Tend.
    ! This is a wrapper for backwards consistency, in particular for the
    ! ND58 diagnostics.
    IF ( Input_Opt%ITS_A_CH4_SIM .OR.            &
       ( IDTCH4 > 0 .and. Input_Opt%LCH4EMIS ) ) THEN
       CALL EMISSCH4( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
       IF ( RC /= GIGC_SUCCESS ) RETURN 
    ENDIF

    ! For UCX, use Seb's routines for stratospheric species for now.
#if defined( UCX )
    IF ( Input_Opt%LBASICEMIS ) THEN
       CALL EMISS_BASIC( am_I_Root, Input_Opt, State_Met, State_Chm )
    ENDIF
#endif

    ! For mercury, use old emissions code for now
    IF ( Input_Opt%ITS_A_MERCURY_SIM ) THEN
       CALL EMISSMERCURY ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
    ENDIF

    ! Return w/ success
    RC = GIGC_SUCCESS

  END SUBROUTINE EMISSIONS_RUN
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: EMISSIONS_FINAL
!
! !DESCRIPTION: Subroutine EMISSIONS\_FINAL calls the HEMCO - GEOS-Chem
! interface finalization routines.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE EMISSIONS_FINAL( am_I_Root )
!
! !USES:
!
    USE HCOI_GC_MAIN_MOD, ONLY : HCOI_GC_FINAL
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! root CPU?
!
! !REVISION HISTORY: 
!  27 Aug 2014 - C. Keller    - Initial version 
!EOP
!------------------------------------------------------------------------------
!BOC
 
    !=================================================================
    ! EMISSIONS_FINAL begins here!
    !=================================================================

    CALL HCOI_GC_Final( am_I_Root )

  END SUBROUTINE EMISSIONS_FINAL
!EOC
END MODULE EMISSIONS_MOD
