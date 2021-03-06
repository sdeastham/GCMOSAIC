!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: paranox_util_mod.F
!
! !DESCRIPTION: Module ParaNOx\_Util\_Mod contains subroutines for reading 
! and interpolating look up tables necessary for the PARANOX (PARAmeterization
! of emitted NOX) ship plume model developed by G.C.M. Vinken.  We follow
! the implementation of PARANOX into the GEOS-Chem CTM.
!\\
!\\
! !INTERFACE: 
!
      MODULE PARANOX_UTIL_MOD
!
! !USES:
!
      USE HCO_ERROR_MOD
      USE inquireMod, ONLY : findFreeLUN

      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: READ_PARANOX_LUT
      PUBLIC  :: INTERPOLATE_LUT2
!
! !MODULE VARIABLES
!
      ! fracnox                 = look up table for fraction of NOx remaining
      !                           for ship emissions (gvinken, 6/6/10)
      ! intope                  = look up table for integrated Ozone Production
      !                           Efficiency for ship emiss (gvinken, 6/6/10)
      REAL*4 ::   fracnox(4,4,4,12,12,4,5)
      REAL*4 ::   intope(4,4,4,12,12,4,5)
!
! !REMARKS
!  References:
!  ============================================================================
!  (1 ) Vinken, G.C.M., Boersma, K.F., Jacob, D.J., and Meijer, E.W.:
!       Accounting for non-linear chemistry of ship plumes in the GEOS-Chem
!       global chemistry transport model, Atmos. Chem. Phys., 11, 11707-11722,
!       doi:10.5194/acp-11-11707-2011, 2011.
!
! !REVISION HISTORY:
!  06 Feb 2012 - M. Payer    - Initial version
!  01 Mar 2012 - R. Yantosca - Use updated GET_LOCALTIME from time_mod.F
!  03 Aug 2012 - R. Yantosca - Move calls to findFreeLUN out of DEVEL block
!  23 Aug 2013 - C. Keller   - Renamed from paranox_mod to paranox_util_mod
!  28 Jul 2014 - C. Keller   - Broke dependencies from GEOS-Chem code: Now
!                              use HEMCO routines for error handling etc.
!EOP
!------------------------------------------------------------------------------
!BOC
      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Read_ParaNOx_Lut
!
! !DESCRIPTION: Subroutine ReaD\_ParaNOx\_Lut reads look up tables for use in
!  the PARANOX ship plume model (G.C.M. Vinken)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_PARANOX_LUT( FracNOx_FILE, IntOPE_FILE, RC )
!
! !INPUT PARAMETERS:
!
      CHARACTER(LEN=*), INTENT(IN)    :: FracNOx_FILE
      CHARACTER(LEN=*), INTENT(IN)    :: IntOPE_FILE
      INTEGER,          INTENT(INOUT) :: RC
!
! !REVISION HISTORY:
!  06 Feb 2012 - M. Payer    - Initial version modified from code provided by
!                              G.C.M. Vinken
!  01 Aug 2012 - R. Yantosca - Add reference to findFreeLUN from inqure_mod.F90
!  03 Aug 2012 - R. Yantosca - Move calls to findFreeLUN out of DEVEL block
!  24 Jun 2014 - R. Yantosca - Now pass FILENAME via the argument list

!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: IOS
      INTEGER            :: IU_FILE
      CHARACTER(LEN=255) :: MSG
      CHARACTER(LEN= 63), PARAMETER ::
     &   LOC = 'READ_PARANOX_LUT (paranox_util_mod.F)'

      !=================================================================
      ! READ_PARANOX_LUT begins here
      !=================================================================

      !=================================================================
      ! Read look up table for fraction of NOx remaining for ship
      ! emissions [unitless]
      !=================================================================

      ! Echo info
      WRITE( MSG, 100 ) TRIM( FracNOx_FILE )
 100  FORMAT( 'READ_PARANOX_LUT: Reading ', a )
      CALL HCO_MSG( MSG )
 
      ! Find a free file LUN
      IU_FILE = findFreeLUN()

      ! Open file to read  
      OPEN( IU_FILE, FILE=TRIM( FracNOx_FILE ), FORM="binary", 
     &      IOSTAT=IOS )
      IF ( IOS /= 0 ) THEN 
         CALL HCO_ERROR( 'read_paranox_lut:1', RC, THISLOC=LOC )
         RETURN
      ENDIF

      ! Read file
      READ( IU_FILE, IOSTAT=IOS ) FRACNOX
      IF ( IOS /= 0 ) THEN 
         CALL HCO_ERROR( 'read_paranox_lut:2', RC, THISLOC=LOC )
         RETURN
      ENDIF

      ! Close file
      CLOSE( IU_FILE )

      !=================================================================
      ! Read look up table for integrated Ozone Production Efficiency
      ! for ship emissions [molec O3 produced / molec NOx lost]
      !=================================================================

      ! Echo info
      WRITE( MSG, 101 ) TRIM( IntOPE_FILE )
 101  FORMAT( 'READ_PARANOX_LUT: Reading ', a )
      CALL HCO_MSG(MSG)

      ! Find a free file LUN
      IU_FILE = findFreeLUN()

      ! Open file to read
      OPEN( IU_FILE, FILE=TRIM( IntOPE_FILE ), FORM="BINARY", 
     &      IOSTAT=IOS )
      IF ( IOS /= 0 ) THEN 
         CALL HCO_ERROR( 'read_paranox_lut:3', RC, THISLOC=LOC )
         RETURN
      ENDIF

      ! Read file
      READ( IU_FILE, IOSTAT=IOS ) INTOPE
      IF ( IOS /= 0 ) THEN 
         CALL HCO_ERROR( 'read_paranox_lut:4', RC, THISLOC=LOC )
         RETURN
      ENDIF

      ! Close file
      CLOSE( IU_FILE )

      ! Return w/ success
      RC = HCO_SUCCESS

      END SUBROUTINE READ_PARANOX_LUT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Interpolate_Iut2
!
! !DESCRIPTION:  Subroutine Interpolate\_Lut2 returns FracNOx or IntOPE from
! the lookup tables (G.C.M. Vinken, KNMI, June 2010)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INTERPOLATE_LUT2( HcoState, 
     &                             I, J, o3, no, no2, dens, JO1D, JNO2,
     &                             TS, SUNCOSmid5, SUNCOSmid,
     &                             fraction_nox, int_ope, RC )
!
! !USES:
!
      USE HCO_STATE_MOD,   ONLY : HCO_STATE
!
! !INPUT PARAMETERS: 
!
      TYPE(HCO_STATE), POINTER     :: HcoState
      INTEGER,         INTENT(IN)  :: I, J
      REAL*8,          INTENT(IN)  :: o3, no, no2, dens, JNO2, JO1D
      REAL*8,          INTENT(IN)  :: TS, SUNCOSmid5, SUNCOSmid 
!
! !OUTPUT PARAMETERS: 
!
      REAL*4,          INTENT(OUT) :: fraction_nox, int_ope
!
! !INPUT/OUTPUT PARAMETERS: 
!
      INTEGER,          INTENT(INOUT) :: RC 
!
! !REVISION HISTORY:
!     Jun 2010 - G.C.M. Vinken - Initial version
!  21 Feb 2011 - G.C.M. Vinken - Updated for NOx in LUT
!  06 Feb 2012 - M. Payer      - Moved from emissions_mod.F to paranox_mod.F;
!                                Added ProTeX headers
!  15 Feb 2012 - M. Payer      - Add error trap to ensure 0 < fracnox < 1.
!  09 Nov 2012 - M. Payer      - Replaced all met field arrays with State_Met
!                                derived type object
!  28 Nov 2012 - R. Yantosca   - Replace SUNCOS_MID w/ State_Met%SUNCOSmid
!  28 Nov 2012 - R. Yantosca   - Replace SUNCOS_MID_5hr w/ State_Met%SUNCOSmid5
!  17 Jun 2013 - R. Yantosca   - Bug fix: declare all REAL variables with
!                                REAL*4 in order to avoid numerical precision
!                                errors when compiling with OMP=yes.
!  29 May 2013 - R. Yantosca   - Avoid numerical instability when computing
!                                the ratio JO1D/JNO2.  J-values go to zero
!                                at night, which can lead to a div-by-zero.
!  28 Jul 2014 - C. Keller     - Remove calls to GC routines.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

      !=======================================================================
      ! temp   : model temperature
      ! jno2   : J(NO2) value
      ! cao3   : concentration O3 in ambient air
      ! alfa0  : solar zenith angle 5 hours ago
      ! alfa5  : solar zenith angle at this time
      ! jo1d   : ratio J(O1D)/J(NO2)
      ! canox  : concentration NOx in ambient air 
      ! 
      ! o3     : incoming o3 concentration
      ! no     : incoming no
      ! no2    : incoming no2
      ! dens   : incoming air density
      !=======================================================================

      INTEGER                    :: IJLOOP
      INTEGER, PARAMETER         :: ntemp  = 4
      INTEGER, PARAMETER         :: njno2  = 4
      INTEGER, PARAMETER         :: ncao3  = 4
      INTEGER, PARAMETER         :: nalfa0 = 12
      INTEGER, PARAMETER         :: nalfa5 = 12
      INTEGER, PARAMETER         :: njo1d  = 4
      INTEGER, PARAMETER         :: ncanox = 5 
 
      REAL*4,  DIMENSION(ntemp)  :: templev
      REAL*4,  DIMENSION(njno2)  :: jno2lev
      REAL*4,  DIMENSION(ncao3)  :: cao3lev
      REAL*4,  DIMENSION(nalfa0) :: alfa0lev
      REAL*4,  DIMENSION(nalfa5) :: alfa5lev      
      REAL*4,  DIMENSION(njo1d)  :: jo1dlev
      REAL*4,  DIMENSION(ncanox) :: canoxlev

      ! Temporary variable storage
      REAL*4                     :: temp_tmp,  jno2_tmp,  cao3_tmp
      REAL*4                     :: alfa0_tmp, alfa5_tmp, jo1d_tmp
      REAL*4                     :: canox_tmp
      
      ! Interpolation parameters
      REAL*4,  DIMENSION(2)      :: xtemp,     xjno2,     xcao3, xalfa0
      REAL*4,  DIMENSION(2)      :: xalfa5,    xjo1d,     xcanox
      
      ! For loops
      INTEGER                    :: itemp,     ijno2,     icao3, ialfa0
      INTEGER                    :: ialfa5,    ijo1d,     icanox
      INTEGER                    :: i0,i1,i2,i3,i4,i5,i6,i7
 
      ! array contain temp, jno2, cao3, alfa_0, alfa_5, jo1d, canox
      REAL*4,  DIMENSION(7)      :: var_array

      CHARACTER(LEN=255)         :: MSG
      CHARACTER(LEN= 63), PARAMETER ::
     &   LOC = 'INTERPOLATE_LUT2 (paranox_util_mod.F)'

      !=================================================================
      ! INTERPOLATE_LUT2 begins here!
      !=================================================================

      ! Set the levels that were chosen in the look up table
      templev  = (/ 275.,  280.,   285.,   310.         /)
      jno2lev  = (/ 5.e-4, 0.0025, 0.0050, 0.012        /)
      cao3lev  = (/   5.,   20.,    35.,    75.         /)
      alfa0lev = (/ -90.,  -60.,   -45.,   -30., 
     $              -15.,    0.,    15.,    30.,
     $               45.,   60.,    75.,    90.         /)
      alfa5lev = (/ -90.,  -60.,   -45.,   -30., 
     $              -15.,    0.,    15.,    30.,
     $               45.,   60.,    75.,    90.         /)
      jo1dlev  = (/ 5.e-4, 0.0015, 0.0025, 0.0055       /)
      canoxlev = (/  10.,  200.,   1000.,  2000., 6000. /) 

!      PRINT*,"Temperature levels are: ",templev
!      PRINT*,"This is grid cell: ",I,J
      
      ! Temperature
!      PRINT*,"Temperature here is: ",State_Met%TS(I,J)
!      PRINT*,"USA: ",State_Met%TS(32,64)

      ! Tracer concentrations in v/v
!      PRINT*,"[O3] is: ",STT(I,J,1,IDTO3)/ State_Met%AD(I,J,1) * TCVV(IDTO3)
!      PRINT*,"[CO] is: ",STT(I,J,1,IDTCO)/ State_Met%AD(I,J,1) * TCVV(IDTCO)
!      PRINT*,"IDTO3 is: ", IDTO3
!      PRINT*,"IDO3 is: ", IDO3
!      PRINT*,"In USA: ",STT(32,64,1,IDTO3)/State_Met%AD(32,64,1) * TCVV(IDTO3)
      
      ! SOLAR ZENITH ANGLES IN DEGREES
!      IJLOOP = ( (J-1) * IIPAR ) + I
!      PRINT*,"Local Time: ",GET_LOCALTIME(I,1,1)
!      PRINT*,"Solar Zenith Angle at this location: ", 
!     $            ASIND(State_Met%SUNCOSmid(I,J))
!      IJLOOP = ( (64-1) * IIPAR ) + 32
!      PRINT*,"Local USA time: ", GET_LOCALTIME(32,1,1)
!      PRINT*,"Solar Zenith Angle at USA: ", 
!     &  ASIND(State_Met%SUNCOSmid(I,J))
!      PRINT*,"Solar Zenith Angle at USA - 5: ",
!     &  ASIND(State_Met%SUNCOSmid5(I,J))
      
      ! Set the variables
      IJLOOP       = ( (J-1) * HcoState%NX ) + I
      var_array(1) = TS                               ! Temperature
      var_array(2) = JNO2                             ! J(NO2), 1/s
      var_array(3) = o3 / dens * 1.E9                 ! [O3] in ppbv
      var_array(4) = ASIND(SUNCOSmid5)                ! alfa0
      var_array(5) = ASIND(SUNCOSmid)                 ! alfa5
      !--------------------------------------------------------------------
      ! Bug fix!  Make sure JNO2 is not zero so that it is safe to divide
      ! This was found by compiling with FPE=yes. (bmy, 8/23/13)
      !var_array(6) = JO1D / JNO2                     ! J(O1D)/J(NO2)
      !--------------------------------------------------------------------
      ! Avoid call to SAFE_DIV to make this module independent from GC
      ! (ckeller, 7/28/14)
      !var_array(6) = SAFE_DIV( JO1D, JNO2, 0d0 )      ! J(O1D)/J(NO2)
      !--------------------------------------------------------------------
      IF ( JNO2 > 0d0 ) THEN
         var_array(6) = JO1D / JNO2
      ELSE
         var_array(6) = 0d0
      ENDIF
      var_array(7) = (no + no2) / dens * 1.E12        ! [NOx] in pptv
      
      ! prevent NaN when jvalues are 0.
      IF (JNO2 .eq. 0.) var_array(6) = 0. 

      ! First some error checking
      ! ########### MAYBE CHECK HERE FOR NEGATIVE VALUES?##########
 
      !
      ! Determine reference index ( itemp,  ijno2, icao3, ialfa0,
      !                             ialfa5, ialfa, ijo1d, icaco )
      !
      !========================================================================
      ! Find smallest temperature reference level (i) for which actual
      ! temperature is smaller, then do
      !
      ! x(1) = ( temperature_level(i+1) - actual temperature   )
      !        ------------------------------------------------- 
      !        ( temperature_level(i+1) - temperature_level(i) )
      !
      ! then x(2) = 1.0 - x(1)
      !
      !========================================================================

      !---------------------
      ! Temperature:
      !---------------------
      temp_tmp = var_array(1)

      ! If temperature larger than largest in LUT, assign largest temp
      IF ( var_array(1) > templev( ntemp ) ) temp_tmp = templev(ntemp)

      ! If temp smaller, assign smallest temp level
      IF ( var_array(1) < templev(1) )       temp_tmp = templev(1)

      DO i0=1,ntemp-1
         itemp = i0
         IF( templev( itemp+1 ) > temp_tmp ) EXIT 
      END DO

      xtemp(1) = ( templev( itemp+1 ) - temp_tmp       ) /
     $           ( templev( itemp+1 ) - templev( itemp ) )
      xtemp(2) = 1.0 - xtemp(1)       

      !---------------------
      ! J(NO2):
      !---------------------
      jno2_tmp = var_array(2)

      ! If larger than largest in LUT, assign largest level values
      IF ( var_array(2) > jno2lev( njno2 ) ) jno2_tmp = jno2lev(njno2)

      ! If smaller, assign smallest level value
      IF ( var_array(2) < jno2lev(1) )       jno2_tmp = jno2lev(1)

      DO i0=1,njno2-1
         ijno2 = i0
         IF( jno2lev( ijno2+1 ) > jno2_tmp ) EXIT 
      END DO

      xjno2(1) = ( jno2lev( ijno2+1 ) - jno2_tmp       ) /
     $           ( jno2lev( ijno2+1 ) - jno2lev( ijno2 ) )
      xjno2(2) = 1.0 - xjno2(1)       

      !---------------------
      ! [O3]:
      !---------------------
      cao3_tmp = var_array(3)

      ! If larger than largest in LUT, assign largest level values
      IF ( var_array(3) > cao3lev( ncao3 ) ) cao3_tmp = cao3lev(ncao3)

      ! If smaller, assign smallest level value
      IF ( var_array(3) < cao3lev(1) )       cao3_tmp = cao3lev(1)

      DO i0=1,ncao3-1
         icao3 = i0
         IF( cao3lev( icao3+1 ) > cao3_tmp ) EXIT 
      END DO

      xcao3(1) = ( cao3lev( icao3+1 ) - cao3_tmp       ) /
     $           ( cao3lev( icao3+1 ) - cao3lev( icao3 ) )
      xcao3(2) = 1.0 - xcao3(1)       

      !---------------------
      ! alfa0:
      !---------------------
      alfa0_tmp = var_array(4)

      ! If larger than largest in LUT, assign largest level values
      IF ( var_array(4) > alfa0lev( nalfa0 ) ) alfa0_tmp = 
     $                                         alfa0lev(nalfa0)

      ! If smaller, assign smallest level value
      IF ( var_array(4) < alfa0lev(1) )        alfa0_tmp = alfa0lev(1)

      DO i0=1,nalfa0-1
         ialfa0 = i0
         IF( alfa0lev( ialfa0+1 ) > alfa0_tmp ) EXIT 
      END DO

      xalfa0(1) = ( alfa0lev( ialfa0+1 ) - alfa0_tmp        ) /
     $            ( alfa0lev( ialfa0+1 ) - alfa0lev( ialfa0 ) )
      xalfa0(2) = 1.0 - xalfa0(1)       

      !---------------------
      ! alfa5:
      !---------------------
      alfa5_tmp = var_array(5)

      ! If larger than largest in LUT, assign largest level values
      IF ( var_array(5) > alfa5lev( nalfa5 ) ) alfa5_tmp = 
     $                                         alfa5lev(nalfa5)

      ! If smaller, assign smallest level value
      IF ( var_array(5) < alfa5lev(1) )        alfa5_tmp = alfa5lev(1)

      DO i0=1,nalfa5-1
         ialfa5 = i0
         IF( alfa5lev( ialfa5+1 ) > alfa5_tmp ) EXIT 
      END DO

      xalfa5(1) = ( alfa5lev( ialfa5+1 ) - alfa5_tmp        ) /
     $            ( alfa5lev( ialfa5+1 ) - alfa5lev( ialfa5 ) )
      xalfa5(2) = 1.0 - xalfa5(1)       

      !---------------------
      ! jo1d:
      !---------------------
      jo1d_tmp = var_array(6)

      ! If larger than largest in LUT, assign largest level values
      IF ( var_array(6) > jo1dlev( njo1d ) ) jo1d_tmp = jo1dlev(njo1d)

      ! If smaller, assign smallest level value
      IF ( var_array(6) < jo1dlev(1) )       jo1d_tmp = jo1dlev(1)

      DO i0=1,njo1d-1
         ijo1d = i0
         IF( jo1dlev( ijo1d+1 ) > jo1d_tmp ) EXIT 
      END DO

      xjo1d(1) = ( jo1dlev( ijo1d+1 ) - jo1d_tmp       ) /
     $           ( jo1dlev( ijo1d+1 ) - jo1dlev( ijo1d ) )
      xjo1d(2) = 1.0 - xjo1d(1)       

      !---------------------
      ! [NOx]:
      !---------------------
      canox_tmp = var_array(7)

      ! If larger than largest in LUT, assign largest level values
      IF ( var_array(7) > canoxlev( ncanox ) ) canox_tmp = 
     $                                         canoxlev(ncanox)

      ! If smaller, assign smallest level value
      IF ( var_array(7) < canoxlev(1) )        canox_tmp = canoxlev(1)

      DO i0=1,ncanox-1
         icanox = i0
         IF( canoxlev( icanox+1 ) > canox_tmp ) EXIT 
      END DO

      xcanox(1) = ( canoxlev( icanox+1 ) - canox_tmp        ) /
     $            ( canoxlev( icanox+1 ) - canoxlev( icanox ) )
      xcanox(2) = 1.0 - xcanox(1)       

!      PRINT*,"The i-values are:",    itemp, ijno2, icao3, ialfa0,
!     $                               ialfa5, ijo1d, icanox
!      PRINT*,"Variables are: ",      var_array
!      PRINT*,"For testing, xtemp: ", xtemp
      
      !======================
      ! Linear interpolation
      !======================

      fraction_nox = 0.0
      int_ope      = 0.0

      DO i1=1,2
      DO i2=1,2
      DO i3=1,2
      DO i4=1,2
      DO i5=1,2
      DO i6=1,2
      DO i7=1,2 

         !IF ENCOUNTER -999 IN THE LUT PRINT ERROR!!       
         IF ( ( fracnox( itemp+i1-1,  ijno2+i2-1,  icao3+i3-1,
     $                   ialfa0+i4-1, ialfa5+i5-1, ijo1d+i6-1,
     $                   icanox+i7-1  ) < 0. )     .or.
     $        ( fracnox( itemp+i1-1,  ijno2+i2-1,  icao3+i3-1,
     $                   ialfa0+i4-1, ialfa5+i5-1, ijo1d+i6-1, 
     $                   icanox+i7-1  ) > 1. ) )   THEN

            WRITE(MSG,*), 'INTERPOLATE_LUT2: fracnox = ,',
     $          fracnox( itemp+i1-1,  ijno2+i2-1,  icao3+i3-1,
     $                   ialfa0+i4-1, ialfa5+i5-1, ijo1d+i6-1,
     $                   icanox+i7-1  )
            CALL HCO_MSG(MSG)

            MSG = 'LUT error: Fracnox should be between 0 and 1!'
            CALL HCO_ERROR( MSG, RC, THISLOC=LOC ) 
            RETURN 
         ENDIF

         ! fracnox is the array with the actual lut data
         fraction_nox = fraction_nox + xtemp(i1)    * xjno2(i2)  *
     $                  xcao3(i3)    * xalfa0(i4)   * xalfa5(i5) *
     $                  xjo1d(i6)    * xcanox(i7)   *
     $                  fracnox(       itemp+i1-1,    ijno2+i2-1,
     $                  icao3+i3-1,    ialfa0+i4-1,   ialfa5+i5-1,
     $                  ijo1d+i6-1,    icanox+i7-1    )

         ! intope is the array with the actual lut data
         int_ope =      int_ope      + xtemp(i1)    * xjno2(i2)  * 
     $                  xcao3(i3)    * xalfa0(i4)   * xalfa5(i5) *
     $                  xjo1d(i6)    * xcanox(i7)   * 
     $                  intope(        itemp+i1-1,    ijno2+i2-1,
     $                  icao3+i3-1,    ialfa0+i4-1,   ialfa5+i5-1,
     $                  ijo1d+i6-1,    icanox+i7-1    )

      END DO
      END DO
      END DO
      END DO
      END DO
      END DO
      END DO
      
!      IF ((I .eq. 108) .and. (J .eq. 49)) THEN
!         PRINT*,"-----INTERPOLATE_LUT2, for 108,49-----"
!         PRINT*,"Fraction_nox and int_OPE: ", fraction_nox,
!     &                                        int_ope
!         PRINT*,"Jvalues are: ",              JNO2, JO1D
!         PRINT*,"Vars are: ",                 var_array
!         PRINT*,"[O3] in interpolate_lut: ",  var_array(3)
!         PRINT*,"[NOx] in interpolate_lut: ", var_array(7)
!         PRINT*,"J(O1D)/J(NO2) : ",           var_array(6)
!         PRINT*,"The i-values are:",          itemp,  ijno2,  icao3,
!     $                                        ialfa0, ialfa5, ijo1d,
!     $                                        icanox
!         PRINT*,"Interpolation parameters: ", xtemp,  xjno2,  xcao3,
!     $                                        xalfa0, xalfa5, xjo1d,
!     $                                        xcanox 
!         PRINT*,"---------------------------------"
!      ENDIF
!       
!      IF ((I .eq. 73) .and. (J .eq. 76)) THEN
!         PRINT*,"-----INTERPOLATE_LUT2, for 73,76-----"
!         PRINT*,"Fraction_nox and int_OPE: ", fraction_nox,
!     &                                        int_ope
!         PRINT*,"Jvalues are: ",              JNO2, JO1D
!         PRINT*,"Vars are: ",                 var_array
!         PRINT*,"[O3] in interpolate_lut: ",  var_array(3)
!         PRINT*,"[NOx] in interpolate_lut: ", var_array(7)
!         PRINT*,"J(O1D)/J(NO2) : ",           var_array(6)
!         PRINT*,"The i-values are:",          itemp,  ijno2,  icao3,
!     $                                        ialfa0, ialfa5, ijo1d,
!     $                                        icanox
!         PRINT*,"Interpolation parameters: ", xtemp,  xjno2,  xcao3,
!     $                                        xalfa0, xalfa5, xjo1d,
!     $                                        xcanox 
!         PRINT*,"------------------------------------"
!      ENDIF
     
      ! Return w/ success
      RC = HCO_SUCCESS
 
      END SUBROUTINE INTERPOLATE_LUT2
!EOC
      END MODULE PARANOX_UTIL_MOD
