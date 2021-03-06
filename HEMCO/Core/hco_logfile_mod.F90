!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hco_logfile_mod
!
! !DESCRIPTION: Module HCO\_LOGFILE\_MOD contains some wrapper routines to 
! write data into the HEMCO logfile.
!\\
!\\
! !INTERFACE: 
!
MODULE HCO_LOGFILE_MOD 
! 
! !USES:
!
  USE HCO_ERROR_MOD

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: HCO_Spec2Log
  PUBLIC  :: HCO_PrintList
  PUBLIC  :: HCO_PrintDataCont
!
! !REVISION HISTORY:
!  27 May 2014 - C. Keller   - Initialization
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
! !IROUTINE: hco_spec2log
!
! !DESCRIPTION: Subroutine HCO\_Spec2Log writes information of a species
! to the logfile.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_Spec2Log( am_I_Root, HcoState, ID )
!
! !USES:
!
    USE HCO_STATE_MOD, ONLY : HCO_State
!
!
! !INPUT PARAMETER
!
    LOGICAL,          INTENT(IN)     :: am_I_Root  ! Root CPU
    TYPE(HCO_State),  POINTER        :: HcoState   ! HEMCO state object
    INTEGER,          INTENT(IN)     :: ID         ! HEMCO species ID
!
! !REVISION HISTORY: 
!  27 May 2014 - C. Keller   - Initialization
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL             :: verb
    CHARACTER(LEN=255)  :: MSG 

    !=================================================================
    ! HCO_Spec2Log begins here 
    !=================================================================

    ! Verbose mode?
    verb = HCO_VERBOSE_CHECK() .AND. am_I_Root

    MSG = 'Species ' // TRIM(HcoState%Spc(ID)%SpcName)
    CALL HCO_MSG(MSG)
    IF ( verb ) THEN
       write(MSG,*) '--> HcoID         : ', HcoState%Spc(ID)%HcoID
       CALL HCO_MSG(MSG)
       write(MSG,*) '--> ModID         : ', HcoState%Spc(ID)%ModID
       CALL HCO_MSG(MSG)
       write(MSG,*) '--> MW (g/mol)    : ', HcoState%Spc(ID)%MW_g
       CALL HCO_MSG(MSG)
       write(MSG,*) '--> emitted MW    : ', HcoState%Spc(ID)%EmMW_g
       CALL HCO_MSG(MSG)
       write(MSG,*) '--> Molecule ratio: ', HcoState%Spc(ID)%MolecRatio
       CALL HCO_MSG(MSG)
       write(MSG,*) '--> Henry constant: ', HcoState%Spc(ID)%HenryK0
       CALL HCO_MSG(MSG)
       write(MSG,*) '--> Henry temp.   : ', HcoState%Spc(ID)%HenryCR
       CALL HCO_MSG(MSG)
       write(MSG,*) '--> Henry pKA     : ', HcoState%Spc(ID)%HenryPKA
       CALL HCO_MSG(MSG)
    ENDIF    

  END SUBROUTINE HCO_Spec2Log
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: HCO_PrintList
!
! !DESCRIPTION: Subroutine HCO\_PrintList displays the content of List. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE HCO_PrintList ( List, Verbose )
!
! !USES:
!
      USE HCO_DATACONT_MOD,     ONLY : ListCont
!
! !INPUT ARGUMENTS:
!
      TYPE(ListCont), POINTER    :: List
      LOGICAL,        INTENT(IN) :: Verbose
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !ARGUMENTS:
!
      TYPE(ListCont), POINTER   :: TmpLct => NULL()
      CHARACTER(LEN=255)        :: MSG 

      ! ================================================================
      ! HCO_PrintList begins here
      ! ================================================================

      ! Point to first element
      TmpLct => List
      DO WHILE ( ASSOCIATED(TmpLct) ) 
         IF ( ASSOCIATED(TmpLct%Dct) ) THEN
            CALL HCO_PrintDataCont(TmpLct%Dct,Verbose)
         ENDIF
         TmpLct => TmpLct%NextCont
      ENDDO        

      TmpLct => NULL()

      END SUBROUTINE HCO_PrintList
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: HCO_PrintDataCont
!
! !DESCRIPTION: Subroutine HCO\_PrintDataCont displays the content of the
! data container Dct. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE HCO_PrintDataCont ( Dct, Verbose )
!
! !USES
!
      USE HCO_DATACONT_MOD,     ONLY : DataCont
!
! !INPUT ARGUMENTS:
!
      TYPE(DataCont), POINTER    :: Dct
      LOGICAL,        INTENT(IN) :: Verbose
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !ARGUMENTS:
!
      CHARACTER(LEN=255) :: MSG 
      INTEGER            :: nx, ny, nz, nt      
      REAL(sp)           :: sm

      ! ================================================================
      ! HCO_PrintDataCont begins here
      ! ================================================================

      sm = 0.0_sp 
      nx = 0 
      ny = 0
      nz = 0
      nt = Dct%Dta%nt
      IF ( nt > 0 ) THEN
         IF ( Dct%Dta%spaceDim<=2 ) THEN
            IF ( ASSOCIATED(Dct%Dta%V2) ) THEN
               nx = SIZE(Dct%Dta%V2(1)%Val,1)
               ny = SIZE(Dct%Dta%V2(1)%Val,2)
               sm = SUM(Dct%Dta%V2(1)%Val)
            ENDIF
         ELSE
            IF ( ASSOCIATED(Dct%Dta%V3) ) THEN
               nx = SIZE(Dct%Dta%V3(1)%Val,1)
               ny = SIZE(Dct%Dta%V3(1)%Val,2)
               nz = SIZE(Dct%Dta%V3(1)%Val,3)
               sm = SUM(Dct%Dta%V3(1)%Val)
            ENDIF
         ENDIF
      ENDIF

      ! Always print name 
      MSG = 'Container ' // TRIM(Dct%cName)
      CALL HCO_MSG(MSG)

      ! Eventually add details
      IF ( verbose ) THEN

         ! General information
         write(MSG,*) '   -->Data type       : ', Dct%DctType
         CALL HCO_MSG(MSG)
         write(MSG,*) '   -->Container ID    : ', Dct%cID
         CALL HCO_MSG(MSG)
         write(MSG,*) '   -->Target ID       : ', Dct%targetID
         CALL HCO_MSG(MSG)
         write(MSG,*) '   -->File data home?   ', Dct%DtaHome
         CALL HCO_MSG(MSG)
         write(MSG,*) '   -->Source file     : ', TRIM(Dct%Dta%ncFile)
         CALL HCO_MSG(MSG)
         write(MSG,*) '   -->ncRead?           ', Dct%Dta%ncRead
         CALL HCO_MSG(MSG)
         write(MSG,*) '   -->Shared data file? ', Dct%Dta%DoShare
         CALL HCO_MSG(MSG)
         IF ( Dct%Dta%ncRead ) THEN
            write(MSG,*) '   -->Source parameter: ', TRIM(Dct%Dta%ncPara)
            CALL HCO_MSG(MSG)
            write(MSG,*) '   -->Year range      : ', Dct%Dta%ncYrs 
            CALL HCO_MSG(MSG)
            write(MSG,*) '   -->Month range     : ', Dct%Dta%ncMts 
            CALL HCO_MSG(MSG)
            write(MSG,*) '   -->Day range       : ', Dct%Dta%ncDys 
            CALL HCO_MSG(MSG)
            write(MSG,*) '   -->Hour range      : ', Dct%Dta%ncHrs 
            CALL HCO_MSG(MSG)
            write(MSG,*) '   -->SpaceDim        : ', Dct%Dta%SpaceDim
            CALL HCO_MSG(MSG)
         ENDIF
         IF ( NZ > 0 ) THEN
            write(MSG,*) '   -->Array dimension : ', nx,ny,nz
         ELSE
            write(MSG,*) '   -->Array dimension : ', nx,ny
         ENDIF
         CALL HCO_MSG(MSG)
         write(MSG,*) '   -->Array sum       : ', sm 
         CALL HCO_MSG(MSG)
         write(MSG,*) '   -->Time dimension  : ', nt 
         CALL HCO_MSG(MSG)
         write(MSG,*) '   -->Delta t[h]      : ', Dct%Dta%DeltaT
         CALL HCO_MSG(MSG)
         IF ( ASSOCIATED(Dct%Dta%tIDx) ) THEN
            write(MSG,*) '   -->Tempres         : ', &
               TRIM(Dct%Dta%tIDx%TempRes)
            CALL HCO_MSG(MSG)
         ENDIF
         write(MSG,*) '   -->OrigUnit        : ',TRIM(Dct%Dta%OrigUnit)
         CALL HCO_MSG(MSG)
         write(MSG,*) '   -->Concentration?    ', Dct%Dta%IsConc
         CALL HCO_MSG(MSG)
         write(MSG,*) '   -->Coverage        : ', Dct%Dta%Cover
         CALL HCO_MSG(MSG)

         ! For base emissions
         IF ( Dct%DctType==1 ) THEN
            write(MSG,*) '   -->Extension Nr    : ', Dct%ExtNr
            CALL HCO_MSG(MSG)
            write(MSG,*) '   -->Species name    : ',TRIM(Dct%SpcName)
            CALL HCO_MSG(MSG)
            write(MSG,*) '   -->HEMCO species ID: ', Dct%HcoID
            CALL HCO_MSG(MSG)
            write(MSG,*) '   -->Category        : ', Dct%Cat
            CALL HCO_MSG(MSG)
            write(MSG,*) '   -->Hierarchy       : ', Dct%Hier
         CALL HCO_MSG(MSG)

         ! For scale factors
         ELSEIF ( Dct%DctType>1 ) THEN
            write(MSG,*) '   -->Scal ID         : ', Dct%ScalID
            CALL HCO_MSG(MSG)
            write(MSG,*) '   -->Operator        : ', Dct%Oper
            CALL HCO_MSG(MSG)
         ENDIF
      ENDIF

      END SUBROUTINE HCO_PrintDataCont
!EOC
END MODULE HCO_LOGFILE_MOD
