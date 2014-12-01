# 1 "main/mosaic_support.f90"
# 1 "<built-in>"
# 1 "<command line>"
# 1 "main/mosaic_support.f90"

module mosaic_support
  !Purpose: This module contains subroutines which have codes which depend upon
  ! the host code (CAM, WRF etc.). #defines are used to seprate codes
  ! which depends on the host code






  implicit none
  private

  public:: mosaic_warn_mess
  public:: mosaic_err_mess


contains

  subroutine mosaic_warn_mess(message)
    !Purpose: Print out the warning messages from Mosaic code

    character(len=*), intent(in) :: message

    !Local variables
    character(len=16), parameter :: warn_str = 'MOSAIC WARNING: '






    write(*,*)warn_str,message




  end subroutine mosaic_warn_mess


  subroutine mosaic_err_mess(message)
    !Purpose
    character(len=*), intent(in) :: message

    !Local variables
    character(len=14), parameter :: err_str = 'MOSAIC ERROR: '
    character(len=500) :: str_to_prnt

    write(str_to_prnt,*)err_str,message






    write(*,*)(trim(adjustl(str_to_prnt)))
    stop

  end subroutine mosaic_err_mess



end module mosaic_support
