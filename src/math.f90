module math
  !===============================================================================!
  !                         __  __       _   _                                    !
  !                        |  \/  | __ _| |_| |__  ___                            !
  !                        | |\/| |/ _` | __| '_ \/ __|                           !
  !                        | |  | | (_| | |_| | | \__ \                           !
  !                        |_|  |_|\__,_|\__|_| |_|___/                           !
  ! This modules contains some general useful functions that are not part of      !
  ! Fortran intrinsics.                                                           !
  ! ------------------------------------------------------------------------------!
  ! Modules used                                                                  !
  ! standard Fortran environment: constants                                       !
  !===============================================================================!
  use constants, only : dp

contains

  logical function math_isclose(a,b,atol,rtol)
    implicit none
    complex(kind=dp) :: a,b
    real(kind=dp),optional :: atol,rtol

    real(kind=dp) :: l_atol, l_rtol

    l_atol = 1e-6_dp
    l_rtol = 1e-5_dp

    if(present(atol)) l_atol=atol
    if(present(rtol)) l_rtol=rtol

    math_isclose = abs(a-b)<= (l_atol + l_rtol*abs(a))
    return
  end function math_isclose
end module math
