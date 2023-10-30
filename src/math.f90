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
  private

  !---------------------------------------------------------------------------!
  !                  Interfaces for overloaded routines                       !
  !---------------------------------------------------------------------------!
  interface math_isclose
     module procedure math_isclose_cmplx
     module procedure math_aisclose_cmplx
     module procedure math_isclose_real
     module procedure math_aisclose_real
  end interface math_isclose

  !---------------------------------------------------------------------------!
  !                       P u b l i c   R o u t i n e s                       !
  !---------------------------------------------------------------------------!
  public :: math_isclose
  public :: math_get_vec_angle
contains

  logical function math_isclose_cmplx(a,b,atol,rtol)
    !============================================================!
    ! Checks if two complex values are equal within              !
    ! a given tolerance. The formula used to check this is       !
    ! abs(a-b) <= (atol + rtol*abs(b))                           !
    ! Note this comparison is not symmetric and it is assumed    !
    ! that b is the reference value.                             !
    !------------------------------------------------------------!
    ! Arguments                                                  !
    ! a(in) :: the value to compare                              !
    ! b(in) :: the reference value                               !
    ! atol  :: absolute tolerance (Default: 1e-6)                !
    ! btol  :: relative tolerance (Default: 1e-5)                !
    !============================================================!
    implicit none
    complex(kind=dp) :: a,b
    real(kind=dp),optional :: atol,rtol

    real(kind=dp) :: l_atol, l_rtol

    l_atol = 1e-6_dp
    l_rtol = 1e-5_dp
    if(present(atol)) l_atol=atol
    if(present(rtol)) l_rtol=rtol

    math_isclose_cmplx = abs(a-b)<= (l_atol + l_rtol*abs(b))
    return
  end function math_isclose_cmplx

  logical function math_aisclose_cmplx(a,b,atol,rtol)
    !============================================================!
    ! Vectorised version math_isclose_val.                       !
    ! Compares if two complex arraysare element-wise equal       !
    ! within tolerance.                                          !
    ! NOTE: Array must be flattened before calling routine       !
    ! if dim>1                                                   !
    !------------------------------------------------------------!
    ! Arguments                                                  !
    ! a(in) :: the array to compare                              !
    ! b(in) :: the reference array                               !
    ! atol  :: absolute tolerance (Default: 1e-6)                !
    ! btol  :: relative tolerance (Default: 1e-5)                !
    !------------------------------------------------------------!
    ! Necesary Conditions                                        !
    ! a and b have the same size and must be flattened.          !
    !============================================================!
    implicit none
    complex(kind=dp) :: a(:),b(:)
    real(kind=dp),optional :: atol,rtol

    real(kind=dp) :: l_atol, l_rtol
    integer :: i

    l_atol = 1e-6_dp
    l_rtol = 1e-5_dp
    if(present(atol)) l_atol=atol
    if(present(rtol)) l_rtol=rtol

    ! Assume true unless told otherwise
    math_aisclose_cmplx=.true.
    do i=1,size(a)
       if (.not.math_isclose_cmplx(a(i),b(i),atol=l_atol,rtol=l_rtol)) then
          math_aisclose_cmplx=.false.
          ! Can return here as there is no need to check the rest if not equal.
          return
       end if
    end do
  end function math_aisclose_cmplx

  logical function math_isclose_real(a,b,atol,rtol)
    !============================================================!
    ! Checks if two real values are equal within                 !
    ! a given tolerance. The formula used to check this is       !
    ! abs(a-b) <= (atol + rtol*abs(b))                           !
    ! Note this comparison is not symmetric and it is assumed    !
    ! that b is the reference value.                             !
    !------------------------------------------------------------!
    ! Arguments                                                  !
    ! a(in) :: the value to compare                              !
    ! b(in) :: the reference value                               !
    ! atol  :: absolute tolerance (Default: 1e-6)                !
    ! btol  :: relative tolerance (Default: 1e-5)                !
    !============================================================!
    implicit none
    real(kind=dp) :: a,b
    real(kind=dp),optional :: atol,rtol

    real(kind=dp) :: l_atol, l_rtol

    l_atol = 1e-6_dp
    l_rtol = 1e-5_dp
    if(present(atol)) l_atol=atol
    if(present(rtol)) l_rtol=rtol

    math_isclose_real = abs(a-b)<= (l_atol + l_rtol*abs(b))
    return
  end function math_isclose_real

  logical function math_aisclose_real(a,b,atol,rtol)
    !============================================================!
    ! Vectorised version math_isclose_real_val.                  !
    ! Compares if two real arrays are element-wise equal         !
    ! within tolerance.                                          !
    ! NOTE: Array must be flattened before calling routine       !
    ! if dim>1                                                   !
    !------------------------------------------------------------!
    ! Arguments                                                  !
    ! a(in) :: the array to compare                              !
    ! b(in) :: the reference array                               !
    ! atol  :: absolute tolerance (Default: 1e-6)                !
    ! btol  :: relative tolerance (Default: 1e-5)                !
    !------------------------------------------------------------!
    ! Necesary Conditions                                        !
    ! a and b have the same size and must be flattened.          !
    !============================================================!
    implicit none
    real(kind=dp) :: a(:),b(:)
    real(kind=dp),optional :: atol,rtol

    real(kind=dp) :: l_atol, l_rtol
    integer :: i

    l_atol = 1e-6_dp
    l_rtol = 1e-5_dp
    if(present(atol)) l_atol=atol
    if(present(rtol)) l_rtol=rtol

    ! Assume true unless told otherwise
    math_aisclose_real=.true.
    do i=1,size(a)
       if (.not.math_isclose_real(a(i),b(i),atol=l_atol,rtol=l_rtol)) then
          math_aisclose_real=.false.
          ! Can return here as there is no need to check the rest if not equal.
          return
       end if
    end do
  end function math_aisclose_real

  function math_get_vec_angle(a,b) result(theta)
    !============================================================!
    ! Calculate the angle theta between two vectors in degrees   !
    ! noting that a dot b = |a||b|cos(theta)                     !
    !------------------------------------------------------------!
    ! Arguments                                                  !
    ! a(in) :: the first vector                                  !
    ! b(in) :: the second vector                                 !
    !------------------------------------------------------------!
    ! Necesary Conditions                                        !
    ! a and b must be arrays with dimension(3).                  !
    !============================================================!
    use constants, only : pi

    implicit none
    real(kind=dp),dimension(:) :: a,b
    real(kind=dp) :: theta
    real(kind=dp) :: norm_a,norm_b

    ! Calculate the norm of the two vectors
    norm_a = sqrt(sum(a**2))
    norm_b = sqrt(sum(b**2))

    theta = acos(dot_product(a,b)/(norm_a*norm_b))
    theta = theta * 180.0_dp/pi
  end function math_get_vec_angle
end module math
