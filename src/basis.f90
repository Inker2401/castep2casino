module basis
  !===============================================================================!
  !                           ____            _                                   !
  !                          | __ )  __ _ ___(_)___                               !
  !                          |  _ \ / _` / __| / __|                              !
  !                          | |_) | (_| \__ \ \__ \                              !
  !                          |____/ \__,_|___/_|___/                              !
  !                                                                               !
  ! This module contains routines related to the management of the grid           !
  ! All  routines related to FFTs are located here as well.                       !
  ! ------------------------------------------------------------------------------!
  ! Modules used                                                                  !
  ! Constants, IO                                                                 !
  ! ------------------------------------------------------------------------------!
  ! Public variables:                                                             !
  ! basis_data: derived type for the user's basis                                 !
  ! castep_data: public 'instance' of basis_data for CASTEP grid                  !
  !===============================================================================!
  use constants,only : dp,int_dp
  private

  !---------------------------------------------------------------------------!
  !                       P u b l i c   V a r i a b l e s                     !
  !---------------------------------------------------------------------------!
  type,public :: basis_data
     integer                      :: ngx,ngy,ngz          ! total number of grid points along x,y,z
     integer                      :: total_grid_points    ! total number of points in grid
     integer(kind=int_dp)         :: plan_for             ! forward FFT plan
     integer(kind=int_dp)         :: plan_back            ! backward FFT plan
  end type basis_data

  type(basis_data),public,save::castep_basis              ! the CASTEP grid

  !---------------------------------------------------------------------------!
  !                       P u b l i c   R o u t i n e s                       !
  !---------------------------------------------------------------------------!

  interface basis_shift
     module procedure basis_shift_cmplx
     module procedure basis_shift_real
  end interface basis_shift

  public :: basis_initialise
  public :: basis_deallocate
  public :: basis_real_to_recip
  public :: basis_recip_to_real
  public :: basis_shift

  !---------------------------------------------------------------------------!
  !                     P r i v a t e   V a r i a b l e s                     !
  !---------------------------------------------------------------------------!
  integer,parameter :: FFTW_FORWARD=-1
  integer,parameter :: FFTW_BACKWARD=1
  integer,parameter :: FFTW_ESTIMATE=0

contains

  subroutine basis_initialise(ngx,ngy,ngz)
    !============================================================!
    ! This routine initialises the grid used for CASTEP as well  !
    ! and initialises all FFT routines by creating the necessary !
    ! plans.                                                     !
    ! This routine should be called BEFORE using any other       !
    ! routine in this module.                                    !
    !------------------------------------------------------------!
    ! Arguments                                                  !
    ! ngx,ngy,ngz(in) :: number of grid points along each        !
    !                    grid dimension                          !
    !------------------------------------------------------------!
    ! Necessary Conditions                                       !
    ! None                                                       !
    !============================================================!
    implicit none
    integer,intent(in) :: ngx,ngy,ngz
    complex(kind=dp),allocatable :: grid(:,:,:)
    integer :: stat

    ! Set the correct attributes
    castep_basis%ngx=ngx
    castep_basis%ngy=ngy
    castep_basis%ngz=ngz
    castep_basis%total_grid_points = ngx*ngy*ngz

    ! Allocate the grid for basis
    allocate(grid(ngx,ngy,ngz),stat=stat)
    if(stat/=0) error stop 'basis_initialise: Failed to allocate grid.'

    ! Create the FFT plans
    call dfftw_plan_dft_3d(castep_basis%plan_for,ngx,ngy,ngz,grid,grid,FFTW_FORWARD,FFTW_ESTIMATE)
    call dfftw_plan_dft_3d(castep_basis%plan_back,ngx,ngy,ngz,grid,grid,FFTW_BACKWARD,FFTW_ESTIMATE)

  end subroutine basis_initialise

  subroutine basis_deallocate()
    !============================================================!
    ! This routine deallocates the CASTEP basis and destroys     !
    ! any FFT plans. This should only be called at the end of    !
    ! the program, or at the very least after the FFTs are done. !
    !============================================================!

    implicit none

    ! Set all other variables to 0
    castep_basis%ngx = 0
    castep_basis%ngy = 0
    castep_basis%ngz = 0
    castep_basis%total_grid_points = 0

    ! Destroy FFT plans
    call dfftw_destroy_plan(castep_basis%plan_for)
    call dfftw_destroy_plan(castep_basis%plan_back)

  end subroutine basis_deallocate

  subroutine basis_real_to_recip(grid)
    !============================================================!
    ! This routines transforms a quantity on a grid from real    !
    ! space to reciprocal space.                                 !
    ! Note that we must divide by the total number of grid points!
    ! on the transformed grid as FFT multiples by the grid points!
    ! when doing each transform.                                 !
    !------------------------------------------------------------!
    ! Arguments                                                  !
    ! grid(inout) :: the grid to transform to reciprocal space.  !
    !------------------------------------------------------------!
    ! Necessary Conditions                                       !
    ! basis_initialise must have been called.                    !
    !============================================================!
    implicit  none
    complex(kind=dp),intent(inout) :: grid(:,:,:)

    call dfftw_execute_dft(castep_basis%plan_for,grid,grid)

    ! Normalise reciprocal space grid
    grid = grid/castep_basis%total_grid_points
  end subroutine basis_real_to_recip

  subroutine basis_recip_to_real(grid)
    !============================================================!
    ! This routines transforms a quantity on a grid from         !
    ! reciprocal space to reciprocal space.                      !
    !------------------------------------------------------------!
    ! Arguments                                                  !
    ! grid(inout) :: the grid to transform to reciprocal space.  !
    !------------------------------------------------------------!
    ! Necessary Conditions                                       !
    ! basis_initialise must have been called.                    !
    !============================================================!
    implicit  none
    complex(kind=dp),intent(inout) :: grid(:,:,:)

    call dfftw_execute_dft(castep_basis%plan_back,grid,grid)
  end subroutine basis_recip_to_real

  subroutine basis_shift_cmplx(grid)
    !============================================================!
    ! This routine shifts a real space grid by an amount         !
    ! specified in fractional coordinates.                       !
    ! This routine is for real space grids/fields that are       !
    ! complex-valued.                                            !
    !------------------------------------------------------------!
    ! Arguments                                                  !
    ! grid(inout) :: the grid to be shifted                      !
    !------------------------------------------------------------!
    ! Necessary Conditions                                       !
    ! basis_initialise must have been called.                    !
    !============================================================!
    use latt,only  : user_params

    implicit none
    complex(kind=dp),intent(inout) :: grid(:,:,:)
    integer :: xshift,yshift,zshift

    xshift = nint(castep_basis%ngx*user_params%shift_frac(1))
    yshift = nint(castep_basis%ngy*user_params%shift_frac(2))
    zshift = nint(castep_basis%ngz*user_params%shift_frac(3))

    ! write(*,*) castep_basis%ngx*user_params%shift_frac(1),castep_basis%ngy*user_params%shift_frac(2),&
    !      castep_basis%ngz*user_params%shift_frac(3)
    ! write(*,*) xshift,yshift,zshift

    grid = cshift(grid,xshift,1)
    grid = cshift(grid,yshift,2)
    grid = cshift(grid,zshift,3)
  end subroutine basis_shift_cmplx

  subroutine basis_shift_real(grid)
    !============================================================!
    ! This routine shifts a real space grid by an amount         !
    ! specified in fractional coordinates.                       !
    ! This routine is for real space grids/fields that are       !
    ! real-valued.                                               !
    !------------------------------------------------------------!
    ! Arguments                                                  !
    ! grid(inout) :: the grid to be shifted                      !
    !------------------------------------------------------------!
    ! Necessary Conditions                                       !
    ! basis_initialise must have been called.                    !
    !============================================================!
    use latt,only  : user_params

    implicit none
    real(kind=dp),intent(inout) :: grid(:,:,:)
    integer :: xshift,yshift,zshift

    xshift = nint(castep_basis%ngx*user_params%shift_frac(1))
    yshift = nint(castep_basis%ngy*user_params%shift_frac(2))
    zshift = nint(castep_basis%ngz*user_params%shift_frac(3))

    ! write(*,*) castep_basis%ngx*user_params%shift_frac(1),castep_basis%ngy*user_params%shift_frac(2),&
    !      castep_basis%ngz*user_params%shift_frac(3)
    ! write(*,*) xshift,yshift,zshift

    grid = cshift(grid,xshift,1)
    grid = cshift(grid,yshift,2)
    grid = cshift(grid,zshift,3)
  end subroutine basis_shift_real
end module basis
