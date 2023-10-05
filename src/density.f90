module density
  !===============================================================================!
  !                      ____                 _ _                                 !
  !                     |  _ \  ___ _ __  ___(_) |_ _   _                         !
  !                     | | | |/ _ \ '_ \/ __| | __| | | |                        !
  !                     | |_| |  __/ | | \__ \ | |_| |_| |                        !
  !                     |____/ \___|_| |_|___/_|\__|\__, |                        !
  !                                                 |___/                         !
  ! This module contains the routines for dealing with the properties of densities!
  ! in a CASTEP format                                                            !
  ! Note that basis_initialise should be called before using any routine within   !
  ! this module.                                                                  !
  ! ------------------------------------------------------------------------------!
  ! Modules used                                                                  !
  ! Constants,IO,Basis                                                            !
  ! ------------------------------------------------------------------------------!
  ! Public variables:                                                             !
  ! elec_density : stores the properties of an electron density                   !
  !===============================================================================!
  use constants, only : dp
  use basis, only     : castep_basis
  implicit none

  private

  !---------------------------------------------------------------------------!
  !                       P u b l i c   V a r i a b l e s                     !
  !---------------------------------------------------------------------------!

  type,public :: elec_density
     ! All grid values are stored times the number of grid points.
     ! E.g. charge density norm is such that the sum(charge)/total_grid_points
     ! is equal to the number of electrons.
     real(kind=dp),allocatable      :: real_charge(:,:,:)  ! real charge density(nx, ny, nz)
     complex(kind=dp),allocatable   :: charge(:,:,:)       ! complex charge density(nx, ny, nz)
     logical                        :: have_cmplx_den      ! have a complex density
   contains
     procedure :: norm => density_norm
  end type elec_density

  !---------------------------------------------------------------------------!
  !                       P u b l i c   R o u t i n e s                       !
  !---------------------------------------------------------------------------!
  public :: density_allocate
  public :: density_deallocate
  public :: density_complex_to_real
  public :: density_real_to_complex
  public :: density_zero
  public :: density_recip_to_real
  public :: density_write
contains

  subroutine density_allocate(den,cmplx_den)
    !============================================================!
    ! Allocates an electron density based for a given grid       !
    !------------------------------------------------------------!
    ! Arguments                                                  !
    ! den(out)       :: density to allocate                      !
    ! cmplx_den(in)  :: is the density complex (default: true)   !
    !------------------------------------------------------------!
    ! Modules used                                               !
    ! Constants, Basis                                           !
    !------------------------------------------------------------!
    ! Necessary Conditions                                       !
    ! castep_basis should be initialised, i.e. basis_initialise  !
    ! should be called prior to this routine.                    !
    !============================================================!
    use constants, only : cmplx_0
    use basis, only : basis_data

    implicit none

    ! Arguments
    type(elec_density),intent(out)       :: den
    logical, optional,intent(in)         :: cmplx_den

    logical :: l_cmplx_den
    integer :: stat

    ! Decide if density is real or complex(assume complex unless told otherwise)
    l_cmplx_den=.true.
    if(present(cmplx_den)) l_cmplx_den = cmplx_den

    ! Set the appropriate flag in density
    den%have_cmplx_den = l_cmplx_den

    ! Allocate charge density and initialise with zeroes.
    if (den%have_cmplx_den) then
       allocate(den%charge(castep_basis%ngx,castep_basis%ngy,castep_basis%ngz),stat=stat)
       if(stat/=0) error stop 'density_allocate: Failed to allocate complex charge density'
       den%charge=cmplx_0
    else ! real charge density
       allocate(den%real_charge(castep_basis%ngx,castep_basis%ngy,castep_basis%ngz),stat=stat)
       if(stat/=0) error stop 'density_allocate: Failed to allocate real charge density'
       den%charge=0.0_dp
    end if
  end subroutine density_allocate

  subroutine  density_deallocate(den)
    !============================================================!
    ! Deallocates all the arrays in an electron density and      !
    ! 'resets' all other variables.                              !
    !------------------------------------------------------------!
    ! Arguments                                                  !
    ! den(inout)       :: density to allocate                    !
    !============================================================!
    implicit none
    type(elec_density),intent(inout) :: den
    integer :: stat

    ! Deallocate real charge, if allocated
    if(allocated(den%real_charge))then
       deallocate(den%real_charge, stat=stat)
       if(stat/=0) error stop 'density_deallocate: Failed to deallocate real charge density'
    end if

    ! Deallocate complex charge, if allocated
    if(allocated(den%charge))then
       deallocate(den%charge, stat=stat)
       if(stat/=0) error stop 'density_deallocate: Failed to deallocate complex charge density'
    end if

    ! Reset complex flag to true
    den%have_cmplx_den=.true.
  end subroutine density_deallocate

  subroutine density_complex_to_real(den)
    !============================================================!
    ! Converts a density with a complex charge to a density      !
    ! with real charge values.                                   !
    ! NOTE: This is done by discarding the imaginary part        !
    ! of the complex number (and NOT by taking the modulus)      !
    !------------------------------------------------------------!
    ! Arguments                                                  !
    ! den(inout) :: the density to be modified                   !
    !------------------------------------------------------------!
    ! Modules used                                               !
    ! None                                                       !
    !------------------------------------------------------------!
    ! Necessary Conditions                                       !
    ! den contains a valid density                               !
    !============================================================!
    implicit none

    type(elec_density),intent(inout) :: den
    integer :: stat

    ! Check if density is already real
    if (den%have_cmplx_den .eqv. .false.) then
       ! Density is real so nothing to do
       return
    else
       ! Allocate real charge
       allocate(den%real_charge(castep_basis%ngx,castep_basis%ngy,castep_basis%ngz),stat=stat)
       if (stat/=0) error stop 'density_complex_to_real: Failed to allocate real charge.'

       ! Now copy charge over
       den%real_charge=real(den%charge,dp)

       ! Deallocate complex charge and set appropriate flag
       den%have_cmplx_den=.false.
       deallocate(den%charge,stat=stat)
       if (stat/=0) error stop 'density_complex_to_real: Failed to deallocate complex charge.'
    end if
  end subroutine density_complex_to_real

  subroutine density_real_to_complex(den)
    !============================================================!
    ! Converts a density with a real charge to a density         !
    ! with complex charge values.                                !
    ! NOTE: This is done by setting the imaginary part of the    !
    ! complex number to 0                                        !
    !------------------------------------------------------------!
    ! Arguments                                                  !
    ! den(inout) :: the density to be modified                   !
    !------------------------------------------------------------!
    ! Modules used                                               !
    ! None                                                       !
    !------------------------------------------------------------!
    ! Necessary Conditions                                       !
    ! den contains a valid density                               !
    !============================================================!
    implicit none

    type(elec_density),intent(inout) :: den

    integer :: stat

    ! Check if density is already complex
    if (den%have_cmplx_den .eqv. .true.) then
       ! Density is complex so nothing to do
       return
    else
       ! Allocate complex charge
       allocate(den%charge(castep_basis%ngx,castep_basis%ngy,castep_basis%ngz),stat=stat)
       if (stat/=0) error stop 'density_real_to_complex: Failed to allocate complex charge.'

       ! Now copy charge over
       den%charge=cmplx(den%real_charge,0.0_dp,dp)

       ! Deallocate real charge and set appropriate flag
       den%have_cmplx_den=.true.
       deallocate(den%real_charge,stat=stat)
       if (stat/=0) error stop 'density_real_to_complex: Failed to deallocate real charge.'
    end if
  end subroutine density_real_to_complex

  ! TODO Remove redundant routine - initialisation is already done in density_allocate
  subroutine density_zero(den)
    !============================================================!
    ! Initialises the charge density to zero everywhere on a grid!
    !------------------------------------------------------------!
    ! Arguments                                                  !
    ! den(inout) :: the density to be initialised with zeros     !
    !------------------------------------------------------------!
    ! Modules used                                               !
    ! Constants                                                  !
    !------------------------------------------------------------!
    ! Necessary Conditions                                       !
    ! Density should be allocated                                !
    !============================================================!
    use constants,only : cmplx_0
    implicit none
    type(elec_density),intent(inout) :: den
    if (den%have_cmplx_den) then
       den%charge = cmplx_0
    else
       den%real_charge = 0.0_dp
    end if
  end subroutine density_zero

  subroutine density_recip_to_real(recip_den,real_den)
    use basis,only : castep_basis,basis_recip_to_real

    implicit none

    type(elec_density),intent(inout)          :: recip_den
    type(elec_density),intent(inout),optional :: real_den

    complex(kind=dp),allocatable :: real_grid(:,:,:)

    integer :: stat

    allocate(real_grid(castep_basis%ngx,castep_basis%ngy,castep_basis%ngz),stat=stat)
    if(stat/=0) error stop 'density_recip_to_real: Failed to allocate real space grid.'

    ! Initialise real space grid for charge density and initialise with reciprocal space grid
    if (recip_den%have_cmplx_den) then
       real_grid = recip_den%charge
    else
       real_grid = cmplx(recip_den%real_charge,0.0_dp,dp)
    end if

    ! Perform the Fourier transform to get into real space.
    call basis_recip_to_real(real_grid)

    ! Store results of Fourier transform appropriately.
    if (present(real_den)) then
       if (real_den%have_cmplx_den) then
          real_den%charge = real_grid
       else
          real_den%real_charge = real(real_grid,dp)
       endif
    else
       if (recip_den%have_cmplx_den) then
          recip_den%charge = real_grid
       else
          recip_den%real_charge = real(real_grid,dp)
       endif
    end if

    ! Clean up
    deallocate(real_grid,stat=stat)
    if(stat/=0) error stop 'density_recip_to_real: Failed to deallocate real space grid.'
  end subroutine density_recip_to_real

  ! TODO density_real_to_recip

  subroutine density_write(den,fmt)
    !============================================================!
    ! Write the density on a grid in a format that CASTEP        !
    ! understands.                                               !
    ! The format (fmt) used is as follows                        !
    ! 'R' - real density                                         !
    ! 'C' - complex density, this writes only the real part of   !
    !       of the density                                       !
    ! 'A' - absolute value, for a complex density this takes     !
    !       takes the absolute value                             !
    !------------------------------------------------------------!
    ! Arguments                                                  !
    ! den(inout) :: the density to be modified                   !
    !------------------------------------------------------------!
    ! Modules used                                               !
    ! None                                                       !
    !------------------------------------------------------------!
    ! Known issues                                               !
    ! Lattice parameters are hard-coded rather than calculated   !
    ! afresh.                                                    !
    !------------------------------------------------------------!
    ! Necessary Conditions                                       !
    ! den contains a valid density                               !
    ! castep_basis should still be initialised                   !
    !============================================================!
    use basis, only : castep_basis
    use latt,  only : platt,den_fmt_file
    use basis, only : basis_data
    use io,    only : stderr,stdout

    implicit none

    ! Arguments
    type(elec_density),intent(in) :: den         ! density to write to file
    character(len=1), intent(in),optional :: fmt ! 'R' - real, 'C' - complex , 'A' - absolute value/complex modulus

    ! Local variables
    integer :: den_unit ! file unit to write density to
    character(len=1) :: l_fmt
    character(len=80):: iomsg
    integer :: i,ix,iy,iz
    integer :: iostat

    ! Determine format to use from density and then override if needed
    if (den%have_cmplx_den) then
       l_fmt = 'C'
    else
       l_fmt = 'R'
    end if
    if(present(fmt)) l_fmt=fmt

    ! Open the file
    open(file=trim(den_fmt_file),newunit=den_unit,status='REPLACE',action='WRITE',&
         iostat=iostat,iomsg=iomsg)
    if (iostat/=0) then
       write(stderr,'(A7,A)') 'ERROR: ',trim(iomsg)
       error stop 'density_read: Failed to open density file.'
    end if

    write(stdout,'(A32,A)') ' Writing real space density to: ',trim(den_fmt_file)

    ! Write the CASTEP file header
    write(den_unit,'(A12)') 'BEGIN header'
    write(den_unit,*) ''
    write(den_unit,'(A)') 'Real Lattice(Bohr)               Lattice parameters(Bohr)    Cell Angles'

    ! TODO Get lattice parameters from real lattice
    write(den_unit,100) (platt(1,i),i=1,3), 5.13157067_dp/sqrt(2.0_dp), 60.0_dp
    write(den_unit,200) (platt(2,i),i=1,3), 5.13157067_dp/sqrt(2.0_dp), 60.0_dp
    write(den_unit,300) (platt(3,i),i=1,3), 5.13157067_dp/sqrt(2.0_dp), 60.0_dp
    write(den_unit,*) ''
100 format(3f14.7,5x,'a =',f12.6,2x,'alpha =',f12.6)
200 format(3f14.7,5x,'b =',f12.6,2x,'beta  =',f12.6)
300 format(3f14.7,5x,'c =',f12.6,2x,'gamma =',f12.6)

    ! Write spin information and grid information
    write(den_unit,'(A58)') '1   F                        ! nspins, non-collinear spin'
    write(den_unit,'(3(i4,2x),T30,a)') castep_basis%ngx,castep_basis%ngy,castep_basis%ngz, &
         '! fine FFT grid along <a,b,c>'
    select case(l_fmt)
    case('A')
       write(den_unit,'(A)')'END header: data is "<a b c> abs(charge) in units of electrons/grid_point * number of grid_points'
    case('R')
       write(den_unit,'(A)')'END header: data is "<a b c> charge in units of electrons/grid_point * number of grid_points'
    case('C')
       write(den_unit,'(A)')'END header: data is "<a b c> real(charge) imag(charge) &
            & in units of electrons/grid_point * number of grid_points'
    case default
       error stop 'density_write: Unknown format specified.'
    end select
    write(den_unit,*) ''

    select case(l_fmt)
    case('A')
       write(stdout,'(A17,A)') '  Output format: ', 'absolute value'
       do ix=1,castep_basis%ngx
          do iy=1,castep_basis%ngy
             do iz=1,castep_basis%ngz
                if(den%have_cmplx_den)then
                   write(den_unit,'(3I6,f20.12)') ix,iy,iz, &
                        abs(den%charge(ix,iy,iz))
                else
                   write(den_unit,'(3I6,f20.12)') ix,iy,iz, &
                        abs(den%real_charge(ix,iy,iz))
                end if
             end do
          end do
       end do

    case('R')
       write(stdout,'(A17,A)') '  Output format: ', 'real'
       do ix=1,castep_basis%ngx
          do iy=1,castep_basis%ngy
             do iz=1,castep_basis%ngz
                if(den%have_cmplx_den)then
                   write(den_unit,'(3I6,f20.12)') ix,iy,iz, &
                        real(den%charge(ix,iy,iz),dp)
                else
                   write(den_unit,'(3I6,f20.12)') ix,iy,iz, &
                        den%real_charge(ix,iy,iz)
                end if
             end do
          end do
       end do

    case('C')
       write(stdout,'(A17,A)') '  Output format: ', 'complex'
       do ix=1,castep_basis%ngx
          do iy=1,castep_basis%ngy
             do iz=1,castep_basis%ngz
                if(den%have_cmplx_den)then
                   write(den_unit,'(3I6,2f20.12)') ix,iy,iz, &
                        den%charge(ix,iy,iz)
                else
                   write(den_unit,'(3I6,2f20.12)') ix,iy,iz, &
                        cmplx(den%real_charge(ix,iy,iz),0.0_dp,dp)
                end if
             end do
          end do
       end do
    end select

    close(den_unit,iostat=iostat)
    if(iostat/=0) error stop 'density_write: Failed to close formatted density file.'

  end subroutine density_write

  function density_norm(den) result(norm)
    !============================================================!
    ! Calculates the norm (number of electrons) of a charge      !
    ! density                                                    !
    !------------------------------------------------------------!
    ! Necessary Conditions                                       !
    ! density passed must be allocated and in real space         !
    !============================================================!
    use basis,only : castep_basis
    implicit none
    class(elec_density), intent(in) :: den
    real(kind=dp) :: norm

    ! Integrate/sum over all space
    if (den%have_cmplx_den) then
       norm = sum(abs(den%charge))
    else
       norm = sum(den%real_charge)
    end if

    ! Divide by the number of grid points for correct normalisation
    norm = norm/castep_basis%total_grid_points
  end function density_norm
end module density
