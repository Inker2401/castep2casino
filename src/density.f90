! Copyright (C) 2023 Visagan Ravindran
! This file is part of castep2casino
! castep2casino is free software: you can redistribute it and/or modify it under the terms of
! the GNU General Public License as published by the Free Software Foundation,
! either version 3 of the License, or (at your option) any later version.
!
! castep2casino is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
! without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
! See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with castep2casino.
! If not, see <https://www.gnu.org/licenses/>.
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
     real(kind=dp),allocatable      :: real_spin(:,:,:)    ! real spin density(nx,ny,nz)
     complex(kind=dp),allocatable   :: charge(:,:,:)       ! complex charge density(nx, ny, nz)
     complex(kind=dp),allocatable   :: spin(:,:,:)         ! complex spin density(nx,ny,nz)
     logical                        :: have_cmplx_den      ! have a complex density
     integer                        :: nspins              ! number of spin components
                                                           ! 1 - spin degenerate
                                                           ! 2 - spin polarised
   contains
     procedure :: norm => density_norm
     procedure :: mag_moment => density_mag_moment
  end type elec_density

  !---------------------------------------------------------------------------!
  !                       P u b l i c   R o u t i n e s                       !
  !---------------------------------------------------------------------------!
  public :: density_allocate
  public :: density_deallocate
  public :: density_copy
  public :: density_complex_to_real
  public :: density_real_to_complex
  public :: density_recip_to_real
  public :: density_write
  public :: density_shift
  public :: density_to_spin_ch

contains

  subroutine density_allocate(den,nspins,cmplx_den)
    !============================================================!
    ! Allocates an electron density based for a given grid       !
    !------------------------------------------------------------!
    ! Arguments                                                  !
    ! den(out)       :: density to allocate                      !
    ! nspins(in)     :: number of spin components                !
    !                   1 - spin degenerate                      !
    !                   2 - spin polarised                       !
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
    integer, intent(in)                  :: nspins
    logical, optional,intent(in)         :: cmplx_den

    logical :: l_cmplx_den
    integer :: stat

    ! Decide if density is real or complex(assume complex unless told otherwise)
    l_cmplx_den=.true.
    if(present(cmplx_den)) l_cmplx_den = cmplx_den

    ! Get number of spin components
    select case(nspins)
    case (1,2)
       den%nspins = nspins
    case default
       error stop 'density_allocate: nspins should be 1 or 2.'
    end select

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
       den%real_charge=0.0_dp
    end if

    ! Allocate spin and initialise with zeroes.
    if ( den%nspins==2) then
       if (den%have_cmplx_den) then
          allocate(den%spin(castep_basis%ngx,castep_basis%ngy,castep_basis%ngz),stat=stat)
          if(stat/=0) error stop 'density_allocate: Failed to allocate complex spin density'
          den%spin=cmplx_0
       else ! real spin density
          allocate(den%real_spin(castep_basis%ngx,castep_basis%ngy,castep_basis%ngz),stat=stat)
          if(stat/=0) error stop 'density_allocate: Failed to allocate real spin density'
          den%real_spin=0.0_dp
       end if
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

    ! Deallocate real spin, if allocated
    if(allocated(den%real_spin))then
       deallocate(den%real_spin, stat=stat)
       if(stat/=0) error stop 'density_deallocate: Failed to deallocate real spin density'
    end if

    ! Deallocate complex spin, if allocated
    if(allocated(den%spin))then
       deallocate(den%spin, stat=stat)
       if(stat/=0) error stop 'density_deallocate: Failed to deallocate complex spin density'
    end if

    ! Reset complex flag to true
    den%have_cmplx_den=.true.

    ! Reset spins to 1
    den%nspins = 1
  end subroutine density_deallocate

  subroutine density_copy(den1, den2, do_alloc,force_copy)
    !============================================================!
    ! Copies the density 1 into density 2                        !
    ! This routine assumes that den2 is already allocated unless !
    ! do_alloc is set to true.                                   !
    ! Note when copying from complex to real density,            !
    ! an error will be raised unless force_copy is set to true   !
    !------------------------------------------------------------!
    ! Arguments                                                  !
    ! den1(in) :: input density                                  !
    ! den2(inout) :: output density                              !
    ! do_alloc :: allocate den2 with same allocation as den1     !
    !             (default : False)                              !
    ! force_copy:: force a copy of a complex density into real   !
    !              density (default : False)                     !
    !------------------------------------------------------------!
    ! Modules used                                               !
    ! None                                                       !
    !------------------------------------------------------------!
    ! Necessary Conditions                                       !
    ! den1 contains a valid density                              !
    ! den2 is allocated (if do_alloc is false)                   !
    !============================================================!
    implicit none
    type(elec_density), intent(in)    :: den1    ! input density
    type(elec_density), intent(inout) :: den2    ! output density
    logical, optional,  intent(in) :: do_alloc   ! Do allocation of den2 (Default : False)
    logical, optional,  intent(in) :: force_copy ! Do allocation of den2 (Default : False)


    logical :: l_alloc, l_force

    l_alloc = .false. ; l_force = .false.
    if (present(do_alloc)) l_alloc=do_alloc
    if (present(force_copy)) l_force=force_copy

    ! Check if we need to allocate the density
    if (l_alloc) then
       call density_allocate(den2, den1%nspins, den1%have_cmplx_den)
    end if
    if (den1%nspins/=den2%nspins) error stop 'density_copy: Densities do not have the same number of spins'

    ! Now copy density over
    if (den1%have_cmplx_den) then
       if (den2%have_cmplx_den) then
          den2%charge = den1%charge
          if(den1%nspins==2) den2%spin = den1%spin
       else
          if (l_force) then
             den2%real_charge = real(den1%charge, dp)
             if(den1%nspins==2) den2%real_spin = real(den1%spin, dp)
          else
             error stop 'density_copy: den1 is complex but den2 is real'
          end if
       end if

    else ! den1 is real

       if (den2%have_cmplx_den) then
          den2%charge = cmplx(den1%real_charge, 0.0_dp, dp)
          if(den1%nspins==2) den2%spin = cmplx(den1%real_spin, 0.0_dp, dp)
       else
          den2%real_charge = den1%real_charge
          if(den1%nspins==2) den2%real_spin = den1%real_spin
       end if
    end if
  end subroutine density_copy

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

       ! Do the same with spin densities
       if (den%nspins==2) then
          allocate(den%real_spin(castep_basis%ngx,castep_basis%ngy,castep_basis%ngz),stat=stat)
          if (stat/=0) error stop 'density_complex_to_real: Failed to allocate real spin.'

          ! Now copy spin over
          den%real_spin=real(den%spin,dp)

          ! Deallocate complex spin and set appropriate flag
          den%have_cmplx_den=.false.
          deallocate(den%spin,stat=stat)
          if (stat/=0) error stop 'density_complex_to_real: Failed to deallocate complex spin.'
       end if
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

       if (den%nspins==2) then
          ! Allocate complex spin
          allocate(den%spin(castep_basis%ngx,castep_basis%ngy,castep_basis%ngz),stat=stat)
          if (stat/=0) error stop 'density_real_to_complex: Failed to allocate complex spin.'

          ! Now copy spin over
          den%spin=cmplx(den%real_spin,0.0_dp,dp)

          ! Deallocate real spin and set appropriate flag
          den%have_cmplx_den=.true.
          deallocate(den%real_spin,stat=stat)
          if (stat/=0) error stop 'density_real_to_complex: Failed to deallocate real spin.'
       end if
    end if
  end subroutine density_real_to_complex

  subroutine density_recip_to_real(recip_den,real_den)
    !============================================================!
    ! Transforms a density in reciprocal space to a density      !
    ! in real space.                                             !
    ! The real space density overwrites the existing density     !
    ! data in recip_den unless real_den is specified in which    !
    ! case it is placed there.                                   !
    !------------------------------------------------------------!
    ! Arguments                                                  !
    ! recip_den(inout) : the reciprocal space density            !
    ! real_den(inout),optional : the density in whcih we should  !
    !                            store the FFT of recip_den      !
    !------------------------------------------------------------!
    ! Modules used                                               !
    ! Basie                                                      !
    !------------------------------------------------------------!
    ! Necessary conditions                                       !
    ! If real_den is supplied, it should be allocated            !
    !============================================================!
    use basis,only : castep_basis,basis_recip_to_real

    implicit none

    type(elec_density),intent(inout)          :: recip_den
    type(elec_density),intent(inout),optional :: real_den

    complex(kind=dp),allocatable :: real_grid(:,:,:)

    integer :: stat

    allocate(real_grid(castep_basis%ngx,castep_basis%ngy,castep_basis%ngz),stat=stat)
    if(stat/=0) error stop 'density_recip_to_real: Failed to allocate real space grid.'

    ! Sanity check that the real and reciprocal space densities have the same number of spins
    if (present(real_den)) then
       if(real_den%nspins/=recip_den%nspins) &
            error stop 'density_real_to_recip: Real and reciprocal space densities do not have same nspins!'
    end if

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

    ! Now do the same thing with spins if needed
    if (recip_den%nspins==2) then
       ! Initialise real space grid for spin density and initialise with reciprocal space grid
       if (recip_den%have_cmplx_den) then
          real_grid = recip_den%spin
       else
          real_grid = cmplx(recip_den%real_spin,0.0_dp,dp)
       end if

       ! Perform the Fourier transform to get into real space.
       call basis_recip_to_real(real_grid)

       ! Store results of Fourier transform appropriately.
       if (present(real_den)) then
          if (real_den%have_cmplx_den) then
             real_den%spin = real_grid
          else
             real_den%real_spin = real(real_grid,dp)
          endif
       else
          if (recip_den%have_cmplx_den) then
             recip_den%spin = real_grid
          else
             recip_den%real_spin = real(real_grid,dp)
          endif
       end if
    end if

    ! Clean up
    deallocate(real_grid,stat=stat)
    if(stat/=0) error stop 'density_recip_to_real: Failed to deallocate real space grid.'
  end subroutine density_recip_to_real

  subroutine density_write(den,den_file,fmt)
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
    ! den_file(in) :: file to write the density to               !
    ! fmt(in), optional :: format to use when writing            !
    !------------------------------------------------------------!
    ! Modules used                                               !
    ! None                                                       !
    !------------------------------------------------------------!
    ! Necessary Conditions                                       !
    ! den contains a valid density                               !
    ! castep_basis should still be initialised                   !
    !============================================================!
    use basis, only : castep_basis
    use latt,  only : user_params
    use basis, only : basis_data
    use io,    only : stderr,stdout

    implicit none

    ! Arguments
    type(elec_density),intent(in) :: den         ! density to write to file
    character(len=*), intent(in)          :: den_file ! file to write density to
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
    open(file=trim(den_file),newunit=den_unit,status='REPLACE',action='WRITE',&
         iostat=iostat,iomsg=iomsg)
    if (iostat/=0) then
       write(stderr,'(A7,A)') 'ERROR: ',trim(iomsg)
       error stop 'density_read: Failed to open density file.'
    end if

    ! Write the CASTEP file header
    write(den_unit,'(A12)') 'BEGIN header'
    write(den_unit,*) ''
    write(den_unit,'(15x,A)') 'Real Lattice(Bohr)               Lattice parameters(Bohr)    Cell Angles'

    write(den_unit,100) (user_params%platt(1,i),i=1,3), user_params%lat_consts(1), user_params%lat_angles(1)
    write(den_unit,200) (user_params%platt(2,i),i=1,3), user_params%lat_consts(2), user_params%lat_angles(2)
    write(den_unit,300) (user_params%platt(3,i),i=1,3), user_params%lat_consts(3), user_params%lat_angles(3)
    write(den_unit,*) ''
100 format(3f14.7,5x,'a =',f12.6,2x,'alpha =',f12.6)
200 format(3f14.7,5x,'b =',f12.6,2x,'beta  =',f12.6)
300 format(3f14.7,5x,'c =',f12.6,2x,'gamma =',f12.6)

    ! Write spin information and grid information
    write(den_unit,'(1x,I1,3x,A53)') den%nspins, 'F                        ! nspins, non-collinear spin'
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
                   if(den%nspins==2) then
                      write(den_unit,'(3I6,2f20.12)') ix,iy,iz, &
                           abs(den%charge(ix,iy,iz)), abs(den%spin(ix,iy,iz))
                   else
                      write(den_unit,'(3I6,f20.12)') ix,iy,iz, &
                           abs(den%charge(ix,iy,iz))
                   end if
                else
                   if(den%nspins==2) then
                      write(den_unit,'(3I6,2f20.12)') ix,iy,iz, &
                           abs(den%real_charge(ix,iy,iz)), abs(den%real_spin(ix,iy,iz))
                   else
                      write(den_unit,'(3I6,f20.12)') ix,iy,iz, &
                           abs(den%real_charge(ix,iy,iz))
                   end if
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
                   if(den%nspins==2) then
                      write(den_unit,'(3I6,2f20.12)') ix,iy,iz, &
                           real(den%charge(ix,iy,iz),dp), real(den%spin(ix,iy,iz),dp)
                   else
                      write(den_unit,'(3I6,f20.12)') ix,iy,iz, &
                           real(den%charge(ix,iy,iz),dp)
                   end if
                else
                   if(den%nspins==2) then
                      write(den_unit,'(3I6,2f20.12)') ix,iy,iz, &
                           den%real_charge(ix,iy,iz), den%real_spin(ix,iy,iz)
                   else
                      write(den_unit,'(3I6,f20.12)') ix,iy,iz, &
                           den%real_charge(ix,iy,iz)
                   end if
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
                   if(den%nspins==2) then
                      write(den_unit,'(3I6,2f20.12)') ix,iy,iz, &
                           den%charge(ix,iy,iz), den%spin(ix,iy,iz)
                   else
                      write(den_unit,'(3I6,f20.12)') ix,iy,iz, &
                           den%charge(ix,iy,iz)
                   end if
                else
                   if(den%nspins==2) then
                      write(den_unit,'(3I6,2f20.12)') ix,iy,iz, &
                           cmplx(den%real_charge(ix,iy,iz),0.0_dp,dp), cmplx(den%real_spin(ix,iy,iz),0.0_dp,dp)
                   else
                      write(den_unit,'(3I6,f20.12)') ix,iy,iz, &
                           cmplx(den%real_charge(ix,iy,iz),0.0_dp,dp)
                   end if
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
       norm = abs(sum(den%charge))
    else
       norm = sum(den%real_charge)
    end if

    ! Divide by the number of grid points for correct normalisation
    norm = norm/castep_basis%total_grid_points
  end function density_norm

  function density_mag_moment(den,absval) result(mag_moment)
    !============================================================!
    ! Calculates the (abs) magnetic moment for a spin density    !
    ! in atomic units                                            !
    !------------------------------------------------------------!
    ! Necessary Conditions                                       !
    ! density passed must be allocated and in real space         !
    ! den%nspins=2, i.e. must have a spin density                !
    !============================================================!
    use basis,only : castep_basis
    implicit none
    class(elec_density), intent(in) :: den
    logical,intent(in),optional     :: absval

    real(kind=dp) :: mag_moment
    logical :: l_abs
    l_abs = .false.
    if(present(absval)) l_abs = absval

    ! Integrate/sum over all space
    if (den%have_cmplx_den) then
       if (l_abs) then
          mag_moment = sum(abs(real(den%spin,dp)))
       else
          mag_moment = sum(real(den%spin,dp))
       end if
    else
       if (l_abs) then
          mag_moment = abs(sum(den%real_spin))
       else
          mag_moment = sum(den%real_spin)
       end if
    end if

    ! Divide by the number of grid points for correct normalisation
    mag_moment = mag_moment/castep_basis%total_grid_points
  end function density_mag_moment

  subroutine density_shift(den)
    !============================================================!
    ! Shift the density on a grid by an amount specified in      !
    ! fractional coordinates.                                    !
    !------------------------------------------------------------!
    ! Arguments                                                  !
    ! den(inout) :: the density to be shifted                    !
    !------------------------------------------------------------!
    ! Necessary Conditions                                       !
    ! density passed must be allocated and in real space         !
    !============================================================!
    use io,only    : stdout
    use basis,only : basis_shift
    use latt,only  : user_params

    implicit none

    type(elec_density),intent(inout) :: den

    ! Perform a shift of the charge density grid but only if requested.
    if (user_params%shift_grid) then
       write(stdout,'(/,A64)') ' WARNING: Shifted real space grid - DO NOT USE FOR CALCULATION!'
       write(stdout,'(A20,3F10.5,/)') '  Real space shift: ',user_params%shift_frac

       if (den%have_cmplx_den) then
          call basis_shift(den%charge)
          if(den%nspins==2) call basis_shift(den%spin)
       else ! real density
          call basis_shift(den%real_charge)
          if(den%nspins==2) call basis_shift(den%real_spin)
       endif

    end if

  end subroutine density_shift

  subroutine density_to_spin_ch(den1,den2)
    !============================================================!
    ! Turn charge and spin density into the density for each     !
    ! spin channel, i.e. rho^up and rho^down                     !
    ! These are stored in the charge and spin arrays of the      !
    ! electron_density type respectively.                        !
    !------------------------------------------------------------!
    ! Arguments                                                  !
    ! den1(in)  :: charge and spin density                       !
    ! den2(out) :: density for each spin channel                 !
    !============================================================!
    implicit none
    type(elec_density), intent(in)  :: den1
    type(elec_density), intent(out) :: den2

    ! If no spin, then we have nothing to do
    if (den1%nspins/=2) return

    if (den1%have_cmplx_den) then
       call density_allocate(den2,2,.true.)
    else
       call density_allocate(den2,2,.false.)
    end if

    if (den1%have_cmplx_den) then
       den2%charge = (den1%charge + den1%spin)/2.0_dp
       den2%spin =(den1%charge - den1%spin)/2.0_dp
    else
       den2%real_charge = (den1%real_charge + den1%real_spin)/2.0_dp
       den2%real_spin =(den1%real_charge - den1%real_spin)/2.0_dp
    end if

  end subroutine density_to_spin_ch
end module density
