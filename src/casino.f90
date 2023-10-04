module casino
  !===============================================================================!
  !                      ____          _                                          !
  !                     / ___|__ _ ___(_)_ __   ___                               !
  !                    | |   / _` / __| | '_ \ / _ \                              !
  !                    | |__| (_| \__ \ | | | | (_) |                             !
  !                     \____\__,_|___/_|_| |_|\___/                              !
  ! This module handles all processing tasks related to the CASINO density,       !
  ! in reciprocal space. The main bit involves the reading of G-vectors           !
  ! but also checking for inversion symmetry etc.                                 !
  !-------------------------------------------------------------------------------!
  ! Modules used                                                                  !
  ! See below                                                                     !
  !-------------------------------------------------------------------------------!
  ! Public variables:                                                             !
  ! ngvec - the number of G-vectors in the file                                   !
  ! gvecs - the G-vectors of the non-zero Fourier components                      !
  ! rho_gs - the Fourier components of the density in reciprocal space            !
  !===============================================================================!
  use constants,only : dp
  use io,       only : stdout

  private
  !---------------------------------------------------------------------------!
  !                       P u b l i c   V a r i a b l e s                     !
  !---------------------------------------------------------------------------!
  integer,public,save                      :: ngvec       ! number of G-vectors
  real(kind=dp),public,save,allocatable    :: gvecs(:,:)  ! Gvector, component
  complex(kind=dp),public,save,allocatable :: rho_gs(:)   ! Fourier components of the density

  !---------------------------------------------------------------------------!
  !                     P r i v a t e   V a r i a b l e s                     !
  !---------------------------------------------------------------------------!

  !---------------------------------------------------------------------------!
  !                       P u b l i c   R o u t i n e s                       !
  !---------------------------------------------------------------------------!
  public :: casino_read
  public :: casino_to_castep
  public :: casino_check_symmetry

contains

  subroutine casino_read(filename)
    !============================================================!
    ! Reads the CASINO file and obtains the number of G-vectors  !
    ! the G-vectors themselves and the Fourier components of the !
    ! density.                                                   !
    ! This routine should be called before calling any other     !
    ! routines in this module.                                   !
    !------------------------------------------------------------!
    ! Arguments                                                  !
    ! filename(in) : the name of the CASINO file                 !
    !------------------------------------------------------------!
    ! Modules used                                               !
    ! IO                                                         !
    !------------------------------------------------------------!
    ! Necessary Conditions                                       !
    ! filename must point to a FORMATTED file                    !
    !============================================================!
    use io,only : io_skip_header

    implicit none

    character(len=*),intent(in) :: filename ! CASINO file name

    ! File headers to skip past in CASINO file
    character(len=20),parameter :: ngvec_header='Number of G vectors:'
    character(len=25),parameter :: gvec_header='G vectors (Hartree a.u.):'
    character(len=38),parameter :: rho_g_header='Fourier coefficients of density set 1:'

    integer :: unit,iostat,stat
    integer :: i
    character(len=60) :: iomsg

    ! Open the file
    open(file=trim(filename),newunit=unit,status='OLD',action='READ',&
         iostat=iostat,iomsg=iomsg)
    if(iostat/=0) then
       write(*,'(A7,A)') 'ERROR: ',trim(iomsg)
       stop 'Could not open CASINO file'
    end if

    ! Obtain the number of G-vectors
    call io_skip_header(unit,ngvec_header)
    read(unit,*,iostat=iostat) ngvec
    if(iostat/=0) error stop 'casino_read: Failed to read number of G-vectors'
    write(stdout,'(A22,I7)') ' Number of G-vectors: ', ngvec

    ! Obtain the G-vectors
    allocate(gvecs(ngvec,3),stat=stat)
    if(stat/=0) error stop 'casino_read: Failed to allocate G-vectors array.'
    call io_skip_header(unit,gvec_header)
    do i=1,ngvec
       read(unit,*,iostat=iostat) gvecs(i,:)
       if(iostat/=0) error stop 'casino_read: Failed to read G-vectors'
    end do
    ! DEBUG G-vector read
    ! write(*,*) 'G-vectors'
    ! write(*,*) gvecs(1,:)
    ! write(*,*) gvecs(2,:)
    ! write(*,*) gvecs(ngvec-1,:)
    ! write(*,*) gvecs(ngvec,:)

    ! Obtain the Fourier components
    allocate(rho_gs(ngvec),stat=stat)
    if(stat/=0) error stop 'casino_read: Failed to allocate Fourier components array.'
    call io_skip_header(unit,rho_g_header)
    do i=1,ngvec
       read(unit,*,iostat=iostat) rho_gs(i)
       if(iostat/=0) error stop 'casino_read: Failed to read Fourier components'
    end do

    ! DEBUG Fourier componets
    ! write(*,*) 'Fourier components'
    ! write(*,*) rho_gs(1)
    ! write(*,*) rho_gs(2)
    ! write(*,*) rho_gs(ngvec-2)
    ! write(*,*) rho_gs(ngvec-1)
    ! write(*,*) rho_gs(ngvec)

    close(unit,iostat=iostat,status='KEEP')
    if(iostat/=0) error stop 'casino_read: Failed to close CASINO file.'

  end subroutine casino_read

  subroutine casino_to_castep(den)
    !============================================================!
    ! This routine takes the reciprocal density in CASINO format !
    ! and turns it into a CASTEP formatted density               !
    ! NB - still in reciprocal space.                            !
    !------------------------------------------------------------!
    ! Arguments                                                  !
    ! den(out) : the reciprocal space density on a grid          !
    !------------------------------------------------------------!
    ! Modules used                                               !
    ! IO,Basis,Density                                           !
    !------------------------------------------------------------!
    ! Necessary Conditions                                       !
    ! The primitive lattice vectors, platt, from lattice module  !
    ! must have been initialised via a call to latt_read         !
    ! casino_read must be called before calling this routine.    !
    !============================================================!
    use io,only : stderr
    use basis, only : castep_basis
    use density,only: electron_density,density_allocate

    implicit none
    type(electron_density),intent(out) :: den  ! the density in reciprocal space

    integer,allocatable :: gvec_int(:,:) ! G-vecs as multiples of recip. latt. vectors
    integer :: min_x,min_y,min_z         ! minimum grid size required
    logical :: grid_ok                   ! is the CASTEP grid large enough

    ! Turn the G-vectors into integer multiples of primitive lattice vectors
    call casino_G_to_int(gvec_int)

    ! Check if CASTEP grid can contain all the G-vectors
    min_x = 2*maxval(gvec_int(:,1))+1
    min_y = 2*maxval(gvec_int(:,2))+1
    min_z = 2*maxval(gvec_int(:,3))+1
    write(stdout,'(A31,3I4)') ' Minimum CASTEP grid required: ', min_x, min_y, min_z

    ! Assume the CASTEP grid is large enough and then check.
    ! Continue checking all dimensions before then stopping execution.
    grid_ok=.true.
    if (castep_basis%ngx<=min_z) then
       grid_ok=.false.
       write(stderr,*) 'ERROR: x-dimension of CASTEP grid too small'
    end if
    if (castep_basis%ngy<=min_z) then
       grid_ok=.false.
       write(stderr,*) 'ERROR: y-dimension of CASTEP grid too small'
    end if
    if (castep_basis%ngz<=min_z) then
       grid_ok=.false.
       write(stderr,*) 'ERROR: z-dimension of CASTEP grid too small'
    end if
    if (.not.grid_ok) then
       stop
    end if

    ! Check for symmetries in the density just to make sure it's been read in correctly.
    call casino_check_symmetry(gvec_int)

    ! Now we read the density for real
    ! TODO - Seeing as its real for inversion symmetry, we can reduce memory costs and gain speed for FFT using the corresponding real routines.
    call density_allocate(den)
    call casino_read_recip_grid(den,gvec_int)

  end subroutine casino_to_castep

  subroutine casino_G_to_int(gvec_int)
    !============================================================!
    ! This routine turns the G-vectors into integer multiples of !
    ! the primitive reciprocal lattice vectors. This is done     !
    ! by usage of the orthogonality condition between the        !
    ! primitive real and reciprocal lattice vectors              !
    !          G_i R_j = 2pi * delta_ij                          !
    !------------------------------------------------------------!
    ! Arguments                                                  !
    ! gvec_int(out) : G-vectors as integer multiples of recip    !
    !                 lattice vectors.                           !
    !------------------------------------------------------------!
    ! Modules used                                               !
    ! Constants,latt                                             !
    !------------------------------------------------------------!
    ! Necessary Conditions                                       !
    ! The primitive lattice vectors, platt, from lattice module  !
    ! must have been initialised via a call to latt_read         !
    ! casino_read must be called before calling this routine.    !
    !============================================================!
    use constants, only : pi
    use latt, only      : platt

    implicit none
    integer,allocatable,intent(out) :: gvec_int(:,:)


    integer :: i
    integer :: stat

    allocate(gvec_int(ngvec,3),stat=stat)
    if(stat/=0) error stop 'casino_G_to_int: Failed to allocate integer G-grid'

    ! Multiply each G-vector component by corresponding reciprocal lattice vector
    ! making sure to multiply by pre-factor of 1/2*pi
    do i=1,ngvec
       gvec_int(i,:) = nint(matmul(platt,gvecs(i,:))/(2.0_dp*pi))
    end do

    ! DEBUG int grid
    ! do i=1,ngvec
    !    write(99,'(3f10.4)') real(gvec_int(i,:),dp)
    ! end do

    ! Check symmetric grid
    if (maxval(gvec_int(:,1))/=-minval(gvec_int(:,1))) then
       write(stdout,'(A54)') ' WARNING: G-vectors do not appear symmetrised along x'
    elseif (maxval(gvec_int(:,2))/=-minval(gvec_int(:,2))) then
       write(stdout,'(A54)') ' WARNING: G-vectors do not appear symmetrised along y'
    elseif (maxval(gvec_int(:,3))/=-minval(gvec_int(:,3))) then
       write(stdout,'(A54)') ' WARNING: G-vectors do not appear symmetrised along z'
    end if
  end subroutine casino_G_to_int

  subroutine casino_check_symmetry(gvec_int,have_inv_sym)
    !============================================================!
    ! This routine checks the CASINO density to ensure that we   !
    ! have the correct expected symmetries.                      !
    !------------------------------------------------------------!
    ! Arguments                                                  !
    ! gvec_int(in) : G-vectors as integer multiples of recip     !
    !                 lattice vectors.                           !
    ! have_inv_sym(out) : Does the density have inversion        !
    !                     symmetry                               !
    !------------------------------------------------------------!
    ! Modules used                                               !
    ! Constants,Basis                                            !
    !------------------------------------------------------------!
    ! Necessary Conditions                                       !
    ! basis_initialised should have been called.                 !
    ! casino_read must be called before calling this routine.    !
    !============================================================!
    use constants,only : cmplx_0
    use basis,only : castep_basis

    implicit none
    integer,intent(in) :: gvec_int(:,:)
    logical,intent(out),optional :: have_inv_sym    ! Does density have inversion symmetry

    complex(kind=dp),allocatable :: den_grid(:,:,:) ! Temporary complex grid for symmetries

    integer :: nx,ny,nz    ! max symmetric grid length along x,y,z
    logical :: l_inv_sym   ! local copy of inv_sym
    integer :: i           ! loop counters
    integer :: stat

    ! Find out max grid length needed for number of points
    nx = (castep_basis%ngx-1)/2
    ny = (castep_basis%ngy-1)/2
    nz = (castep_basis%ngy-1)/2

    ! This is not the order we will use to actually store the density but it is easier to check for symmetries this way
    allocate(den_grid(-nx:nx,-ny:ny,-nz:nz),stat=stat)
    if(stat/=0) error stop 'casino_check_symmetry: Failed to allocate grid'

    ! Pad with zeros and then overwrite with non-zero Fourier components
    den_grid = cmplx_0
    do i=1,ngvec
       den_grid(gvec_int(i,1),gvec_int(i,2),gvec_int(i,3)) = rho_gs(i)
    end do

    ! First check for complex conjugate symmetry
    if(.not.casino_check_cmplx_conj(den_grid,nx,ny,nz)) then
       error stop 'Reciprocal space density does not appear to have complex conjugate symmetry'
    ! else
    !    write(stdout,*) ' Reciprocal density has complex conjugate symmetry. <-- SYMMETRY_CHECK'
    end if

    l_inv_sym = casino_check_inver_sym(den_grid,nx,ny,nz)
    if(.not.l_inv_sym) then
       write(stdout,'(A61)') ' Reciprocal space density does not have inversion symmetry   '
    else
       write(stdout,'(A61)') ' SYMMETRY DETECTED: Reciprocal density has inversion symmetry'
    end if
    write(stdout,*) ''

    if(present(have_inv_sym))have_inv_sym=l_inv_sym

    ! Clean up
    deallocate(den_grid)

  contains

    logical function casino_check_cmplx_conj(den_grid,nx,ny,nz)
      !============================================================!
      ! This function checks for complex conjugate symmetry of the !
      ! the density. Since the real space density is a real        !
      ! function, in reciprocal space, the Fourier components      !
      ! satisfy                                                    !
      !       rho_G = rho_-G^*                                     !
      !------------------------------------------------------------!
      ! Arguments                                                  !
      ! den_grid(in) : The density on a grid with symmetric bounds !
      ! nx,ny,nz(in) : upper bound of den_grid                     !
      !------------------------------------------------------------!
      ! Modules used                                               !
      ! Constants,Basis                                            !
      !------------------------------------------------------------!
      ! Necessary Conditions                                       !
      ! den_grid must have array bounds that are symmetric, i.e.   !
      ! if there are a TOTAL of ngx points along x, the bounds for !
      ! x must run from -x0 to x0 where x0=(ngx-1)/2               !
      !============================================================!
      use math, only : math_isclose

      implicit none
      complex(kind=dp),allocatable :: den_grid(:,:,:) ! density stored on grid with symmetric bounds
      integer :: nx,ny,nz

      integer :: ix,iy,iz

      ! Assume true until told otherwise
      casino_check_cmplx_conj=.true.

      ! Loop over all grid points.
      ! Since we are just checking the symmetry, we just need to loop over half of them
      ! as they should be the same!
      do ix=0,nx
         do iy=0,ny
            do iz=0,nz
               if(.not. math_isclose(den_grid(ix,iy,iz),den_grid(-ix,-iy,-iz)) ) then
                  casino_check_cmplx_conj=.false.
                  return
               end if
            end do
         end do
      end do
    end function casino_check_cmplx_conj

    logical function casino_check_inver_sym(den_grid,nx,ny,nz)
      !============================================================!
      ! This function checks for inversion symmetry of the         !
      ! the density where the Fourier components satisfy           !
      !       rho_G = rho_G^*                                      !
      ! i.e. the components are real.                              !
      !------------------------------------------------------------!
      ! Arguments                                                  !
      ! den_grid(in) : The density on a grid with symmetric bounds !
      ! nx,ny,nz(in) : upper bound of den_grid                     !
      !------------------------------------------------------------!
      ! Modules used                                               !
      ! Constants,Basis                                            !
      !------------------------------------------------------------!
      ! Necessary Conditions                                       !
      ! den_grid must have array bounds that are symmetric, i.e.   !
      ! if there are a TOTAL of ngx points along x, the bounds for !
      ! x must run from -x0 to x0 where x0=(ngx-1)/2               !
      !============================================================!
      use math, only : math_isclose
      implicit none
      complex(kind=dp),allocatable :: den_grid(:,:,:) ! density stored on grid with symmetric bounds
      integer :: nx,ny,nz

      integer :: ix,iy,iz

      ! Assume true until told otherwise
      casino_check_inver_sym=.true.

      ! Loop over all grid points.
      ! Since we are just checking the symmetry, we just need to loop over half of them
      ! as they should be the same!
      do ix=-nx,nx
         do iy=-ny,ny
            do iz=-nz,nz
               if(.not. math_isclose(den_grid(ix,iy,iz),conjg(den_grid(ix,iy,iz))) ) then
                  casino_check_inver_sym=.false.
                  return
               end if
            end do
         end do
      end do
    end function casino_check_inver_sym

  end subroutine casino_check_symmetry

  subroutine casino_read_recip_grid(den,gvec_int)
    !============================================================!
    ! This routine takes the reciprocal density in CASINO format !
    ! and turns it into a CASTEP formatted density               !
    ! NB - still in reciprocal space.                            !
    !------------------------------------------------------------!
    ! Arguments                                                  !
    ! den(out) : the reciprocal space density on a grid          !
    ! gvec_int(in) : G-vectors as integer multiples of           !
    !            reciprocal lattice vectors.                     !
    !------------------------------------------------------------!
    ! Modules used                                               !
    ! Basis,Density                                              !
    !------------------------------------------------------------!
    ! Necessary Conditions                                       !
    ! gvec_int should contain grid locations assuming a symmetric!
    ! grid, i.e. from -nx to nx, rather than 1 to 2*nx+1         !
    ! We deal with this slight difference of convention in this  !
    ! routine!!!!                                                !
    !============================================================!
    use density,only : electron_density,density_zero
    use basis,only   : castep_basis

    implicit none

    type(electron_density),intent(inout) :: den
    integer,intent(in)    :: gvec_int(:,:)

    integer :: ix,iy,iz
    integer :: i

    ! Zero the density
    call density_zero(den)

    ! Now read in non-zero Fourier components
    do i=1,ngvec
       ! Currently we have a symmetric grid dimensions from -nx to nx for a total of ngx=2*nx+1 points.
       ! We want this to start from 1 to ngx.
       ix = gvec_int(i,1) + 1
       iy = gvec_int(i,2) + 1
       iz = gvec_int(i,3) + 1
       if(ix<1) ix = castep_basis%ngx+ix
       if(iy<1) iy = castep_basis%ngy+iy
       if(iz<1) iz = castep_basis%ngz+iz
       if (den%have_cmplx_den) then
          den%charge(ix,iy,iz) = rho_gs(i)
       else
          den%real_charge(ix,iy,iz) = real(rho_gs(i),dp)
       end if
    end do
  end subroutine casino_read_recip_grid
end module casino
