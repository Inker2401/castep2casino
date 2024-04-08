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
  integer,public,save                      :: ngvec         ! number of G-vectors
  integer,public,save                      :: nspins        ! number of spins in density

  !---------------------------------------------------------------------------!
  !                     P r i v a t e   V a r i a b l e s                     !
  !---------------------------------------------------------------------------!
  real(kind=dp),save,allocatable           :: gvecs(:,:)    ! Gvector, component
  complex(kind=dp),save,allocatable        :: charge_gs(:)  ! Fourier components of total charge density
                                                            ! i.e. rho^up + rho^down
  complex(kind=dp),save,allocatable        :: spin_gs(:)    ! Fourier components of total charge density
                                                            ! i.e. rho^up - rho^down

  !---------------------------------------------------------------------------!
  !                       P u b l i c   R o u t i n e s                       !
  !---------------------------------------------------------------------------!
  public :: casino_read
  public :: casino_to_castep
  public :: casino_check_symmetry
  public :: casino_truncate_gs

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
    ! Key internal variables                                     !
    ! rho_gs : - Fourier component for the density of each spin  !
    !            channel, i.e rho^up and rho^down                !
    !            These are normalised to the number of electrons !
    !            for each spin                                   !
    !------------------------------------------------------------!
    ! Parent variables used                                      !
    ! ngvec: number of G-vectors                                 !
    ! charge_gs : Fourier components of total charge density     !
    !             rho = rho^up + rho^down                        !
    ! spin_gs   : Fourier components of total spin density       !
    !             s = rho^up - rho^down                          !
    !------------------------------------------------------------!
    ! Necessary Conditions                                       !
    ! filename must point to a FORMATTED file                    !
    !============================================================!
    use io,only : io_skip_header
    use latt, only : user_params

    implicit none

    character(len=*),intent(in) :: filename ! CASINO file name

    ! File headers to skip past in CASINO file
    character(len=20),parameter :: ngvec_header='Number of G vectors:'
    character(len=25),parameter :: gvec_header='G vectors (Hartree a.u.):'
    character(len=27),parameter :: nspins_header='Number of charge densities:'
    character(len=38),parameter :: rho_g_header_up='Fourier coefficients of density set 1:'
    character(len=38),parameter :: rho_g_header_down='Fourier coefficients of density set 2:'

    ! For writing G-vectors to a file
    character(len=15),parameter :: untrun_gfile='untrun_rhoG.txt'
    integer :: unit,iostat,stat
    integer :: i
    character(len=60) :: iomsg

    complex(kind=dp),allocatable :: rho_gs(:,:)   ! Fourier components of the density for each spin channel, spin,
                                                  ! i.e. rho^up rho^down

    ! Open the file
    open(file=trim(filename),newunit=unit,status='OLD',action='READ',&
         iostat=iostat,iomsg=iomsg)
    if(iostat/=0) then
       write(*,'(A7,A)') 'ERROR: ',trim(iomsg)
       stop 'Could not open CASINO file'
    end if

    ! Obtain the number of G-vectors
    call io_skip_header(unit,ngvec_header,case_sensitive=.true.) ! Make case sensitive just in case file format changes...
    read(unit,*,iostat=iostat) ngvec
    if(iostat/=0) error stop 'casino_read: Failed to read number of G-vectors'
    write(stdout,'(A22,I7)') ' Number of G-vectors: ', ngvec

    ! Obtain the G-vectors
    allocate(gvecs(ngvec,3),stat=stat)
    if(stat/=0) error stop 'casino_read: Failed to allocate G-vectors array.'
    call io_skip_header(unit,gvec_header,case_sensitive=.true.)
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

    ! Obtain the number of spins 05/12/2023
    call io_skip_header(unit,nspins_header,case_sensitive=.true.)
    read(unit,*,iostat=iostat) nspins
    if(iostat/=0) error stop 'casino_read: Failed to read number of spins.'
    select case(nspins)
    case(1,2)
       write(stdout,'(A28,I1)') ' Number of spin components: ', nspins
    case default
       stop 'casino_read: Error in CASINO file. Should have only 1 or 2 spin components!'
    end select

    ! Obtain the Fourier components...
    allocate(rho_gs(ngvec,nspins),stat=stat)
    if(stat/=0) error stop 'casino_read: Failed to allocate Fourier components array.'
    call io_skip_header(unit,rho_g_header_up,case_sensitive=.true.)
    do i=1,ngvec
       read(unit,*,iostat=iostat) rho_gs(i,1)
       if(iostat/=0) error stop 'casino_read: Failed to read Fourier components'
    end do

    ! ... remembering to do the second spin if we have it
    if (nspins==2) then
       call io_skip_header(unit,rho_g_header_down,case_sensitive=.true.)
       do i=1,ngvec
          read(unit,*,iostat=iostat) rho_gs(i,2)
          if(iostat/=0) error stop 'casino_read: Failed to read Fourier components for second spin component'
       end do
    end if

    ! DEBUG Fourier componets
    ! write(*,*) 'Fourier components'
    ! write(*,*) rho_gs(1,1)
    ! write(*,*) rho_gs(2,1)
    ! write(*,*) rho_gs(ngvec-2,1)
    ! write(*,*) rho_gs(ngvec-1,1)
    ! write(*,*) rho_gs(ngvec,1)

    ! Now turn the densities into total charge densities and spin densities
    ! NOTE that the density for each spin channel is normalised to the number of electrons for that channel
    allocate(charge_gs(ngvec),stat=stat)
    if(stat/=0) error stop 'casino_read: Failed to allocate total charge density'
    if (nspins==2) then
       charge_gs = rho_gs(:,1) + rho_gs(:,2)
       allocate(spin_gs(ngvec),stat=stat)
       if(stat/=0) error stop 'casino_read: Failed to allocate spin density'
       spin_gs = rho_gs(:,1) - rho_gs(:,2)
    else
       charge_gs = rho_gs(:,1)
    end if

    ! We have the total charge density and spin density so we don't need the ones for each spin channel anymore
    deallocate(rho_gs,stat=stat)
    if(stat/=0) error stop 'casino_read: Failed to deallocate rho_gs'

    if(.not. user_params%no_write_untrun_g) then
       ! Write the Fourier components out as an absolute G-vector out to a file
       write(stdout,'(A56,A,/)') ' Writing out planewave energy of Fourier components to: ',trim(untrun_gfile)
       call casino_plot_Gs(charge_gs,untrun_gfile)
    end if

    close(unit,iostat=iostat,status='KEEP')
    if(iostat/=0) error stop 'casino_read: Failed to close CASINO file.'

  end subroutine casino_read

  subroutine casino_plot_Gs(den_gs,g_filename)
    !============================================================!
    ! Write the Fourier components of the density as a function  !
    ! of the planewave energy |G|^2/2  to a file                 !
    !------------------------------------------------------------!
    ! Arguments                                                  !
    ! den_gs(in)     : the Fourier components of charge/spin     !
    !                  density                                   !
    ! g_filename(in) : filename containing the G-vectors         !
    !------------------------------------------------------------!
    ! Modules used                                               !
    ! IO                                                         !
    !------------------------------------------------------------!
    ! Parent module variables used                               !
    ! ngvec     : number of G-vectors                            !
    !                                                            !
    ! gvecs  : the G-vectors at which the charge density is      !
    !          non-zero in reciprocal space                      !
    !------------------------------------------------------------!
    ! Necessary conditions                                       !
    ! casino_read should be called before calling this routine   !
    !============================================================!
    implicit none

    ! Arguments
    complex(kind=dp),intent(in) :: den_gs(:)
    character(len=*),intent(in) :: g_filename
    integer :: g_unit
    character(len=90) :: iomsg

    integer :: i
    integer :: iostat

    open(file=trim(g_filename),newunit=g_unit,status='UNKNOWN',action='WRITE',iostat=iostat,iomsg=iomsg)
    if(iostat/=0) then
       write(stdout,*) trim(iomsg)
       stop 'casino_plot_Gs: Failed to open file'
    end if

    ! Write a header so we know what we have
    write(g_unit,'(A30)') '#      |G|^2/2(Ha)     |rho_G|'
    ! Now write the planewave cutoffs |G|^2/2 and the Fourier components of the density.
    do i=1,ngvec
       ! NOTE: Here we take the absolute value as the density coefficients can in general be complex if we don't have inversion symmetry.
       write(g_unit,'(F18.10,ES20.11)') sum(gvecs(i,:)**2.0_dp)/2.0_dp, abs(den_gs(i))
    end do

    close(g_unit,iostat=iostat)
    if(iostat/=0) error stop 'casino_plot_Gs: Failed to close file.'
  end subroutine casino_plot_Gs

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
    ! Parent module variables used                               !
    ! charge_gs, spin_gs : Fourier components of charge/spin     !
    !                      densities                             !
    !------------------------------------------------------------!
    ! Necessary Conditions                                       !
    ! The user's input parameters must be initialised via call   !
    ! to latt_read must be called before calling this routine.   !
    !============================================================!
    use io,     only : stderr
    use basis,  only : castep_basis
    use density,only : elec_density

    implicit none
    type(elec_density),intent(out) :: den  ! the density in reciprocal space

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
    call casino_check_symmetry(charge_gs,gvec_int,label='C')
    if(nspins==2) call casino_check_symmetry(spin_gs,gvec_int,label='S')
    write(stdout,*) ''

    ! Now we store the non-zero reciprocal space density components onto a grid
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
    ! Parent module variables used                               !
    ! gvecs : the G-vectors in reciprocal space in units of      !
    !         inverse Bohr                                       !
    !------------------------------------------------------------!
    ! Necessary Conditions                                       !
    ! The primitive lattice vectors, platt, from lattice module  !
    ! must have been initialised via a call to latt_read         !
    ! casino_read must be called before calling this routine.    !
    !============================================================!
    use constants, only : pi
    use latt,      only : user_params

    implicit none
    integer,allocatable,intent(out) :: gvec_int(:,:)


    integer :: i
    integer :: stat

    allocate(gvec_int(ngvec,3),stat=stat)
    if(stat/=0) error stop 'casino_G_to_int: Failed to allocate integer G-grid'

    ! Multiply each G-vector component by corresponding reciprocal lattice vector
    ! making sure to multiply by pre-factor of 1/2*pi
    do i=1,ngvec
       gvec_int(i,:) = nint(matmul(user_params%platt,gvecs(i,:))/(2.0_dp*pi))
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

  subroutine casino_check_symmetry(dens_gs, gvec_int,have_inv_sym, label)
    !============================================================!
    ! This routine checks the CASINO density to ensure that we   !
    ! have the correct expected symmetries.                      !
    !------------------------------------------------------------!
    ! Arguments                                                  !
    ! dens_gs(in)  : the Fourier components of the density       !
    !                (either total charge/spin density)          !
    ! gvec_int(in) : G-vectors as integer multiples of recip     !
    !                 lattice vectors.                           !
    ! have_inv_sym(out), optional : Does the density             !
    !                     have inversion symmetry                !
    ! label(in,optional): label for density                      !
    !                     C-charge                               !
    !                     S-spin                                 !
    !------------------------------------------------------------!
    ! Modules used                                               !
    ! Constants,Basis                                            !
    !------------------------------------------------------------!
    ! Necessary Conditions                                       !
    ! basis_initialised should have been called.                 !
    ! casino_read must be called before calling this routine.    !
    !============================================================!
    use constants,only : cmplx_0
    use basis,    only : castep_basis

    implicit none
    complex(kind=dp),intent(in)  :: dens_gs(:)
    integer,intent(in)           :: gvec_int(:,:)
    logical,intent(out),optional :: have_inv_sym    ! Does density have inversion symmetry
    character          ,optional :: label           ! label to use when writing output

    character :: l_label
    character(6) :: tmpstr

    complex(kind=dp),allocatable :: den_grid(:,:,:) ! Temporary complex grid for symmetries
    integer   :: nx,ny,nz    ! max symmetric grid length along x,y,z
    logical   :: l_inv_sym   ! local copy of inv_sym
    integer   :: i           ! loop counters
    integer   :: stat

    l_label = 'C'

    ! Set the label to use when writing out the density
    if (present(label)) l_label = label
    select case (l_label)
    case('C')
       tmpstr = 'charge'
    case('S')
       tmpstr = 'spin'
    end select

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
       den_grid(gvec_int(i,1),gvec_int(i,2),gvec_int(i,3)) = dens_gs(i)
    end do

    ! First check for complex conjugate symmetry - this is REQUIRED or we have done something wrong!
    if(.not.casino_check_cmplx_conj(den_grid,nx,ny,nz)) then
       error stop 'Reciprocal space density does not appear to have complex conjugate symmetry'
    end if

    ! Now check for inversion symmetry - not necessary unless it's actually there!
    l_inv_sym = casino_check_inver_sym(den_grid,nx,ny,nz)
    if(.not.l_inv_sym) then
       write(stdout,'(A)') ' Reciprocal space '//trim(tmpstr)//' density does not have inversion symmetry   '
    else
       write(stdout,'(A)') ' SYMMETRY DETECTED: Reciprocal space '//trim(tmpstr)//' density has inversion symmetry'
    end if

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
               if(.not. math_isclose(den_grid(ix,iy,iz), conjg(den_grid(-ix,-iy,-iz))) ) then
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
    !                                                            !
    ! We deal with this slight difference of convention in this  !
    ! routine!!!!                                                !
    !============================================================!
    use density,only : elec_density, density_allocate
    use basis,  only : castep_basis

    implicit none

    type(elec_density),intent(inout) :: den
    integer,intent(in)    :: gvec_int(:,:)

    integer :: ix,iy,iz
    integer :: i

    ! TODO - Seeing as its real for inversion symmetry, we can reduce memory costs and gain speed for FFT using the corresponding real routines.
    ! Allocate and initialise densities with zero
    call density_allocate(den,nspins)

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
          den%charge(ix,iy,iz) = charge_gs(i)
          if(nspins==2) den%spin(ix,iy,iz) = spin_gs(i)
       else
          den%real_charge(ix,iy,iz) = real(charge_gs(i),dp)
          if(nspins==2) den%real_spin(ix,iy,iz) = real(spin_gs(i),dp)
       end if
    end do
  end subroutine casino_read_recip_grid

  subroutine casino_truncate_gs(cutoff,write_gs)
    !============================================================!
    ! Sets all Fourier components of the charge density outside  !
    ! the sphere of radius |G_cut| in reciprocal space to zero.  !
    !------------------------------------------------------------!
    ! Arguments                                                  !
    ! cutoff(in) : cutoff radius to truncate the G-vectors at    !
    ! write_gs(in,optional) : write out the truncated G-vectors  !
    !                         to a file.  (Default: False)       !
    !------------------------------------------------------------!
    ! Modules used                                               !
    ! constants, io                                              !
    !------------------------------------------------------------!
    ! Necessary conditions                                       !
    ! casino_read should be called BEFORE this routine so that   !
    ! charge_gs and spin_gs are allocated and initialised.       !
    !============================================================!

    use constants, only : cmplx_0, eV
    use io,        only : stderr

    implicit none
    ! Arguments
    real(dp),intent(in)           :: cutoff    ! Planewave cutoff |G_cut|^2/2 after which Fourier components should be set to zero (Hartrees)
    logical, intent(in), optional :: write_gs  ! write out truncated G-vectors to a file

    ! Local Variables
    character(len=*), parameter :: g_file='trun_rho_g.txt'
    logical :: l_write_gs

    real(dp):: g_ke               ! Planewave cutoff of a given G-vector
    real(dp):: ke_max             ! Maximum planewave cutoff of the given set of G-vectors (i.e. what the calculation used)
    integer :: ndel               ! number of Fourier components we removed

    integer :: ng

    l_write_gs = .false.
    if (present(write_gs)) l_write_gs=write_gs

    ! Check that we have the data in the right format allocated already.
    if(.not.allocated(charge_gs)) error stop 'casino_truncate_gs: Reciprocal charge Fourier components not allocated/stored.'
    if(.not.allocated(spin_gs).and.nspins==2) &
         error stop 'casino_truncate_gs: Reciprocal spin Fourier components not allocated/stored.'

    if(nspins==2) then
       write(stdout,'(/,A)') ' Truncating Fourier components for charge density and spin density.'
    else
       write(stdout,'(/,A)') ' Truncating Fourier components for charge density.'
    end if
    write(stdout,100) ' Applying cutoff ', cutoff, cutoff*eV
    ke_max = 0.0_dp
    ndel = 0
    do ng=1,ngvec
       ! Calculate the kinetic energy associated with each G-vector
       g_ke = sum(gvecs(ng,:)**2.0_dp)/2.0_dp
       ! Note the current maximum cutoff we have.
       ! When we do reach the actual maximum, this is the cutoff of the calculation.
       if(g_ke>ke_max) ke_max = g_ke

       if (g_ke>cutoff) then
          ! If the G-vector has a radius (energy) outside the cutoff radius, set that Fourier component to zero
          ndel = ndel + 1
          charge_gs(ng) = cmplx_0
          if (nspins==2) spin_gs(ng) = cmplx_0
       end if
    end do ! ng

    ! Now check that the cutoff specified by the user is smaller than the calculation's cutoff,
    ! otherwise, nothing happens...
    write(stdout,100) ' Calculation cutoff ', ke_max, ke_max*eV
    if (ke_max<=cutoff) then
       write(stderr,'(A62,/,A13)') ' Your planewave cutoff is larger than that of the calculation!',' No effect...'
    else
       write(stdout,101) ndel, real(ndel,dp)/real(ngvec,dp)
    end if
100 format(A,T30,'= ',F10.2,' Hartrees= ',F14.2, ' eV')
101 format(' Number of G-vectors removed: ',I6,/,T30,'= ',F8.4, ' of all G-vectors')

    ! Write out truncated G-vectors to file
    if (l_write_gs) then
       write(stdout,'(A46,A,/)') ' Writing out truncated Fourier components to: ',trim(g_file)
       call casino_plot_Gs(charge_gs,g_file)
    end if
  end subroutine casino_truncate_gs
end module casino
