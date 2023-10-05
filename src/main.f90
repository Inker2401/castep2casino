program castep2casino
  !=================================================================================!
  ! This program takes a CASINO density and turns it into a format that CASTEP can  !
  ! understand. The following are the main steps that happen:                       !
  !                                                                                 !
  ! 1. Get the size of the user's CASTEP grid and then initialise the grid and FFT  !
  !    plan. This step must be done BEFORE any further input arrays                 !
  !    that will be Fourier transformed are initialised.                            !
  !                                                                                 !
  ! 2. Obtain the primitive lattice vectors from the user's cell file               !
  !                                                                                 !
  ! 3. Read the CASINO file and obtain the G-vectors and Fourier components.        !
  !                                                                                 !
  ! 4. Store the G-vectors onto reciprocal space grid taking into account the       !
  !    the array ordering that FFTW expects.                                        !
  !                                                                                 !
  ! 5. Perform the FFT to get the density to real space before writing to CASTEP    !
  !---------------------------------------------------------------------------------!
  ! Important variables                                                             !
  !---------------------------------------------------------------------------------!
  ! Written by : Visagan Ravindran, Stewart Clark                                   !
  !              Department of Physics, University of Durham                        !
  !---------------------------------------------------------------------------------!
  use io,       only : file_maxpath, stdout
  use latt,     only : latt_read,ngx,ngy,ngz
  use basis,    only : basis_initialise,basis_deallocate
  use casino,   only : casino_read,casino_to_castep
  use density,  only : elec_density,density_allocate,density_deallocate,density_recip_to_real,density_write

  implicit none

  character(len=file_maxpath) :: casino_file, latt_file ! CASINO file, lattice geometry (user options) file
  type(elec_density)          :: recip_den,real_den     ! the density on a grid in reciprocal space, real space

  ! Get arguments from the command line
  call get_input_files(casino_file,latt_file)

  ! Open the user options (latt_geom) file and obtain user parameters like grid etc.
  call latt_read(trim(latt_file))

  ! Initialise FFT grids
  call basis_initialise(ngx,ngy,ngz)

  ! Read CASINO file
  call casino_read(trim(casino_file))

  ! Store CASINO density to CASTEP grid - still in reciprocal space
  call casino_to_castep(recip_den)

  ! Perform FFT from reciprocal space to real space
  ! TODO - if we know it is going to be real, e.g. inversion symmetry, should allocate real instead
  call density_allocate(real_den)
  call density_recip_to_real(recip_den,real_den)
  ! Done with reciprocal space density, can deallocate it now.
  call density_deallocate(recip_den)

  ! Sanity checks on real space density
  call check_real_den(real_den)
  write(stdout, '(A21,F12.5)') ' Number of electrons:', real_den%norm()
  if (real_den%have_cmplx_den) then
     write(stdout,'(A18,3x,F18.10)') ' Maximum density: ',maxval(abs(real_den%charge))
     write(stdout,'(A18,3x,F18.10)') ' Minimum density: ',minval(abs(real_den%charge))
  else
     write(stdout,'(A18,3x,F18.10)') ' Maximum density: ',maxval(real_den%real_charge)
     write(stdout,'(A18,3x,F18.10)') ' Minimum density: ',minval(real_den%real_charge)
  end if

  ! Write real space density in CASTEP format.
  ! TODO - rather than specifying format manually,density should always be real so we should only need to allocate it as a real type - see previous TODO
  call density_write(real_den,fmt='R')

  ! Clean up
  call density_deallocate(real_den)
  ! Destroy FFT plan and deallocate FFT grids.
  call basis_deallocate()

contains

  subroutine get_input_files(casino_file,latt_file)
    !============================================================!
    ! Gets the file names for the input files from the           !
    ! command line.                                              !
    ! If a help flag is specified, print_help is called          !
    ! and a help message is displayed.                           !
    ! If only a single argument is passed on the command line,   !
    ! the user's latt geom file is assumed to have the same      !
    ! seedname as the CASINO file.                               !
    !============================================================!
    use io, only : io_strip_extension
    implicit none
    character(len=file_maxpath),intent(out) :: casino_file, latt_file ! CASINO file, lattice geometry (user options) file

    if (command_argument_count() == 0) then
       ! Print help if program is called by itself.
       call print_help()
    elseif (command_argument_count() == 1) then
       ! Assume user lattice geometry file is CASINO seedname with appropriate file extension.

       call get_command_argument(1,casino_file)
       ! Check if help flags specified.
       if (trim(casino_file)=='-h' .or. trim(casino_file).eq.'--help') call print_help()
       ! Remove file extension
       latt_file = trim(io_strip_extension(trim(casino_file)))//'.latt'
    else
       call get_command_argument(1,casino_file)
       ! Check if help flags specified.
       if (trim(casino_file)=='-h' .or. trim(casino_file).eq.'--help') call print_help()
       call get_command_argument(2,latt_file)
    end if

    write(stdout,'(A32,A)') ' Reading CASINO File:           ', trim(casino_file)
    write(stdout,'(A32,A)') ' Reading user parameters from:  ', trim(latt_file)
  end subroutine get_input_files

  subroutine print_help()
    !============================================================!
    ! Outputs a help message to STDOUT and then terminates the   !
    ! program.                                                   !
    !============================================================!
    implicit none
    write(stdout,'(A60)') 'Usage:                                                      '
    write(stdout,'(A60)') ' castep2casino.e [-h/--help] <casino_file> [lat_geom_file]  '
    write(stdout,'(A60)') ' -h, --help   : print this help message                     '
    write(stdout,'(A60)') repeat(' ',60)
    write(stdout,'(A60)') repeat(' ',60)
    write(stdout,'(A60)') 'Lat Geom File keywords:                                     '
    write(stdout,'(A60)') ' castep_grid  : specify number of grid points for CASTEP    '
    write(stdout,'(A60)') ' prim_lat_vec : primitive lattice vectors                   '
    write(stdout,'(A60)') '                specified as block (see README.md)          '
    write(stdout,'(A60)') ' latt_bohr    : specify lattice vectors in Bohr,            '
    write(stdout,'(A60)') '                otherwise will use angstroms instead.       '
    write(stdout,'(A60)') ' output_file  : name of output file to use                  '
    write(stdout,'(A60)') '             Default - seedname from <casino_file>.den_fmt  '
    stop
  end subroutine print_help

  subroutine check_real_den(real_den)
    !============================================================!
    ! Checks if the real space density is real, i.e. has no      !
    ! imaginary components. Due to some noise from the FFT, if   !
    ! it was allocated as a complex density, there might be      !
    ! very small imaginary components of the order of 10**-18.   !
    ! Therefore, check is done using tolerance so that they are  !
    ! close to 0, not strictly equal to 0.                       !
    !------------------------------------------------------------!
    ! Arguments                                                  !
    ! real_den(in) :: the real space density on a grid           !
    !------------------------------------------------------------!
    ! Necessary conditions:                                      !
    ! Should go without saying but this assumes the density      !
    ! actually is in real space, i.e. call density_real_to_recip !
    ! first before using this!!!                                 !
    !============================================================!
    use constants, only : dp
    use math,only  : math_isclose
    use basis,only : castep_basis

    implicit none

    type(elec_density),intent(in) :: real_den

    real(kind=dp),allocatable :: flat_grid(:) ! flattened real space grid containing imaginary parts of density values
    real(kind=dp),allocatable    :: zero_arr(:)  ! array of zeros

    integer :: stat

    if (real_den%have_cmplx_den) then
       allocate(flat_grid(castep_basis%total_grid_points),stat=stat)
       if(stat/=0) error stop 'check_real_den: Failed to allocated flattened grid.'
       allocate(zero_arr((castep_basis%total_grid_points)))
       zero_arr = 0.0_dp

       ! Flatten grid and store store imaginary parts only since that's what we want to check
       flat_grid = aimag(reshape(real_den%charge, (/castep_basis%total_grid_points/)))

       ! Check no non-zero imaginary parts in density (within tolerance!)
       if (.not. math_isclose(flat_grid, zero_arr)) then
          write(stdout,'(A53)') ' WARNING: Real space density appears to be complex...'
       ! else
       !    write(stdout,*) 'Real space density is fine...'
       end if

       deallocate(flat_grid,stat=stat)
       if(stat/=0) error stop 'check_real_den: Failed to deallocate flattened grid.'
       deallocate(zero_arr)
    else
       ! Nothing to check as it is already real
       return
    end if
  end subroutine check_real_den
end program castep2casino
