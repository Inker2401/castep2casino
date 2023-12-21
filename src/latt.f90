module latt
  !===============================================================================!
  !                       _          _   _   _                                    !
  !                      | |    __ _| |_| |_(_) ___ ___                           !
  !                      | |   / _` | __| __| |/ __/ _ \                          !
  !                      | |__| (_| | |_| |_| | (_|  __/                          !
  !                      |_____\__,_|\__|\__|_|\___\___|                          !
  ! This module handles dealing with the contents of the lattice_geom file        !
  ! effectively user input and importantly the reading of the primitive lattice   !
  ! vectors.                                                                      !
  !-------------------------------------------------------------------------------!
  ! Modules used                                                                  !
  ! IO, Constants                                                                 !
  !-------------------------------------------------------------------------------!
  ! Public variables:                                                             !
  ! user_params :: the public instance of the params type containing user         !
  !                user input parameters.                                         !
  !===============================================================================!

  use constants, only : dp
  use io, only        : file_maxpath

  private
  !---------------------------------------------------------------------------!
  !                       P u b l i c   V a r i a b l e s                     !
  !---------------------------------------------------------------------------!
  type params
     real(kind=dp),dimension(3,3)   :: platt               ! lattice vector, component
     real(kind=dp),dimension(3)     :: lat_consts          ! lattice constants (in Bohr) a,b,c
     real(kind=dp),dimension(3)     :: lat_angles          ! lattice angles (in degrees) alpha,beta,gamma
     character(len=file_maxpath)    :: den_fmt_file        ! name of formatted density file
     integer                        :: ngx,ngy,ngz         ! CASTEP grid size
     logical                        :: shift_grid          ! Does user want to shift real space grid?
     real(kind=dp),allocatable      :: shift_frac(:)       ! the shift applied to the user's grid in fractional coordinates along x,y,z
     logical                        :: write_mathematica   ! write data to file as a Mathematica list
  end type params

  type(params),public,save :: user_params   ! public instance of params type
  !---------------------------------------------------------------------------!
  !                       P u b l i c   R o u t i n e s                       !
  !---------------------------------------------------------------------------!

  public :: latt_read

  !---------------------------------------------------------------------------!
  !                     P r i v a t e   V a r i a b l e s                     !
  !---------------------------------------------------------------------------!
contains
  subroutine latt_read(filename)
    !============================================================!
    ! Reads the latt_geom file and gets the user's input         !
    ! parameters.                                                !
    ! First we read the size of the user's FFT grid and then     !
    ! read the primitive lattice vectors                         !
    !------------------------------------------------------------!
    ! Arguments                                                  !
    ! filename(in) : the name of the latt_geom file              !
    !------------------------------------------------------------!
    ! Modules used                                               !
    ! constants, IO                                              !
    !------------------------------------------------------------!
    ! Necessary Conditions                                       !
    ! filename must point to a FORMATTED file                    !
    !============================================================!
    use constants, only  : bohr_radius
    use io, only         : file_maxpath,stdout,stderr,io_file_present,io_file_code,io_read_block,&
                           io_strip_extension
    implicit none

    ! File variables
    character(len=*),intent(in) :: filename
    character(len=file_maxpath),allocatable :: tmpc(:)

    integer :: i
    character(len=60) :: iomsg
    integer :: unit,iostat,stat

    ! Open the file
    open(file=trim(filename),newunit=unit,status='OLD',action='READ',&
         iostat=iostat,iomsg=iomsg)
    if(iostat/=0) then
       write(stderr,'(A7,A)') 'ERROR: ',trim(iomsg)
       stop
    end if

    ! First read the size of the user's grid
    allocate(tmpc(1))
    if (io_file_present(unit,'castep_grid')) then
       tmpc = trim(io_file_code(unit,'castep_grid',whole=.true.))
       read(tmpc,*,iostat=iostat) user_params%ngx,user_params%ngy,user_params%ngz
       if(iostat/=0) error stop 'latt_read: Failed to read CASTEP grid size'
       write(stdout,'(A20,10x,3I4)') ' CASTEP grid size:  ',user_params%ngx,user_params%ngy,user_params%ngz
    else
       stop 'ERROR: Did not specify CASTEP grid size'
    end if
    deallocate(tmpc)

    ! Check if the grid contains an odd number of grid points
    call latt_check_grid(user_params%ngx,user_params%ngy,user_params%ngz,required=.false.)

    ! Obtain the user's primitive lattice vectors
    call io_read_block(unit,'prim_latt_cart',tmpc)
    if (size(tmpc)/=3) stop 'ERROR: prim_latt_cart block incorrectly specified.'
    do i=1,3
       read(tmpc(i),*,iostat=iostat) user_params%platt(i,1), user_params%platt(i,2), user_params%platt(i,3)
       if(iostat/=0) error stop 'latt_read: Failed to read prim_latt_cart block'
    end do
    deallocate(tmpc)

    ! Now convert the lattice vectors to Bohr (if required)
    if (.not.io_file_present(unit,'unit_bohr')) then
       user_params%platt = user_params%platt/bohr_radius
    end if

    write(stdout,'(15x,"Real Lattice (",A1,")", 40x, "Real Lattice (",A4,")")') 'A', 'Bohr'
    do i=1,3
       write(stdout,100) user_params%platt(i,:)*bohr_radius, user_params%platt(i,:)
    end do
100 format(3f16.10,5x,3f16.10)

    ! Calculate the lattice constants - NB output is handled within latt_calc_const_and_angles
    call latt_calc_const_and_angles()

    ! Get density output file
    user_params%den_fmt_file = trim(io_strip_extension(trim(filename)))//'.den_fmt'
    if (io_file_present(unit,'output_file')) then
       user_params%den_fmt_file = trim(io_file_code(unit,'output_file',code_keep_case=.true.))
    end if

    ! Check if user wants to perform a shift of grid for visualisation purposes
    user_params%shift_grid=.false.
    if (io_file_present(unit,'shift_grid')) then
       allocate(tmpc(1))
       allocate(user_params%shift_frac(3),stat=stat)
       if(stat/=0) error stop 'Failed to allocate shift_frac'

       ! Get shift in fractional coordinates
       tmpc = trim(io_file_code(unit,'shift_grid',whole=.true.))
       read(tmpc,*,iostat=iostat) user_params%shift_frac(1),user_params%shift_frac(2),user_params%shift_frac(3)
       deallocate(tmpc)

       ! Set flag to true
       user_params%shift_grid=.true.
       write(stdout,'(A18,3F13.6)') ' Real Space Shift:', user_params%shift_frac
    end if
    ! write(*,*) 'shift_user_grid: ', user_params%shift_grid

    ! Write real space density as Mathematica list - WRITE_MATHEMATICA
    user_params%write_mathematica=io_file_present(unit,'write_mathematica')

    write(stdout,*) ''
    ! Finished reading parameters
    close(unit,iostat=iostat)
    if(iostat/=0) error stop 'casino_read: Failed to close lattice geometry file.'
  end subroutine latt_read

  subroutine latt_check_grid(ngx,ngy,ngz,required)
    !============================================================!
    ! Checks if the user specified an odd number of points along !
    ! each CASTEP grid dimension. The reason is because the      !
    ! the Fourier components are contained within a SPHERE       !
    ! of reciprocal space of radius 2G.                          !
    ! This means the CASTEP grid should effectively run along    !
    ! e.g. -x0 to x0 which contains 2x0+1 grid points.           !
    !------------------------------------------------------------!
    ! Arguments                                                  !
    ! ngx,ngy,ngz(in) :: number of points along each grid dim    !
    ! required(in) :: Do we halt execution and return an error   !
    !                 or just return a warning                   !
    !                 (Default: True - return error and stop)    !
    !------------------------------------------------------------!
    ! Modules used                                               !
    ! IO                                                         !
    !============================================================!
    use io, only : stderr
    implicit none
    integer,intent(in)          :: ngx,ngy,ngz  ! user's CASTEP grid size along x, y and z
    logical,intent(in),optional :: required     ! do we raise this as an error(true) just a warning(false)
                                                ! (Default: true)
    logical :: bad_vals    ! are the number of grid points sufficient
    logical :: l_require
    l_require = .true.
    if(present(required))l_require=required

    ! Check if grids contain odd number of grid points
    ! (so that x runs from say -x to x contains 2x+1 points which is manifestly odd).
    bad_vals=.false.
    if(mod(ngx,2)==0) bad_vals=.true.
    if(mod(ngy,2)==0) bad_vals=.true.
    if(mod(ngz,2)==0) bad_vals=.true.

    ! Now warning message or stop program
    if (bad_vals) then
       if(l_require) then
          write(stderr,'(A50)') ' ERROR: CASTEP grid dimensions should be odd.     '
          stop

       else
          write(stderr,'(A50)') ' WARNING: CASTEP grid dimensions should be odd.   '
       end if
    end if
  end subroutine latt_check_grid

  subroutine latt_calc_const_and_angles(silent)
    !============================================================!
    ! Calculates the lattice angles and lattice constants from   !
    ! lattice vectors.                                           !
    !------------------------------------------------------------!
    ! Arguments                                                  !
    !------------------------------------------------------------!
    ! silent(in,optional) :: do not output calculated results    !
    !                        (Default : False)                   !
    !------------------------------------------------------------!
    ! Modules used                                               !
    ! math,constants,stdout                                      !
    !------------------------------------------------------------!
    ! Necessary Conditions                                       !
    ! The primitive lattice vectors contained in                 !
    ! user_params%platt must contain the valid set of            !
    ! lattice vectors (in Bohr!!) before this routine is called. !
    !============================================================!
    use math, only : math_get_vec_angle
    use constants, only : bohr_radius
    use io, only : stdout

    implicit none
    logical,intent(in),optional :: silent
    logical :: l_silent

    integer :: dir
    character(len=1), dimension(3), parameter :: side_labels=(/ 'a','b','c' /)
    character(len=5), dimension(3), parameter :: angle_labels=(/ 'alpha','beta ','gamma' /)

    l_silent=.false.
    if(present(silent))l_silent=silent

    ! Norm of lattice vector gives the lattice constant
    do dir=1,3
       user_params%lat_consts(dir) = sqrt(sum(user_params%platt(dir,:)**2))
    end do

    ! Now get angles
    user_params%lat_angles(1) = math_get_vec_angle(user_params%platt(2,:),user_params%platt(3,:))
    user_params%lat_angles(2) = math_get_vec_angle(user_params%platt(1,:),user_params%platt(3,:))
    user_params%lat_angles(3) = math_get_vec_angle(user_params%platt(1,:),user_params%platt(2,:))

    if (.not.l_silent) then
       write(stdout,*) ''
       write(stdout,'(35x,A18)') 'Lattice Parameters'
       do dir=1,3
          write(stdout,101) side_labels(dir), user_params%lat_consts(dir)*bohr_radius, user_params%lat_consts(dir), &
               angle_labels(dir), user_params%lat_angles(dir)
       end do
       write(stdout,*) ''
    end if
101 format(10x,A1,' = ',f14.8,' A  = ',f14.8,' Bohr, ',A5,' = ',f12.6)
  end subroutine latt_calc_const_and_angles
end module latt
