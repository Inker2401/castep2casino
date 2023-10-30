module io
  !===============================================================================!
  !                                 __ ___                                        !
  !                                |_ _/ _ \                                      !
  !                                 | | | | |                                     !
  !                                 | | |_| |                                     !
  !                                |___\___/                                      !
  ! This module contains some additional IO utilities to supplement those         !
  ! to supplement those that are already available in Fortran.                    !
  ! As a rule of thumb, all the routines here are case-sensitive unless           !
  ! explicitly stated otherwise.                                                  !
  ! ------------------------------------------------------------------------------!
  ! Modules used                                                                  !
  ! standard Fortran environment: iso_fortran_env                                 !
  ! ------------------------------------------------------------------------------!
  ! Public variables:                                                             !
  ! stdout: the unit to write to standard output                                  !
  ! stderr: the unit to write error messgaes to                                   !
  ! stderr: EOF error code                                                        !
  ! file_maxpath: max line length and path to a file                              !
  !===============================================================================!
  use iso_fortran_env, only : output_unit, error_unit, iostat_end

  implicit none

  !---------------------------------------------------------------------------!
  !                       P u b l i c   V a r i a b l e s                     !
  !---------------------------------------------------------------------------!
  ! Define IO constants
  integer,public,parameter :: stdout=output_unit
  integer,public,parameter :: stderr=error_unit
  integer,public,parameter :: EOF=iostat_end       ! status code for having reached end-of-file
  integer,public,parameter :: file_maxpath=256     ! maximum length of a line in a file / path to file

contains

  function io_strip_extension(filename) result(seedname)
    !============================================================!
    ! This routine removes the filename from a file and returns  !
    ! the seedname.                                              !
    ! For example: my_file.txt causes this routine to            !
    !              return my_file                                !
    !------------------------------------------------------------!
    ! Arguments                                                  !
    ! filename(in) : full filename                               !
    ! result(out)  : filename with extension removed             !
    !============================================================!
    implicit none
    character(len=*) :: filename            ! full filename
    character(len=file_maxpath) :: seedname ! filename without file extension
    integer :: extpos  ! extension position in filename

    ! Find position of last full stop and assume that's the file extension.
    extpos = scan(trim(filename),'.',BACK=.true.)

    seedname = filename(:extpos-1)
  end function io_strip_extension

  function io_lowercase(str) result(lower_str)
    !============================================================!
    ! This routine takes a string and turns everything into      !
    ! lower case                                                 !
    !------------------------------------------------------------!
    ! Arguments                                                  !
    ! str(in) :: the string to turn into lower case              !
    !============================================================!
    implicit none
    character(len=*),intent(in) :: str
    character(len=len(str))     :: lower_str

    integer :: ia, iz, ic, ishift
    integer :: i

    ! Get indexes in for lower case characters in character set
    ia=ichar('A')
    iz=ichar('Z')
    ishift = ichar('a') - ia

    do i=1,len(str)
       ic = ichar(str(i:i))

       if((ic>=ia .and. ic<=iz)) then
          ! If upper case, convert to lower case
          lower_str(i:i) = char(ic+ishift)
       else
          lower_str(i:i) = str(i:i)
       end if
    end do

    ! Check for non ASCII characters (note each is 2 bytes)
    do i=1,len(lower_str)-1
       select case(lower_str(i:i+1))
       ! German umlauts and scharfes S
       case('Ä')
          lower_str(i:i+1)='ä'
       case('Ö')
          lower_str(i:i+1)='ö'
       case('Ü')
          lower_str(i:i+1)='ü'
       case('ẞ')
          lower_str(i:i+1)='ß'
       ! Acute accents
       case('Á')
          lower_str(i:i+1)='á'
       case('É')
          lower_str(i:i+1)='é'
       case('Í')
          lower_str(i:i+1)='í'
       case('Ó')
          lower_str(i:i+1)='ó'
       case('Ú')
          lower_str(i:i+1)='ú'
       ! Grave accents
       case('À')
          lower_str(i:i+1)='à'
       case('È')
          lower_str(i:i+1)='è'
       case('Ì')
          lower_str(i:i+1)='ì'
       case('Ò')
          lower_str(i:i+1)='ò'
       case('Ù')
          lower_str(i:i+1)='ù'
       ! Cedilla
       case('Ç')
          lower_str(i:i+1)='ç'
       ! Circumflex accents
       case('Â')
          lower_str(i:i+1)='â'
       case('Ê')
          lower_str(i:i+1)='ê'
       case('Î')
          lower_str(i:i+1)='î'
       case('Ô')
          lower_str(i:i+1)='ô'
       case('Û')
          lower_str(i:i+1)='û'
       ! Spanish eñe
       case('Ñ')
          lower_str(i:i+1)='ñ'
       end select
   end do
  end function io_lowercase

  subroutine io_skip_header(unit,header)
    !============================================================!
    ! This routine skips past the header or to a point in a file !
    ! based on a search string given in the variable 'header'    !
    ! NOTE that the search is case-insensitive.                  !
    !------------------------------------------------------------!
    ! Arguments                                                  !
    ! unit (in)  :: unit for file we want to skip past           !
    ! header(in) :: search string we want to skip past           !
    !------------------------------------------------------------!
    ! Necessary Conditions                                       !
    ! File must be OPEN and FORMATTED                            !
    !============================================================!
    implicit none
    ! Arguments
    integer, intent(in)          :: unit   ! unit for file we want to read
    character(len=*), intent(in) :: header ! string we want to skip past in the file

    ! Local variables
    character(len=file_maxpath) :: filename
    logical :: opened
    character(len=file_maxpath) :: line

    integer :: iostat

    ! Check if the file is open
    inquire(unit=unit, opened=opened, name=filename)
    if (.not. opened) then
       write(stderr,*) trim(filename), ' is not open'
       error stop 'io_skip_header: File was not opened.'
    end if

    ! Rewind the file
    rewind(unit=unit, iostat=iostat)
    if(iostat/=0) error stop 'io_skip_header: Failed to rewind file.'

    ! Skip file contents to point in header
    do while(index(trim(adjustl(line)), trim(adjustl(header))) == 0)
       read(unit,'(A)',iostat=iostat) line
       if (iostat/=0) error stop 'io_skip_header: Unable to read line in file'
    end do
  end subroutine io_skip_header

  subroutine io_read_block(unit,blockname,blockdata,startblock,endblock)
    !============================================================!
    ! This routine reads the text contained within a block       !
    ! delimited by start_block and end_block.                    !
    ! For example,                                               !
    ! start_block blockname                                      !
    ! ...                                                        !
    ! block contents                                             !
    ! ...                                                        !
    ! end_block blockname                                        !
    ! If the block was unable to be found, an error is raised    !
    ! and program execution is halted.                           !
    ! NOTE that this routine is CASE-SENSITIVE.                  !
    !------------------------------------------------------------!
    ! Arguments                                                  !
    ! unit (in) : file unit to read block                        !
    ! blockname (in) : name of block to read                     !
    ! blockdata (out): contents of block                         !
    ! startblock (in): starting delimter for block               !
    !                  (default: '%block' )                      !
    ! endblock (in): ending delimter for block                   !
    !                  (default: '%endblock' )                   !
    !------------------------------------------------------------!
    ! Necessary Conditions                                       !
    ! File must be OPEN and FORMATTED                            !
    !============================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: unit ! data to read block from
    character(len=*), intent(in) :: blockname ! name of block
    character(len=file_maxpath), allocatable, intent(out) :: blockdata(:) ! contents of block(line)
    character(len=file_maxpath), intent(in), optional :: startblock ! start delimiter of block
    character(len=file_maxpath), intent(in), optional :: endblock   ! end delimiter of block

    ! Local Variables
    character(len=file_maxpath) :: filename
    character(len=file_maxpath) :: local_startblock, local_endblock

    integer :: nlines ! number of lines inside the block
    integer :: n      ! line counter

    character(len=file_maxpath) :: line, str1, str2
    logical :: have_startblock, have_endblock, opened

    integer :: iostat, stat
    local_startblock='%block'
    if(present(startblock)) local_startblock=startblock
    local_endblock='%endblock'
    if(present(endblock)) local_endblock=endblock

    ! Check if the file is open
    inquire(unit=unit, opened=opened, name=filename)
    if (.not. opened) then
       error stop 'io_read_block: File was not opened.'
    end if

    ! Rewind the file to the beginning.
    rewind(unit=unit,iostat=iostat)
    if(iostat/=0) error stop 'io_read_block: Unable to perform initial rewind of file.'

    have_startblock=.false.
    have_endblock=.false.
    nlines = 0

    ! Read the block
    do
       read(unit,'(A)',iostat=stat) line
       if (stat==EOF) error stop 'io_read_block: Reached end of file but could not find block.'

       read(line, *, iostat=iostat) str1, str2
       if(iostat==0) then
          ! Check if we are in a block
          if(io_is_block_delim(line,local_startblock,blockname)) then
             have_startblock = .true.
             cycle
          elseif(io_is_block_delim(line,local_endblock,blockname)) then
             have_endblock = .true.
             exit
          endif
       end if

       ! Start counting lines while we are inside the block
       if (have_startblock .and. .not. have_endblock) then
          nlines = nlines + 1
       end if
    end do

    ! Check if we reached end of file without encountering the block
    if (.not. have_startblock) then
       write(stderr,*) 'ERROR: Block ', trim(blockname), ' in ' , trim(filename), 'is should be a block'
       stop
    end if

    ! Check that it was closed
    if (.not. have_endblock) then
       write(stderr,*) 'ERROR: Block ', trim(blockname), ' in ' , trim(filename), 'is should be closed'
       stop
    end if
    ! write(*,*) 'Number of lines in block', nlines

    ! Rewind the file and go back to the block's position
    rewind(unit=unit,iostat=iostat)
    if(iostat/=0) error stop 'io_read_block: Unable to rewind file to block positon'
    have_startblock=.false.
    do
       read(unit,'(A)',iostat=stat) line
       read(line, *, iostat=iostat) str1, str2
       if(iostat==0) then
          if(io_is_block_delim(line,local_startblock,blockname)) then
             have_startblock = .true.
             exit
          endif
       end if
    end do

    ! Allocate the array and now store the block contents
    allocate(blockdata(nlines), stat=stat)
    do n=1,nlines
       read(unit,'(A)',iostat=iostat) blockdata(n)
       if(iostat/=0) error stop 'Unable to store block data'
    end do
    if(stat/=0) error stop 'io_read_block: Unable to allocate block data'

  end subroutine io_read_block

  logical function io_is_block_delim(line, delimiter, blockname)
    !============================================================!
    ! This checks if a line in a file indicates a block delimiter!
    ! i.e. is the line DELIMITER BLOCKNAME                       !
    ! NOTE: This is a CASE-SENSITIVE routine                     !
    !------------------------------------------------------------!
    ! line :: the current line in the file (or any string)       !
    ! delimiter :: the string that marks the line as a block,    !
    !              see above                                     !
    ! blockname :: the name of the block, see above              !
    !------------------------------------------------------------!
    ! Known issues                                               !
    ! None though this function does not need to be called with  !
    ! trim in the arguments as they are trimmed in here.         !
    !============================================================!
    implicit none

    ! Arguments
    character(len=*), intent(in) :: line       ! current line in file
    character(len=*), intent(in) :: delimiter  ! block delimiter (e.g. %block)
    character(len=*), intent(in) :: blockname  ! name of block
    character(len=file_maxpath) :: str1, str2
    integer :: iostat

    read(line, *, iostat=iostat) str1, str2

    io_is_block_delim = .false.
    if(iostat==0) then
       if (trim(str1) .eq. trim(delimiter) .and. trim(str2) .eq. trim(blockname)) then
          io_is_block_delim = .true.
          return
       endif
    end if

  end function io_is_block_delim

  logical function io_str_present(str, key)
    !============================================================!
    ! Checks if a string contains a substring/key                !
    ! NOTE: This is a CASE-SENSITIVE routine                     !
    !------------------------------------------------------------!
    ! str :: the string to search                                !
    ! key :: the substring to search for within str              !
    !------------------------------------------------------------!
    ! Known issues                                               !
    ! None though this function does not need to be called with  !
    ! trim in the arguments as they are trimmed in here.         !
    !============================================================!
    implicit none
    ! Arguments
    character(len=*) :: str   ! The string to search
    character(len=*) :: key   ! The substring to search for within the main string

    ! Is the key is contained within str?
    io_str_present = index(trim(str), trim(key))>0
  end function io_str_present

  character(len=file_maxpath) function io_str_code(str, key, whole)
    !============================================================!
    ! Returns the value string contains a key                    !
    ! If we have a string like                                   !
    ! 'foo bar key : val                                         !
    ! This routine will return val                               !
    ! NOTE: This is a CASE-SENSITIVE routine                     !
    ! In addition, anycolons (:) and equal signs (=) between     !
    ! the key and its value are ignored as are whitespaces.      !
    !------------------------------------------------------------!
    ! str :: the string to search                                !
    ! key :: the substring to search for within str              !
    ! whole :: return the whole line following the key           !
    !          default: False                                    !
    !------------------------------------------------------------!
    ! Known issues                                               !
    ! None though this function does not need to be called with  !
    ! trim in the arguments as they are trimmed in here.         !
    !============================================================!
    implicit none
    ! Arguments
    character(len=*) :: str     ! The string to search
    character(len=*) :: key     ! The substring to search for within the main string
    logical,optional :: whole   ! Return the whole line in the string (Default: False)

    ! Local variables
    logical :: l_whole
    integer :: pos        ! postion of key within the string
    logical :: keypresent ! is key present in the string?

    integer :: iostat

    l_whole = .false.
    if (present(whole)) l_whole = whole

    ! Check if the string is present to begin with and if so, obtain the position at which it starts.
    pos = index(trim(str), trim(key))
    keypresent = pos>0

    ! Trim out the white spaces in the key so that we can use the pos to advance to just after the str.
    pos = pos + len_trim(key)

    if (keypresent) then
       ! Advance past any separation characters
       do while (pos < len(str) .and. scan(str(pos:pos),' =:') > 0)
          pos = pos + 1
       end do
       if (l_whole) then
          read(str(pos:),'(A)',iostat=iostat) io_str_code
       else
          read(str(pos:),*,iostat=iostat) io_str_code
       end if
       if (iostat/=0) io_str_code = ''
    else
       io_str_code = ''
    end if

  end function io_str_code

  logical function io_file_present(unit, key)
    !============================================================!
    ! This is effectively a wrapper to io_str_present except     !
    ! it searches a file instead.                                !
    !------------------------------------------------------------!
    ! unit :: the file to search through                         !
    ! key :: the substring to search for within str              !
    !------------------------------------------------------------!
    ! Known issues                                               !
    ! None though this function does not need to be called with  !
    ! trim in the arguments as they are trimmed within           !
    ! io_str_present.                                            !
    !                                                            !
    ! Additionally, the file must be opened                      !
    !------------------------------------------------------------!
    ! Necessary Conditions                                       !
    ! File must be OPEN and FORMATTED                            !
    !============================================================!
    implicit none
    ! Arguments
    integer :: unit         ! file unit to search for substring
    character(len=*) :: key ! substring to search for within the file

    ! Local variables
    character(len=file_maxpath) :: line  ! current line in file
    integer :: have_eof                  ! Have we reached an End of File

    logical :: opened
    character(len=file_maxpath) :: filename
    integer :: iostat

    inquire(unit=unit, opened=opened, name=filename)
    if (.not. opened) then
       error stop 'io_file_present: File was not opened.'
    end if

    ! Rewind the file
    rewind(unit, iostat=iostat)
    if (iostat/=0) error stop 'io_file_present: Failed to rewind file'

    io_file_present = .false.
    do
       read(unit,'(A)',iostat=have_eof) line
       if (have_eof==eof) return
       if (io_file_present.eqv..true.) return
       io_file_present = io_str_present(line, key)
    end do

  end function io_file_present

  character(len=file_maxpath) function io_file_code(unit,key,whole)
    !============================================================!
    ! This is effectively a wrapper to io_str_present except     !
    ! it searches a file instead.                                !
    !------------------------------------------------------------!
    ! unit :: the file to search through                         !
    ! key :: the substring to search for within str              !
    ! whole :: return the whole line following the key           !
    !          default: False                                    !
    !------------------------------------------------------------!
    ! Known issues                                               !
    ! None though this function does not need to be called with  !
    ! trim in the arguments as they are trimmed within           !
    ! io_str_present.                                            !
    !                                                            !
    ! Additionally, the file must be opened                      !
    !------------------------------------------------------------!
    ! Necessary Conditions                                       !
    ! File must be OPEN and FORMATTED                            !
    !============================================================!

    implicit none
    ! Arguments
    integer :: unit         ! file unit to search for substring
    character(len=*) :: key ! substring to search for within the file
    logical,optional :: whole   ! Return the whole line in the string (Default: False)

    ! Local variables
    logical :: l_whole
    character(len=file_maxpath) :: line  ! current line in file
    integer :: have_eof                  ! Have we reached an End of File

    logical :: opened
    character(len=file_maxpath) :: filename
    integer :: iostat

    l_whole = .false.
    if(present(whole)) l_whole = whole

    inquire(unit=unit, opened=opened, name=filename)
    if (.not. opened) then
       error stop 'io_file_present: File was not opened.'
    end if

    ! Rewind the file
    rewind(unit, iostat=iostat)
    if (iostat/=0) error stop 'io_file_present: Failed to rewind file'

    io_file_code = ''
    do
       read(unit,'(A)',iostat=have_eof) line
       if(have_eof==eof) return
       io_file_code = io_str_code(line,key,whole=l_whole)
       if (trim(io_file_code)/='') return
    end do
  end function io_file_code
end module io
