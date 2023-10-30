program io_test_suite
  use io

  implicit none

  character(len=300) :: line
  character(len=file_maxpath), allocatable :: blockdata(:)
  integer :: i

  print *, 'io_strip_extension: ', io_strip_extension('file.txt')
  print *, 'io_lowercase: ', io_lowercase('THIS IS A STRing')

  open(unit=90,file='io_test.txt',status='UNKNOWN')

  call io_skip_header(90,'END HEADER',case_sensitive=.true.)
  read(90,'(A)') line
  print *, 'skip_header case-sensitive:   ', trim(line)

  call io_skip_header(90,'end header')
  read(90,'(A)') line
  print *, 'skip_header case-insensitive: ', trim(line)

  print *, ''
  print *, 'io_is_block_delim(case-sensitive)   ', io_is_block_delim('%BLOCK MY_BLOCK','%BLOCK','MY_BLOCK',case_sensitive=.true.) &
       , io_is_block_delim('%BLOCK MY_BLOCK','%block','MY_BLOCK',case_sensitive=.true.)
  print *, 'io_is_block_delim(case-insensitive) ', io_is_block_delim('%BLOCK MY_BLOCK','%BLOCK','MY_BLOCK',case_sensitive=.false.) &
       , io_is_block_delim('%BLOCK MY_BLOCK','%BLOCK','MY_BLOCK',case_sensitive=.false.)

  print *, ''
  call io_read_block(90,'MY_BLOCK',blockdata,case_sensitive=.false.)
  print *, 'io_read_block insensitive ', (trim(blockdata(i))//' ', i=1,size(blockdata))
  call io_read_block(90,'MY_BLOCK',blockdata,'%BLOCK', '%ENDBLOCK', case_sensitive=.true.)
  print *, 'io_read_block sensitive   ', (trim(blockdata(i))//' ', i=1,size(blockdata))

  print *, ''
  print *, 'io_str_present case-insensitive ', &
       io_str_present('This is a string containing MY_CODE','my_code',case_sensitive=.false.)
  print *, 'io_str_present case-sensitive   ', &
       io_str_present('This is a string containing MY_CODE','MY_CODE',case_sensitive=.true.)
  print *, 'io_str_code case-insensitive ', trim(io_str_code('MY_CODE : val','my_CODE', case_sensitive=.false.))
  print *, 'io_str_code case-sensitive   ', trim(io_str_code('MY_CODE : VAL','MY_CODE', case_sensitive=.true.))

  print *, 'File codes'
  line = io_file_code(90,'GRID_SIZE',whole=.true.)
  print *, 'Grid size', trim(io_file_code(90,'GRID_SIZE',whole=.true.))
  print *, 'no of particles ', trim(io_file_code(90,'NPARTICLES',whole=.true.))
end program io_test_suite
