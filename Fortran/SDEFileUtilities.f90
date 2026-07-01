MODULE SDEFileUtilities

  USE iso_fortran_env, ONLY : REAL64
  IMPLICIT NONE

  INTEGER, PARAMETER, PRIVATE :: maxbins = 11000 ! Prevent infinite loop if delim not found
  INTEGER, PARAMETER, PRIVATE :: max_buf = 2**20 ! Needed to read the entire cross section

 CONTAINS

  ! Read a single line from 'unit' and convert to an array
  ! Uses a persistent buffer
  ! pass dealloc = .TRUE. to deallocate this to free memory
  SUBROUTINE readLineOfReals(unit, row, err, delim, dealloc)
    INTEGER, INTENT(IN) :: unit
    REAL(KIND=REAL64), ALLOCATABLE, DIMENSION(:), INTENT(INOUT) :: row
    INTEGER, INTENT(INOUT) :: err
    LOGICAL, INTENT(IN), OPTIONAL :: dealloc
    CHARACTER(LEN=1), VALUE, OPTIONAL :: delim
    CHARACTER(LEN=:), ALLOCATABLE, SAVE :: buffer
    CHARACTER(LEN=30) :: fmt
    INTEGER :: sz
   
    IF(PRESENT(dealloc)) THEN
      IF(dealloc .AND. ALLOCATED(buffer)) DEALLOCATE(buffer)
      RETURN
    END IF

    WRITE(fmt, *) max_buf
    fmt = TRIM("(A"//ADJUSTL(fmt))//")"
    IF(.NOT. ALLOCATED(buffer)) ALLOCATE(CHARACTER(LEN=max_buf)::buffer)
 
    ! Read a whole line
    READ(unit, fmt, SIZE=sz, IOSTAT=err, ADVANCE='NO') buffer
    IF(err == -1) RETURN ! END OF FILE, return to caller!
    IF(err == -2 .AND. sz >= LEN(buffer)) ERROR STOP "Line buffer size "//fmt//"too small for file. Increase max_buf and try again"
    row = lineToArray(buffer)

 END SUBROUTINE

  FUNCTION lineToArray(line, delim_in) RESULT(row)
    CHARACTER(LEN=*), INTENT(IN) :: line
    REAL(KIND=REAL64), ALLOCATABLE, DIMENSION(:) :: row
    CHARACTER(LEN=1), VALUE, OPTIONAL :: delim_in
    CHARACTER(LEN=1) :: delim
    INTEGER :: err, j, st, st_old, ind

    IF(PRESENT(delim_in)) THEN
        delim = delim_in
    ELSE
        delim = " "
    END IF

    st = 1
    st_old = 1
    ind = 1
    ! Find out how many substrings there are
    DO j = 1, maxbins
      IF(st >= LEN(TRIM(line))) EXIT
      ind = SCAN(line(st:), delim)
      st = st + ind
    END DO
    IF(j >= maxbins) ERROR STOP "More values in this line than we can handle"
    ALLOCATE(row(j-1))
    st = 1
    st_old = 1
    ind = 1
    READ(line, *, IOSTAT=err) row
    IF(err /= 0) ERROR STOP "Unknown error converting line: "//line
  END FUNCTION

END MODULE