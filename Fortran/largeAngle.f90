MODULE largeAngle
    USE iso_fortran_env, only: real64 ! REPLACE WITH kinds!
    IMPLICIT NONE

    CHARACTER(LEN=50), PARAMETER :: datadir = "../Splines/" ! Temporary

    TYPE crossSectionFromData
        LOGICAL :: ready = .FALSE. ! Debug/development 
        REAL(KIND=REAL64), DIMENSION(:), ALLOCATABLE :: energies, values
    END TYPE

    ! Associative arrays, index will be material number
    TYPE(crossSectionFromData), DIMENSION(:), ALLOCATABLE :: lA_crossSections

    CONTAINS

    SUBROUTINE defineCrossSections(names)
        ! Read data for all necessary materials
        CHARACTER(LEN=30), DIMENSION(:) :: names
        CHARACTER(LEN=30) :: file
        INTEGER :: i

        ALLOCATE(lA_crossSections(SIZE(names)))
        DO i = 1, SIZE(names)
            file = ADJUSTL(TRIM(names(i)))//"_ne_rate.txt"
            CALL fillCrossSections(file, lA_crossSections(i))
        END DO
    END SUBROUTINE

    SUBROUTINE fillCrossSections(file, crossSec)
        CHARACTER(LEN=30), INTENT(IN) :: file
        TYPE(crossSectionFromData), INTENT(INOUT) :: crossSec
        CHARACTER(LEN=50) :: fullpath
        REAL(KIND=REAL64), DIMENSION(:), ALLOCATABLE :: energies, values
        INTEGER :: unit, err, ct

        fullpath = TRIM(datadir)//ADJUSTL(TRIM(file))
        OPEN(newunit=unit, FILE=fullpath, ACTION="READ", IOSTAT=err)

        IF(err /= 0) THEN
            PRINT*, "Error opening File "//TRIM(fullpath)
            ERROR STOP
        END IF
        !Suggest starting each file with the count - its a lot easier
        READ(unit, *) ct
        ! Read the pre-prepared data files by material name/number ?
        ALLOCATE(energies(ct), values(ct))
        ! Header row
        READ(unit, *) energies
        READ(unit, *) values
        PRINT*, energies, values

        CALL MOVE_ALLOC(energies, crossSec%energies)
        CALL MOVE_ALLOC(values, crossSec%values)
        crossSec%ready = .TRUE.
        
        CLOSE(unit)

    END SUBROUTINE

END MODULE