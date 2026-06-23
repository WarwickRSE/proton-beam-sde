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
        INTEGER :: i

        ALLOCATE(lA_crossSections(SIZE(names)))
        DO i = 1, SIZE(names)
            CALL fillCrossSections(names(i), i)
        END DO
    END SUBROUTINE

    SUBROUTINE fillCrossSections(name, mat_num)
        CHARACTER(LEN=30), INTENT(IN) :: name
        INTEGER, INTENT(IN) :: mat_num
        CHARACTER(LEN=50) :: path, fullpath
        REAL(KIND=REAL64), DIMENSION(:), ALLOCATABLE :: energies, values
        INTEGER :: unit, err, ct

        path = "carbon_ne_rate.txt"
        fullpath = TRIM(datadir)//ADJUSTL(TRIM(path))
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

        CALL MOVE_ALLOC(energies, lA_crossSections(mat_num)%energies)
        CALL MOVE_ALLOC(values, lA_crossSections(mat_num)%values)
        
        CLOSE(unit)

    END SUBROUTINE

END MODULE