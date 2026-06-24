MODULE dataDefinedScattering
    USE iso_fortran_env, only: real64 ! REPLACE WITH kinds!
    IMPLICIT NONE

    CHARACTER(LEN=50), PARAMETER :: datadir = "../Splines/" ! Temporary

    TYPE crossSection1D
        LOGICAL :: ready = .FALSE. ! Debug/development 
        REAL(KIND=REAL64), DIMENSION(:), ALLOCATABLE :: energies, values
    END TYPE

    TYPE cdfRow
        REAL(KIND=REAL64), DIMENSION(:), ALLOCATABLE :: angles, values
    END TYPE
    TYPE crossSection2D
        LOGICAL :: ready = .FALSE. ! Debug/development 
        REAL(KIND=REAL64), DIMENSION(:), ALLOCATABLE :: energies
        TYPE(cdfRow), DIMENSION(:), ALLOCATABLE :: cdf
    END TYPE

    TYPE crossSection3D
        LOGICAL :: ready = .FALSE. ! Debug/development 
        REAL(KIND=REAL64), DIMENSION(:), ALLOCATABLE :: energies, values
    END TYPE

    ! Associative arrays, index will be ATOMIC number
    TYPE(crossSection1D), DIMENSION(:), ALLOCATABLE :: NE_crossSections, RU_crossSections
    TYPE(crossSection2D), DIMENSION(:), ALLOCATABLE :: RU_angle_cdf

    INTERFACE fillCrossSections
      MODULE PROCEDURE fillCrossSections1D
      MODULE PROCEDURE fillCrossSections1DRuth
      MODULE PROCEDURE fillCrossSections2DRuth
    END INTERFACE

    CONTAINS

    !> \brief Populate the cross section data from files
    SUBROUTINE defineCrossSections(names, ru_cutoff)
        ! Read data for all necessary materials
        CHARACTER(LEN=30), DIMENSION(:), INTENT(IN) :: names
        REAL(KIND=REAL64) :: ru_cutoff
        CHARACTER(LEN=50) :: file
        INTEGER :: i

        ALLOCATE(NE_crossSections(SIZE(names)))
        DO i = 1, SIZE(names)
            file = ADJUSTL(TRIM(names(i)))//"_ne_rate.txt"
            CALL fillCrossSections(file, NE_crossSections(i))
        END DO
!        ALLOCATE(RU_crossSections(SIZE(names)))
!        DO i = 1, SIZE(names)
!            file = ADJUSTL(TRIM(names(i)))//"_el_ruth_cross_sec.txt"
!            CALL fillCrossSections(file, RU_crossSections(i), ru_cutoff)
!        END DO
        ALLOCATE(RU_angle_cdf(SIZE(names)))
        DO i = 1, SIZE(names)
            file = ADJUSTL(TRIM(names(i)))//"_el_ruth_cross_sec.txt"
            CALL fillCrossSections(file, RU_angle_cdf(i), ru_cutoff)
        END DO

    END SUBROUTINE


    !> \brief Helper - read a 1-D section
    SUBROUTINE fillCrossSections1D(file, crossSec)
        CHARACTER(LEN=50), INTENT(IN) :: file
        TYPE(crossSection1D), INTENT(INOUT) :: crossSec
        CHARACTER(LEN=80) :: fullpath
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

    !> \brief Helper - read a 1-D Rutherford section
    SUBROUTINE fillCrossSections1DRuth(file, crossSec, cutoff)
        ! Stash the needful when reading the full 2D and dont re-read the files for this
        CHARACTER(LEN=50), INTENT(IN) :: file
        TYPE(crossSection1D), INTENT(INOUT) :: crossSec
        CHARACTER(LEN=80) :: fullpath
        REAL(KIND=REAL64) :: cutoff
    END SUBROUTINE

    !> \brief Helper - read the 2-D Rutherford section
    SUBROUTINE fillCrossSections2DRuth(file, crossSec, cutoff)
        CHARACTER(LEN=50), INTENT(IN) :: file
        TYPE(crossSection2D), INTENT(INOUT) :: crossSec
        CHARACTER(LEN=80) :: fullpath
        REAL(KIND=REAL64) :: cutoff
        REAL(KIND=REAL64) :: interp
        REAL(KIND=REAL64), DIMENSION(:), ALLOCATABLE :: energies, tmp
        INTEGER :: i, unit, err, ct, a_ct, f_ct

        fullpath = TRIM(datadir)//ADJUSTL(TRIM(file))
        OPEN(newunit=unit, FILE=fullpath, ACTION="READ", IOSTAT=err)

        IF(err /= 0) THEN
            PRINT*, "Error opening File "//TRIM(fullpath)
            ERROR STOP
        END IF
        !Suggest starting each file with the (line length) count - its a lot easier
        READ(unit, *) ct, a_ct
        ! Read the pre-prepared data files by material name/number ?
        ALLOCATE(energies(ct), crossSec%cdf(ct), tmp(a_ct))

        ! Header row
        READ(unit, *) energies

        DO i = 1, ct
            PRINT*, i
            ! For each row:
            ! Read the angles
            READ(unit, *, IOSTAT=err) tmp
            IF(err /= 0) ERROR STOP "Missing Energy Value in File "//TRIM(fullpath)
            ! Find cutoff index
            ! TODO - check for off-by-one
            f_ct = MINLOC(tmp, DIM=1, MASK = (tmp < cutoff))
            ! Move the angles array
            crossSec%cdf(i)%angles = tmp(1:f_ct)
            crossSec%cdf(i)%angles(f_ct) = cutoff ! Force last angle to cutoff
            ! Allocate and read the cdf row including one value past the cutoff
            ALLOCATE(crossSec%cdf(i)%values(f_ct))
            READ(unit, *) crossSec%cdf(i)%values

            ! Correct the cdf value at the last angle (currently just past the cutoff, interpolate back)
            interp = (cutoff - tmp(f_ct - 1)) / (tmp(f_ct) - tmp(f_ct -1))
            crossSec%cdf(i)%values(f_ct) = crossSec%cdf(i)%values(f_ct) * interp + (1.0_REAL64 - interp) * crossSec%cdf(i)%values(f_ct - 1)
        END DO
        PRINT*, energies
        DO i = 1, ct
            print*, minVAL(crossSec%cdf(i)%angles), maxval(crossSec%cdf(i)%angles)
        end do

        CALL MOVE_ALLOC(energies, crossSec%energies)
        crossSec%ready = .TRUE.
        
        CLOSE(unit)

    END SUBROUTINE



END MODULE