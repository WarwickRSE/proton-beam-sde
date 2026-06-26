MODULE dataDefinedScattering
    USE iso_fortran_env, only: real64 ! REPLACE WITH kinds!
    IMPLICIT NONE

    CHARACTER(LEN=50), PARAMETER :: datadir = "../Splines/" ! Temporary
    REAL(KIND=REAL64), PARAMETER :: pi = 3.14159265_REAL64

    TYPE randomGen
    ! This exists just to flag where random numbers are needed here
    END TYPE
    TYPE randomVal
      REAL(KIND=REAL64) :: v
    ! DItto
    END TYPE

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
    CHARACTER(LEN=30), DIMENSION(:), ALLOCATABLE :: atom_names
    TYPE(crossSection1D), DIMENSION(:), ALLOCATABLE :: NE_crossSections, RU_crossSections
    TYPE(crossSection2D), DIMENSION(:), ALLOCATABLE :: RU_angle_cdf

    INTERFACE fillCrossSections
      MODULE PROCEDURE fillCrossSections1D
      MODULE PROCEDURE fillCrossSectionsRuth
      MODULE PROCEDURE fillCrossSectionsHydrogen
    END INTERFACE

    INTERFACE evaluate
      MODULE PROCEDURE evaluate1DCrossSection
    END INTERFACE

    CONTAINS

    !> \brief Populate the cross section data from files
    SUBROUTINE defineCrossSections(names, ru_cutoff, bs_cutoff)
        ! Read data for all necessary materials
        CHARACTER(LEN=30), DIMENSION(:), INTENT(IN) :: names
        REAL(KIND=REAL64), INTENT(IN) :: ru_cutoff, bs_cutoff
        CHARACTER(LEN=50) :: file
        INTEGER :: i

        atom_names = names

        ALLOCATE(NE_crossSections(SIZE(names)))
        DO i = 1, SIZE(names)
            IF(TRIM(names(i)) /= 'hydrogen' ) THEN
              file = ADJUSTL(TRIM(names(i)))//"_ne_rate.txt"
              CALL fillCrossSections(file, NE_crossSections(i))
            END IF
        END DO
        ALLOCATE(RU_crossSections(SIZE(names)))
        ALLOCATE(RU_angle_cdf(SIZE(names)))
        DO i = 1, SIZE(names)
            file = ADJUSTL(TRIM(names(i)))//"_el_ruth_cross_sec.txt"
            IF(TRIM(names(i)) /= 'hydrogen' ) THEN
              CALL fillCrossSections(file, RU_angle_cdf(i), RU_crossSections(i), ru_cutoff)
            ELSE
              CALL fillCrossSections(file, RU_angle_cdf(i), RU_crossSections(i), ru_cutoff, bs_cutoff)
            END IF
        END DO

    END SUBROUTINE

    FUNCTION getCrossSectionIndex(name) RESULT(ind)
        CHARACTER(LEN=30) :: name
        INTEGER :: ind

        ind = MINLOC(atom_names, DIM=1, MASK=(atom_names == name))
    END FUNCTION

    FUNCTION getNECrossSection(name) RESULT(X)
        CHARACTER(LEN=30) :: name
        TYPE(crossSection1D) :: X
        INTEGER :: ind

        ind = MINLOC(atom_names, DIM=1, MASK=(atom_names == name))

        X = NE_crossSections(ind)
    END FUNCTION
    FUNCTION getRUCrossSection(name) RESULT(X)
        CHARACTER(LEN=30) :: name
        TYPE(crossSection1D) :: X
        INTEGER :: ind

        ind = MINLOC(atom_names, DIM=1, MASK=(atom_names == name))

        X = RU_crossSections(ind)
    END FUNCTION
    FUNCTION getRU2DCrossSection(name) RESULT(X)
        CHARACTER(LEN=30) :: name
        TYPE(crossSection2D) :: X
        INTEGER :: ind

        ind = MINLOC(atom_names, DIM=1, MASK=(atom_names == name))

        X = RU_angle_cdf(ind)
    END FUNCTION

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

        CALL MOVE_ALLOC(energies, crossSec%energies)
        CALL MOVE_ALLOC(values, crossSec%values)
        crossSec%ready = .TRUE.
        
        CLOSE(unit)

    END SUBROUTINE

    !> \brief Helper - read the Rutherford sections - both the 2D and the 1D in a single read
    SUBROUTINE fillCrossSectionsRuth(file, crossSec, crossSec1D, cutoff)
        CHARACTER(LEN=50), INTENT(IN) :: file
        TYPE(crossSection2D), INTENT(INOUT) :: crossSec
        TYPE(crossSection1D), INTENT(INOUT) :: crossSec1D
        CHARACTER(LEN=80) :: fullpath
        REAL(KIND=REAL64), INTENT(IN) :: cutoff
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

        ! Allocation
        crossSec1D%energies = energies
        ALLOCATE(crossSec1D%values(SIZE(energies)))
        CALL MOVE_ALLOC(energies, crossSec%energies)

        DO i = 1, ct
            ! For each row:
            ! Read the angles
            READ(unit, *, IOSTAT=err) tmp
            IF(err /= 0) ERROR STOP "Missing Energy Value in File "//TRIM(fullpath)
            ! Find cutoff index
            ! TODO - check for off-by-one
            f_ct = MINLOC(tmp, DIM=1, MASK=(tmp >= cutoff)) + 1
            ! Move the angles array
            crossSec%cdf(i)%angles = tmp(1:f_ct)
            crossSec%cdf(i)%angles(f_ct) = cutoff ! Force last angle to cutoff
            ! Allocate and read the cdf row including one value past the cutoff
            ALLOCATE(crossSec%cdf(i)%values(f_ct))
            READ(unit, *) crossSec%cdf(i)%values

            ! Correct the cdf value at the last angle (currently just past the cutoff, interpolate back)
            ! TODO double check interpolation
            interp = (cutoff - tmp(f_ct - 1)) / (tmp(f_ct) - tmp(f_ct -1))
            crossSec%cdf(i)%values(f_ct) = crossSec%cdf(i)%values(f_ct) * interp + (1.0_REAL64 - interp) * crossSec%cdf(i)%values(f_ct-1)

            ! Making a copy for the 1D X-section
            crossSec1D%values(i) = crossSec%cdf(i)%values(f_ct)
            ! Re-normalise CDF so that last value is 1
            crossSec%cdf(i)%values = crossSec%cdf(i)%values / crossSec%cdf(i)%values(f_ct) 
        END DO
        
        !DO i = 1, ct
        !    print*, minVAL(crossSec%cdf(i)%angles), maxval(crossSec%cdf(i)%angles)
        !end do

        crossSec%ready = .TRUE.
        crossSec1D%ready = .TRUE.

        CLOSE(unit)

    END SUBROUTINE

    PURE FUNCTION hydrogen_cm_to_lab(ang, En) RESULT(out)
        REAL(KIND=REAL64), VALUE :: ang, En
        REAL(KIND=REAL64) :: out
        REAL(KIND=REAL64), PARAMETER :: mp = 938.346
        REAL(KIND=REAL64) :: p, u, g, e, v_ratio
    
        ang = pi - ang
        p = SQRT(en * (En + 2 * mp))
        u = p / (En + 2 * mp)
        g = 1 / SQRT(1 - u * u)
        e = En + mp
        v_ratio = u * (e - u * p) / (p - u * e)
        if (ABS(g * (cos(ang) + v_ratio)) == 0) THEN
            out = PI / 2
        else
           out = atan(sin(ang) / (g * (cos(ang) + v_ratio)))
        END IF
        if (out < 0) THEN
        out = out +  PI
        END IF
    END FUNCTION

    !> \brief Helper - read the Rutherford sections - both the 2D and the 1D in a single read - but for Hydrogen
    SUBROUTINE fillCrossSectionsHydrogen(file, crossSec, crossSec1D, cutoff, bs_cutoff)
        CHARACTER(LEN=50), INTENT(IN) :: file
        TYPE(crossSection2D), INTENT(INOUT) :: crossSec
        TYPE(crossSection1D), INTENT(INOUT) :: crossSec1D
        CHARACTER(LEN=80) :: fullpath
        REAL(KIND=REAL64), INTENT(IN) :: cutoff, bs_cutoff
        REAL(KIND=REAL64) :: interp, lab_ang_cutoff
        REAL(KIND=REAL64), DIMENSION(:), ALLOCATABLE :: energies, tmp
        INTEGER :: i, unit, err, ct, a_ct, f_ct, b_ct

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

        ! Allocation
        crossSec1D%energies = energies
        ALLOCATE(crossSec1D%values(SIZE(energies)))
        CALL MOVE_ALLOC(energies, crossSec%energies)

        DO i = 1, ct
            ! For each row:
            ! Transform the cutoff
            lab_ang_cutoff = hydrogen_cm_to_lab(bs_cutoff, crossSec%energies(i))
            ! Read the angles
            READ(unit, *, IOSTAT=err) tmp
            IF(err /= 0) ERROR STOP "Missing Energy Value in File "//TRIM(fullpath)
            ! Find cutoff index
            ! TODO - check for off-by-one
            f_ct = MINLOC(tmp, DIM=1, MASK=(tmp > cutoff))
            b_ct = MINLOC(tmp, DIM=1, MASK=(tmp > lab_ang_cutoff))
            IF(b_ct > f_ct) ERROR STOP "I don't think the cutoffs can be this way round"
            ! Move the angles array
            crossSec%cdf(i)%angles = tmp(b_ct:f_ct)
   !         PRINT*, b_ct, f_ct, f_ct-b_ct, SIZE(tmp(b_ct:f_ct))
            crossSec%cdf(i)%angles(f_ct-b_ct+1) = cutoff ! Force last angle to cutoff
            ! Allocate and read the cdf row including one value past the cutoff
            !ALLOCATE(crossSec%cdf(i)%values(f_ct-b_ct))
            READ(unit, *) tmp(1:f_ct)
            crossSec%cdf(i)%values = tmp(b_ct:f_ct)

            ! f_ct is the size from here on
            f_ct = f_ct - b_ct + 1
            ! Correct the cdf value at the last angle (currently just past the cutoff, interpolate back)
            interp = (cutoff - crossSec%cdf(i)%angles(f_ct - 1)) / (crossSec%cdf(i)%angles(f_ct) - crossSec%cdf(i)%angles(f_ct -1))
            crossSec%cdf(i)%values(f_ct) = crossSec%cdf(i)%values(f_ct) * interp + (1.0_REAL64 - interp) * crossSec%cdf(i)%values(f_ct-1)

            ! Making a copy for the 1D X-section
            crossSec1D%values(i) = crossSec%cdf(i)%values(f_ct) - crossSec%cdf(i)%values(1)
            ! Re-normalise CDF so that last value is 1
            crossSec%cdf(i)%values = crossSec%cdf(i)%values / crossSec%cdf(i)%values(f_ct) 
        END DO
        
        !DO i = 1, ct
        !    print*, minVAL(crossSec%cdf(i)%angles), maxval(crossSec%cdf(i)%angles)
        !end do

        crossSec%ready = .TRUE.
        crossSec1D%ready = .TRUE.

        CLOSE(unit)

    END SUBROUTINE



    PURE FUNCTION evaluate1DCrossSection(crossSection, energy) RESULT(val)
        TYPE(crossSection1D), INTENT(IN) :: crossSection
        REAL(KIND=REAL64), INTENT(IN) :: energy
        REAL(KIND=REAL64) :: val
        REAL(KIND=REAL64), PARAMETER :: tol = 1.0d-7
        INTEGER :: ct, ind
    
        IF(crossSection%ready) THEN
            ct = SIZE(crossSection%energies)
            IF(energy < crossSection%energies(1)) THEN
                val = crossSection%values(1) ! Best guess - the lowest value
            ELSE IF(energy > crossSection%energies(ct)) THEN
                val = crossSection%values(ct)
            ELSE
                ! Location of first energy which exceeds target
                ind = MINLOC(crossSection%energies, DIM=1, MASK=(crossSection%energies > energy))
                ! TODO can ind == 1 here?
                ! Interpolate if energies are not too close together
                ! TODO - better to soften the division ?
                IF(crossSection%energies(ind) - crossSection%energies(ind-1) > tol) THEN
                  val = ((crossSection%energies(ind) - energy) * crossSection%values(ind-1) + &
                      (energy - crossSection%energies(ind-1)) * crossSection%values(ind))/ &
                      (crossSection%energies(ind) - crossSection%energies(ind-1))
                ELSE
                  val = crossSection%values(ind)
                END IF
            END IF
        ELSE
            ERROR STOP "Trying to evaluate an unpopulated or empty cross section"
        ENDIF
    
    END FUNCTION

    PURE FUNCTION sampleAtEnergy(cdf, u) RESULT(val)
        TYPE(cdfRow), INTENT(IN) :: cdf
        REAL(KIND=REAL64), INTENT(IN) :: u
        REAL(KIND=REAL64) :: val, diff
        REAL(KIND=REAL64), PARAMETER :: tol = 1.0d-7
        INTEGER :: ct, ind
 
        ct = SIZE(cdf%angles)
        IF(u <= cdf%values(1)) THEN
            val = cdf%angles(1)
        ELSE IF(u >= cdf%values(ct)) THEN
            val = cdf%angles(ct)
        ELSE
            ! Location of first value which exceeds target
            ind = MINLOC(cdf%values, DIM=1, MASK=(cdf%values > u))
            ! Interpolate if values are not too close together
            ! NOTE: interpolate the angle based on the cdf spacing
            ! TODO - better to soften the division ?
            IF(cdf%values(ind) - cdf%values(ind-1) > tol) THEN
              diff = (u - cdf%values(ind-1)) / &
                (cdf%values(ind) - cdf%values(ind-1))
              val = cdf%angles(ind) * diff + &
                cdf%angles(ind-1) * (1.0_REAL64 - diff)
            ELSE
              val = cdf%angles(ind)
            END IF
        END IF
    END FUNCTION

    PURE FUNCTION sampleAngleFromSection(crossSection, energy, val) RESULT(angle)
        TYPE(crossSection2D), INTENT(IN) :: crossSection
        REAL(KIND=REAL64), INTENT(IN) :: energy
        TYPE(randomVal), VALUE :: val
        REAL(KIND=REAL64), PARAMETER :: tol = 1.0d-7
        REAL(KIND=REAL64) :: angle, tmp_angle, diff
        INTEGER :: sz, ind

        sz = SIZE(crossSection%energies)
        IF(energy <= crossSection%energies(1)) THEN
            ! Use lowest energy strand
            angle = sampleAtEnergy(crossSection%cdf(1), val%v)
        ELSE IF(energy >= crossSection%energies(sz)) THEN
            ! Highest energy strand
            angle = sampleAtEnergy(crossSection%cdf(sz), val%v)
        ELSE
            ind = MINLOC(crossSection%energies, DIM=1, MASK=(crossSection%energies > energy))
            ! Do angle at ind
            angle = sampleAtEnergy(crossSection%cdf(ind), val%v)
            IF((crossSection%energies(ind) - crossSection%energies(ind-1)) > tol) THEN
                ! If needed do one bin lower and interpolate
                tmp_angle = sampleAtEnergy(crossSection%cdf(ind-1), val%v)
                diff = (energy - crossSection%energies(ind-1)) / &
                  (crossSection%energies(ind) - crossSection%energies(ind-1))
                angle = angle * diff + tmp_angle * (1.0_REAL64 - diff)
            END IF
        END IF

    END FUNCTION

END MODULE