PROGRAM crossSectionTests

    USE dataDefinedScattering
    IMPLICIT NONE


    CHARACTER(LEN=30) :: name = "carbon"
    REAL(KIND=REAL64) :: val
    TYPE(crossSection2D) :: xsec

    CALL defineCrossSections([name], 0.04_REAL64)
    xsec = getRU2DCrossSection(name)

    val = sampleAngleFromSection(xsec, 100.0_REAL64, randomVal(0.1_REAL64))
    PRINT*, val
    val = sampleAngleFromSection(xsec, 73.0_REAL64, randomVal(0.1_REAL64))
    PRINT*, val
    val = sampleAngleFromSection(xsec, 5.3_REAL64, randomVal(0.5_REAL64))
    PRINT*, val

    val = sampleAngleFromSection(xsec, 1.0_REAL64, randomVal(0.0_REAL64))
    PRINT*, val
    val = sampleAngleFromSection(xsec, 1.0_REAL64, randomVal(1.0_REAL64))
    PRINT*, val


    val = sampleAngleFromSection(xsec, 160.0_REAL64, randomVal(0.0_REAL64))
    PRINT*, val
    val = sampleAngleFromSection(xsec, 160.0_REAL64, randomVal(1.0_REAL64))
    PRINT*, val


END PROGRAM