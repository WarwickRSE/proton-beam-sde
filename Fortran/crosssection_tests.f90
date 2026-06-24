PROGRAM crossSectionTests

    USE dataDefinedScattering
    IMPLICIT NONE


    CHARACTER(LEN=30) :: name = "carbon"
    REAL(KIND=REAL64) :: a = 12.011_REAL64, val
    INTEGER :: z = 6
    TYPE(crossSection1D) :: xsec

    CALL defineCrossSections([name], 0.04_REAL64)
    xsec = getNECrossSection(name)

    val = evaluate(xsec, 100.0_REAL64)
    PRINT*, val
    val = evaluate(xsec, 73.0_REAL64)
    PRINT*, val
    val = evaluate(xsec, 5.3_REAL64)
    PRINT*, val

    val = evaluate(xsec, 4.0_REAL64)
    PRINT*, val

    val = evaluate(xsec, 160.0_REAL64)
    PRINT*, val

END PROGRAM