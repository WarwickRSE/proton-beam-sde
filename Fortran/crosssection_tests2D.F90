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

    val = sampleAngleFromSection(xsec, 150.0_REAL64, randomVal(0.1_REAL64))
    PRINT*, val
    val = sampleAngleFromSection(xsec, 150.0_REAL64, randomVal(0.67_REAL64))
    PRINT*, val

#ifdef undef
c++
0.349714
0.395738
0.0563128
3.14159
0.04
3.14159
0.04
0.276185
0.0763119

Fortran now

  0.34979011233748014     
  0.39575874285697199     
   5.6362699557183399E-002
   3.1415926535897931     
   4.0000000000000001E-002
   3.1415926535897931     
   4.0000000000000001E-002
  0.27606350020288201     
   7.5981097636682760E-002
#endif

END PROGRAM