PROGRAM crossSectionTests

    USE dataDefinedScattering
    IMPLICIT NONE


    CHARACTER(LEN=30) :: name = "carbon"
    TYPE(crossSection3D) :: xsec
    TYPE(xsecSample) :: val

    CALL defineCrossSections([name], 0.04_REAL64, 0.04_REAL64, "../Splines")
    xsec = get3DCrossSection(name)

    !val = sampleAngleFromSection(xsec, 100.0_REAL64, randomVal(0.1_REAL64))
    val = sampleAtEnergy(xsec%cdf(10), xsec%exit_energy(10), xsec%rvalue(10), 0.1_REAL64)
    PRINT*, val%e, val%r

    val = sampleAngleFromSection(xsec, 73.0_REAL64, randomVal(0.1_REAL64))
     PRINT*, val%e, val%r

    val = sampleAngleFromSection(xsec, 5.3_REAL64, randomVal(0.5_REAL64))
    PRINT*, val%e, val%r

    val = sampleAngleFromSection(xsec, 1.0_REAL64, randomVal(0.0_REAL64))
    PRINT*, val%e, val%r

    val = sampleAngleFromSection(xsec, 1.0_REAL64, randomVal(1.0_REAL64))
    PRINT*, val%e, val%r

    val = sampleAngleFromSection(xsec, 160.0_REAL64, randomVal(0.0_REAL64))
    PRINT*, val%e, val%r

    val = sampleAngleFromSection(xsec, 160.0_REAL64, randomVal(1.0_REAL64))
    PRINT*, val%e, val%r

    val = sampleAngleFromSection(xsec, 150.0_REAL64, randomVal(0.1_REAL64))
    PRINT*, val%e, val%r

    val = sampleAngleFromSection(xsec, 150.0_REAL64, randomVal(0.67_REAL64))
    PRINT*, val%e, val%r

#ifdef undef

carbon
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

  0.34971414586013017     
  0.39573750571348182     
   5.6312822623283409E-002
   3.1415926535897931     
   4.0000000000000001E-002
   3.1415926535897931     
   4.0000000000000001E-002
  0.27618468284961678     
   7.6311936705110564E-002

hydrogen c++

NOTE: lab_ang_cutoff is about 1.55

1.44996
1.53925
1.53162
1.55071
0.04
1.54993
0.04
1.37805
0.540784

Fortran

   1.5502732698094992     
   1.5504138807854488     
   1.5507707812916556     
   1.5507936655375631     
   4.0000000000000001E-002
   1.5500154014037679     
   4.0000000000000001E-002
   1.5500154014037679     
  0.61981994317006173   

#endif

END PROGRAM