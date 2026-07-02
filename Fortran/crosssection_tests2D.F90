PROGRAM crossSectionTests

    USE dataDefinedScattering
    IMPLICIT NONE


    CHARACTER(LEN=30) :: name = "oxygen"
    REAL(KIND=REAL64) :: val
    TYPE(crossSection2D) :: xsec

    CALL defineCrossSections([name], 0.04_REAL64, 0.04_REAL64, "../Splines")
    xsec = getRU2DCrossSection(name)

    val = sampleAngleFromSection(xsec, 100.0_REAL64, 0.1_REAL64)
    PRINT*, val
    val = sampleAngleFromSection(xsec, 73.0_REAL64, 0.1_REAL64)
    PRINT*, val
    val = sampleAngleFromSection(xsec, 5.3_REAL64, 0.5_REAL64)
    PRINT*, val

    val = sampleAngleFromSection(xsec, 1.0_REAL64, 0.0_REAL64)
    PRINT*, val
    val = sampleAngleFromSection(xsec, 1.0_REAL64, 1.0_REAL64)
    PRINT*, val


    val = sampleAngleFromSection(xsec, 160.0_REAL64, 0.0_REAL64)
    PRINT*, val
    val = sampleAngleFromSection(xsec, 160.0_REAL64, 1.0_REAL64)
    PRINT*, val

    val = sampleAngleFromSection(xsec, 150.0_REAL64, 0.1_REAL64)
    PRINT*, val
    val = sampleAngleFromSection(xsec, 150.0_REAL64, 0.67_REAL64)
    PRINT*, val

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

   1.4499630036505826     
   1.5392489413257029     
   1.5316214947367850     
   1.5507131931784290     
   3.9497439654753559E-002
   1.5499317998036497     
   3.9935425037242712E-002
   1.3780458774630606     
  0.54078352971086474  
  

c++



Fortran
  0.32835793450464956     
  0.37487729503721190     
   5.6465170354506974E-002
   3.1415926535897931     
   4.0000000000000001E-002
   3.1415926535897931     
   4.0000000000000001E-002
  0.26027328163412861     
   5.9294522246662894E-002

#endif

END PROGRAM