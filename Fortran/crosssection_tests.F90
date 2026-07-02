PROGRAM crossSectionTests

    USE dataDefinedScattering
    IMPLICIT NONE


    CHARACTER(LEN=30) :: name = "carbon"
    REAL(KIND=REAL64) :: val
    TYPE(crossSection1D) :: xsec

    CALL defineCrossSections([name], 0.04_REAL64, 0.04_REAL64, "../Splines")
    !xsec = getNECrossSection(name)
    xsec = getRUCrossSection(name)

    val = evaluate(xsec, 100.0_REAL64)
    PRINT*, val
    val = evaluate(xsec, 73.0_REAL64)
    PRINT*, val
    val = evaluate(xsec, 5.3_REAL64)
    PRINT*, val

    val = evaluate(xsec, 1.0_REAL64)
    PRINT*, val

    val = evaluate(xsec, 160.0_REAL64)
    PRINT*, val

#ifdef undef
 RU c++

 0.287953
0.501884
51.9379
1466.2
0.134808

Fortran
  0.28795293012259132     
  0.50188429709679916     
   51.937851797007539     
   1466.1996235474701     
  0.13480769369545301 
CORRECTED Fortran
  0.28820004297554463     
  0.50198180187412156     
   52.031999119473987     
   1470.7285436663133     
  0.13458716353250083   


NE c++
0.227
0.2735
0.096396
0
0.22235

Fortran
  0.22699989999999998     
  0.27349983000000000     
   9.6395982999999991E-002
   0.0000000000000000     
  0.22234999999999999   
#endif

END PROGRAM