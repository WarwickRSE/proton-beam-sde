PROGRAM crossSectionTests

    USE dataDefinedScattering
    IMPLICIT NONE


    CHARACTER(LEN=30) :: name =  "carbon"
    TYPE(crossSection3D) :: xsec
    TYPE(xsecSample) :: val

    CALL defineCrossSections([name], 0.04_REAL64, 0.04_REAL64, "../Splines")
    xsec = get3DCrossSection(name)

    !val = sampleAngleFromSection(xsec, 100.0_REAL64, randomVal(0.1_REAL64))
    !val = sampleAtEnergy(xsec%cdf(10), xsec%exit_energy(10), xsec%rvalue(10), 0.1_REAL64)
    val = sampleAngleFromSection(xsec, 100.0_REAL64, randomVal(0.1_REAL64))
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

Fortran carbon
   1.4041869910565683       0.17159733479902800    
   1.4140471563248960       0.24688291745864485     
  0.20756731976943665        0.0000000000000000     
   0.0000000000000000        0.0000000000000000     
  0.10000000000000001        0.0000000000000000     
   0.0000000000000000        4.8160000000000001E-002
   123.15660000000000        1.0000000000000000  
   1.2086428171361061       0.11536204363870320     
   38.104322992441880       0.92511071143752877   

c++
1.40419  0.171597
1.41405  0.246883
0.207567  0
0  0
0.1  0
0  0.04816
123.157  1
1.20864  0.115362
38.1043  0.925111

oxygen Fortran
   1.7773278409592070       0.15590412030261663     
   1.6220764026698675       0.21446108851036197     
   0.0000000000000000        0.0000000000000000     
   0.0000000000000000        0.0000000000000000     
  0.10000000000000001        0.0000000000000000     
   0.0000000000000000        4.4729999999999999E-002
   126.52760000000001        1.0000000000000000  
   1.6867649196392400       0.10683722145087210     
   27.841716680174056       0.73990071622759379 

c++

1.77733 0.155904
1.62208 0.214461
0 0
0 0
0.1 0
0 0.5573
126.528 1
1.68676 0.106837
27.8417 0.739901

#endif

END PROGRAM