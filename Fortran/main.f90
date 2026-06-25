
PROGRAM main
 
    USE MaterialFake
    USE protonEffects
    USE dataDefinedScattering
    IMPLICIT NONE

    TYPE(cp_material) :: testMaterial
    REAL(KIND=REAL64) :: tmp

    ! main exists just to force compilation

    ! SKIPPING hydrogen for now  TODO hydrogen
    ![character(len=30) :: "argon", "calcium", "carbon", "chlorine", "fluorine", "nitrogen", "oxygen", "phosphorus", "potassium", "sodium", "sulfur"]
    !argon, chlorine, potassium have one line too few in the Rutherford files TODO - diagnose or fix
    CALL defineCrossSections([character(len=30) :: "calcium", "carbon", "fluorine", "nitrogen", "oxygen", "phosphorus", "sodium", "sulfur"], 0.04_REAL64)

    ! Start with just a carbon material
    testMaterial%no_nucs = 1
    testMaterial%density = 1.0_REAL64
    testMaterial%mee = 1.0_REAL64
    ALLOCATE(testMaterial%nucs(1))
    testMaterial%nucs(1)%name = 'carbon'
    testMaterial%nucs(1)%num_dens = 1.0_REAL64
    testMaterial%nucs(1)%Z = 6
    testMaterial%nucs(1)%A = 12.011
    testMaterial%nucs(1)%massFraction = 1.0

    ! Store - don't want to look up the sections every time!
    testMaterial%nucs(1)%xsec_ind = getCrossSectionIndex(testMaterial%nucs(1)%name)

    tmp = bethe_bloch_loss(testMaterial, 1.0_REAL64)
    PRINT*, tmp

    tmp = rutherford_and_elastic_rate(testMaterial, 1.0_REAL64)
    PRINT*, tmp

END PROGRAM