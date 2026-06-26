MODULE MaterialFake

  USE iso_fortran_env, only: real64 ! REPLACE WITH kinds!

  TYPE cp_nuclide
    CHARACTER(LEN=50) :: name=' '       ! Nuclide Name
    REAL :: num_dens=0.0                ! Number Density
    INTEGER :: Z    ! Atomic number
    REAL(KIND=REAL64) :: A    ! Atomic weight
    REAL(KIND=REAL64) :: massFraction !

    INTEGER :: xsec_ind
  END TYPE
  TYPE cp_material
    INTEGER :: no_nucs=0                ! Number of Nuclides
    TYPE(cp_nuclide), ALLOCATABLE :: nucs(:)    ! Nuclide data
    REAL(KIND=REAL64) :: density ! Overall density in g/cm^3 expected
    REAL(KIND=REAL64) :: mee ! Mean Excitation Energy
  END TYPE


  CONTAINS

  FUNCTION createCarbon() RESULT(testMaterial)
    TYPE(cp_material) :: testMaterial
    testMaterial%no_nucs = 1
    testMaterial%density = 1.0_REAL64
    testMaterial%mee = 1.0_REAL64
    ALLOCATE(testMaterial%nucs(1))
    testMaterial%nucs(1)%name = 'carbon'
    testMaterial%nucs(1)%num_dens = 1.0_REAL64
    testMaterial%nucs(1)%Z = 6
    testMaterial%nucs(1)%A = 12.011
    testMaterial%nucs(1)%massFraction = 1.0

  END FUNCTION

  FUNCTION createWater() RESULT(testMaterial)
    TYPE(cp_material) :: testMaterial
    testMaterial%no_nucs = 2
    testMaterial%density = 1.0_REAL64
    testMaterial%mee = 1.0_REAL64

    ALLOCATE(testMaterial%nucs(2))
    testMaterial%nucs(1)%name = 'hydrogen'
    testMaterial%nucs(1)%num_dens = 0.33
    testMaterial%nucs(1)%Z = 1
    testMaterial%nucs(1)%A = 1.008
    testMaterial%nucs(1)%massFraction = 0.33

    testMaterial%nucs(2)%name = 'oxygen'
    testMaterial%nucs(2)%num_dens = 0.67
    testMaterial%nucs(2)%Z = 8
    testMaterial%nucs(2)%A = 15.999
    testMaterial%nucs(2)%massFraction = 0.67

  END FUNCTION


END MODULE