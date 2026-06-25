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

END MODULE