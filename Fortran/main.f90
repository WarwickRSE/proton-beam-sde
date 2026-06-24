
PROGRAM main
    
    USE protonEffects
    IMPLICIT NONE

    ! main exists just to force compilation

    ! SKIPPING hydrogen for now  TODO hydrogen
    ![character(len=30) :: "argon", "calcium", "carbon", "chlorine", "fluorine", "nitrogen", "oxygen", "phosphorus", "potassium", "sodium", "sulfur"]
    !argon, chlorine, potassium have one line too few in the Rutherford files TODO - diagnose or fix
    CALL defineCrossSections([character(len=30) :: "calcium", "carbon", "fluorine", "nitrogen", "oxygen", "phosphorus", "sodium", "sulfur"], 0.04_REAL64)

END PROGRAM