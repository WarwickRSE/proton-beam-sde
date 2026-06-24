
PROGRAM main
    
    USE protonEffects
    IMPLICIT NONE

    ! main exists just to force compilation

    CALL defineCrossSections([character(len=30) :: "carbon", "argon"])
END PROGRAM