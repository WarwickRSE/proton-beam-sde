PROGRAM main

    USE randomMod

    TYPE(KissRNGState) :: state
    TYPE(BoxMullerRNGState) :: b_state
    INTEGER :: i


    CALL random_warmup(state, 10)
    DO i = 1, 1000
       WRITE(10, *) random(state)
    END DO

    DO i = 1, 1000
       WRITE(20, *) rejection_sample(state, tophat)
    END DO

    b_state%k_state = state
    DO i = 1, 1000
       WRITE(30, *) random_box_muller(0.5_REAL64, b_state)
    END DO


CONTAINS

FUNCTION tophat(x)
    REAL(KIND=REAL64), INTENT(IN) :: x
    REAL(KIND=REAL64) :: tophat

    IF( x < 0.25 .OR. x > 0.75) THEN
        tophat = 0.5
    ELSE 
        tophat = 1.0
    ENDIF
END FUNCTION

END PROGRAM