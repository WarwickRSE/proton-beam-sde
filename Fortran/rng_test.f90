PROGRAM main

    USE randomMod
!    USE ziggurat

    TYPE(KissRNGState) :: state
    TYPE(BoxMullerRNGState) :: b_state
    !TYPE(zigguratState) :: z_state

    INTEGER :: i
    REAL(KIND=REAL64) :: res


    CALL random_warmup(state, 10)
    DO i = 1, 1000
       WRITE(10, *) random(state)
    END DO

    DO i = 1, 10000
       WRITE(20, *) rejection_sample_beta(state, beta_pdf, 50)
    END DO

    DO i = 1, 1000
        WRITE(30, *) beta_pdf(REAL(i, KIND=REAL64)*0.001, 5)
    END DO
    !DO i = 1, 1000
    !   WRITE(20, *) rejection_sample(state, tophat)
    !END DO

    !b_state%k_state = state
    !DO i = 1, 1000
    !   WRITE(30, *) random_box_muller(0.5_REAL64, b_state)
    !END DO

 !   z_state%k_state = state
 !   CALL z_init(z_state)
 !   DO i = 1, 1
 !       res = z_sample(z_state)
 !   END DO

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