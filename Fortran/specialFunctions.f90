MODULE specialFunctions

  USE iso_fortran_env, only: real64 ! REPLACE WITH kinds!

  IMPLICIT NONE

  CONTAINS
  
  !> \brief Compute the logarithm of the factorial of n
  !> \param n The integer for which to compute the logarithm of the factorial
  !> \return The logarithm of the factorial of n
  PURE FUNCTION log_factorial(n) RESULT(log_fact)
    INTEGER, INTENT(IN) :: n
    REAL(KIND=REAL64) :: log_fact
    INTEGER :: i

    IF (n < 0) THEN
      log_fact = -1.0_REAL64
    ELSE IF (n == 0 .OR. n == 1) THEN
      log_fact = 0.0_REAL64
    ELSE
      log_fact = SUM(LOG(REAL([(i, i=2,n)], KIND=REAL64)))
    END IF
  END FUNCTION

  !> \brief Compute the logarithm of the Pochhammer symbol (a)_x
  !> \param a 
  !> \param x 
  !> \return The logarithm of the Pochhammer symbol (a)_x
  PURE FUNCTION log_pochhammer(a, x) RESULT(log_poch)
    REAL(KIND=REAL64), INTENT(IN) :: a, x
    REAL(KIND=REAL64) :: log_poch

    log_poch = LOG_GAMMA(a + x) - LOG_GAMMA(a)
  END FUNCTION

END MODULE