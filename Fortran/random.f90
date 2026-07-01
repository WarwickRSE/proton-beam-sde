MODULE randomMod

    !WARNING - these are test-suitable RNGS but only of moderate quality. DO NOT rely on these for high quality randomness

    USE iso_fortran_env

    IMPLICIT NONE
    TYPE KissRNGState
      INTEGER :: x = 123456789, y = 362436069, z = 521288629, w = 916191069
    END TYPE

    TYPE :: BoxMullerRNGState
      TYPE(KissRNGState) :: k_state
      LOGICAL :: has_cache
      REAL(KIND=REAL64) :: cached_value
    END TYPE BoxMullerRNGState

    ABSTRACT INTERFACE
     FUNCTION pdf(x) result(res)
       IMPORT REAL64
       real(KIND=REAL64), intent(in) :: x
       real(KIND=REAL64) :: res
     end function pdf
    END INTERFACE
    ABSTRACT INTERFACE
     FUNCTION pdf_2v(x, b) result(res)
       IMPORT REAL64
       real(KIND=REAL64), intent(in) :: x
       INTEGER, INTENT(IN) :: b
       real(KIND=REAL64) :: res
     end function pdf_2v
    END INTERFACE

  CONTAINS 

  ! George Marsaglia - type KISS generator
  ! The  KISS (Keep It Simple Stupid) random number generator. Combines:
  ! (1) The congruential generator x(n)=69069*x(n-1)+1327217885, period 2^32.
  ! (2) A 3-shift shift-register generator, period 2^32-1,
  ! (3) Two 16-bit multiply-with-carry generators,
  !     period 597273182964842497>2^59
  ! Overall period>2^123;  Default seeds x,y,z,w.

  FUNCTION random(state)

    TYPE(KissRNGState), INTENT(INOUT) :: state
    REAL(KIND=REAL64) :: random
    INTEGER :: kiss, a1, b1, a2, b2, a3, b3

    state%x = 69069 * state%x + 1327217885

    a1 = state%y
    b1 = 13
    a2 = IEOR(a1, ISHFT(a1, b1))
    b2 = -17
    a3 = IEOR(a2, ISHFT(a2, b2))
    b3 = 5
    state%y  = IEOR(a3, ISHFT(a3, b3))

    state%z = 18000 * IAND(state%z, 65535) + ISHFT(state%z, - 16)
    state%w = 30903 * IAND(state%w, 65535) + ISHFT(state%w, - 16)

    kiss = state%x + state%y + ISHFT(state%z, 16) + state%w

    random = (REAL(kiss, KIND=REAL64) + 2147483648.0_REAL64) / 4294967296.0_REAL64

  END FUNCTION random

  SUBROUTINE random_warmup(state, seed)

    TYPE(KissRNGState), INTENT(INOUT) :: state
    INTEGER, INTENT(IN) :: seed
    INTEGER :: i
    REAL(KIND=REAL64) :: dummy

    state%x = state%x + seed
    state%y = state%y + seed
    state%z = state%z + seed
    state%w = state%w + seed

    ! 'Warm-up' the generator by cycling through a few times
    DO i = 1, 1000
      dummy = random(state)
    ENDDO

  END SUBROUTINE random_warmup

  ! Basic rejection sampling from PDF defined as follows
  ! range is assumed to be [0,1], pdf is assumed normalised to peak at 1
  ! pdf is assumed to be strictly > 0 else the sampling may loop-out
  FUNCTION rejection_sample(state, pdf_fn) RESULT(val)

    TYPE(KissRNGState), INTENT(INOUT) :: state
    PROCEDURE(pdf) :: pdf_fn
    REAL(KIND=REAL64) :: val
    REAL(KIND=REAL64) :: x, y
    INTEGER :: i
    INTEGER, PARAMETER :: max_it = 1000

    val = 0.0_REAL64
    DO i = 1, max_it
      ! Without a better guess, sample uniform in x and y
      x = random(state)
      y = random(state)
        PRINT*, i, x, y, pdf_fn(x)
      IF(pdf_fn(x) > y) THEN
        val = x
        RETURN
      END IF
    END DO

    IF(i > max_it) ERROR STOP "Failed to find a valid random sample"

  END FUNCTION

  ! Basic rejection sampling from PDF defined as follows
  ! range is assumed to be [0,1], pdf is assumed normalised to peak at 1
  ! pdf is assumed to be strictly > 0 else the sampling may loop-out
  ! alpha is assumed == 1, beta is as given
  FUNCTION rejection_sample_beta(state, pdf_fn, beta) RESULT(val)

    TYPE(KissRNGState), INTENT(INOUT) :: state
    PROCEDURE(pdf_2v) :: pdf_fn
    INTEGER :: beta
    REAL(KIND=REAL64) :: val
    REAL(KIND=REAL64) :: x, y
    INTEGER :: i
    INTEGER, PARAMETER :: max_it = 1000

    val = 0.0_REAL64
    DO i = 1, max_it
      ! Without a better guess, sample uniform in x and y
      x = random(state)
      y = random(state)
        PRINT*, i, x, y, pdf_fn(x, beta)
      IF(pdf_fn(x, beta) > y) THEN
        val = x
        RETURN
      END IF
    END DO

    IF(i > max_it) ERROR STOP "Failed to find a valid random sample"

  END FUNCTION

  ! Polar Box_muller
  ! Generates 2 random values per use, so caches the second for the next call
  FUNCTION random_box_muller(stdev, state) RESULT(val)

    REAL(KIND=REAL64), INTENT(IN) :: stdev
    TYPE(BoxMullerRNGState), INTENT(INOUT) :: state
    REAL(KIND=REAL64) :: val

    REAL(KIND=REAL64) :: rand1, rand2, w
    REAL(KIND=REAL64), PARAMETER :: c_tiny = TINY(1.0_REAL64)

    IF (state%has_cache) THEN
      state%has_cache = .FALSE.
      val = state%cached_value * stdev
    ELSE
      state%has_cache = .TRUE.

      DO
        rand1 = random(state%k_state)
        rand2 = random(state%k_state)

        rand1 = 2.0_REAL64 * rand1 - 1.0_REAL64
        rand2 = 2.0_REAL64 * rand2 - 1.0_REAL64

        w = rand1**2 + rand2**2

        IF (w > c_tiny .AND. w < 1.0_REAL64) EXIT
      END DO

      w = SQRT((-2.0_REAL64 * LOG(w)) / w)

      val = rand1 * w * stdev
      state%cached_value = rand2 * w
    END IF

  END FUNCTION random_box_muller

  PURE FUNCTION beta_fn(a, b)
    INTEGER, INTENT(IN) :: a, b
    REAL(KIND=REAL64) :: beta_fn
    beta_fn = gamma(REAL(a)) * gamma(REAL(b)) / gamma(REAL(a) + REAL(b))
  END FUNCTION
  
  FUNCTION beta_pdf(x, beta)
    REAL(KIND=REAL64), INTENT(IN) :: x
    INTEGER, INTENT(IN) :: beta
    REAL(KIND=REAL64) :: beta_pdf
  
    beta_pdf = gamma(1.0 + beta) / gamma(REAL(beta)) * (1.0 - x)**(beta - 1) / REAL(beta)
  END FUNCTION

  FUNCTION random_beta(state, beta) RESULT(ran)
    TYPE(KissRNGState), INTENT(INOUT) :: state
    REAL(KIND=REAL64) :: ran
    INTEGER, INTENT(IN) :: beta

    !ran = rejection_sample_beta(state, beta_pdf, beta)
    ran = beta_fn(1+beta, 1) * (1.0 + beta) * random(state)
    ran = ran**(1.0/(1.0 + beta))
    ran = 1.0_REAL64 - ran
  END FUNCTION

END MODULE