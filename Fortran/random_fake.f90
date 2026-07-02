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


    REAL(KIND=REAL64), DIMENSION(:), ALLOCATABLE, PRIVATE :: vec
    INTEGER, PRIVATE :: ind

  CONTAINS 

  ! George Marsaglia - type KISS generator
  ! The  KISS (Keep It Simple Stupid) random number generator. Combines:
  ! (1) The congruential generator x(n)=69069*x(n-1)+1327217885, period 2^32.
  ! (2) A 3-shift shift-register generator, period 2^32-1,
  ! (3) Two 16-bit multiply-with-carry generators,
  !     period 597273182964842497>2^59
  ! Overall period>2^123;  Default seeds x,y,z,w.

  SUBROUTINE prime(arr)
    REAL(KIND=REAL64), DIMENSION(:) :: arr
    vec = arr
    ind = 1
  END SUBROUTINE

  FUNCTION yield()
    REAL(KIND=REAL64) :: yield
    IF(ind <= SIZE(vec)) THEN
      yield = vec(ind)
      ind = ind+1
    ELSE
      ERROR STOP "Run of of random numbers, sorry!"
    END IF
  END FUNCTION yield

  FUNCTION random(state)

    TYPE(KissRNGState), INTENT(INOUT) :: state
    REAL(KIND=REAL64) :: random
    random = yield()
  END FUNCTION random

  SUBROUTINE random_warmup(state, seed)

    TYPE(KissRNGState), INTENT(INOUT) :: state
    INTEGER, INTENT(IN) :: seed

  END SUBROUTINE random_warmup

  ! Polar Box_muller
  ! Generates 2 random values per use, so caches the second for the next call
  FUNCTION random_box_muller(stdev, state) RESULT(val)

    REAL(KIND=REAL64), INTENT(IN) :: stdev
    TYPE(BoxMullerRNGState), INTENT(INOUT) :: state
    REAL(KIND=REAL64) :: val
    val = yield()
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
    INTEGER :: beta
    ran = yield()
  END FUNCTION

END MODULE