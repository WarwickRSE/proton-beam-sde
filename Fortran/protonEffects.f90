MODULE protonEffects

    USE iso_fortran_env, only: real64 ! REPLACE WITH kinds!
    USE dataDefinedScattering
    USE MaterialFake

    IMPLICIT NONE

    CONTAINS


    PURE FUNCTION bethe_bloch_loss(material, energy) RESULT(loss)
      TYPE(cp_material), INTENT(IN) :: material
      REAL(KIND=REAL64), INTENT(IN) :: energy
      REAL(KIND=REAL64), PARAMETER :: mecsq = 0.511_REAL64, mpcsq = 938.346_REAL64, param1 = 0.3072_REAL64
      REAL(KIND=REAL64) :: loss
      REAL(KIND=REAL64) :: betasq
      INTEGER :: i

      betasq = (2.0 * mpcsq + energy) * energy /(mpcsq + energy)**2

      loss = 0.0_REAL64
      DO i = 1, material%no_nucs
         loss = loss + material%nucs(i)%massFraction * param1 * material%nucs(i)%Z * material%density * &
                (LOG(2.0_REAL64 * mecsq * betasq / (material%mee * (1.0_REAL64 - betasq))) - betasq) / &
                (betasq * material%nucs(i)%A) ! MeV/cm
      END DO
!    for (unsigned int i = 0; i < at.size(); i++) {
!      ret += x[i] * 0.3072 * at[i].z * density *
!             (log(2 * mecsq * betasq / (I * (1 - betasq))) - betasq) /
!             (betasq * at[i].a); // MeV / cm

    END FUNCTION


END MODULE