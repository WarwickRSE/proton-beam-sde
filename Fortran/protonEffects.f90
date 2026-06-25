MODULE protonEffects

    USE iso_fortran_env, only: real64 ! REPLACE WITH kinds!
    USE dataDefinedScattering
    USE MaterialFake

    IMPLICIT NONE
 
    SAVE
    REAL(KIND=REAL64), PARAMETER :: mecsq = 0.511_REAL64, mpcsq = 938.346_REAL64
    
    CONTAINS

    PURE FUNCTION bethe_bloch_loss(material, energy) RESULT(loss)
      TYPE(cp_material), INTENT(IN) :: material
      REAL(KIND=REAL64), INTENT(IN) :: energy
      REAL(KIND=REAL64), PARAMETER :: param1 = 0.3072_REAL64
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


    PURE FUNCTION rutherford_and_elastic_rate(material, energy) RESULT(rate)
      TYPE(cp_material), INTENT(IN) :: material
      REAL(KIND=REAL64), INTENT(IN) :: energy
      REAL(KIND=REAL64) :: rate, a
      REAL(KIND=REAL64), PARAMETER :: log_avogadro = LOG(6.0) + 23.0 * LOG(10.0), log_barns_to_cmsq = -24.0 * LOG(10.0)
      INTEGER :: i

      a = SUM(material%nucs%massFraction * material%nucs%A)
      rate = 0.0_REAL64
      DO i = 1, material%no_nucs
        rate = rate + material%nucs(i)%massFraction * &
          evaluate(RU_crossSections(material%nucs(i)%xsec_ind), energy)
      END DO
      rate = rate * EXP(LOG(material%density) + log_avogadro - LOG(a) + log_barns_to_cmsq) ! rate per cm
    END FUNCTION

    !for (unsigned int i = 0; i < at.size(); i++) {
    !  a += x[i] * at[i].a; // average molar mass
    !  ret += x[i] * at[i].el_ruth_rate.evaluate(e);
    !}
    !double log_molecule_density =
    !    log(density) + log_avogadro - log(a); // molecules / cm^3
    !ret *= exp(log_barns_to_cmsq + log_molecule_density);

END MODULE