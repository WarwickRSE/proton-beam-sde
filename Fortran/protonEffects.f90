MODULE protonEffects

    USE iso_fortran_env, only: real64 ! REPLACE WITH kinds!
    USE dataDefinedScattering
    USE MaterialFake

    IMPLICIT NONE
 
    SAVE
    REAL(KIND=REAL64), PARAMETER :: mecsq = 0.511_REAL64, mpcsq = 938.346_REAL64
    REAL(KIND=REAL64), PARAMETER :: log_avogadro = LOG(6.0) + 23.0 * LOG(10.0), log_barns_to_cmsq = -24.0 * LOG(10.0)
    REAL(KIND=REAL64), PARAMETER :: alpha = 1.0_REAL64 / 137.0_REAL64
    
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
    END FUNCTION

    PURE FUNCTION nonelastic_rate(material, energy) RESULT(rate)
      TYPE(cp_material), INTENT(IN) :: material
      REAL(KIND=REAL64), INTENT(IN) :: energy
      REAL(KIND=REAL64) :: rate, a
      
      INTEGER :: i

      a = SUM(material%nucs%massFraction * material%nucs%A)
      rate = 0.0_REAL64
      DO i = 1, material%no_nucs
        rate = rate + material%nucs(i)%massFraction * &
          evaluate(NE_crossSections(material%nucs(i)%xsec_ind), energy)
      END DO
      rate = rate * EXP(LOG(material%density) + log_avogadro - LOG(a) + log_barns_to_cmsq) ! rate per cm
    END FUNCTION

    PURE FUNCTION rutherford_and_elastic_rate(material, energy) RESULT(rate)
      TYPE(cp_material), INTENT(IN) :: material
      REAL(KIND=REAL64), INTENT(IN) :: energy
      REAL(KIND=REAL64) :: rate, a
      INTEGER :: i

      a = SUM(material%nucs%massFraction * material%nucs%A)
      rate = 0.0_REAL64
      DO i = 1, material%no_nucs
        rate = rate + material%nucs(i)%massFraction * &
          evaluate(RU_crossSections(material%nucs(i)%xsec_ind), energy)
      END DO
      rate = rate * EXP(LOG(material%density) + log_avogadro - LOG(a) + log_barns_to_cmsq) ! rate per cm
    END FUNCTION


    !> \brief Computes the standard deviation of Moliere scattering
    !> Based on p. 7 and 8 of [1]
    !> References:
    !>  [1] https://doi.org/10.1088/1361-6560/ae5586
    !> \param material The material through which the proton travels
    !> \param energy The energy of the proton
    !> \param time_step The time step for the simulation
    !> \return The standard deviation of Moliere scattering
    PURE FUNCTION moliere_scattering_sd(material, energy, time_step) RESULT(sd)
      TYPE(cp_material), INTENT(IN) :: material
      REAL(KIND=REAL64), INTENT(IN) :: energy, time_step
      REAL(KIND=REAL64) :: chi_a_sq, chi_c_sq, pv_sq, beta_sq, sd, omega, temp1, temp2
      INTEGER :: i

      ! beta^2 see underneath eq. (3) on p. 6 of [1]
      beta_sq = (2.0 * mpcsq + energy) * energy /(mpcsq + energy)**2

      ! (p*beta)^2
      pv_sq = (2.0 * mpcsq + energy) * energy / (mpcsq + energy)
      pv_sq = pv_sq**2

      chi_c_sq = 0.0_REAL64
      chi_a_sq = 0.0_REAL64
      ! chi_c_sq is sum of individual contributions from each nuclide, 
      ! chi_a_sq is a weighted average on a log scale
      DO i = 1, material%no_nucs
        ! Z_i(Z_i+1)/A_i
        temp1 = material%nucs(i)%massFraction * material%nucs(i)%Z * (material%nucs(i)%Z + 1.0_REAL64) / material%nucs(i)%A
        chi_c_sq = chi_c_sq + temp1
        ! (chi_alpha,i)^2, note pv_sq = (p * beta)^2
        temp2 = 2.007E-5_REAL64 * REAL(material%nucs(i)%Z, KIND=REAL64)**(2.0/3.0) * &
          (1.0_REAL64 + 3.34_REAL64 * (material%nucs(i)%Z * alpha)**2 / beta_sq) * beta_sq / pv_sq
        chi_a_sq = chi_a_sq + temp1 * LOG(temp2)
      END DO
      ! normalise and eliminate the log, 
      ! note denominator in log(chi_a_sq) same as chi_c_sq before multiplying by nucleide independent parameters
      chi_a_sq = EXP(chi_a_sq / chi_c_sq)
      ! multiply chi_c_sq by time_step and and parameters independent of the individual nucleides
      chi_c_sq = chi_c_sq * 0.157_REAL64 * time_step * material%density / pv_sq
      omega = chi_c_sq / (chi_a_sq * 2.0_REAL64 * (1.0_REAL64 - 0.98_REAL64)) ! 0.98 - truncation parameter, see p. 8 [1]
      ! standard deviation
      sd = SQRT(chi_c_sq * ((1.0_REAL64 + omega) * LOG(1.0_REAL64 + omega) / omega - 1.0_REAL64) / (1.0_REAL64 + 0.98_REAL64**2))

    END FUNCTION

END MODULE