MODULE protonEffects

    USE iso_fortran_env, only: real64 ! REPLACE WITH kinds!
    USE dataDefinedScattering
    USE MaterialFake
    USE randomMod
    USE specialFunctions

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
    
    !> \brief log(a_k^(m, theta)) defined in Eq. (5) in [3]
    !> References:
    !>  [3] https://doi.org/10.1214/16-AAP1236
    !> \param k
    !> \param m 
    !> \param theta
    !> \return The value of log(a_k^(m, theta))
    PURE FUNCTION compute_log_a_param(k, m, theta) RESULT(log_a_param)
      INTEGER, INTENT(IN) :: k, m
      REAL(KIND=REAL64), INTENT(IN) :: theta
      REAL(KIND=REAL64) :: log_a_param

      log_a_param = LOG(theta + 2 * k - 1) + log_pochhammer(theta + m, REAL(k - 1, KIND=REAL64)) - &
        log_factorial(m) - log_factorial(k - m)
    END FUNCTION

    !> \brief Compute b_k^(t, theta) defined in Proposition 1 in [3]
    !> References:
    !>  [3] https://doi.org/10.1214/16-AAP1236
    !> \param k The number of blocks at time t
    !> \param m The number of blocks at time 0
    !> \param t The time to simulate
    !> \param theta
    !> \return The value of b_k^(t, theta)
    PURE FUNCTION compute_b_param(k, m, t, theta) RESULT(b_param)
      INTEGER, INTENT(IN) :: k, m
      REAL(KIND=REAL64), INTENT(IN) :: t, theta
      REAL(KIND=REAL64) :: b_param

      b_param = 1.0_REAL64
      IF (k > 0) THEN
        b_param = EXP(compute_log_a_param(k, m, theta) - k * (k + theta - 1.0_REAL64) * t / 2.0_REAL64)
      END IF
    END FUNCTION

    !> \brief Compute C_m^(theta, t) as defined in Eq. (7) in [3]
    !> References:
    !>  [3] https://doi.org/10.1214/16-AAP1236
    !> \param m The number of blocks at time t
    !> \param t The time to simulate
    !> \param theta 
    !> \return The value of C_m^(theta, t)
    PURE FUNCTION compute_c_param(m, t, theta) RESULT(c_param)
      INTEGER, INTENT(IN) :: m
      REAL(KIND=REAL64), INTENT(IN) :: t, theta
      REAL(KIND=REAL64) :: b_curr, b_next
      INTEGER :: c_param 
      c_param = 0
      b_curr = compute_b_param(m, m, t, theta)
      b_next = compute_b_param(m + 1, m, t, theta)
      ! c_param irst index such that b_next < b_curr (infimum in Eq. (7) in [3])
      DO 
        If (b_next < b_curr) THEN
          EXIT 
        ELSE
          c_param = c_param + 1
          b_curr = b_next
          b_next = compute_b_param(c_param + m + 1, m, t, theta)
        END IF
      END DO
    END FUNCTION

    !> \brief Simulation of ancestral process of Kingsman's coalescent with mutation
    !> Based on Algorithm 2 in [3]
    !> References:
    !>  [3] https://doi.org/10.1214/16-AAP1236
    !> \param t The time to simulate
    !> \param state The state of the random number generator
    !> \param b_state The state of the Box-Muller random number generator
    !> \return The number of blocks at time t
    FUNCTION number_of_blocks(t, state, b_state) RESULT(n_blocks)
      REAL(KIND=REAL64), INTENT(IN) :: t
      TYPE(KissRNGState), INTENT(INOUT) :: state
      TYPE(BoxMullerRNGState), INTENT(INOUT) :: b_state
      REAL(KIND=REAL64) :: mu, sigma, u, theta, smin, smax, sincrement
      INTEGER :: n_blocks, i
      INTEGER, ALLOCATABLE :: k(:)
      
      IF (t <= 0.07_REAL64) THEN
        ! In text underneath Alg. 2 in [2]:
        ! For t < 0.05 a normal approximation can be used.
        mu = 2/t
        sigma = SQRT(2.0_REAL64 / (3 * t))
        n_blocks = NINT(mu + sigma * random_box_muller(1.0_REAL64, b_state))
      ELSE
        ! Positive uniform random variable U ~ Uniform((0, 1]) 
        DO
          u = random(state)
          IF (u > 0) EXIT
        END DO

        smin = 0
        smax = 0
        sincrement = 0
        n_blocks = 0
        k = [0]
        theta = 1.0_REAL64
        DO 
          k(n_blocks + 1) = CEILING(compute_c_param(n_blocks, t, theta) / 2.0_REAL64) ! line 4 Alg. 2 in [3]
          
          !smin =S_k^-(m), smax = S_k^+(m) (Eq 9 in [3])
          DO i = 0, k(n_blocks + 1)-1
            sincrement = compute_b_param(n_blocks + 2 * i, n_blocks, t, theta) - compute_b_param(n_blocks + 2 * i + 1, n_blocks, t, theta)
            smin = smin + sincrement
            smax = smax + sincrement
          END DO
          sincrement = compute_b_param(n_blocks + 2 * k(n_blocks+1), n_blocks, t, theta)
          smin = smin + sincrement - compute_b_param(n_blocks + 2 * k(n_blocks+1) + 1, n_blocks, t, theta)
          smax = smax + sincrement

          ! lines 5-7 Alg 2 in [3]
          DO
            IF (smin >= u .OR. smax <= u) EXIT

            DO i = 0, n_blocks
              k(i + 1) = k(i + 1) + 1
              sincrement = compute_b_param(i + 2 * k(i + 1), i, t, theta)
              smax = smin + sincrement
              smin = smin + sincrement - compute_b_param(i + 2 * k(i + 1) + 1, i, t, theta)
            END DO
          END DO

          ! lines 9-14 Alg 2 in [3]
          IF (smin > u) THEN
            EXIT
          ELSE
            n_blocks = n_blocks + 1
            k = [ (k(i), i=1, SIZE(k)) , 0]
          END IF
        END DO
      END IF
    END FUNCTION

    !> \brief Simulation of the Wright-Fisher diffusion process
    !> Based on Algorithm 2 in [2] with the followwing parameters:
    !> L = 0, theta_1 = theta_2 = (d - 1)/2 = 1 (line 1 Alg. 1 in [2] for d=3)
    !> References:
    !>  [2] https://doi.org/10.1016/j.spl.2020.108836
    !>  [3] https://doi.org/10.1214/16-AAP1236
    !> \param r The current value of the process
    !> \param state The state of the random number generator
    !> \param b_state The state of the Box-Muller random number generator
    !> \return The value of the process at the next time step
    FUNCTION wright_fisher_diffusion(r, state, b_state) RESULT(y)
      REAL(KIND=REAL64), INTENT(IN) :: r
      TYPE(KissRNGState), INTENT(INOUT) :: state
      TYPE(BoxMullerRNGState), INTENT(INOUT) :: b_state
      REAL(KIND=REAL64) :: y
      INTEGER :: n_blocks

      IF (r > 1E-9_REAL64) THEN
        n_blocks = number_of_blocks(r, state, b_state) ! Using Alg. 2 in [3]
        y = random_beta(state, 1 + n_blocks) ! line 3 in Alg. 2 in [2]
      ELSE
        y = r/2
        y = ABS(random_box_muller(SQRT(r * y * (1.0_REAL64 - y)), b_state))
      END IF
    END FUNCTION

    !> \brief Simulation of the spherical Brownian motion process
    !> Based on Algorithm 1 in [2] 
    !> References:
    !>  [1] https://doi.org/10.1088/1361-6560/ae5586
    !>  [2] https://doi.org/10.1016/j.spl.2020.108836
    !> \param dt The time step for the simulation
    !> \param energy The energy of the proton
    !> \param material The material through which the proton travels
    !> \param direction_in The current direction of the proton in spherical coordinates
    !> \param state The state of the random number generator
    !> \param b_state The state of the Box-Muller random number generator
    !> \return The new direction of the proton in spherical coordinates
    FUNCTION spherical_bm(dt, energy, material, direction_in, state, b_state) RESULT(direction_out)
      REAL(KIND=REAL64), INTENT(IN) :: dt, energy, direction_in(2)
      TYPE(cp_material), INTENT(IN) :: material
      TYPE(KissRNGState), INTENT(INOUT) :: state
      TYPE(BoxMullerRNGState), INTENT(INOUT) :: b_state
      REAL(KIND=REAL64) :: direction_out(2), z(3), u(3), w(3), y, denom, theta
    
      ! Convert to Cartesian coordinates
      z(1) = sin(direction_in(1)) * cos(direction_in(2))
      z(2) = sin(direction_in(1)) * sin(direction_in(2))
      z(3) = cos(direction_in(1))

      y = wright_fisher_diffusion(moliere_scattering_sd(material, energy, dt)**2, state, b_state)
      theta = 2 * PI * random(state)
      
      ! Set up defaults for when z is near (0, 0, 1)
      u = [1.0_REAL64/SQRT(2.0_REAL64), &
          1.0_REAL64/SQRT(2.0_REAL64), &
          0.0_REAL64]

      ! u = (e_3 -z)/|e_3 - z| (line 3, Algorithm 1 [2])
      denom = SQRT(z(1)**2 + z(2)**2 + (z(3) - 1.0_REAL64)**2)
      u = -z / denom
      u(3) = u(3) + 1.0_REAL64/denom

      !Evaluate expression to the right of O(z) in line 4, Algorithm 1 [2]
      z(1) = 2 * SQRT(y * (1.0_REAL64 - y)) * cos(theta)
      z(2) = 2 * SQRT(y * (1.0_REAL64 - y)) * sin(theta)
      z(3) = 1.0_REAL64 - 2 * y

      ! Evaluate O(z)z with O(z)=I-2uu^T (line 3/4 Algorithm 1 [2])
      w = z - 2 * u * (DOT_PRODUCT(u, z))
      
      ! New direction in spherical coordinates
      direction_out(1) = ACOS(w(3))
      direction_out(2) = ATAN2(w(2), w(1))
    END FUNCTION

    !> \brief Energy straggling Eq. (5) in [1]
    !> Depends on particle Lorentz factor (thus energy)
    !> Atomic mass (A), atomic number (Z)
    !> References:
    !>  [1] https://doi.org/10.1088/1361-6560/ae5586
    !> \param energy The energy of the proton
    !> \param material The material through which the proton travels
    PURE FUNCTION energy_straggling_sd(energy, material) RESULT(sd)
      REAL(KIND=REAL64), INTENT(IN) :: energy
      TYPE(cp_material), INTENT(IN) :: material
      REAL(KIND=REAL64) :: sd, betasq, z
      REAL(KIND=REAL64) , PARAMETER :: log_c = LOG(29979245800.0_REAL64) ! log(c) in cm/s
      REAL(KIND=REAL64) , PARAMETER :: log_h_bar = -21 * LOG(10.0_REAL64) + LOG(4.136_REAL64) - LOG(2 * PI) ! MeV * s
      INTEGER :: i

      betasq = (2.0 * mpcsq + energy) * energy /(mpcsq + energy)**2
      z = 0
      DO i = 1, material%no_nucs
        z = z + material%nucs(i)%massFraction * material%nucs(i)%Z / material%nucs(i)%A ! electrons per average molecule
      END DO
      ! log_avogadro + LOG(material%density) is the molecule density in molecules/cm^3
      sd = 4 * PI * z * (1 - betasq/2) / SQRT(1 - betasq) * &
        EXP(2 * (log(alpha) + log_c + log_h_bar) + log_avogadro + LOG(material%density)) ! MeV/cm

      sd = SQRT(sd)
    END FUNCTION

    !> \brief Update the energy of the proton based on Bethe-Bloch and energy straggling
    !> Energy update in Eq. (1) in [1]
    !> References:
    !>  [1] https://doi.org/10.1088/1361-6560/ae5586
    !> \param dt The time step for the simulation
    !> \param energy The current energy of the proton
    !> \param material The material through which the proton travels
    !> \param b_state The state of the Box-Muller random number generator
    !> \return The new energy of the proton
    FUNCTION energy_update(dt, energy, material, b_state) RESULT(energy_out)
      REAL(KIND=REAL64), INTENT(IN) :: dt, energy
      TYPE(cp_material), INTENT(IN) :: material
      TYPE(BoxMullerRNGState), INTENT(INOUT) :: b_state
      REAL(KIND=REAL64) :: loss, energy_out

      loss = bethe_bloch_loss(material, energy) * dt
      energy_out = energy - MIN(MAX(loss + SQRT(dt) * energy_straggling_sd(energy, material) * & 
        random_box_muller(1.0_REAL64, b_state), 0.0_REAL64), 2.0 * loss)

      energy_out = MAX(energy_out, 0.0_REAL64)
    END FUNCTION
END MODULE