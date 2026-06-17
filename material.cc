#include "cross_sections.cc"
#include <cstdlib>
#include <vector>

#ifndef MAT
#define MAT

struct Atom {
  Atom(const double a0, const int z0, const std::string ne_r,
       const std::string el_ruth_cs, const std::string ne_ea,
       const double cutoff)
      : a(a0), z(z0), el_ruth_rate(el_ruth_cs, cutoff), ne_rate(ne_r),
        el_ruth_angle_cdf(el_ruth_cs, cutoff), ne_energy_angle(ne_ea) {}

  // Constructor for zero non-elastic rate for hydrogen
  Atom(const double a0, const int z0, const std::string el_ruth_cs,
       const double cutoff, const double back_cutoff)
      : a(a0), z(z0), el_ruth_rate(el_ruth_cs, cutoff, back_cutoff), ne_rate(),
        el_ruth_angle_cdf(el_ruth_cs, cutoff, back_cutoff), ne_energy_angle() {}

  Atom(const Atom &other)
      : a(other.a), z(other.z), el_ruth_rate(other.el_ruth_rate),
        ne_rate(other.ne_rate), el_ruth_angle_cdf(other.el_ruth_angle_cdf),
        ne_energy_angle(other.ne_energy_angle) {}

  /** Compute separation energies for incident particles (S_a)
   * S_a as in Section 6.2.3.2 on p. 137 [4]
   * References:
   * [4] https://doi.org/10.2172/1425114
   * @return S_a
   */
  double s() const {
    double a_c = a + 1; // mass number compound nuclei
    double n_c = a - z; // neutron number compound nuclei
    double z_c = z + 1; // proton number compound nuclei
    double a_a = a; // mass number target nuclei
    double n_a = a - z; // neutron number target nuclei
    double z_a = z; // proton number target nuclei
    // S_a p. 137 [4]
    double ret = 15.68 * (a_c - a_a) -
                 28.07 * (pow(n_c - z_c, 2) / a_c - pow(n_a - z_a, 2) / a_a) -
                 18.56 * (pow(a_c, 2 / 3) - pow(a_a, 2 / 3)) +
                 33.22 * (pow(n_c - z_c, 2) / pow(a_c, 4.0 / 3) -
                          pow(n_a - z_a, 2) / pow(a_a, 4.0 / 3)) -
                 0.717 * (z_c * z_c / pow(a_c, 1.0 / 3) -
                          z_a * z_a / pow(a_a, 1.0 / 3)) +
                 1.211 * (z_c * z_c / a_c - z_a * z_a / a_a);
    return ret;
  }

  /** Sample from distribution of nonelastic collision
   * 
   * References:
   * [1] https://doi.org/10.1088/1361-6560/ae5586
   * [4] https://doi.org/10.2172/1425114
   * @param e
   * @param alpha
   * @param gen
   * @return 
   */
  void sample_nonelastic_collision(double &e, double &alpha,
                                   gsl_rng *gen) const {
    double out_rvalue, out_energy_cm;
    ne_energy_angle.sample(e, out_rvalue, out_energy_cm, gen);

    // [4] Section 6.2.3.2 after Eq (6.4) on p. 137
    double eps_a = a * e / (a + 1); // enrance channel energy p. 136 with a=AWR_A A-target particle, AWR_a=1 a-incident projectile
    double eps_b = (a + 1) * out_energy_cm / a; // emission channel energy p.136
    double e_a = eps_a + s(); // s = S_a
    double e_b = eps_b + s(); // s = S_b
    double x1 = fmin(e_a, 130) * e_b / e_a;
    double x3 = fmin(e_a, 41) * e_b / e_a;
    double aval = 0.04 * x1 + 1.8 * 1e-6 * pow(x1, 3) + 6.7 * 1e-7 * pow(x3, 4);

    
    double cdfc2 = out_rvalue * cosh(aval) - sinh(aval);
    double cdfc1 = 2 * sinh(aval);
    double u2 = gsl_rng_uniform(gen);
    double z1 = cdfc1 * u2 + cdfc2; // C [1] p.10
    double z2 =
        (z1 + sqrt(pow(z1, 2) - pow(out_rvalue, 2) + 1)) / (out_rvalue + 1); // \mu p.10 [1]
    
    
    double out_angle_cm = log(z2) / aval;
    double out_energy_lab =
        out_energy_cm + e / pow(a + 1, 2) +
        2 * sqrt(out_energy_cm * e) * out_angle_cm / (a + 1);
    double out_angle_lab = sqrt(out_energy_cm / out_energy_lab) * out_angle_cm +
                           sqrt(e / out_energy_lab) / (a + 1);
    e = out_energy_lab;
    alpha = out_angle_lab;
    if (out_energy_cm == 0) {
      e = 0;
      alpha = 1; // If outgoing energy is 0, then out_angle_lab should be 1,
                 // rounding errors allow it to be slightly above 1 which is
                 // invalid.
    }
  }

  const double a;
  const int z;
  CS_1d el_ruth_rate, ne_rate;
  CS_2d el_ruth_angle_cdf;
  CS_3d ne_energy_angle;
};

struct Material {

  Material(std::vector<Atom> &atoms, const std::vector<int> &id,
           const std::vector<double> &x0, const double d0, const double I0)
      : density(d0), I(I0 / 1e6), x(x0), at() {
    for (unsigned int i = 0; i < id.size(); i++) {
      at.push_back(atoms[id[i]]);
    }
  }

  Material() : density(), I(), x(), at() {}

  void read_material(const std::string filename, std::vector<Atom> &atoms) {
    std::ifstream file;
    file.open(filename + ".txt");
    std::string line, token;
    getline(file, line);
    density = atof(line.c_str());
    getline(file, line);
    I = atof(line.c_str()) / 1e6;
    while (getline(file, line)) {
      std::stringstream iss;
      iss << line;
      getline(iss, token, ' ');
      at.push_back(atoms[atoi(token.c_str()) - 1]);
      getline(iss, token, ' ');
      x.push_back(atof(token.c_str()));
    }
    file.close();
    return;
  }

  /**
   * Bethe-bloch ineleastic mean energy loss, Eq 3 (Eq 4 for mixture)
  per unit track length
  * Depends on density (rho), atomic number (Z), atomic Mass (A) and
   mean excitation energy (I)
   * ALSO depends on e, the energy of the Particle
  @param e Energy of particle being scattered
   */
  double bethe_bloch(const double e) const {
    double mecsq = 0.511;   // mass of electron * speed of light squared, MeV
    double mpcsq = 938.346; // mass of proton * speed of light squared, MeV
    double betasq = (2 * mpcsq + e) * e / pow(mpcsq + e, 2);
    double ret = 0;
    for (unsigned int i = 0; i < at.size(); i++) {
      ret += x[i] * 0.3072 * at[i].z * density *
             (log(2 * mecsq * betasq / (I * (1 - betasq))) - betasq) /
             (betasq * at[i].a); // MeV / cm
    }
    return ret;
  }

/** Moliere's small angle elastic
* See P 7 and 8
* chi_c and chi_a and omega are partial factors which are not named
* v is used for w in eq for sigma
* RETURNS the std-dev
* @param e energy of particle being scattered
* @param dt distance? travelled?
* @return sigma-E for this process
*/
  double multiple_scattering_sd(const double e, const double dt) const {
    double mpcsq = 938.346; // mass of proton * speed of light squared, MeV
    double pv = (2 * mpcsq + e) * e / (mpcsq + e);
    double betasq = (2 * mpcsq + e) * e / pow(mpcsq + e, 2);
    double p = pv / sqrt(betasq); // momentum in MeV / c.
    // effective chi_c_sq is just the sum of individual elements
    double chi_c_sq = 0;
    std::vector<double> chi_a_sq_vec(at.size());
    for (unsigned int i = 0; i < at.size(); i++) {
      chi_c_sq += x[i] * at[i].z * (at[i].z + 1.0) / at[i].a;
      chi_a_sq_vec[i] = 2.007e-5 * pow(at[i].z, 2 / 3) *
                        (1 + 3.34 * pow(at[i].z / (137 * sqrt(betasq)), 2)) /
                        (p * p);
    }
    chi_c_sq *= 0.157 * dt * density / (pv * pv);
    // effective chi_a_sq is a weighted average on the log-scale
    double chi_a_sq = 0;
    double denom = 0;
    for (unsigned int i = 0; i < at.size(); i++) {
      chi_a_sq +=
          x[i] * at[i].z * (at[i].z + 1) * log(chi_a_sq_vec[i]) / at[i].a;
      denom += x[i] * at[i].z * (at[i].z + 1.0) / at[i].a;
    }
    chi_a_sq = exp(chi_a_sq / denom);
    double omega = chi_c_sq / chi_a_sq;
    double F = 0.98;
    double v = omega / (2 * (1 - F));
    double ret = sqrt(chi_c_sq * ((1 + v) * log(1 + v) / v - 1) / (1 + F * F));
    return ret;
  }

/** Energy straggling Eq 5
* Depends on particle Lorentz factor (thus energy)
* Atomic mass (A), atomic number (Z), 
* @param e Energy of particle being scattered
* @return The volatility (which relates to SD)
*/
  double energy_straggling_sd(const double e) const {
    double alpha = 1 / 137.0;
    double log_hbar = -21 * log(10) + log(4.136) - log(2 * M_PI); // MeV * s
    double log_c = log(29979245800);                              // cm / s
    double log_avogadro = log(6) + 23 * log(10);
    double mpcsq = 938.346; // mass of proton * speed of light squared, MeV
    double betasq = (2 * mpcsq + e) * e / pow(mpcsq + e, 2);
    double z = 0;
    for (unsigned int i = 0; i < at.size(); i++) {
      z += x[i] * at[i].z / at[i].a; // electrons per average molecule
    }
    double log_molecule_density =
        log(density) + log_avogadro; // molecules / cm^3
    double ret =
        4 * M_PI * z * (1 - betasq / 2) / sqrt(1 - betasq) *
        exp(2 * (log(alpha) + log_hbar + log_c) + log_molecule_density);
    return sqrt(ret);
  }

  /**
   * @brief 
   * 
   * @param e 
   * @return 
   */
  double nonelastic_rate(const double e) const {
    double log_avogadro = log(6) + 23 * log(10); //log (N_A)
    double log_barns_to_cmsq = -24 * log(10); // log (10^{-24})
    double a = 0;
    double ret = 0;
    for (unsigned int i = 0; i < at.size(); i++) {
      a += x[i] * at[i].a; // average molar mass A
      ret += x[i] * at[i].ne_rate.evaluate(e);
    }
    double log_molecule_density =
        log(density) + log_avogadro - log(a); // log(rho N_A 10^{-24}/A) molecules / cm^3
    ret *= exp(log_barns_to_cmsq + log_molecule_density);
    return ret; // rate per cm
  }

  double rutherford_and_elastic_rate(const double e) const {
    double log_avogadro = log(6) + 23 * log(10);
    double log_barns_to_cmsq = -24 * log(10); // log (10^{-24})
    double a = 0;
    double ret = 0;
    for (unsigned int i = 0; i < at.size(); i++) {
      a += x[i] * at[i].a; // average molar mass
      ret += x[i] * at[i].el_ruth_rate.evaluate(e);
    }
    double log_molecule_density =
        log(density) + log_avogadro - log(a); // molecules / cm^3
    ret *= exp(log_barns_to_cmsq + log_molecule_density);
    return ret; // rate per cm
  }
  /**
   * @brief Updates the direction of transport according to scattering angles
   * 
   * @param ang - direction of transport given in spherical coordinates on unit sphere
   * @param alpha - polar scattering angle
   * @param beta - azimutal scattering angle 
   * @return 
   */
  void compute_new_angle(std::vector<double> &ang, const double alpha,
                         const double beta) const {

    // Compute new transport direction 
    // omega_new = cos(alpha)e_r + sin(alpha)(sin(beta) e_theta + cos(beta) e_phi)
    // Where:
    // theta=ang[0], phi=ang[1]
    // e_r = (sin(theta)cos(phi), sin(theta)sin(phi), cos(theta))
    // e_theta = (cos(theta)cos(phi),cos(theta)sin(phi), -sin(theta))
    // e_phi = (-sin(phi), cos(phi), 0)  
    double omega_new1 =
        sin(ang[0]) * cos(ang[1]) * cos(alpha) +
        (cos(ang[0]) * cos(ang[1]) * sin(beta) - sin(ang[1]) * cos(beta)) *
            sin(alpha);
    double omega_new2 =
        sin(ang[0]) * sin(ang[1]) * cos(alpha) +
        (cos(ang[0]) * sin(ang[1]) * sin(beta) + cos(ang[1]) * cos(beta)) *
            sin(alpha);
    double omega_new3 =
        cos(ang[0]) * cos(alpha) - sin(ang[0]) * sin(beta) * sin(alpha);
    double magnitude =
        std::sqrt(omega_new1 * omega_new1 + omega_new2 * omega_new2 +
                  omega_new3 * omega_new3);
    omega_new1 /= magnitude;
    omega_new2 /= magnitude;
    omega_new3 /= magnitude;

    // Convert to spherical coordinates
    ang[0] = acos(omega_new3);
    ang[1] = atan2(omega_new2, omega_new1);
    return;
  }

  void nonelastic_scatter(std::vector<double> &ang, double &e,
                          gsl_rng *gen) const {
    double beta = 2 * M_PI * gsl_rng_uniform(gen);
    double rate = 0;
    // QUERY - is ne_rate.evaluate pure? (note ind is changing) Where does e change in this function?
    for (unsigned int i = 0; i < at.size(); i++) {
      rate += x[i] * at[i].ne_rate.evaluate(e);
    }
    double u = gsl_rng_uniform(gen);
    double ind = 0;
    double tmp = x[ind] * at[ind].ne_rate.evaluate(e) / rate;
    while (tmp < u) {
      ind++;
      tmp += x[ind] * at[ind].ne_rate.evaluate(e) / rate;
    }
    double alpha;
    // ENDF non-elastic scattering, both energy + angle from CM to LAB
    at[ind].sample_nonelastic_collision(e, alpha, gen);
    alpha = acos(alpha);
    compute_new_angle(ang, alpha, beta);
    return;
  }

  void rutherford_elastic_scatter(std::vector<double> &ang, double &e,
                                  gsl_rng *gen) const {
    double beta = 2 * M_PI * gsl_rng_uniform(gen);
    double rate = 0;
    // QUERY - is el_ruth_rate.evaluate pure? Where does e change in this function?
    for (unsigned int i = 0; i < at.size(); i++) {
      rate += x[i] * at[i].el_ruth_rate.evaluate(e);
    }
    double u = gsl_rng_uniform(gen);
    double ind = 0;
    double tmp = x[ind] * at[ind].el_ruth_rate.evaluate(e) / rate;
    while (tmp < u) {
      ind++;
      tmp += x[ind] * at[ind].el_ruth_rate.evaluate(e) / rate;
    }
    double alpha;
    alpha = at[ind].el_ruth_angle_cdf.sample(e, gen);
    compute_new_angle(ang, alpha, beta);
    return;
  }

  double density; // density, g / cm^3
  double I;       // mean excitation energy, MeV
  std::vector<double> x;
  std::vector<Atom> at;
};

#endif
