#include "material.cc"
#include <cmath>
#include <cstdlib>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_gamma.h>
#include <iostream>
#include <vector>

#ifndef PB
#define PB

struct proton_path {

  /**
   * @brief Construct a proton path
   * 
   * Starts with 
   * @param e0 
   * @param dt 
   * @param absorption_e 
   * @param change_points_x 
   * @param change_points_y 
   * @param interval_materials 
   * @param materials 
   */
  proton_path(const double e0, const double dt, const double absorption_e,
              const std::vector<double> &change_points_x,
              const std::vector<double> &change_points_y,
              const std::vector<std::vector<int>> &interval_materials,
              std::vector<Material> &materials)
      : energy(1), s(1), x(1), omega(1), u(3, 0), z(3, 0), w(3, 0) {
    unsigned int n =
        solve_central_ode(e0, dt, absorption_e, change_points_x,
                          change_points_y, interval_materials, materials)
            .size();
    std::vector<double> tmp_x(3, 0);
    std::vector<double> tmp_w(2, 0);
    energy.resize(n, 0);
    s.resize(n, 0);
    x.resize(n, tmp_x);
    omega.resize(n, tmp_w);
  }

  void reset(const double e0, const std::vector<double> x0,
             const std::vector<double> w0) {
    x[0] = x0;
    omega[0] = w0;
    energy[0] = e0;
    s[0] = 0;
    return;
  }

  /**
   * @brief ??
   *
   * Run proton track through specified material setup applying only BetheBloch. If energy gets below absorption_e the track ends
   * @param e0 The initial energy
   * @param dt The track step size
   * @param absorption_e The minimum energy to continue tracking
   * @param change_points_x Co-ordinate location of change to next material
   * @param change_points_y Co-ordinate location of change to next material 
   * @param interval_materials 2-D list of materials in x and y for each region
   * @param materials Map from Material number to actual material
   * @return Vector of energy at each step
   */
  std::vector<double>
  solve_central_ode(const double e0, const double dt, const double absorption_e,
                    const std::vector<double> &change_points_x,
                    const std::vector<double> &change_points_y,
                    const std::vector<std::vector<int>> &interval_materials,
                    std::vector<Material> &materials) {
    std::vector<double> e(1, e0);
    double x = 0;
    int material_index = 1;
    int y_half = 0;
    while (e.back() > absorption_e) {
      while (x >= change_points_x[material_index]) {
        material_index++;
      }
      if (change_points_y[material_index] > 0) {
        y_half = 0;
      } else {
        y_half = 1;
      }
      e.push_back(
          e.back() -
          materials[interval_materials[material_index - 1][y_half]].bethe_bloch(
              e.back()) *
              dt);
      x += dt;
    }
    return e;
  }

  /**
   * @brief Compute log(a_{km}^theta) based on the definition given in Eq (5) in [3]
   * 
   * References: 
   * [3] https://doi.org/10.1214/16-AAP1236
   * @param k
   * @param m
   * @param t
   * @param theta 
   * @return 
   */
  double log_a(const int k, const int m, const double theta) const {
    // gsl_sf_lnpoch(a, x) - logarithm of Pochhammer symbol (log(Gamma(a+x)/Gamma(a)))
    // gsl_sf_lnfact - logarithm of factorial
    double ret = log(theta + 2 * k - 1) + gsl_sf_lnpoch(theta + m, k - 1) -
                 gsl_sf_lnfact(m) - gsl_sf_lnfact(k - m);
    return ret;
  }

  /**
   * @brief Compute b_k^{(t, theta)}(m) given in Proposition above Eq (7) in [3]
   * 
   * References: 
   * [3] https://doi.org/10.1214/16-AAP1236
   * @param k
   * @param m
   * @param t
   * @param theta 
   * @return 
   */
  double b(const int k, const int m, const double t, const double theta) const {
    double ret = 1;
    if (k > 0) {
      ret = exp(log_a(k, m, theta) - k * (k + theta - 1) * t / 2);
    }
    return ret;
  }
  /**
   * @brief Compute C_m^{(t, theta)} given by Eq (7) in [3]
   * 
   * References: 
   * [3] https://doi.org/10.1214/16-AAP1236
   * @param m
   * @param t
   * @param theta 
   * @return 
   */
  int c(const int m, const double t, const double theta) const {
    int i = 0;
    double b_curr = b(m, m, t, theta);
    double b_next = b(m + 1, m, t, theta);
    while (b_next >= b_curr) {
      i++;
      b_curr = b_next;
      b_next = b(i + m + 1, m, t, theta);
    }
    return i;
  }

  /**
   * @brief Simulation of ancestral process of Kingsman's coalescent with mutation
   * 
   * 
   * Based on Algorithm 2 in [3]
   * References: 
   * [2] https://doi.org/10.1016/j.spl.2020.108836
   * [3] https://doi.org/10.1214/16-AAP1236
   * @param t
   * @param gen 
   * @return 
   */
  int number_of_blocks(const double t, gsl_rng *gen) const {
    int m = 0;
    double theta = 1;
    if (t < 0.07) {
      // In text underneath Alg. 2 in [2]:
      // For t < 0.05 a normal approximation can be used.
      double mu = 2 / t;
      double sigma = sqrt(2 / (3 * t));
      m = round(mu + sigma * gsl_ran_gaussian_ziggurat(gen, 1));
    } else {
      std::vector<int> k(1, 0);
      bool proceed = true;
      double u = gsl_rng_uniform_pos(gen); // U ~ Uniform([0, 1])
      double smin = 0, smax = 0, increment = 0;
      while (proceed) {
        k[m] = ceil(c(m, t, theta) / 2.0); // line 4 Alg. 2 in [3]

        // smin =S_k^-(m), smax = S_k^+(m) (Eq 9 in [3])
        for (int i = 0; i < k[m]; i++) {
          increment = b(m + 2 * i, m, t, theta) - b(m + 2 * i + 1, m, t, theta);
          smin += increment;
          smax += increment;
        }
        increment = b(m + 2 * k[m], m, t, theta);
        smin += increment - b(m + 2 * k[m] + 1, m, t, theta);
        smax += increment;

        // lines 5-7 Alg 2 in [3]
        while (smin < u && u < smax) {
          for (int i = 0; i <= m; i++) {
            k[i]++;
            increment = b(i + 2 * k[i], i, t, theta);
            smax = smin + increment;
            smin += increment - b(i + 2 * k[i] + 1, i, t, theta);
          }
        }

        // lines 9-14 Alg 2 in [3]
        if (smin > u) {
          proceed = false;
        } else {
          k.push_back(0);
          m++;
        }
      }
    }
    return m;
  }
  /**
   * @brief Simulation of Wright-Fisher Diffusion
   * 
   * 
   * Based on Algorithm 2 in [2] with the following paramters fixed:
   * x = 0 -> L ~ Binomial(M, x) = 0 (line 1 in Alg. 1 in [2] and line 2 in Alg. 2 in [2])
   * theta_1 = theta_2 = (d - 1)/2 = 1 (line 1 Alg. 1 in [2] for d=3)
   * References: 
   * [2] https://doi.org/10.1016/j.spl.2020.108836
   * [3] https://doi.org/10.1214/16-AAP1236
   * @param r
   * @param gen 
   * @return 
   */
  double wright_fisher(const double r, gsl_rng *gen) const {
    double y;
    if (r > 1e-9) {
      int m = number_of_blocks(r, gen); // Using Alg. 2 [3]
      y = gsl_ran_beta(gen, 1, 1 + m); // line 3 Alg. 2 [2] with theta_1=theta_2 = 1 and L =0
    } else {
      // ? Haven't found a reference for this yet.
      y = r / 2;
      y = fabs(gsl_ran_gaussian_ziggurat(gen, sqrt(r * y * (1 - y))));
    }
    return y;
  }

  /**
   * @brief Spherical Brownian motion evaluation - updates energy, energy difference s, position x and direction of transport omega
   * 
   * Eq (1) in [1]. Energy, position and direction are updated.
   * Based on Algorithm 1 in [2]
   * References: 
   * [1] https://doi.org/10.1088/1361-6560/ae5586
   * [2] https://doi.org/10.1016/j.spl.2020.108836
   * @param dt 
   * @param ix 
   * @param gen 
   * @param mat 
   * @param prev_x_change 
   * @param next_x_change 
   * @param y_change 
   * @return 
   */
  double spherical_bm(const double dt, int &ix, gsl_rng *gen,
                      const Material &mat, const double prev_x_change,
                      const double next_x_change, const double y_change) {

    // Convert spherical coordinates into cartesian coordinates
    z[0] = sin(omega[ix - 1][0]) * cos(omega[ix - 1][1]);
    z[1] = sin(omega[ix - 1][0]) * sin(omega[ix - 1][1]);
    z[2] = cos(omega[ix - 1][0]);

    // Radial component (line 1, Algorithm 1 [2])
    double y = wright_fisher(
        pow(mat.multiple_scattering_sd(energy[ix - 1], dt), 2), gen);
    
    // Angular component uniform on unit sphere (line 2, Algorithm 1 [2])
    double theta = 2 * M_PI * gsl_rng_uniform(gen);

    // Set up defaults for when z is near (0, 0, 1)
    u[0] = 1 / sqrt(2);
    u[1] = 1 / sqrt(2);
    u[2] = 0;

    // u = (e_3 -z)/|e_3 - z| (line 3, Algorithm 1 [2])
    double denom = sqrt(z[0] * z[0] + z[1] * z[1] + (z[2] - 1) * (z[2] - 1));
    if (denom > 1e-10) {
      u[0] = -z[0] / denom;
      u[1] = -z[1] / denom;
      u[2] = (1 - z[2]) / denom;
    }

    // Evaluate expression to the right of O(z) in line 4, Algorithm 1 [2]
    z[0] = 2 * sqrt(y * (1 - y)) * cos(theta);
    z[1] = 2 * sqrt(y * (1 - y)) * sin(theta);
    z[2] = 1 - 2 * y;

    // Evaluate O(z)z with O(z)=I-uu^T (line 3/4 Algorithm 1 [2])
    w[0] = (1 - 2 * u[0] * u[0]) * z[0] - 2 * u[0] * u[1] * z[1] -
           2 * u[0] * u[2] * z[2];
    w[1] = (1 - 2 * u[1] * u[1]) * z[1] - 2 * u[0] * u[1] * z[0] -
           2 * u[1] * u[2] * z[2];
    w[2] = (1 - 2 * u[2] * u[2]) * z[2] - 2 * u[0] * u[2] * z[0] -
           2 * u[1] * u[2] * z[1];

    // Resize data structures
    if (ix == int(omega.size())) {
      omega.resize(2 * omega.size(), omega.back());
      energy.resize(2 * energy.size(), energy.back());
      x.resize(2 * x.size(), x.back());
      s.resize(2 * s.size(), s.back());
    }

    // Convert to spherical coordinates
    omega[ix][0] = acos(w[2]);
    omega[ix][1] = atan2(w[1], w[0]);

    // I'm not entirely sure, but it looks like the following bit of code is
    // computing an average direction between omega[ix-1] and omega[ix]
    // Something like
    // direction_x = int_0^1 sin(v0 + (v1-v0)t)cos(w0 + (w1-w0)t)dt
    // direction_y = int_0^1 sin(v0 + (v1-v0)t)sin(w0 + (w1-w0)t)dt
    // direction_z = int_0^1 cos(v0 + (v1-v0)t)dt
    // If the denominator in the integration is too small, the integral
    // average is simply replaced with the arithmetic mean of
    // the two end points.
    double v0 = omega[ix - 1][0];
    double v1 = omega[ix][0];
    double w0 = omega[ix - 1][1];
    double w1 = omega[ix][1];

    // Check for division by zero in x and y position updates
    double denom_xy = (v0 - v1 + w0 - w1) * (v0 - v1 - w0 + w1);
    double direction_x, direction_y;
    double time_step = dt;
    if (fabs(denom_xy) > 1e-9) {
      direction_x = ((v0 - v1) * (cos(v0) * cos(w0) - cos(v1) * cos(w1)) +
                     (w0 - w1) * (sin(v0) * sin(w0) - sin(v1) * sin(w1))) /
                    denom_xy;
      direction_y = ((w0 - w1) * (cos(w0) * sin(v0) - cos(w1) * sin(v1)) -
                     (v0 - v1) * (cos(v0) * sin(w0) - cos(v1) * sin(w1))) /
                    denom_xy;\


      // ? Ensures timestep is small enought to stay in same material
      // ? Looks a bit like it would all be clearer if direction_x had the opposite sign
      if (direction_x < 0) {
        time_step =
            fmin(time_step, -(next_x_change - x[ix - 1][0]) / direction_x);
      } else {
        time_step =
            fmin(time_step, (x[ix - 1][0] - prev_x_change) / direction_x);
      }
      if ((direction_y > 0 && x[ix - 1][1] < y_change) ||
          (direction_y < 0 && x[ix - 1][1] > y_change)) {
        time_step = fmin(time_step, (y_change - x[ix - 1][1]) / direction_y);
      }

      // Update position x (x and y component)
      x[ix][0] = x[ix - 1][0] - time_step * direction_x;
      x[ix][1] = x[ix - 1][1] + time_step * direction_y;
    } else {
      // Linear approximation when denominator is too small
      direction_x = (sin(v0) * cos(w0) + sin(v1) * cos(w1)) / 2;
      direction_y = (sin(v0) * sin(w0) + sin(v1) * sin(w1)) / 2;

      // ? Ensures timestep is small enought to stay in same material
      // ? Follows the same logic as above
      if (direction_x > 0) {
        time_step =
            fmin(time_step, (next_x_change - x[ix - 1][0]) / direction_x);
      } else {
        time_step =
            fmin(time_step, -(x[ix - 1][0] - prev_x_change) / direction_x);
      }
      if ((direction_y > 0 && x[ix - 1][1] < y_change) ||
          (direction_y < 0 && x[ix - 1][1] > y_change)) {
        time_step = fmin(time_step, (y_change - x[ix - 1][1]) / direction_y);
      }

      // Update position x (x and y component)
      x[ix][0] = x[ix - 1][0] + time_step * direction_x;
      x[ix][1] = x[ix - 1][1] + time_step * direction_y;
    }
    // Z position update
    if (fabs(v0 - v1) > 1e-9) {
      x[ix][2] = x[ix - 1][2] + (sin(v0) - sin(v1)) * time_step / (v0 - v1);
    } else {
      x[ix][2] = x[ix - 1][2] + time_step * (cos(v0) + cos(v1)) / 2;
    }

    // Eq (1) in [1]. Energy update (first line)
    // ? I get the wrong sign for one of the terms when comparing with the paper
    // ? Does not have (1-u_n) term (is this computed inside simulate?)
    double bethe_block_update = mat.bethe_bloch(energy[ix - 1]) * time_step;
    energy[ix] = energy[ix - 1] -
                 fmin(fmax(bethe_block_update +
                               sqrt(time_step) *
                                   mat.energy_straggling_sd(energy[ix - 1]) *
                                   gsl_ran_gaussian_ziggurat(gen, 1),
                           0),
                      2 * bethe_block_update);

    // Ensure engergy is positive
    energy[ix] = fmax(energy[ix], 0);

    // Update energy increment
    s[ix] = energy[ix - 1] - energy[ix];

    // ? Would it be a bit clearer to update this inside the while loop in simulate?
    ix++;
    return time_step;
  }

  int simulate(const double dt, const double absorption_energy,
               const std::vector<double> &change_points_x,
               const std::vector<double> &change_points_y,
               const std::vector<std::vector<int>> &interval_materials,
               const std::vector<Material> &materials, gsl_rng *gen) {
    double nonelastic_jump_rate;
    double rutherford_elastic_jump_rate;
    double alpha;
    double time_step = dt;
    int ix = 1;
    int material_index = 1;
    int y_half = 0;
    if (change_points_y[material_index - 1] < 0) {
      y_half = 1;
    }
    while (energy[ix - 1] > absorption_energy) {
      time_step = spherical_bm(
          dt, ix, gen,
          materials[interval_materials[material_index - 1][y_half]],
          change_points_x[material_index - 1], change_points_x[material_index],
          change_points_y[material_index - 1]);

      // NOTE: ix was incremented inside spherical_bm. 
      if (energy[ix - 1] > absorption_energy) {
        nonelastic_jump_rate =
            materials[interval_materials[material_index - 1][y_half]]
                .nonelastic_rate(energy[ix - 1]);
        rutherford_elastic_jump_rate =
            materials[interval_materials[material_index - 1][y_half]]
                .rutherford_and_elastic_rate(energy[ix - 1]);
        alpha = rutherford_elastic_jump_rate + nonelastic_jump_rate;
        if (gsl_rng_uniform(gen) < 1 - exp(-alpha * time_step)) {
          if (gsl_rng_uniform(gen) < rutherford_elastic_jump_rate / alpha) {
            materials[interval_materials[material_index - 1][y_half]]
                .rutherford_elastic_scatter(omega[ix - 1], energy[ix - 1], gen);
          } else {
            materials[interval_materials[material_index - 1][y_half]]
                .nonelastic_scatter(omega[ix - 1], energy[ix - 1], gen);
          }
        }
      }
      if (fabs(x[ix - 1][0] - change_points_x[material_index]) < 1e-9) {
        material_index++;
        if (x[ix - 1][1] > change_points_y[material_index - 1]) {
          y_half = 1;
        } else {
          y_half = 0;
        }
      } else if (fabs(x[ix - 1][0] - change_points_x[material_index - 1]) <
                 1e-9) {
        material_index--;
        if (x[ix - 1][1] > change_points_y[material_index - 1]) {
          y_half = 1;
        } else {
          y_half = 0;
        }
      } else if (fabs(x[ix - 1][1] - change_points_y[material_index - 1]) <
                 1e-9) {
        y_half = (y_half + 1) % 2;
      }
    }
    return ix;
  }

  std::vector<double> energy, s;
  // x contains the position at each step, omega the direction of tranport
  // omega describes a vector on the unit sphere in spherical coordinates
  // with the corresponding cartesian coordinates at ix given by 
  // (sin(omega[ix][0])cos(omega[ix][1]), sin(omega[ix][0])sin(omega[ix][1]), cos(omega[ix][0]))
  std::vector<std::vector<double>> x, omega;
  // Dummy vectors for spherical BM
  std::vector<double> u, z, w;
};

#endif
