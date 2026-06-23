#include <cmath>
#include <cstdlib>
#include <fstream>
#include <gsl/gsl_rng.h>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#ifndef CS
#define CS

/**
 * Struct for scattering rates
 * 
 * Given by Eq. (6) in [1] (or Eq. (12) for hydrogen) for elastic scattering
 * and directly form data for inelastic scattering
 * 
 * References:
 * [1] https://doi.org/10.1088/1361-6560/ae5586
 */
struct CS_1d {
  /** 
   * @brief Constructor: Retrieve scattering rate \sigma_e for non-elastic scattering from data
   * 
   * The data in the file is expected in the following format:
   *  - First line contains the energy values.
   *  - Second line contains the scattering rates corresponding to the energy values.
   * NOTE: A factor of 10^{-24} * N_A * rho/A needs to be multiplied to get the correct rates.
   *    
   * 
   * @param filename
   * @param cutoff
   */
  CS_1d(const std::string filename) : energy(), rate() {
    std::ifstream file;
    file.open(filename);
    std::string line, token;
    std::stringstream iss;

    // First line contains energy values
    getline(file, line);
    iss << line;
    while (getline(iss, token, ' ')) {
      energy.push_back(atof(token.c_str()));
    }

    // Second line contains scattering rates corresponding to the energy values
    getline(file, line);
    std::stringstream iss2;
    iss2 << line;
    while (getline(iss2, token, ' ')) {
      rate.push_back(atof(token.c_str()));
    }
  }
  /** 
   * @brief Constructor: Retrieve scattering rate \sigma_e for non-elastic scattering from data
   * 
   * The data in the file is expected in the following format:
   *  - First line contains the energy values.
   *  - For each energy value, there are two more lines in the file, 
   *    the first of which contains the scattering angles and the second 
   *    contains the scattering rate (up to a factor) corresponding to each angle
   *    (i.e. the upper limit of the integral in Eq. (7) in [1] is cos(angle)).
   *  - The angles are given in radians in decreasing order between pi and 0.
   * 
   * NOTE: The factor still missing to obtain the rate given in Eq. (7) is 10^{-24} * N_A * rho/A
   *  This factor is added in the 'rutherford_and_elastic_rate' method of the Material.
   *    
   * 
   * @param filename
   * @param cutoff
   */
  CS_1d(const std::string filename, const double cuttoff) : energy(), rate() {
    std::ifstream file;
    file.open(filename);
    std::string line, token;

    // First line contains energy values
    getline(file, line);
    std::stringstream iss;
    iss << line;
    while (getline(iss, token, ' ')) {
      energy.push_back(atof(token.c_str()));
    }
    double tmp_val, tmp_val_old = 0, lin_inter_val = 0;
    int tmp_count, tmp_count_2;
    bool lin_inter_bool;
    while (getline(file, line)) {
      // Find angle closest to cutoff
      tmp_count = 0;
      lin_inter_bool = true;
      std::stringstream iss2;
      iss2 << line;
      while (getline(iss2, token, ' ')) {
        tmp_val = atof(token.c_str()); // current angle
        if (tmp_val > cuttoff) {
          tmp_count++;
          tmp_val_old = tmp_val; // next angle if larger than cutoff
        } else if (lin_inter_bool) {
          // if current angle is smaller than cutoff,
          // precompute factor for linear interpolation of rate
          lin_inter_val = (tmp_val - cuttoff) / (tmp_val - tmp_val_old);
          lin_inter_bool = false; // Only compute lin_inter_val for first angle smaller than cutoff
        }
      }

      // Find rate corresponding to angles closest to cutoff angle and
      // compute rate for cutoff angle by linear inerpolation
      tmp_count_2 = 0;
      lin_inter_bool = true;
      getline(file, line);
      std::stringstream iss3;
      iss3 << line;
      while (getline(iss3, token, ' ')) {
        tmp_val = atof(token.c_str());
        if (tmp_count_2 < tmp_count) {
          tmp_count_2++;
          tmp_val_old = tmp_val;
          // Update value of rate until the correct angle is reached
        } else if (lin_inter_bool) {
          // Linear inerpolation of rate
          rate.push_back(tmp_val * lin_inter_val +
                         (1 - lin_inter_val) * tmp_val_old);
          lin_inter_bool = false;
        }
      }
    }
    file.close();
  }

  /**
   * Convert angle from centre-of-mass frame to laboratory frame for hydrogen
   * 
   * NOTE: This function is repeated as a method in the CS_2d struct.
   * @param ang angle in centre-of-mass frame
   * @param E energy
   * @return angle in laboratory frame
   */
  double hydrogen_cm_to_lab(double ang, const double E) {
    ang = M_PI - ang;
    double mp = 938.346;
    double p = sqrt(E * (E + 2 * mp));
    double u = p / (E + 2 * mp);
    double g = 1 / sqrt(1 - u * u);
    double e = E + mp;
    double v_ratio = u * (e - u * p) / (p - u * e);
    double out;
    if (fabs(g * (cos(ang) + v_ratio)) == 0) {
      out = M_PI / 2;
    } else {
      out = atan(sin(ang) / (g * (cos(ang) + v_ratio)));
    }
    if (out < 0) {
      out += M_PI;
    }
    return out;
  }

  /** 
   * @brief Constructor: Retrieve scattering rate \sigma_e for large angle elastic scattering from data (for hydrogen)
   * 
   * The data in the file is expected in the following format:
   *  - First line contains the energy values.
   *  - For each energy value, there are two more lines in the file, 
   *    the first of which contains the scattering angles and the second 
   *    contains the scattering rate (up to a factor) corresponding to each angle
   *    (i.e. the upper limit of the integral in Eq. (12) in [1] is cos(angle)).
   *  - The angles are given in radians in decreasing order between pi and 0.
   * 
   * NOTE: The factor still missing to obtain the rate given in Eq. (12) is 10^{-24} * N_A * rho/A
   *  This factor is added in the 'rutherford_and_elastic_rate' method of the Material.
   *    
   * 
   * @param filename
   * @param cutoff
   * @param back_cutoff
   */
  CS_1d(const std::string filename, const double cuttoff,
        const double back_cuttoff)
      : energy(), rate() {
    std::ifstream file;
    file.open(filename);
    std::string line, token;

    // First line contains energy values
    getline(file, line);
    std::stringstream iss;
    iss << line;
    while (getline(iss, token, ' ')) {
      energy.push_back(atof(token.c_str()));
    }
    std::vector<double> tmp_vec;
    double lab_ang_cutoff, tmp_val, tmp_val_old = 0, top_rate, bottom_rate,
                                    total_rate, lin_inter_val = 0;
    int tmp_count, tmp_count_2, tmp_count_back, tmp_count_back_2,
        energy_index = 0;
    bool lin_inter_bool;
    while (getline(file, line)) {
      lab_ang_cutoff = hydrogen_cm_to_lab(back_cuttoff, energy[energy_index]);
      energy_index++;
      tmp_count = 0;
      tmp_count_back = 0;
      lin_inter_bool = true;
      std::stringstream iss2;
      iss2 << line;
      // Find angle closest to cutoff
      while (getline(iss2, token, ' ')) {
        // Read current angle
        tmp_val = atof(token.c_str());
        if (tmp_val > lab_ang_cutoff) {
          // Skip all angles larger than back_cutoff
          tmp_count++;
          tmp_count_back++;
        } else if (tmp_val > cuttoff) {
           // next angle if larger than cutoff
          tmp_count++;
          tmp_val_old = tmp_val;
        } else if (lin_inter_bool) {
          // if current angle is smaller than cutoff,
          // precompute factor for linear interpolation of rate
          lin_inter_val = (cuttoff - tmp_val_old) / (tmp_val - tmp_val_old);
          // Only compute lin_inter_val for first angle smaller than cutoff
          lin_inter_bool = false;
        }
      }

      // Find rate corresponding to angles closest to cutoff angle and
      // compute rate for cutoff angle by linear inerpolation
      tmp_count_2 = 0;
      tmp_count_back_2 = 0;
      lin_inter_bool = true;
      getline(file, line);
      tmp_vec.clear();
      std::stringstream iss3;
      iss3 << line;
      while (getline(iss3, token, ' ')) {
        tmp_val = atof(token.c_str());
        if (tmp_count_back_2 < tmp_count_back) {
          // Skip values corresponding to angles larger than back_cutoff
          tmp_count_back_2++;
          tmp_count_2++;
        } else if (tmp_count_2 < tmp_count) {
          // Update value of rate until the correct angle is reached
          tmp_count_2++;
          tmp_vec.push_back(tmp_val);
          tmp_val_old = tmp_val;
        } else if (lin_inter_bool) {
          // Linear inerpolation of rate
          tmp_vec.push_back(tmp_val * lin_inter_val +
                            (1 - lin_inter_val) * tmp_val_old);
          lin_inter_bool = false;
        }
      }
      // Compute total rate by subtracting the rate corresponding to back_cutoff from 
      // the rate corresponding to the cutoff angle to account for delta in lower limit of 
      //integral for hydrogen (see text below Eq. (12) in [1])
      top_rate = tmp_vec.back();
      bottom_rate = tmp_vec.front();
      total_rate = top_rate - bottom_rate;
      rate.push_back(total_rate);
    }
    file.close();
  }

  /**
   * Copy constructor
   */
  CS_1d(const CS_1d &other) : energy(other.energy), rate(other.rate) {}

  /**
   * Default constructor
   */
  CS_1d() : energy(), rate() {}

  /**
   * @brief Linear interpolation of rate from data
   * 
   * @param e energy
   * @return rate
   */
  double evaluate(const double e) const {
    double ret = 0;
    int r;
    double tol = 1e-7;
    if (energy.size() > 0) {
      if (e <= energy[0]) {
        ret = rate[0];
      } else if (e >= energy.back()) {
        ret = rate.back();
      } else {
        r = std::distance(energy.begin(),
                          std::lower_bound(energy.begin(), energy.end(), e));
        if (energy[r] - energy[r - 1] > tol) {
          ret =
              ((energy[r] - e) * rate[r - 1] + (e - energy[r - 1]) * rate[r]) /
              (energy[r] - energy[r - 1]);
        } else {
          ret = energy[r];
        }
      }
    }
    return ret;
  }

  std::vector<double> energy, rate;
};

struct CS_3d {
  /** 
   * @brief Constructor
   * 
   * The data in the file is expected in the following format:
   *  - First line contains the energy values.
   *  - For each energy value, there are three more lines in the file, 
   *    the first of which contains the exit energies, the second contains the corresponding CDF, and the 
   *    third contains the rvalues (pre-compound fraction see p. 136 in [4]) corresponding to each exit energy  
   * 
   * @param filename
   */
  CS_3d(const std::string filename) : energy(), exit_energy(), cdf(), rvalue() {
    std::ifstream file;
    file.open(filename);
    std::string line, token;

    // First line contains energy values
    getline(file, line);
    std::stringstream iss;
    iss << line;
    while (getline(iss, token, ' ')) {
      energy.push_back(atof(token.c_str()));
    }
    std::vector<double> tmp_vec;
    while (getline(file, line)) {
      tmp_vec.clear();
      std::stringstream iss3;
      iss3 << line;
      while (getline(iss3, token, ' ')) {
        tmp_vec.push_back(atof(token.c_str()));
      }
      // Store exit energies corresponding to current energy value
      exit_energy.push_back(tmp_vec);
      getline(file, line);
      tmp_vec.clear();
      std::stringstream iss4;
      iss4 << line;
      while (getline(iss4, token, ' ')) {
        tmp_vec.push_back(atof(token.c_str()));
      }
      // Store CDF corresponding to current energy value
      cdf.push_back(tmp_vec);
      getline(file, line);
      tmp_vec.clear();
      std::stringstream iss5;
      iss5 << line;
      while (getline(iss5, token, ' ')) {
        tmp_vec.push_back(atof(token.c_str()));
      }
      // Store rvalues (pre-compound fraction see p. 136 in [4]) corresponding to current energy value
      rvalue.push_back(tmp_vec);
    }
    file.close();
  }

  /**
   * Copy constructor
   */
  CS_3d(const CS_3d &other)
      : energy(other.energy), exit_energy(other.exit_energy), cdf(other.cdf),
        rvalue(other.rvalue) {}

  /**
   * Default constructor
   */  
  CS_3d() : energy(), exit_energy(), cdf(), rvalue() {}

  /** 
   * Sample from the cross section data for a given energy index
   * @param energy_index The index of the energy value
   * @param u The random number between 0 and 1
   * @param out_energy_cm The sampled exit energy in center-of-mass frame
   * @param out_rvalue The sampled rvalue (pre-compound fraction)
   */
  void sample_from_energy_index(const double energy_index, const double u,
                                double &out_energy_cm,
                                double &out_rvalue) const {
    double diff = 0;
    double tol = 1e-7;
    // Find cdf index for first cdf value larger than u
    int density_index =
        std::distance(cdf[energy_index].begin(),
                      std::lower_bound(cdf[energy_index].begin(),
                                       cdf[energy_index].end(), u));

    // Use enpoints outside data range
    if (density_index == 0) {
      out_energy_cm = exit_energy[energy_index][0];
      out_rvalue = rvalue[energy_index][0];
    } else if (density_index == int(cdf[energy_index].size())) {
      out_energy_cm = exit_energy[energy_index].back();
      out_rvalue = rvalue[energy_index].back();
    } else {
      if (cdf[energy_index][density_index] -
              cdf[energy_index][density_index - 1] >
          tol) {
        // Linear interpolation of exit energy and rvalue
        diff = (u - cdf[energy_index][density_index - 1]) /
               (cdf[energy_index][density_index] -
                cdf[energy_index][density_index - 1]);
        out_energy_cm =
            exit_energy[energy_index][density_index] * diff +
            exit_energy[energy_index][density_index - 1] * (1 - diff);
        out_rvalue = rvalue[energy_index][density_index] * diff +
                     rvalue[energy_index][density_index - 1] * (1 - diff);
      } else {
        out_energy_cm = exit_energy[energy_index][density_index];
        out_rvalue = rvalue[energy_index][density_index];
      }
    }
    return;
  }

  /** 
   * Sample from the cross section data for a given energy
   * @param e The energy value
   * @param r The sampled rvalue (pre-compound fraction)
   * @param out_e_cm The sampled exit energy in center-of-mass frame
   * @param gen The random number generator
   */
  void sample(const double e, double &r, double &out_e_cm, gsl_rng *gen) const {
    // Get index of energy data closest to current energy
    int energy_index = std::distance(
        energy.begin(), std::lower_bound(energy.begin(), energy.end(), e));
    double u = gsl_rng_uniform(gen);
    double out_energy_cm;
    double out_energy_cm_2;
    double out_rvalue;
    double out_rvalue_2;
    double diff;
    double tol = 1e-7;
    // Use enpoints outside data range
    if (energy_index == 0) {
      sample_from_energy_index(0, u, out_energy_cm, out_rvalue);
    } else if (energy_index == int(energy.size())) {
      sample_from_energy_index(energy_index - 1, u, out_energy_cm, out_rvalue);
    } else {
      // Linear interpolation of exit energy and rvalue between two energy indices
      sample_from_energy_index(energy_index, u, out_energy_cm, out_rvalue);
      if (energy[energy_index] - energy[energy_index - 1] > tol) {
        sample_from_energy_index(energy_index - 1, u, out_energy_cm_2,
                                 out_rvalue_2);
        diff = (e - energy[energy_index - 1]) /
               (energy[energy_index] - energy[energy_index - 1]);
        out_energy_cm = out_energy_cm * diff + out_energy_cm_2 * (1 - diff);
        out_rvalue = out_rvalue * diff + out_rvalue_2 * (1 - diff);
      }
    }
    r = out_rvalue;
    out_e_cm = out_energy_cm;
    return;
  }

  std::vector<double> energy;
  std::vector<std::vector<double>> exit_energy, cdf, rvalue;
};
/**
 * Struct for sampling exit angles from data for large angle elastic scattering
 * 
 * CDF given in Eq. (9) in [1] 
 * (with lower limit of integral -1+delta for hydrogen, see text below Eq. (12))
 * 
 * References:
 * [1] https://doi.org/10.1088/1361-6560/ae5586
 */
struct CS_2d {
  /**
   * Constructor
   * 
   * Reading in data from file
   * @param filename
   * @param cuttoff
   */
  CS_2d(const std::string filename, const double cuttoff)
      : energy(), exit_angle(), cdf() {
    std::ifstream file;
    file.open(filename);
    std::string line, token;

    // First line contains energy values
    getline(file, line);
    std::stringstream iss;
    iss << line;
    while (getline(iss, token, ' ')) {
      energy.push_back(atof(token.c_str()));
    }
    std::vector<double> tmp_vec;
    double tmp_val, tmp_val_old = 0, total_rate, lin_inter_val = 0;
    int tmp_count, tmp_count_2;
    bool lin_inter_bool;
    while (getline(file, line)) {
      tmp_vec.clear();
      tmp_count = 0;
      lin_inter_bool = true;
      std::stringstream iss2;
      iss2 << line;
      while (getline(iss2, token, ' ')) {
        tmp_val = atof(token.c_str());
        if (tmp_val > cuttoff) {
          // Get all angles larger than cutoff
          tmp_count++;
          tmp_vec.push_back(tmp_val);
          tmp_val_old = tmp_val;
        } else if (lin_inter_bool) {
          // precompute factor for linear interpolation of CDF
          lin_inter_val = (cuttoff - tmp_val_old) / (tmp_val - tmp_val_old);
          // Only compute lin_inter_val for first angle smaller than cutoff
          lin_inter_bool = false;
          // Add cutoff to angles
          tmp_vec.push_back(cuttoff);
        }
      }
      exit_angle.push_back(tmp_vec);
      tmp_count_2 = 0;
      lin_inter_bool = true;
      getline(file, line);
      tmp_vec.clear();
      std::stringstream iss3;
      iss3 << line;
      while (getline(iss3, token, ' ')) {
        tmp_val = atof(token.c_str());
        if (tmp_count_2 < tmp_count) {
          // Get all numerators in Eq. (9) in [1] corresponding to angles larger than cutoff
          tmp_count_2++;
          tmp_vec.push_back(tmp_val);
          tmp_val_old = tmp_val;
        } else if (lin_inter_bool) {
          // Linear interpolation for cutoff angle
          tmp_vec.push_back(tmp_val * lin_inter_val +
                            (1 - lin_inter_val) * tmp_val_old);
          lin_inter_bool = false;
        }
      }
      // Normalise \Pi_e: total_rate = denominator in Eq.(9) in [1]
      total_rate = tmp_vec.back();
      for (double &i : tmp_vec) {
        i /= total_rate;
      }
      cdf.push_back(tmp_vec);
    }
    file.close();
  }
  /**
   * Convert angle from centre-of-mass frame to laboratory frame for hydrogen
   * 
   * NOTE: This function is repeated as a method in the CS_1d struct.
   * @param ang angle in centre-of-mass frame
   * @param E energy
   * @return angle in laboratory frame
   */
  double hydrogen_cm_to_lab(double ang, const double E) {
    ang = M_PI - ang;
    double mp = 938.346;
    double p = sqrt(E * (E + 2 * mp));
    double u = p / (E + 2 * mp);
    double g = 1 / sqrt(1 - u * u);
    double e = E + mp;
    double v_ratio = u * (e - u * p) / (p - u * e);
    double out;
    if (fabs(g * (cos(ang) + v_ratio)) == 0) {
      out = M_PI / 2;
    } else {
      out = atan(sin(ang) / (g * (cos(ang) + v_ratio)));
    }
    if (out < 0) {
      out += M_PI;
    }
    return out;
  }
  /**
   * Constructor
   * 
   * Reading in data from file for hydrogen
   * @param filename
   * @param cuttoff
   * @param back_cutoff
   */
  CS_2d(const std::string filename, const double cuttoff,
        const double back_cuttoff)
      : energy(), exit_angle(), cdf() {
    std::ifstream file;
    file.open(filename);
    std::string line, token;

    // First line contains energy values
    getline(file, line);
    std::stringstream iss;
    iss << line;
    while (getline(iss, token, ' ')) {
      energy.push_back(atof(token.c_str()));
    }
    std::vector<double> tmp_vec;
    double lab_ang_cutoff, tmp_val, tmp_val_old = 0, top_rate, bottom_rate,
                                    total_rate, lin_inter_val = 0;
    int tmp_count, tmp_count_2, tmp_count_back, tmp_count_back_2,
        energy_index = 0;
    bool lin_inter_bool;
    while (getline(file, line)) {
      // Compute lab angle corresponding to back_cutoff in CM frame
      lab_ang_cutoff = hydrogen_cm_to_lab(back_cuttoff, energy[energy_index]);
      energy_index++;
      tmp_vec.clear();
      tmp_count = 0;
      tmp_count_back = 0;
      lin_inter_bool = true;
      std::stringstream iss2;
      iss2 << line;
      while (getline(iss2, token, ' ')) {
        tmp_val = atof(token.c_str());
        if (tmp_val > lab_ang_cutoff) {
          // Skip all angles larger than back_cutoff
          tmp_count++;
          tmp_count_back++;
        } else if (tmp_val > cuttoff) {
          // Get all angles larger than cutoff and smaller than back_cutoff
          tmp_count++;
          tmp_vec.push_back(tmp_val);
          tmp_val_old = tmp_val;
        } else if (lin_inter_bool) {
          // precompute factor for linear interpolation of CDF
          lin_inter_val = (cuttoff - tmp_val_old) / (tmp_val - tmp_val_old);
          // Only compute lin_inter_val for first angle smaller than cutoff
          lin_inter_bool = false;
          // Add cutoff to angle
          tmp_vec.push_back(cuttoff);
        }
      }
      exit_angle.push_back(tmp_vec);
      tmp_count_2 = 0;
      tmp_count_back_2 = 0;
      lin_inter_bool = true;
      getline(file, line);
      tmp_vec.clear();
      std::stringstream iss3;
      iss3 << line;
      while (getline(iss3, token, ' ')) {
        tmp_val = atof(token.c_str());
        if (tmp_count_back_2 < tmp_count_back) {
          // Skip values corresponding to angles larger than back_cutoff
          tmp_count_back_2++;
          tmp_count_2++;
        } else if (tmp_count_2 < tmp_count) {
          // Get all numerators in Eq. (9) in [1] corresponding to angles larger than cutoff
          tmp_count_2++;
          tmp_vec.push_back(tmp_val);
          tmp_val_old = tmp_val;
        } else if (lin_inter_bool) {
          // Linear interpolation for cutoff angle
          tmp_vec.push_back(tmp_val * lin_inter_val +
                            (1 - lin_inter_val) * tmp_val_old);
          lin_inter_bool = false;
        }
      }
      // Normalise \Pi_e: total_rate = denominator in Eq.(9) in [1],
      // accounting for delta in lower limit of integral for hydrogen (see text below Eq. (12) in [1])
      top_rate = tmp_vec.back();
      bottom_rate = tmp_vec.front();
      total_rate = top_rate - bottom_rate;
      for (double &i : tmp_vec) {
        i = (i - bottom_rate) / total_rate;
      }
      cdf.push_back(tmp_vec);
    }
    file.close();
  }
  /**
   * Copy constructor
   */
  CS_2d(const CS_2d &other)
      : energy(other.energy), exit_angle(other.exit_angle), cdf(other.cdf) {}

  /**
   * Default constructor
   */
  CS_2d() : energy(), exit_angle(), cdf() {}

  /**
   * Sample outgoing scattering angles from data 
   * 
   * @param energy_index - index for a specific energy data point
   * @param u - random number for sampling
   */
  double sample_from_energy_index(const double energy_index,
                                  const double u) const {
    double ret = 0;
    double diff = 0;
    double tol = 1e-7;
    // Find cdf index for first cdf value larger than u
    int density_index =
        std::distance(cdf[energy_index].begin(),
                      std::lower_bound(cdf[energy_index].begin(),
                                       cdf[energy_index].end(), u));

    // Use enpoints if outside data range
    if (density_index == 0) {
      ret = exit_angle[energy_index][0];
    } else if (density_index == int(cdf[energy_index].size())) {
      ret = exit_angle[energy_index].back();
    } else {
      // Linear interpolation of angle from data
      if (cdf[energy_index][density_index] -
              cdf[energy_index][density_index - 1] >
          tol) {
        diff = (u - cdf[energy_index][density_index - 1]) /
               (cdf[energy_index][density_index] -
                cdf[energy_index][density_index - 1]);
        ret = exit_angle[energy_index][density_index] * diff +
              exit_angle[energy_index][density_index - 1] * (1 - diff);
      } else {
        ret = exit_angle[energy_index][density_index];
      }
    }
    return ret;
  }

  /**
   * Sample outgoing scattering angles from data 
   * 
   * @param e - energy
   * @param gen
   */
  double sample(const double e, gsl_rng *gen) const {
    // Find closest energy in data that is larger than e
    int energy_index = std::distance(
        energy.begin(), std::lower_bound(energy.begin(), energy.end(), e));
    double u = gsl_rng_uniform(gen);
    double out_angle_cm;
    double out_angle_cm_2;
    double diff;
    double tol = 1e-7;
    // If energy outside data range use closest point in data
    if (energy_index == 0) {
      out_angle_cm = sample_from_energy_index(0, u);
    } else if (energy_index == int(energy.size())) {
      out_angle_cm = sample_from_energy_index(energy_index - 1, u);
    } else {
      // Otherwise linear interpolation from sampled angles
      out_angle_cm = sample_from_energy_index(energy_index, u);
      if (energy[energy_index] - energy[energy_index - 1] > tol) {
        out_angle_cm_2 = sample_from_energy_index(energy_index - 1, u);
        diff = (e - energy[energy_index - 1]) /
               (energy[energy_index] - energy[energy_index - 1]);
        out_angle_cm = out_angle_cm * diff + out_angle_cm_2 * (1 - diff);
      }
    }
    return out_angle_cm;
  }

  std::vector<double> energy; // energy data
  std::vector<std::vector<double>> exit_angle, cdf; // angle data and corresponding cdf
};

#endif
