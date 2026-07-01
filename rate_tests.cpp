
#include "material.h"
#include "cross_sections.h"

int main(int argc, char ** argv){

/*
  double a = 1.008;
  int z = 1;
  std::string name = "hydrogen";

  Atom myAtom(a, z,
               "./Splines/" + name + "_el_ruth_cross_sec.txt",
               0.04, 0.04);

  double a2 = 15.99;
  int z2 = 8;
  std::string name2 = "oxygen";
  Atom myAtom2(a2, z2, "./Splines/" + name2 + "_ne_rate.txt",
               "./Splines/" + name2 + "_el_ruth_cross_sec.txt",
               "./Splines/" + name2 + "_ne_energyangle_cdf.txt",
               0.04);
 
  std::vector<Atom> myAt({myAtom, myAtom2});
  Material myMat(myAt, {0, 1}, {0.67, 0.33}, 1.0, 1.0);
*/
  double a2 = 12.011;
  int z2 = 6;
  std::string name2 = "carbon";
  Atom myAtom2(a2, z2, "./Splines/" + name2 + "_ne_rate.txt",
               "./Splines/" + name2 + "_el_ruth_cross_sec.txt",
               "./Splines/" + name2 + "_ne_energyangle_cdf.txt",
               0.04);
 
  std::vector<Atom> myAt({myAtom2});
  Material myMat(myAt, {0}, {1.0}, 1.0, 1.0e6);

//Material(std::vector<Atom> &atoms, const std::vector<int> &id,
  //         const std::vector<double> &x0, const double d0, const double I0)

  gsl_rng myGen;
  myGen.prime({0.1, 0.1, 0.5, 0.0, 1.0, 0.0, 1.0, 0.1, 0.67});

  std::cout<<myMat.bethe_bloch(1.0)<<std::endl;
  std::cout<<myMat.rutherford_and_elastic_rate(1.0)<<std::endl;
  std::cout<<myMat.multiple_scattering_sd(1.0, 1.0)<<std::endl;
 
}
