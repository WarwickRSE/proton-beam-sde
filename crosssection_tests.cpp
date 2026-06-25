
#include "material.h"
#include "cross_sections.h"

int main(int argc, char ** argv){

  double a = 12.011;
  int z = 6;
   std::string name = "carbon";

  //double a =30.974;
  //int z = 15;
  //std::string name = "phosphorus";


  Atom myAtom(a, z, "./Splines/" + name + "_ne_rate.txt",
               "./Splines/" + name + "_el_ruth_cross_sec.txt",
               "./Splines/" + name + "_ne_energyangle_cdf.txt",
               0.04);

  /*
  std::cout<<myAtom.ne_rate.evaluate(100.0)<<std::endl;
  std::cout<<myAtom.ne_rate.evaluate(73.0)<<std::endl;
  std::cout<<myAtom.ne_rate.evaluate(5.3)<<std::endl;
  std::cout<<myAtom.ne_rate.evaluate(1.0)<<std::endl;
  std::cout<<myAtom.ne_rate.evaluate(160.0)<<std::endl;
*/
/*  std::cout<<myAtom.el_ruth_rate.evaluate(100.0)<<std::endl;
  std::cout<<myAtom.el_ruth_rate.evaluate(73.0)<<std::endl;
  std::cout<<myAtom.el_ruth_rate.evaluate(5.3)<<std::endl;
  std::cout<<myAtom.el_ruth_rate.evaluate(1.0)<<std::endl;
  std::cout<<myAtom.el_ruth_rate.evaluate(160.0)<<std::endl;
*/

  gsl_rng myGen;
  myGen.prime({0.1, 0.1, 0.5, 0.0, 1.0, 0.0, 1.0, 0.0});

  std::cout<<myAtom.el_ruth_angle_cdf.sample(100.0, &myGen)<<std::endl;
  std::cout<<myAtom.el_ruth_angle_cdf.sample(73.0, &myGen)<<std::endl;
  std::cout<<myAtom.el_ruth_angle_cdf.sample(5.3, &myGen)<<std::endl;

  std::cout<<myAtom.el_ruth_angle_cdf.sample(1.0, &myGen)<<std::endl;
  std::cout<<myAtom.el_ruth_angle_cdf.sample(1.0, &myGen)<<std::endl;

  std::cout<<myAtom.el_ruth_angle_cdf.sample(160.0, &myGen)<<std::endl;
  std::cout<<myAtom.el_ruth_angle_cdf.sample(160.0, &myGen)<<std::endl;

}
