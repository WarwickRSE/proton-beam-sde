#ifndef __gsl_mock__
#define __gsl_mock__

#include <vector>
#include <iostream>
// Mock for rng (basic)
struct gsl_fake{
  std::vector<double> values;
  size_t last = 0;

  void prime(std::vector<double> vals){values = vals; last = 0;}
  void clear(){values.clear();}
  double yield(){return values.at(last++);} //at will throw if out of values
};

using gsl_rng = gsl_fake;

double gsl_rng_uniform(gsl_fake * gen){return gen->yield();};

#endif
