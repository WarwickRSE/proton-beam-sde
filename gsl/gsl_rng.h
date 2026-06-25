#ifndef __gsl_mock__
#define __gsl_mock__

// Mock for rng (basic)
struct gsl_fake{
  double value;
};

using gsl_rng = gsl_fake;

double gsl_rng_uniform(gsl_fake * gen){return gen->value;};

#endif
