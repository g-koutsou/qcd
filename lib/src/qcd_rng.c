#include <gsl/gsl_rng.h>
#include <mpi.h>
#include <qcd.h>

void
qcd_rng_init(unsigned int seed, int myid)
{
  rng_state.gsl_rng_state = gsl_rng_alloc(gsl_rng_ranlux389);
  unsigned int process_seed = seed;

  srand(seed);
  for (int i=0; i<myid; i++)
    process_seed = rand();

  gsl_rng_set(rng_state.gsl_rng_state, process_seed);
  return;
}

void
qcd_rng_finalize()
{
  gsl_rng_free(rng_state.gsl_rng_state);
  return;
}

unsigned long int
qcd_get_rand_4()
{
  unsigned long int rand = gsl_rng_uniform_int(rng_state.gsl_rng_state, 4);
  return rand;
}

