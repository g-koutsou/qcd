#ifndef _QCD_RNG_H
#define _QCD_RNG_H 1

void qcd_rng_init(unsigned int, int);
void qcd_rng_finalize();
unsigned long int qcd_get_rand_4();

#endif /* _QCD_RNG_H */
