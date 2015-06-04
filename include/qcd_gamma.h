/* qcd_gamma.h
 * collection of gamma matrices
 * ETMC conventions
 *
 * Tomasz Korzec 2009
 **************************************/

#ifndef H_QCD_GAMMA
#define H_QCD_GAMMA 1


extern const qcd_uint_2 qcd_EPS[6][3];
extern const qcd_int_2 qcd_SGN_EPS[6];
extern const qcd_complex_16 qcd_ONE[4][4];
extern const qcd_complex_16 qcd_GAMMA[8][4][4];
extern const qcd_complex_16 qcd_CGAMMA[6][4][4];
extern const qcd_complex_16 qcd_BAR_CGAMMA[6][4][4];
extern const qcd_complex_16 qcd_G5GAMMA[8][4][4];
extern const qcd_complex_16 qcd_BAR_G5GAMMA[8][4][4];
extern const qcd_complex_16 qcd_ONE_PLUS_GAMMA[6][4][4];
extern const qcd_complex_16 qcd_ONE_MINUS_GAMMA[6][4][4];
extern const qcd_complex_16 qcd_C[4][4];
extern const qcd_complex_16 qcd_BAR_C[4][4];

/* prototypes */
void qcd_gamma5Vector(qcd_vector *v);  
void qcd_gamma5Propagator(qcd_propagator *p);
void qcd_Propagatorgamma5(qcd_propagator *p);

#endif
