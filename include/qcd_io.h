/* qcd_io.h
 *
 * header file for qcd_io.c
 *
 * Tomasz Korzec 2009
 *******************************************/
 
#ifndef H_QCD_IO
#define H_QCD_IO 1 

/* constants */
#define qcd_GF_LIME 0 
#define qcd_PROP_CMI 0
#define qcd_PROP_HMC 1
#define qcd_PROP_HMCV 2
#define qcd_PROP_LIME 3
#define qcd_SOURCE_LIME 4

 
/* prototypes */ 
int qcd_isBigEndian(void);
void qcd_swap_8(double *Rd, size_t N);
void qcd_swap_4(float *Rd, size_t N);
char* qcd_getParams(char* fname,int *len);
char* qcd_getParam(char token[],char* params,int len);
int qcd_getGaugeLime(char *fname, qcd_gaugeField *u);
int qcd_getGaugeField(char *fname, int type, qcd_gaugeField *u);
int qcd_getLimeMessage(char *fname, qcd_geometry *geo, char **message_ptr);
int qcd_writeGaugeLime(char *fname, qcd_gaugeField *u, char *message);
int qcd_writeGaugeField(char *fname, int type, qcd_gaugeField *u, ...);
int qcd_getPropagatorCMI(char *fname, qcd_propagator *p);
int qcd_getPropagatorHMC(char *fname, qcd_propagator *p);
int qcd_getPropagatorLime(char *fname, qcd_propagator *p);
int qcd_getPropagator(char *fname, int type, qcd_propagator *p);
int qcd_getVectorCMI(char *fname, qcd_uint_2 nu, qcd_uint_2 c2, qcd_vector *v);
int qcd_getVectorLime(char *fname, qcd_vector *v);
int qcd_getVectorHMC(char *fname, qcd_uint_2 nu, qcd_uint_2 c2, qcd_vector *v);
int qcd_getVectorHMCV(char *fname, qcd_uint_2 nu, qcd_uint_2 c2, qcd_vector *v);
int qcd_getVector(char *fname, int type, qcd_uint_2 nu, qcd_uint_2 c2, qcd_vector *v);
int qcd_writePropagatorCMI(char *fname, qcd_propagator *p);
int qcd_writePropagator(char *fname, int type, qcd_propagator *p);
int qcd_writePropagatorASCII(char *fname_basis, qcd_propagator *prop);
int qcd_writeVectorCMI(char *fname, qcd_uint_2 nu, qcd_uint_2 c2, qcd_vector *v);
int qcd_writeVectorLime(char *fname, int type, qcd_vector *v);
int qcd_writeVector(char *fname, int type, qcd_uint_2 nu, qcd_uint_2 c2, qcd_vector *v);
#endif
