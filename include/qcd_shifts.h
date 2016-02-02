/* qcd_shifts.h
 *
 * header file for qcd_shifts.c
 *
 *******************************************/
 
#ifndef H_QCD_SHIFTS
#define H_QCD_SHIFTS 1 
 
/* prototypes */ 
int qcd_shiftGaugeP(qcd_gaugeField *out, qcd_gaugeField *in, int dir);
int qcd_shiftGaugeM(qcd_gaugeField *out, qcd_gaugeField *in, int dir);

int qcd_shiftGaugeTransformationP(qcd_gaugeTransformation *out, qcd_gaugeTransformation *in, int dir);
int qcd_shiftGaugeTransformationM(qcd_gaugeTransformation *out, qcd_gaugeTransformation *in, int dir);

int qcd_shiftPropagatorP(qcd_propagator *out, qcd_propagator *in, int dir);
int qcd_shiftPropagatorM(qcd_propagator *out, qcd_propagator *in, int dir);

#endif
