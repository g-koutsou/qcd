/* qcd_communication.h
 *
 * header file for qcd_communication.c
 *
 * Tomasz Korzec 2009
 *******************************************/
 
#ifndef H_QCD_COMMUNICATION
#define H_QCD_COMMUNICATION 1 
 
/* prototypes */ 
int qcd_communicateVectorPMGaugeP(qcd_vector *v, qcd_gaugeField *u);
int qcd_communicateVectorPM(qcd_vector *v);
int qcd_communicateGaugeP(qcd_gaugeField *u);
int qcd_communicateGaugePInverse(qcd_gaugeField *u);
int qcd_communicateGaugeM(qcd_gaugeField *u);
int qcd_communicateGaugePM(qcd_gaugeField *u);
int qcd_communicateTransformationM(qcd_gaugeTransformation *u);
int qcd_communicatePropagatorP(qcd_propagator *p);
int qcd_communicatePropagatorPM(qcd_propagator *p);
void qcd_waitall(qcd_geometry *geo);

#endif
