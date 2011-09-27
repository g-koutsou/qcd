/* qcd_gaugeFixing.h
 *
 * header file for qcd_gaugeFixing.c
 *
 * Tomasz Korzec 2009
 *******************************************/
 
#ifndef H_QCD_GAUGEFIXING
#define H_QCD_GAUGEFIXING 1

#define qcd_LANDAU_GAUGE_FIXING_PRECISION 1.0e-12
#define qcd_LANDAU_GAUGE_FIXING_MAX_ITER 10000

void qcd_gluonField(qcd_gaugeField *a, qcd_gaugeField *u);
int qcd_landauGauge(qcd_gaugeField *landauu, qcd_gaugeField *u, qcd_real_8 overparam);

#endif 
