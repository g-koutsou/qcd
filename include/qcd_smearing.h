/* qcd_smearing.h
 *
 * header file for qcd_smearing.c
 *
 * Tomasz Korzec 2009
 *******************************************/
 
#ifndef H_QCD_SMEARING
#define H_QCD_SMEARING 1

#define qcd_SUM_UP_HOPP(t,a,b)  \
   ({ t[0]  += a[0] + b[0];  \
      t[1]  += a[1] + b[1];  \
      t[2]  += a[2] + b[2];  \
      t[3]  += a[3] + b[3];  \
      t[4]  += a[4] + b[4];  \
      t[5]  += a[5] + b[5];  \
      t[6]  += a[6] + b[6];  \
      t[7]  += a[7] + b[7];  \
      t[8]  += a[8] + b[8];  \
      t[9]  += a[9] + b[9];  \
      t[10] += a[10]+ b[10]; \
      t[11] += a[11]+ b[11]; \
      t[12] += a[12]+ b[12]; \
      t[13] += a[13]+ b[13]; \
      t[14] += a[14]+ b[14]; \
      t[15] += a[15]+ b[15]; \
      t[16] += a[16]+ b[16]; \
      t[17] += a[17]+ b[17]; \
      t[18] += a[18]+ b[18]; \
      t[19] += a[19]+ b[19]; \
      t[20] += a[20]+ b[20]; \
      t[21] += a[21]+ b[21]; \
      t[22] += a[22]+ b[22]; \
      t[23] += a[23]+ b[23]; \
   })      

/* prototypes */
int qcd_gaussIteration3d(qcd_vector *v, qcd_gaugeField *u, qcd_real_8 alpha, qcd_uint_4 t);
int qcd_gaussIteration3dAll(qcd_vector *v, qcd_gaugeField *u, qcd_real_8 alpha, qcd_uint_2 gaugeCom);
int qcd_apeSmear3d(qcd_gaugeField *apeu, qcd_gaugeField *u, qcd_real_8 alpha);

#endif 
