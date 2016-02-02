/* qcd_observables.c
 * 
 * computes simple observables, e.g. plaquette
 *
 * Tomasz Korzec 2009
 ***************************************************/
 
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <qcd.h>
 
qcd_real_8
_qcd_calculatePlaquette(qcd_gaugeField *u)
{
  qcd_complex_16 plaq[3][3],tmp[3][3];
  qcd_real_8 meanplaq = 0;
  qcd_real_8 result = 0;
  qcd_uint_4 l;
  qcd_uint_2 mu,nu;
  
  qcd_communicateGaugeM(u);
  qcd_waitall(u->geo);
  
  for(l=0; l<u->geo->lV; l++)
    for(mu=0; mu<3; mu++)
      for(nu=mu+1; nu<4; nu++)
	{
	  qcd_MUL3x3(plaq,u->D[l][mu],u->D[u->geo->plus[l][mu]][nu]);
	  qcd_MULADJOINT3x3(tmp,plaq,u->D[u->geo->plus[l][nu]][mu]);
	  qcd_MULADJOINT3x3(plaq,tmp,u->D[l][nu]);
	  meanplaq += qcd_SU3TRACER(plaq);
	}
  MPI_Allreduce(&meanplaq, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return result/(u->geo->V * 3.0 * 6.0); 
}

static inline qcd_real_8
plaq_pp(qcd_gaugeField *u_00)
{
  qcd_gaugeField *u_p0 = malloc(sizeof(qcd_gaugeField));
  qcd_gaugeField *u_0p = malloc(sizeof(qcd_gaugeField));
  qcd_initGaugeField(u_p0, u_00->geo);
  qcd_initGaugeField(u_0p, u_00->geo);

  qcd_real_8 p = 0.;
  for(int mu=0; mu<3; mu++) {
    qcd_shiftGaugeM(u_p0, u_00, mu);    
    for(int nu=mu+1; nu<4; nu++) {
      qcd_shiftGaugeM(u_0p, u_00, nu);
      qcd_real_8 p_mu_nu = 0;
      for(int l=0; l<u_00->geo->lV; l++) {
	qcd_complex_16 aux0[3][3], aux1[3][3];
	qcd_MUL3x3(aux0, u_00->D[l][mu], u_p0->D[l][nu]);
	qcd_MULADJOINT3x3(aux1, aux0, u_0p->D[l][mu]);
	qcd_MULADJOINT3x3(aux0, aux1, u_00->D[l][nu]);
	p_mu_nu += qcd_SU3TRACER(aux0);
      }
      //printf("%d %d %f\n", mu, nu, p_mu_nu/(u_00->geo->V * 3.0));
      p += p_mu_nu;
    }
  }
  qcd_destroyGaugeField(u_0p);
  qcd_destroyGaugeField(u_p0);

  free(u_p0);
  free(u_0p);
  
  return p;
}

static inline qcd_real_8
plaq_pm(qcd_gaugeField *u_00)
{
  qcd_gaugeField *u_0m = malloc(sizeof(qcd_gaugeField));
  qcd_gaugeField *u_pm = malloc(sizeof(qcd_gaugeField));
  qcd_initGaugeField(u_0m, u_00->geo);
  qcd_initGaugeField(u_pm, u_00->geo);

  qcd_real_8 p = 0.;
  for(int mu=0; mu<3; mu++) {
    for(int nu=mu+1; nu<4; nu++) {
      qcd_shiftGaugeP(u_0m, u_00, nu);
      qcd_shiftGaugeM(u_pm, u_0m, mu);          
      for(int l=0; l<u_00->geo->lV; l++) {
	qcd_complex_16 aux0[3][3], aux1[3][3];
	qcd_MULADJOINT3x3(aux0, u_00->D[l][mu], u_pm->D[l][nu]);
	qcd_MULADJOINT3x3(aux1, aux0, u_0m->D[l][mu]);
	qcd_MUL3x3(aux0, aux1, u_0m->D[l][nu]);
	p += qcd_SU3TRACER(aux0);
      }
    }
  }
  qcd_destroyGaugeField(u_pm);
  qcd_destroyGaugeField(u_0m);

  free(u_0m);
  free(u_pm);
  
  return p;
}

static inline qcd_real_8
plaq_mp(qcd_gaugeField *u_00)
{
  qcd_gaugeField *u_m0 = malloc(sizeof(qcd_gaugeField));
  qcd_gaugeField *u_mp = malloc(sizeof(qcd_gaugeField));
  qcd_initGaugeField(u_m0, u_00->geo);
  qcd_initGaugeField(u_mp, u_00->geo);

  qcd_real_8 p = 0.;
  for(int mu=0; mu<3; mu++) {
    qcd_shiftGaugeP(u_m0, u_00, mu);    
    for(int nu=mu+1; nu<4; nu++) {
      qcd_shiftGaugeM(u_mp, u_m0, nu);          
      for(int l=0; l<u_00->geo->lV; l++) {
	qcd_complex_16 aux0[3][3], aux1[3][3];
	qcd_MULADJOINT3x3(aux0, u_00->D[l][nu], u_mp->D[l][mu]);
	qcd_MULADJOINT3x3(aux1, aux0, u_m0->D[l][nu]);
	qcd_MUL3x3(aux0, aux1, u_m0->D[l][mu]);
	p += qcd_SU3TRACER(aux0);
      }
    }
  }
  qcd_destroyGaugeField(u_m0);
  qcd_destroyGaugeField(u_mp);

  free(u_m0);
  free(u_mp);
  
  return p;
}

static inline qcd_real_8
plaq_mm(qcd_gaugeField *u_00)
{
  qcd_gaugeField *u_m0 = malloc(sizeof(qcd_gaugeField));
  qcd_gaugeField *u_0m = malloc(sizeof(qcd_gaugeField));
  qcd_gaugeField *u_mm = malloc(sizeof(qcd_gaugeField));
  qcd_initGaugeField(u_m0, u_00->geo);
  qcd_initGaugeField(u_0m, u_00->geo);
  qcd_initGaugeField(u_mm, u_00->geo);

  qcd_real_8 p = 0.;
  for(int mu=0; mu<3; mu++) {
    qcd_shiftGaugeP(u_m0, u_00, mu);    
    for(int nu=mu+1; nu<4; nu++) {
      qcd_shiftGaugeP(u_mm, u_m0, nu);          
      qcd_shiftGaugeP(u_0m, u_00, nu);
      for(int l=0; l<u_00->geo->lV; l++) {
	qcd_complex_16 aux0[3][3], aux1[3][3];
	qcd_MULADJOINT3x3(aux0, u_0m->D[l][nu], u_m0->D[l][mu]);
	qcd_MULADJOINT3x3(aux1, aux0, u_mm->D[l][nu]);
	qcd_MUL3x3(aux0, aux1, u_mm->D[l][mu]);
	p += qcd_SU3TRACER(aux0);
      }
    }
  }
  qcd_destroyGaugeField(u_m0);
  qcd_destroyGaugeField(u_0m);
  qcd_destroyGaugeField(u_mm);

  free(u_m0);
  free(u_0m);
  free(u_mm);  
  return p;
}

qcd_real_8
qcd_calculatePlaquette(qcd_gaugeField *u_00)
{
  qcd_real_8 meanplaq = plaq_pp(u_00);
  qcd_real_8 result = 0.;
  MPI_Allreduce(&meanplaq, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  return result/(u_00->geo->V * 3.0 * 6.0); 
}
