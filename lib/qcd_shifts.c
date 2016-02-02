/* qcd_shift.c
 * 
 * shift lattice fields in given direction
 * 
 **********************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <qcd.h>

int
qcd_shiftGaugeP(qcd_gaugeField *out, qcd_gaugeField *in, int dir)
{
  if(!out->initialized) {
    fprintf(stderr, "Error in %s, gauge fields must be properly initialized\n", __func__);
    exit(1);
  }
  
  qcd_communicateGaugePdir(in, dir);
  qcd_waitall(in->geo);
  
  for(int l=0; l<in->geo->lV; l++) {
    int l_minus_dir = in->geo->minus[l][dir];
    for(int mu=0; mu<4; mu++)
      for(int c0=0; c0<3; c0++)
	for(int c1=0; c1<3; c1++) {
	  out->D[l][mu][c0][c1].re = in->D[l_minus_dir][mu][c0][c1].re;
	  out->D[l][mu][c0][c1].im = in->D[l_minus_dir][mu][c0][c1].im;
	}
  }  
  return 0;
}

int
qcd_shiftGaugeM(qcd_gaugeField *out, qcd_gaugeField *in, int dir)
{
  if(!out->initialized) {
    fprintf(stderr, "Error in %s, gauge fields must be properly initialized\n", __func__);
    exit(1);
  }
  
  qcd_communicateGaugeMdir(in, dir);
  qcd_waitall(in->geo);
  
  for(int l=0; l<in->geo->lV; l++) {
    int l_plus_dir = in->geo->plus[l][dir];
    for(int mu=0; mu<4; mu++)
      for(int c0=0; c0<3; c0++)
	for(int c1=0; c1<3; c1++) {
	  out->D[l][mu][c0][c1].re = in->D[l_plus_dir][mu][c0][c1].re;
	  out->D[l][mu][c0][c1].im = in->D[l_plus_dir][mu][c0][c1].im;
	}
  }  
  return 0;
}

int
qcd_shiftGaugeTransformationP(qcd_gaugeTransformation *out, qcd_gaugeTransformation *in, int dir)
{
  if(!out->initialized) {
    fprintf(stderr, "Error in %s, gauge fields must be properly initialized\n", __func__);
    exit(1);
  }
  
  qcd_communicateGaugeTransformationPdir(in, dir);
  qcd_waitall(in->geo);
  
  for(int l=0; l<in->geo->lV; l++) {
    int l_minus_dir = in->geo->minus[l][dir];
    for(int c0=0; c0<3; c0++)
      for(int c1=0; c1<3; c1++) {
	out->D[l][c0][c1].re = in->D[l_minus_dir][c0][c1].re;
	out->D[l][c0][c1].im = in->D[l_minus_dir][c0][c1].im;
      }
  }  
  return 0;
}

int
qcd_shiftGaugeTransformationM(qcd_gaugeTransformation *out, qcd_gaugeTransformation *in, int dir)
{
  if(!out->initialized) {
    fprintf(stderr, "Error in %s, gauge fields must be properly initialized\n", __func__);
    exit(1);
  }
  
  qcd_communicateGaugeTransformationMdir(in, dir);
  qcd_waitall(in->geo);
  
  for(int l=0; l<in->geo->lV; l++) {
    int l_plus_dir = in->geo->plus[l][dir];
    for(int c0=0; c0<3; c0++)
      for(int c1=0; c1<3; c1++) {
	out->D[l][c0][c1].re = in->D[l_plus_dir][c0][c1].re;
	out->D[l][c0][c1].im = in->D[l_plus_dir][c0][c1].im;
      }
  }  
  return 0;
}

int
qcd_shiftPropagatorP(qcd_propagator *out, qcd_propagator *in, int dir)
{
  if(!out->initialized) {
    fprintf(stderr, "Error in %s, propagator fields must be properly initialized\n", __func__);
    exit(1);
  }
  
  qcd_communicatePropagatorPdir(in, dir);
  qcd_waitall(in->geo);
  
  for(int l=0; l<in->geo->lV; l++) {
    int l_minus_dir = in->geo->minus[l][dir];
    for(int mu=0; mu<4; mu++)
      for(int nu=0; nu<4; nu++)
	for(int c0=0; c0<3; c0++)
	  for(int c1=0; c1<3; c1++) {
	    out->D[l][mu][nu][c0][c1].re = in->D[l_minus_dir][mu][nu][c0][c1].re;
	    out->D[l][mu][nu][c0][c1].im = in->D[l_minus_dir][mu][nu][c0][c1].im;
	  }
  }
  return 0;
}

int
qcd_shiftPropagatorM(qcd_propagator *out, qcd_propagator *in, int dir)
{
  if(!out->initialized) {
    fprintf(stderr, "Error in %s, propagator fields must be properly initialized\n", __func__);
    exit(1);
  }
  
  qcd_communicatePropagatorMdir(in, dir);
  qcd_waitall(in->geo);
  
  for(int l=0; l<in->geo->lV; l++) {
    int l_plus_dir = in->geo->plus[l][dir];
    for(int mu=0; mu<4; mu++)
      for(int nu=0; nu<4; nu++)
	for(int c0=0; c0<3; c0++)
	  for(int c1=0; c1<3; c1++) {
	    out->D[l][mu][nu][c0][c1].re = in->D[l_plus_dir][mu][nu][c0][c1].re;
	    out->D[l][mu][nu][c0][c1].im = in->D[l_plus_dir][mu][nu][c0][c1].im;
	  }
  }  
  return 0;
}
