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

 
qcd_real_8 qcd_calculatePlaquette(qcd_gaugeField *u)
{
   qcd_communicateGaugeM(u);
   qcd_waitall(u->geo);
   qcd_complex_16 plaq[3][3],tmp[3][3];
   qcd_real_8 meanplaq=0;
   qcd_real_8 result=0;
   qcd_uint_4 l;
   qcd_uint_2 mu,nu;
   
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
   return(result/(u->geo->V * 3.0 * 6.0)); // volume * N_color * number of mu-nu combinations
}//end qcd_calculatePlaquette 
