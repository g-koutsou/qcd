/* qcd_smearing.c
 *
 * gauss-smearing of fermion fields
 *
 * Tomasz Korzec 2008
 **********************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <qcd.h>
 

 

/* perform 1 iteration of gaussian smearing 
 * on time-slice t with
 * parameter alpha.
 *
 * v-> (v+alpha*Hv) / (1+6*alpha)
 */
int qcd_gaussIteration3d(qcd_vector *v, qcd_gaugeField *u, qcd_real_8 alpha, qcd_uint_4 t)
{
   qcd_uint_8 i,j;
   qcd_uint_2 c1,mu,nu,b; 
   qcd_uint_2 bl1[4] = {2,2,1,1};
   qcd_uint_2 bl2[4] = {3,3,3,2};
   qcd_uint_4 x,y,z,b0,b1,b2,b3,tt=0;
   qcd_complex_16 tmp[3];
   qcd_real_8 nrm = 1.0/(1.0 + alpha*6.0);
   qcd_vector v2;
   qcd_real_8 *uu;
   qcd_real_8 *psi;
   qcd_real_8 upsi[24];
   qcd_real_8 udaggerpsi[24];
   qcd_real_8 *total;

   if(!v->initialized)
   {
      fprintf(stderr,"Error in qcd_gaussIteration3d! Vector mus be properly initialized\n");
      return(1);
   }
   if(qcd_initVector(&v2, v->geo))
   {
     fprintf(stderr,"process %i: Error in qcd_gaussIteration3d! Could not initialize vector", v->geo->myid);
     return(1);
   }   

   //start communication (4d comm is in principle too much, but hey...
   qcd_communicateVectorPMGaugeP(v,u);

   //smear inner points:   
   qcd_zeroVector(&v2);

   if(v->geo->lL[1]>2 && v->geo->lL[2]>2 && v->geo->lL[3]>2 
      && t>=v->geo->Pos[0]*v->geo->lL[0] && t<(v->geo->Pos[0]+1)*v->geo->lL[0])
   { 
      tt=t-v->geo->Pos[0]*v->geo->lL[0];
      for(z=1;z<v->geo->lL[3]-1;z++)
      for(y=1;y<v->geo->lL[2]-1;y++)
      for(x=1;x<v->geo->lL[1]-1;x++)
      {
         i=qcd_LEXIC(tt,x,y,z,v->geo->lL);
         for(nu=1;nu<4;nu++)
         {
             uu  = (qcd_real_8*) &(u->D[i][nu][0][0].re);
             psi = (qcd_real_8*) &(v->D[v->geo->plus[i][nu]][0][0].re);
             qcd_APPLY_U(uu,upsi,psi);

             uu  = (qcd_real_8*) &(u->D[v->geo->minus[i][nu]][nu][0][0].re);
             psi = (qcd_real_8*) &(v->D[v->geo->minus[i][nu]][0][0].re); 
             qcd_APPLY_U_DAGGER(uu,udaggerpsi,psi);

             total = (qcd_real_8*) &(v2.D[i][0][0].re);
             qcd_SUM_UP_HOPP(total,upsi,udaggerpsi);
         }
      }//end inner-point loop
   }//end inner-points condition

   qcd_waitall(v->geo);

   //now boundary points
   if(t>=v->geo->Pos[0]*v->geo->lL[0] && t<(v->geo->Pos[0]+1)*v->geo->lL[0])
   for(j=0; j<v->geo->edge0Points; j++)
   {
      i=v->geo->edge0[j]*v->geo->lL[0]+tt; // works only with present lexic
      for(nu=1;nu<4;nu++)
      {
         uu  = (qcd_real_8*) &(u->D[i][nu][0][0].re);
         psi = (qcd_real_8*) &(v->D[v->geo->plus[i][nu]][0][0].re);
         qcd_APPLY_U(uu,upsi,psi);

         uu  = (qcd_real_8*) &(u->D[v->geo->minus[i][nu]][nu][0][0].re);
         psi = (qcd_real_8*) &(v->D[v->geo->minus[i][nu]][0][0].re); 
         qcd_APPLY_U_DAGGER(uu,udaggerpsi,psi);

         total = (qcd_real_8*) &(v2.D[i][0][0].re);
         qcd_SUM_UP_HOPP(total,upsi,udaggerpsi);
      }           
   }//end boundaries loop                  
   
   qcd_scaleVector3d(&v2,alpha,t);
   qcd_addVector(v,v,&v2);
   qcd_scaleVector3d(v,nrm,t); 

   qcd_destroyVector(&v2);
   return 0;
} 


/* perform 1 iteration of gaussian smearing 
 * on all time-slices with
 * parameter alpha.
 *
 * v-> (v+alpha*Hv) / (1+6*alpha)
 */
int qcd_gaussIteration3dAll(qcd_vector *v, qcd_gaugeField *u, qcd_real_8 alpha, qcd_uint_2 gaugeCom)
{
   qcd_uint_8 i,j;
   qcd_uint_2 c1,mu,nu,b; 
   qcd_uint_2 bl1[4] = {2,2,1,1};
   qcd_uint_2 bl2[4] = {3,3,3,2};
   qcd_uint_4 x,y,z,b0,b1,b2,b3,t,tt;
   qcd_complex_16 tmp[3];
   qcd_real_8 nrm = 1.0/(1.0 + alpha*6.0);
   qcd_vector v2;
   qcd_real_8 *uu;
   qcd_real_8 *psi;
   qcd_real_8 upsi[24];
   qcd_real_8 udaggerpsi[24];
   qcd_real_8 *total;

   if(!v->initialized)
   {
      fprintf(stderr,"Error in qcd_gaussIteration3dAll! Vector mus be properly initialized\n");
      return(1);
   }
   if(qcd_initVector(&v2, v->geo))
   {
     fprintf(stderr,"process %i: Error in qcd_gaussIteration3dAll! Could not initialize vector", v->geo->myid);
     return(1);
   }   

   if(gaugeCom)
      qcd_communicateVectorPMGaugeP(v,u);
   else
      qcd_communicateVectorPM(v);

   //smear inner points:   
   qcd_zeroVector(&v2);

   if(v->geo->lL[1]>2 && v->geo->lL[2]>2 && v->geo->lL[3]>2)
   {
      for(z=1;z<v->geo->lL[3]-1;z++)
      for(y=1;y<v->geo->lL[2]-1;y++)
      for(x=1;x<v->geo->lL[1]-1;x++)
      for(tt=0;tt<v->geo->lL[0];tt++)
      {
         i=qcd_LEXIC(tt,x,y,z,v->geo->lL);
         for(nu=1;nu<4;nu++)
         {
             uu  = (qcd_real_8*) &(u->D[i][nu][0][0].re);
             psi = (qcd_real_8*) &(v->D[v->geo->plus[i][nu]][0][0].re);
             qcd_APPLY_U(uu,upsi,psi);

             uu  = (qcd_real_8*) &(u->D[v->geo->minus[i][nu]][nu][0][0].re);
             psi = (qcd_real_8*) &(v->D[v->geo->minus[i][nu]][0][0].re); 
             qcd_APPLY_U_DAGGER(uu,udaggerpsi,psi);

             total = (qcd_real_8*) &(v2.D[i][0][0].re);
             qcd_SUM_UP_HOPP(total,upsi,udaggerpsi);
         }
      }//end inner-point loop
   }//end inner-points condition
   qcd_waitall(v->geo);

   //now boundary points
   for(j=0; j<v->geo->edge0Points; j++)
   for(tt=0; tt<v->geo->lL[0]; tt++)
   {
      i=v->geo->edge0[j]*v->geo->lL[0]+tt; // works only with present lexic
      for(nu=1;nu<4;nu++)
      {
         uu  = (qcd_real_8*) &(u->D[i][nu][0][0].re);
         psi = (qcd_real_8*) &(v->D[v->geo->plus[i][nu]][0][0].re);
         qcd_APPLY_U(uu,upsi,psi);

         uu  = (qcd_real_8*) &(u->D[v->geo->minus[i][nu]][nu][0][0].re);
         psi = (qcd_real_8*) &(v->D[v->geo->minus[i][nu]][0][0].re); 
         qcd_APPLY_U_DAGGER(uu,udaggerpsi,psi);

         total = (qcd_real_8*) &(v2.D[i][0][0].re);
         qcd_SUM_UP_HOPP(total,upsi,udaggerpsi);
      }           
   }//end boundaries loop                  

   qcd_scaleVector(&v2,alpha);
   qcd_addVector(v,v,&v2);
   qcd_scaleVector(v,nrm); 

   qcd_destroyVector(&v2);
   return 0;
}//end qcd_gaussIteration3dAll




/* perform 1 iteration of 3d APE-smearing 
 * with parameter alpha.
 *
 * u -> SU3-projection( u + alpha * sum spatial staples)
 */
int qcd_apeSmear3d(qcd_gaugeField *apeu, qcd_gaugeField *u, qcd_real_8 alpha)
{
   qcd_propagator edge;
   qcd_complex_16 stapleForward[3][3];
   qcd_complex_16 stapleBackward[3][3];
   qcd_uint_2 mu,nu,i=0,c1,c2;
   qcd_uint_4 l;
   qcd_complex_16 tmp[3][3];

   qcd_initPropagator(&edge,u->geo); // store edges in a propagator-structure. 

   /* since staples need next-to-nearest neighbors like U(x+mu-nu), this is done in 2 steps
      a) communicate U & calculate edges U_mu(x)U_nu(x+mu)
      b) communicate edges and put them together to staples.
   */
   
   qcd_communicateGaugePM(u);
   qcd_zeroGaugeField(apeu);   
   qcd_waitall(u->geo);

   for(mu=1; mu<4; mu++)
   for(nu=1;nu<4; nu++)
   if(mu!=nu)   
   for(l=0;l<u->geo->lV;l++)
   {
      qcd_MUL3x3(edge.D[l][mu][nu], u->D[l][mu], u->D[u->geo->plus[l][mu]][nu]);
   }
   
   qcd_communicatePropagatorP(&edge);

   //the forward staple doesn't need the edges
   for(mu=1; mu<4; mu++)
   for(nu=1;nu<4; nu++)
   if(mu!=nu)
   for(l=0;l<u->geo->lV;l++)
   {
      qcd_MUL3x3(tmp, u->D[l][nu],u->D[u->geo->plus[l][nu]][mu]);
      qcd_MULADJOINT3x3(stapleForward, tmp, u->D[u->geo->plus[l][mu]][nu]);
      for(c1=0;c1<3;c1++)
      for(c2=0;c2<3;c2++)
         apeu->D[l][mu][c1][c2] = qcd_CADD(apeu->D[l][mu][c1][c2],stapleForward[c1][c2]);
   }   
   
   qcd_waitall(u->geo);
   
   for(mu=1; mu<4; mu++)
   for(nu=1;nu<4; nu++)
   if(mu!=nu)
   for(l=0;l<u->geo->lV;l++)
   {
      qcd_ADJOINTMUL3x3(stapleBackward, u->D[u->geo->minus[l][nu]][nu], edge.D[u->geo->minus[l][nu]][mu][nu]);
      for(c1=0;c1<3;c1++)
      for(c2=0;c2<3;c2++)
         apeu->D[l][mu][c1][c2] = qcd_CADD(apeu->D[l][mu][c1][c2],stapleBackward[c1][c2]);
   }
   
   
   qcd_scaleGaugeField(apeu,alpha);
   qcd_addGaugeField(apeu,u,apeu);
   
   qcd_projectSU33d(apeu);
   
   qcd_destroyPropagator(&edge);
   return(0);
}//end qcd_apeSmear3d
