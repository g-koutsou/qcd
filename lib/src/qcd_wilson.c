/* qcd_wilson.c
 * 
 * collection of Wilson Dirac operators
 *
 * Tomasz Korzec 2008
 **********************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <qcd.h>
 

   

/* apply (1+gamma_mu) on a 4x3 dirac-color-vector */
void qcd_onePlusGamma(qcd_complex_16 result[][3], qcd_complex_16 vec[][3], qcd_uint_2 mu)
{
   qcd_uint_2 alpha,beta,c;
   for(alpha=0; alpha<4; alpha++)
   for(c=0; c<3; c++)
   {
      result[alpha][c]=(qcd_complex_16) {0,0};
      for(beta=0; beta<4; beta++)
      {
         result[alpha][c] = qcd_CADD(result[alpha][c],qcd_CMUL(qcd_ONE_PLUS_GAMMA[mu][alpha][beta],vec[beta][c]));
      }
   }
}//end qcd_onePlusGamma

/* apply (1-gamma_mu) on a 4x3 dirac-color-vector */
void qcd_oneMinusGamma(qcd_complex_16 result[][3], qcd_complex_16 vec[][3], qcd_uint_2 mu)
{
   qcd_uint_2 alpha,beta,c;
   for(alpha=0; alpha<4; alpha++)
   for(c=0; c<3; c++)
   {
      result[alpha][c]=(qcd_complex_16) {0,0};
      for(beta=0; beta<4; beta++)
      {
         result[alpha][c] = qcd_CADD(result[alpha][c],qcd_CMUL(qcd_ONE_MINUS_GAMMA[mu][alpha][beta],vec[beta][c]));
      }
   }
}//end qcd_oneMinusGamma

/* apply gauge field on a 4x3 dirac-color-vector */
void qcd_applyU(qcd_complex_16 result[][3],qcd_complex_16 u[][3],qcd_complex_16 vec[][3])
{
   qcd_uint_2 alpha,c1,c2;
   for(alpha=0; alpha<4; alpha++)
   for(c1=0; c1<3; c1++)
   {
      result[alpha][c1] = (qcd_complex_16) {0,0};
      for(c2=0; c2<3; c2++)
      {
         result[alpha][c1] = qcd_CADD(result[alpha][c1],qcd_CMUL(u[c1][c2],vec[alpha][c2]));
      }
   }
}//end qcd_applyU   

/* apply gauge field on a 4x3 dirac-color-vector */
void qcd_applyUDagger(qcd_complex_16 result[][3], qcd_complex_16 u[][3],qcd_complex_16 vec[][3])
{
   qcd_uint_2 alpha,c1,c2;
   for(alpha=0; alpha<4; alpha++)
   for(c1=0; c1<3; c1++)
   {
      result[alpha][c1] = (qcd_complex_16) {0,0};
      for(c2=0; c2<3; c2++)
      {
         result[alpha][c1] = qcd_CADD(result[alpha][c1],qcd_CMUL(qcd_CONJ(u[c2][c1]),vec[alpha][c2]));
      }
   }
}//end qcd_applyUDagger
   
   
/* apply (1*(m+4) * i*gamma_5*tm) on a 4x3 dirac-color-vector */
void qcd_massTerms(qcd_complex_16 result[][3], qcd_complex_16 vec[][3], qcd_real_8 mp4, qcd_real_8 tm)
{
   qcd_uint_2 alpha, beta, c;
   for(alpha=0; alpha<4; alpha++)
   for(c=0; c<3; c++)
   {
      result[alpha][c] = qcd_CSCALE(vec[alpha][c],mp4);
      for(beta=0; beta<4; beta++)
      {
         result[alpha][c] = qcd_CADD(result[alpha][c],qcd_CMUL(qcd_CMUL(qcd_GAMMA[5][alpha][beta],vec[beta][c]),((qcd_complex_16) {0,tm})));
      }  
   }   
}//end qcd_massTerms
   
   
/* multiply 4x3 dirac-color-vector with complex phase phi */
void qcd_phaseMul(qcd_complex_16 vec[][3], qcd_complex_16 phi)
{
   qcd_uint_2 alpha,c;
   for(alpha=0; alpha<4; alpha++)
   for(c=0; c<3; c++)
   {
      vec[alpha][c] = qcd_CMUL(vec[alpha][c],phi);
   }
}//end qcd_phaseMul
   
   

/* apply the Wilson TM dirac operator on vector v 
 * u  gauge field
 * m  mass parameter
 * tm twisted mass parameter 
 *
 * slightly optimized version, but still 100% architecture independent
 */
int qcd_applyWilsonTMOp(qcd_vector *Dv, qcd_vector *v, qcd_gaugeField *u, qcd_real_8 m, qcd_real_8 tm)
{
   qcd_uint_2 b;

   qcd_uint_4 t,x,y,z;
   qcd_uint_8 l,j;
   qcd_uint_4 b0,b1,b2,b3;
   qcd_uint_2 bl0[4] = {1,0,0,0}; 
   qcd_uint_2 bl1[4] = {2,2,1,1};
   qcd_uint_2 bl2[4] = {3,3,3,2};
   qcd_uint_2 isPeriodic[4];
   
   qcd_real_8 epg0[12], emg0[12], epg1[12], emg1[12]; // vectors to store (1+-gamma_mu)psi
   qcd_real_8 epg2[12], emg2[12], epg3[12], emg3[12];
   qcd_real_8 uepg0[24], uemg0[24], uepg1[24], uemg1[24];
   qcd_real_8 uepg2[24], uemg2[24], uepg3[24], uemg3[24];
   qcd_real_8 diag[24];
   qcd_real_8 *psi;
   qcd_real_8 mplus4=4+m;
   qcd_real_8 *uu;
   qcd_real_8 tmp;
   
   qcd_complex_16 pphases[4];
   qcd_complex_16 mphases[4];
   
   if(!v->initialized)
   {
      fprintf(stderr,"Error in qcd_ApplyWilsonTMOp! Vector must be properly initialized\n");
      return(1);
   }
   if(!u->initialized)
   {
      fprintf(stderr,"Error in qcd_ApplyWilsonTMOp! Gauge field must be properly initialized\n");
      return(1);
   }
   if(v->geo->initialized!=1 ||
      u->geo->initialized!=1 ||
      v->geo->L[0] != u->geo->L[0] ||
      v->geo->L[1] != u->geo->L[1] ||
      v->geo->L[2] != u->geo->L[2] ||
      v->geo->L[3] != u->geo->L[3] ||
      v->geo->lL[0] != u->geo->lL[0] ||
      v->geo->lL[1] != u->geo->lL[1] ||
      v->geo->lL[2] != u->geo->lL[2] ||
      v->geo->lL[3] != u->geo->lL[3] ||
      v->geo->Pos[0] != u->geo->Pos[0] ||
      v->geo->Pos[1] != u->geo->Pos[1] ||
      v->geo->Pos[2] != u->geo->Pos[2] ||
      v->geo->Pos[3] != u->geo->Pos[3])
   {
      fprintf(stderr,"Error in qcd_ApplyWilsonTMOp! Geometry not initialized or missmatched\n");
      return(1);   
   }       
   
   //if(u->geo->myid==0) printf("Wilson operator: starting communication\n");
   qcd_communicateVectorPMGaugeP(v,u);   
      
      
   for(b=0; b<4; b++)
   {
      if(u->geo->theta[b] != 0)
      {
         pphases[b] = (qcd_complex_16) {cos(u->geo->theta[b]/((double) u->geo->L[b])),  sin(u->geo->theta[b]/((double)u->geo->L[b]))};
         mphases[b] = (qcd_complex_16) {cos(u->geo->theta[b]/((double) u->geo->L[b])), -sin(u->geo->theta[b]/((double)u->geo->L[b]))};
         isPeriodic[b] = 0;
      }
      else
      {
         pphases[b] = (qcd_complex_16) {1,0};
         mphases[b] = (qcd_complex_16) {1,0};
         isPeriodic[b] = 1;
      }   
   }
   
   //if(u->geo->myid==0) printf("Wilson operator: starting inner-points loop\n");
   //calculate interior points
   if((u->geo->lL[0] > 2) && (u->geo->lL[1] > 2) && (u->geo->lL[2] > 2) && (u->geo->lL[3] > 2))
   {
      for(z=1; z<u->geo->lL[3]-1; z++)
      for(y=1; y<u->geo->lL[2]-1; y++)
      for(x=1; x<u->geo->lL[1]-1; x++)
      for(t=1; t<u->geo->lL[0]-1; t++) 
      {
         l = qcd_LEXIC(t,x,y,z,u->geo->lL);         
         
         psi = (qcd_real_8*) &(v->D[v->geo->minus[l][0]][0][0].re);         
         qcd_ONE_PLUS_GAMMA0_HALFSPINOR(epg0,psi);
         psi = (qcd_real_8*) &(v->D[v->geo->minus[l][1]][0][0].re);         
         qcd_ONE_PLUS_GAMMA1_HALFSPINOR(epg1,psi);
         psi = (qcd_real_8*) &(v->D[v->geo->minus[l][2]][0][0].re);         
         qcd_ONE_PLUS_GAMMA2_HALFSPINOR(epg2,psi);
         psi = (qcd_real_8*) &(v->D[v->geo->minus[l][3]][0][0].re);         
         qcd_ONE_PLUS_GAMMA3_HALFSPINOR(epg3,psi);
         psi = (qcd_real_8*) &(v->D[v->geo->plus[l][0]][0][0].re);         
         qcd_ONE_MINUS_GAMMA0_HALFSPINOR(emg0,psi);
         psi = (qcd_real_8*) &(v->D[v->geo->plus[l][1]][0][0].re);         
         qcd_ONE_MINUS_GAMMA1_HALFSPINOR(emg1,psi);
         psi = (qcd_real_8*) &(v->D[v->geo->plus[l][2]][0][0].re);         
         qcd_ONE_MINUS_GAMMA2_HALFSPINOR(emg2,psi);
         psi = (qcd_real_8*) &(v->D[v->geo->plus[l][3]][0][0].re);         
         qcd_ONE_MINUS_GAMMA3_HALFSPINOR(emg3,psi); 

         uu = (qcd_real_8*) &(u->D[l][0][0][0].re);
         qcd_APPLY_U_HALFSPINOR(uu,uemg0,emg0);
         uu = (qcd_real_8*) &(u->D[l][1][0][0].re);
         qcd_APPLY_U_HALFSPINOR(uu,uemg1,emg1);
         uu = (qcd_real_8*) &(u->D[l][2][0][0].re);
         qcd_APPLY_U_HALFSPINOR(uu,uemg2,emg2);
         uu = (qcd_real_8*) &(u->D[l][3][0][0].re);
         qcd_APPLY_U_HALFSPINOR(uu,uemg3,emg3);

         uu = (qcd_real_8*) &(u->D[u->geo->minus[l][0]][0][0][0].re);
         qcd_APPLY_U_DAGGER_HALFSPINOR(uu,uepg0,epg0);
         uu = (qcd_real_8*) &(u->D[u->geo->minus[l][1]][1][0][0].re);
         qcd_APPLY_U_DAGGER_HALFSPINOR(uu,uepg1,epg1);
         uu = (qcd_real_8*) &(u->D[u->geo->minus[l][2]][2][0][0].re);
         qcd_APPLY_U_DAGGER_HALFSPINOR(uu,uepg2,epg2);
         uu = (qcd_real_8*) &(u->D[u->geo->minus[l][3]][3][0][0].re);
         qcd_APPLY_U_DAGGER_HALFSPINOR(uu,uepg3,epg3);                           
         
         psi = (qcd_real_8*) &(v->D[l][0][0].re);
         qcd_MASSTERMS(diag,psi,mplus4,tm);
         
         psi = (qcd_real_8*) &(Dv->D[l][0][0].re);
         if(!isPeriodic[0])
         {         
            qcd_PHASE_MUL_HALFSPINOR(uepg0,tmp,mphases[0].re,mphases[0].im);
            qcd_PHASE_MUL_HALFSPINOR(uemg0,tmp,pphases[0].re,pphases[0].im); 
         }
         if(!isPeriodic[1])
         {         
            qcd_PHASE_MUL_HALFSPINOR(uepg1,tmp,mphases[1].re,mphases[1].im);
            qcd_PHASE_MUL_HALFSPINOR(uemg1,tmp,pphases[1].re,pphases[1].im); 
         }
         if(!isPeriodic[2])
         {         
            qcd_PHASE_MUL_HALFSPINOR(uepg2,tmp,mphases[2].re,mphases[2].im);
            qcd_PHASE_MUL_HALFSPINOR(uemg2,tmp,pphases[2].re,pphases[2].im); 
         }
         if(!isPeriodic[3])
         {         
            qcd_PHASE_MUL_HALFSPINOR(uepg3,tmp,mphases[3].re,mphases[3].im);
            qcd_PHASE_MUL_HALFSPINOR(uemg3,tmp,pphases[3].re,pphases[3].im); 
         }
         qcd_SUM_OP_HALFSPINORS(psi,diag,uepg0,uemg0,uepg1,uemg1,uepg2,uemg2,uepg3,uemg3);
      }
   }
   
   //if(u->geo->myid==0) printf("Wilson operator: starting communication finalization\n");
   qcd_waitall(v->geo);
   
   
   //if(u->geo->myid==0) printf("Wilson operator: starting boundaries loop\n");
   //calculate boundary points   
   for(j=0; j<v->geo->edgePoints; j++)
   {
      l=v->geo->edge[j];  
      psi = (qcd_real_8*) &(v->D[v->geo->minus[l][0]][0][0].re);         
      qcd_ONE_PLUS_GAMMA0_HALFSPINOR(epg0,psi);
      psi = (qcd_real_8*) &(v->D[v->geo->minus[l][1]][0][0].re);         
      qcd_ONE_PLUS_GAMMA1_HALFSPINOR(epg1,psi);
      psi = (qcd_real_8*) &(v->D[v->geo->minus[l][2]][0][0].re);         
      qcd_ONE_PLUS_GAMMA2_HALFSPINOR(epg2,psi);
      psi = (qcd_real_8*) &(v->D[v->geo->minus[l][3]][0][0].re);         
      qcd_ONE_PLUS_GAMMA3_HALFSPINOR(epg3,psi);
      psi = (qcd_real_8*) &(v->D[v->geo->plus[l][0]][0][0].re);         
      qcd_ONE_MINUS_GAMMA0_HALFSPINOR(emg0,psi);
      psi = (qcd_real_8*) &(v->D[v->geo->plus[l][1]][0][0].re);         
      qcd_ONE_MINUS_GAMMA1_HALFSPINOR(emg1,psi);
      psi = (qcd_real_8*) &(v->D[v->geo->plus[l][2]][0][0].re);         
      qcd_ONE_MINUS_GAMMA2_HALFSPINOR(emg2,psi);
      psi = (qcd_real_8*) &(v->D[v->geo->plus[l][3]][0][0].re);         
      qcd_ONE_MINUS_GAMMA3_HALFSPINOR(emg3,psi); 

      uu = (qcd_real_8*) &(u->D[l][0][0][0].re);
      qcd_APPLY_U_HALFSPINOR(uu,uemg0,emg0);
      uu = (qcd_real_8*) &(u->D[l][1][0][0].re);
      qcd_APPLY_U_HALFSPINOR(uu,uemg1,emg1);
      uu = (qcd_real_8*) &(u->D[l][2][0][0].re);
      qcd_APPLY_U_HALFSPINOR(uu,uemg2,emg2);
      uu = (qcd_real_8*) &(u->D[l][3][0][0].re);
      qcd_APPLY_U_HALFSPINOR(uu,uemg3,emg3);

      uu = (qcd_real_8*) &(u->D[u->geo->minus[l][0]][0][0][0].re);
      qcd_APPLY_U_DAGGER_HALFSPINOR(uu,uepg0,epg0);
      uu = (qcd_real_8*) &(u->D[u->geo->minus[l][1]][1][0][0].re);
      qcd_APPLY_U_DAGGER_HALFSPINOR(uu,uepg1,epg1);
      uu = (qcd_real_8*) &(u->D[u->geo->minus[l][2]][2][0][0].re);
      qcd_APPLY_U_DAGGER_HALFSPINOR(uu,uepg2,epg2);
      uu = (qcd_real_8*) &(u->D[u->geo->minus[l][3]][3][0][0].re);
      qcd_APPLY_U_DAGGER_HALFSPINOR(uu,uepg3,epg3);                           

      psi = (qcd_real_8*) &(v->D[l][0][0].re);
      qcd_MASSTERMS(diag,psi,mplus4,tm);

      psi = (qcd_real_8*) &(Dv->D[l][0][0].re);
      if(!isPeriodic[0])
      {         
         qcd_PHASE_MUL_HALFSPINOR(uepg0,tmp,mphases[0].re,mphases[0].im);
         qcd_PHASE_MUL_HALFSPINOR(uemg0,tmp,pphases[0].re,pphases[0].im); 
      }
      if(!isPeriodic[1])
      {         
         qcd_PHASE_MUL_HALFSPINOR(uepg1,tmp,mphases[1].re,mphases[1].im);
         qcd_PHASE_MUL_HALFSPINOR(uemg1,tmp,pphases[1].re,pphases[1].im); 
      }
      if(!isPeriodic[2])
      {         
         qcd_PHASE_MUL_HALFSPINOR(uepg2,tmp,mphases[2].re,mphases[2].im);
         qcd_PHASE_MUL_HALFSPINOR(uemg2,tmp,pphases[2].re,pphases[2].im); 
      }
      if(!isPeriodic[3])
      {         
         qcd_PHASE_MUL_HALFSPINOR(uepg3,tmp,mphases[3].re,mphases[3].im);
         qcd_PHASE_MUL_HALFSPINOR(uemg3,tmp,pphases[3].re,pphases[3].im); 
      }
      qcd_SUM_OP_HALFSPINORS(psi,diag,uepg0,uemg0,uepg1,uemg1,uepg2,uemg2,uepg3,uemg3);
   
   }//end loop over boundary points
      
   
   return(0);
}//end applyWilsonTMOp



   
   
   
   
/* apply the Wilson TM dirac operator on vector v 
 * u  gauge field
 * m  mass parameter
 * tm twisted mass parameter 
 *
 * version without pre-processor scripts
 */
int qcd_applyWilsonTMOpNoScripts(qcd_vector *Dv, qcd_vector *v, qcd_gaugeField *u, qcd_real_8 m, qcd_real_8 tm)
{
   qcd_uint_2 b;

   qcd_uint_4 t,x,y,z;
   qcd_uint_8 l,j;
   qcd_uint_4 b0,b1,b2,b3;
   qcd_uint_2 bl0[4] = {1,0,0,0}; 
   qcd_uint_2 bl1[4] = {2,2,1,1};
   qcd_uint_2 bl2[4] = {3,3,3,2};
   qcd_uint_2 isPeriodic[4],di,co;
   
   qcd_complex_16 epg0[4][3], emg0[4][3], epg1[4][3], emg1[4][3]; // vectors to store (1+-gamma_mu)psi
   qcd_complex_16 epg2[4][3], emg2[4][3], epg3[4][3], emg3[4][3];
   qcd_complex_16 uepg0[4][3], uemg0[4][3], uepg1[4][3], uemg1[4][3];
   qcd_complex_16 uepg2[4][3], uemg2[4][3], uepg3[4][3], uemg3[4][3];
   qcd_complex_16 diag[4][3];
   qcd_real_8 mplus4=4+m;
   qcd_real_8 *uu;
   qcd_real_8 tmp;
   
   qcd_complex_16 pphases[4];
   qcd_complex_16 mphases[4];
   
   if(!v->initialized)
   {
      fprintf(stderr,"Error in qcd_ApplyWilsonTMOp! Vector must be properly initialized\n");
      return(1);
   }
   if(!u->initialized)
   {
      fprintf(stderr,"Error in qcd_ApplyWilsonTMOp! Gauge field must be properly initialized\n");
      return(1);
   }
   if(v->geo->initialized!=1 ||
      u->geo->initialized!=1 ||
      v->geo->L[0] != u->geo->L[0] ||
      v->geo->L[1] != u->geo->L[1] ||
      v->geo->L[2] != u->geo->L[2] ||
      v->geo->L[3] != u->geo->L[3] ||
      v->geo->lL[0] != u->geo->lL[0] ||
      v->geo->lL[1] != u->geo->lL[1] ||
      v->geo->lL[2] != u->geo->lL[2] ||
      v->geo->lL[3] != u->geo->lL[3] ||
      v->geo->Pos[0] != u->geo->Pos[0] ||
      v->geo->Pos[1] != u->geo->Pos[1] ||
      v->geo->Pos[2] != u->geo->Pos[2] ||
      v->geo->Pos[3] != u->geo->Pos[3])
   {
      fprintf(stderr,"Error in qcd_ApplyWilsonTMOp! Geometry not initialized or missmatched\n");
      return(1);   
   }       
   
   //if(u->geo->myid==0) printf("Wilson operator: starting communication\n");
   qcd_communicateVectorPMGaugeP(v,u);   
      
      
   for(b=0; b<4; b++)
   {
      if(u->geo->theta[b] != 0)
      {
         pphases[b] = (qcd_complex_16) {cos(u->geo->theta[b]/((double) u->geo->L[b])),  sin(u->geo->theta[b]/((double)u->geo->L[b]))};
         mphases[b] = (qcd_complex_16) {cos(u->geo->theta[b]/((double) u->geo->L[b])), -sin(u->geo->theta[b]/((double)u->geo->L[b]))};
         isPeriodic[b] = 0;
      }
      else
      {
         pphases[b] = (qcd_complex_16) {1,0};
         mphases[b] = (qcd_complex_16) {1,0};
         isPeriodic[b] = 1;
      }   
   }
   
   //if(u->geo->myid==0) printf("Wilson operator: starting inner-points loop\n");
   //calculate interior points
   if((u->geo->lL[0] > 2) && (u->geo->lL[1] > 2) && (u->geo->lL[2] > 2) && (u->geo->lL[3] > 2))
   {
      for(z=1; z<u->geo->lL[3]-1; z++)
      for(y=1; y<u->geo->lL[2]-1; y++)
      for(x=1; x<u->geo->lL[1]-1; x++)
      for(t=1; t<u->geo->lL[0]-1; t++)
      {
         l = qcd_LEXIC(t,x,y,z,u->geo->lL);         
      
         qcd_onePlusGamma(epg0,v->D[v->geo->minus[l][0]],0);
         qcd_onePlusGamma(epg1,v->D[v->geo->minus[l][1]],1);
         qcd_onePlusGamma(epg2,v->D[v->geo->minus[l][2]],2);
         qcd_onePlusGamma(epg3,v->D[v->geo->minus[l][3]],3);
         qcd_oneMinusGamma(emg0,v->D[v->geo->plus[l][0]],0);
         qcd_oneMinusGamma(emg1,v->D[v->geo->plus[l][1]],1);
         qcd_oneMinusGamma(emg2,v->D[v->geo->plus[l][2]],2);
         qcd_oneMinusGamma(emg3,v->D[v->geo->plus[l][3]],3);

         qcd_applyU(uemg0,u->D[l][0],emg0);
         qcd_applyU(uemg1,u->D[l][1],emg1);
         qcd_applyU(uemg2,u->D[l][2],emg2);
         qcd_applyU(uemg3,u->D[l][3],emg3);

         qcd_applyUDagger(uepg0,u->D[u->geo->minus[l][0]][0],epg0);
         qcd_applyUDagger(uepg1,u->D[u->geo->minus[l][1]][1],epg1);
         qcd_applyUDagger(uepg2,u->D[u->geo->minus[l][2]][2],epg2);
         qcd_applyUDagger(uepg3,u->D[u->geo->minus[l][3]][3],epg3);
         
         qcd_massTerms(diag,v->D[l],mplus4,tm);
                 
         if(!isPeriodic[0])
         {         
            qcd_phaseMul(uepg0,mphases[0]);
            qcd_phaseMul(uemg0,pphases[0]); 
         }
         if(!isPeriodic[1])
         {         
            qcd_phaseMul(uepg1,mphases[1]);
            qcd_phaseMul(uemg1,pphases[1]); 
         }
         if(!isPeriodic[2])
         {         
            qcd_phaseMul(uepg2,mphases[2]);
            qcd_phaseMul(uemg2,pphases[2]); 
         }
         if(!isPeriodic[3])
         {         
            qcd_phaseMul(uepg3,mphases[3]);
            qcd_phaseMul(uemg3,pphases[3]); 
         }
         for(di=0; di<4; di++)
         for(co=0; co<3; co++)
         {
            Dv->D[l][di][co] = qcd_CSUB(diag[di][co], qcd_CSCALE(qcd_CADD(uepg0[di][co],
                                                                 qcd_CADD(uemg0[di][co],
                                                                 qcd_CADD(uepg1[di][co],
                                                                 qcd_CADD(uemg1[di][co],
                                                                 qcd_CADD(uepg2[di][co],
                                                                 qcd_CADD(uemg2[di][co],
                                                                 qcd_CADD(uepg3[di][co],uemg3[di][co]))))))),0.5));
         }   
      }
   }
   
   //if(u->geo->myid==0) printf("Wilson operator: starting communication finalization\n");
   qcd_waitall(v->geo);
   
   
   //if(u->geo->myid==0) printf("Wilson operator: starting boundaries loop\n");
   //calculate boundary points
   for(j=0; j<v->geo->edgePoints; j++)
   {
      l=v->geo->edge[j];  
      qcd_onePlusGamma(epg0,v->D[v->geo->minus[l][0]],0);
      qcd_onePlusGamma(epg1,v->D[v->geo->minus[l][1]],1);
      qcd_onePlusGamma(epg2,v->D[v->geo->minus[l][2]],2);
      qcd_onePlusGamma(epg3,v->D[v->geo->minus[l][3]],3);
      qcd_oneMinusGamma(emg0,v->D[v->geo->plus[l][0]],0);
      qcd_oneMinusGamma(emg1,v->D[v->geo->plus[l][1]],1);
      qcd_oneMinusGamma(emg2,v->D[v->geo->plus[l][2]],2);
      qcd_oneMinusGamma(emg3,v->D[v->geo->plus[l][3]],3);

      qcd_applyU(uemg0,u->D[l][0],emg0);
      qcd_applyU(uemg1,u->D[l][1],emg1);
      qcd_applyU(uemg2,u->D[l][2],emg2);
      qcd_applyU(uemg3,u->D[l][3],emg3);

      qcd_applyUDagger(uepg0,u->D[u->geo->minus[l][0]][0],epg0);
      qcd_applyUDagger(uepg1,u->D[u->geo->minus[l][1]][1],epg1);
      qcd_applyUDagger(uepg2,u->D[u->geo->minus[l][2]][2],epg2);
      qcd_applyUDagger(uepg3,u->D[u->geo->minus[l][3]][3],epg3);

      qcd_massTerms(diag,v->D[l],mplus4,tm);

      if(!isPeriodic[0])
      {         
         qcd_phaseMul(uepg0,mphases[0]);
         qcd_phaseMul(uemg0,pphases[0]); 
      }
      if(!isPeriodic[1])
      {         
         qcd_phaseMul(uepg1,mphases[1]);
         qcd_phaseMul(uemg1,pphases[1]); 
      }
      if(!isPeriodic[2])
      {         
         qcd_phaseMul(uepg2,mphases[2]);
         qcd_phaseMul(uemg2,pphases[2]); 
      }
      if(!isPeriodic[3])
      {         
         qcd_phaseMul(uepg3,mphases[3]);
         qcd_phaseMul(uemg3,pphases[3]); 
      }
      for(di=0; di<4; di++)
      for(co=0; co<3; co++)
      {
         Dv->D[l][di][co] = qcd_CSUB(diag[di][co], qcd_CSCALE(qcd_CADD(uepg0[di][co],
                                                              qcd_CADD(uemg0[di][co],
                                                              qcd_CADD(uepg1[di][co],
                                                              qcd_CADD(uemg1[di][co],
                                                              qcd_CADD(uepg2[di][co],
                                                              qcd_CADD(uemg2[di][co],
                                                              qcd_CADD(uepg3[di][co],uemg3[di][co]))))))),0.5));
      }
   }//end loop over boundary points
   
   
   
   return(0);
}//end applyWilsonTMOpNoScripts 








/* apply the Wilson TM dirac operator on vector v 
 * u  gauge field
 * m  mass parameter
 * tm twisted mass parameter 
 */
int qcd_applyWilsonTMOpNoOptim(qcd_vector *Dv, qcd_vector *v, qcd_gaugeField *u, qcd_real_8 m, qcd_real_8 tm)
{
   qcd_uint_2 b;

   qcd_uint_4 t,x,y,z;
   qcd_uint_8 l,j;
   qcd_uint_4 b0,b1,b2,b3;
   qcd_uint_2 bl0[4] = {1,0,0,0}; 
   qcd_uint_2 bl1[4] = {2,2,1,1};
   qcd_uint_2 bl2[4] = {3,3,3,2};
   qcd_uint_2 isPeriodic[4];
   
   qcd_real_8 epg0[24], emg0[24], epg1[24], emg1[24]; // vectors to store (1+-gamma_mu)psi
   qcd_real_8 epg2[24], emg2[24], epg3[24], emg3[24];
   qcd_real_8 uepg0[24], uemg0[24], uepg1[24], uemg1[24];
   qcd_real_8 uepg2[24], uemg2[24], uepg3[24], uemg3[24];
   qcd_real_8 diag[24];
   qcd_real_8 *psi;
   qcd_real_8 mplus4=4+m;
   qcd_real_8 *uu;
   qcd_real_8 tmp;
   
   qcd_complex_16 pphases[4];
   qcd_complex_16 mphases[4];
   
   if(!v->initialized)
   {
      fprintf(stderr,"Error in qcd_ApplyWilsonTMOp! Vector must be properly initialized\n");
      return(1);
   }
   if(!u->initialized)
   {
      fprintf(stderr,"Error in qcd_ApplyWilsonTMOp! Gauge field must be properly initialized\n");
      return(1);
   }
   if(v->geo->initialized!=1 ||
      u->geo->initialized!=1 ||
      v->geo->L[0] != u->geo->L[0] ||
      v->geo->L[1] != u->geo->L[1] ||
      v->geo->L[2] != u->geo->L[2] ||
      v->geo->L[3] != u->geo->L[3] ||
      v->geo->lL[0] != u->geo->lL[0] ||
      v->geo->lL[1] != u->geo->lL[1] ||
      v->geo->lL[2] != u->geo->lL[2] ||
      v->geo->lL[3] != u->geo->lL[3] ||
      v->geo->Pos[0] != u->geo->Pos[0] ||
      v->geo->Pos[1] != u->geo->Pos[1] ||
      v->geo->Pos[2] != u->geo->Pos[2] ||
      v->geo->Pos[3] != u->geo->Pos[3])
   {
      fprintf(stderr,"Error in qcd_ApplyWilsonTMOp! Geometry not initialized or missmatched\n");
      return(1);   
   }       
   
   //if(u->geo->myid==0) printf("Wilson operator: starting communication\n");
   qcd_communicateVectorPMGaugeP(v,u);   
      
      
   for(b=0; b<4; b++)
   {
      if(u->geo->theta[b] != 0)
      {
         pphases[b] = (qcd_complex_16) {cos(u->geo->theta[b]/((double) u->geo->L[b])),  sin(u->geo->theta[b]/((double)u->geo->L[b]))};
         mphases[b] = (qcd_complex_16) {cos(u->geo->theta[b]/((double) u->geo->L[b])), -sin(u->geo->theta[b]/((double)u->geo->L[b]))};
         isPeriodic[b] = 0;
      }
      else
      {
         pphases[b] = (qcd_complex_16) {1,0};
         mphases[b] = (qcd_complex_16) {1,0};
         isPeriodic[b] = 1;
      }   
   }
   
   //if(u->geo->myid==0) printf("Wilson operator: starting inner-points loop\n");
   //calculate interior points
   if((u->geo->lL[0] > 2) && (u->geo->lL[1] > 2) && (u->geo->lL[2] > 2) && (u->geo->lL[3] > 2))
   {
      for(z=1; z<u->geo->lL[3]-1; z++)
      for(y=1; y<u->geo->lL[2]-1; y++)
      for(x=1; x<u->geo->lL[1]-1; x++)
      for(t=1; t<u->geo->lL[0]-1; t++) 
      {
         l = qcd_LEXIC(t,x,y,z,u->geo->lL);         
         
         psi = (qcd_real_8*) &(v->D[v->geo->minus[l][0]][0][0].re);         
         qcd_ONE_PLUS_GAMMA0(epg0,psi);
         psi = (qcd_real_8*) &(v->D[v->geo->minus[l][1]][0][0].re);         
         qcd_ONE_PLUS_GAMMA1(epg1,psi);
         psi = (qcd_real_8*) &(v->D[v->geo->minus[l][2]][0][0].re);         
         qcd_ONE_PLUS_GAMMA2(epg2,psi);
         psi = (qcd_real_8*) &(v->D[v->geo->minus[l][3]][0][0].re);         
         qcd_ONE_PLUS_GAMMA3(epg3,psi);
         psi = (qcd_real_8*) &(v->D[v->geo->plus[l][0]][0][0].re);         
         qcd_ONE_MINUS_GAMMA0(emg0,psi);
         psi = (qcd_real_8*) &(v->D[v->geo->plus[l][1]][0][0].re);         
         qcd_ONE_MINUS_GAMMA1(emg1,psi);
         psi = (qcd_real_8*) &(v->D[v->geo->plus[l][2]][0][0].re);         
         qcd_ONE_MINUS_GAMMA2(emg2,psi);
         psi = (qcd_real_8*) &(v->D[v->geo->plus[l][3]][0][0].re);         
         qcd_ONE_MINUS_GAMMA3(emg3,psi); 

         uu = (qcd_real_8*) &(u->D[l][0][0][0].re);
         qcd_APPLY_U(uu,uemg0,emg0);
         uu = (qcd_real_8*) &(u->D[l][1][0][0].re);
         qcd_APPLY_U(uu,uemg1,emg1);
         uu = (qcd_real_8*) &(u->D[l][2][0][0].re);
         qcd_APPLY_U(uu,uemg2,emg2);
         uu = (qcd_real_8*) &(u->D[l][3][0][0].re);
         qcd_APPLY_U(uu,uemg3,emg3);

         uu = (qcd_real_8*) &(u->D[u->geo->minus[l][0]][0][0][0].re);
         qcd_APPLY_U_DAGGER(uu,uepg0,epg0);
         uu = (qcd_real_8*) &(u->D[u->geo->minus[l][1]][1][0][0].re);
         qcd_APPLY_U_DAGGER(uu,uepg1,epg1);
         uu = (qcd_real_8*) &(u->D[u->geo->minus[l][2]][2][0][0].re);
         qcd_APPLY_U_DAGGER(uu,uepg2,epg2);
         uu = (qcd_real_8*) &(u->D[u->geo->minus[l][3]][3][0][0].re);
         qcd_APPLY_U_DAGGER(uu,uepg3,epg3);                           
         
         psi = (qcd_real_8*) &(v->D[l][0][0].re);
         qcd_MASSTERMS(diag,psi,mplus4,tm);
         
         psi = (qcd_real_8*) &(Dv->D[l][0][0].re);
         if(!isPeriodic[0])
         {         
            qcd_PHASE_MUL(uepg0,tmp,mphases[0].re,mphases[0].im);
            qcd_PHASE_MUL(uemg0,tmp,pphases[0].re,pphases[0].im); 
         }
         if(!isPeriodic[1])
         {         
            qcd_PHASE_MUL(uepg1,tmp,mphases[1].re,mphases[1].im);
            qcd_PHASE_MUL(uemg1,tmp,pphases[1].re,pphases[1].im); 
         }
         if(!isPeriodic[2])
         {         
            qcd_PHASE_MUL(uepg2,tmp,mphases[2].re,mphases[2].im);
            qcd_PHASE_MUL(uemg2,tmp,pphases[2].re,pphases[2].im); 
         }
         if(!isPeriodic[3])
         {         
            qcd_PHASE_MUL(uepg3,tmp,mphases[3].re,mphases[3].im);
            qcd_PHASE_MUL(uemg3,tmp,pphases[3].re,pphases[3].im); 
         }
         qcd_SUM_OP(psi,diag,uepg0,uemg0,uepg1,uemg1,uepg2,uemg2,uepg3,uemg3);
      }
   }
   
   //if(u->geo->myid==0) printf("Wilson operator: starting communication finalization\n");
   qcd_waitall(v->geo);
   
   
   //if(u->geo->myid==0) printf("Wilson operator: starting boundaries loop\n");
   //calculate boundary points   
   for(j=0; j<v->geo->edgePoints; j++)
   {
      l=v->geo->edge[j];  
      psi = (qcd_real_8*) &(v->D[v->geo->minus[l][0]][0][0].re);         
      qcd_ONE_PLUS_GAMMA0(epg0,psi);
      psi = (qcd_real_8*) &(v->D[v->geo->minus[l][1]][0][0].re);         
      qcd_ONE_PLUS_GAMMA1(epg1,psi);
      psi = (qcd_real_8*) &(v->D[v->geo->minus[l][2]][0][0].re);         
      qcd_ONE_PLUS_GAMMA2(epg2,psi);
      psi = (qcd_real_8*) &(v->D[v->geo->minus[l][3]][0][0].re);         
      qcd_ONE_PLUS_GAMMA3(epg3,psi);
      psi = (qcd_real_8*) &(v->D[v->geo->plus[l][0]][0][0].re);         
      qcd_ONE_MINUS_GAMMA0(emg0,psi);
      psi = (qcd_real_8*) &(v->D[v->geo->plus[l][1]][0][0].re);         
      qcd_ONE_MINUS_GAMMA1(emg1,psi);
      psi = (qcd_real_8*) &(v->D[v->geo->plus[l][2]][0][0].re);         
      qcd_ONE_MINUS_GAMMA2(emg2,psi);
      psi = (qcd_real_8*) &(v->D[v->geo->plus[l][3]][0][0].re);         
      qcd_ONE_MINUS_GAMMA3(emg3,psi); 

      uu = (qcd_real_8*) &(u->D[l][0][0][0].re);
      qcd_APPLY_U(uu,uemg0,emg0);
      uu = (qcd_real_8*) &(u->D[l][1][0][0].re);
      qcd_APPLY_U(uu,uemg1,emg1);
      uu = (qcd_real_8*) &(u->D[l][2][0][0].re);
      qcd_APPLY_U(uu,uemg2,emg2);
      uu = (qcd_real_8*) &(u->D[l][3][0][0].re);
      qcd_APPLY_U(uu,uemg3,emg3);

      uu = (qcd_real_8*) &(u->D[u->geo->minus[l][0]][0][0][0].re);
      qcd_APPLY_U_DAGGER(uu,uepg0,epg0);
      uu = (qcd_real_8*) &(u->D[u->geo->minus[l][1]][1][0][0].re);
      qcd_APPLY_U_DAGGER(uu,uepg1,epg1);
      uu = (qcd_real_8*) &(u->D[u->geo->minus[l][2]][2][0][0].re);
      qcd_APPLY_U_DAGGER(uu,uepg2,epg2);
      uu = (qcd_real_8*) &(u->D[u->geo->minus[l][3]][3][0][0].re);
      qcd_APPLY_U_DAGGER(uu,uepg3,epg3);                           

      psi = (qcd_real_8*) &(v->D[l][0][0].re);
      qcd_MASSTERMS(diag,psi,mplus4,tm);

      psi = (qcd_real_8*) &(Dv->D[l][0][0].re);
      if(!isPeriodic[0])
      {         
         qcd_PHASE_MUL(uepg0,tmp,mphases[0].re,mphases[0].im);
         qcd_PHASE_MUL(uemg0,tmp,pphases[0].re,pphases[0].im); 
      }
      if(!isPeriodic[1])
      {         
         qcd_PHASE_MUL(uepg1,tmp,mphases[1].re,mphases[1].im);
         qcd_PHASE_MUL(uemg1,tmp,pphases[1].re,pphases[1].im); 
      }
      if(!isPeriodic[2])
      {         
         qcd_PHASE_MUL(uepg2,tmp,mphases[2].re,mphases[2].im);
         qcd_PHASE_MUL(uemg2,tmp,pphases[2].re,pphases[2].im); 
      }
      if(!isPeriodic[3])
      {         
         qcd_PHASE_MUL(uepg3,tmp,mphases[3].re,mphases[3].im);
         qcd_PHASE_MUL(uemg3,tmp,pphases[3].re,pphases[3].im); 
      }
      qcd_SUM_OP(psi,diag,uepg0,uemg0,uepg1,uemg1,uepg2,uemg2,uepg3,uemg3);                 
   }//end loop over boundary points
   
   
   
   return(0);
}//end applyWilsonTMOpNoOptim 


int qcd_applyQTMOp(qcd_vector *Dv, qcd_vector *v, qcd_gaugeField *u, qcd_real_8 m, qcd_real_8 tm)
{
   /* applies the operator -i*gamma_5*D to vector v, where D is the usual WilsonTM operator 
    * works only in basis with gamma_5 = diag(-1,-1,1,1) */
   
   qcd_uint_8 i;
   qcd_uint_2 mu, col;
   
   qcd_applyWilsonTMOpNoOptim(Dv, v, u,  m, tm);
   
   for(i=0; i<Dv->geo->lV; i++)
   for(col=0; col<3; col++)
   {
      Dv->D[i][0][col] = (qcd_complex_16) {-Dv->D[i][0][col].im, Dv->D[i][0][col].re};
      Dv->D[i][1][col] = (qcd_complex_16) {-Dv->D[i][1][col].im, Dv->D[i][1][col].re};
      Dv->D[i][2][col] = (qcd_complex_16) {Dv->D[i][2][col].im, -Dv->D[i][2][col].re};
      Dv->D[i][3][col] = (qcd_complex_16) {Dv->D[i][3][col].im, -Dv->D[i][3][col].re};
   }
   return(0);
}//end qcd_applyQTMOp
