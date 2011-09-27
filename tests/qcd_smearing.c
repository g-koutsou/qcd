/* qcd_smearing.c
 *
 * gauss-smearing of fermion fields
 *
 * Tomasz Korzec 2008
 **********************************************/
 
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
   qcd_uint_4 x,y,z,b0,b1,b2,b3,tt;
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
      fprintf(stderr,"process %i: Error in qcd_gaussIteration3d! Could not initialize vector");
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
