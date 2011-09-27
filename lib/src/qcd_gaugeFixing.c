/* qcd_gaugeFixing.c
 *
 * collection of gauge-fixing routines
 *
 * Tomasz Korzec 2009
 **********************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <qcd.h>


/* returns the lattice gluon field a, calculated from
   SU(3) gauge-field u
   a = 1/(2i) [u-u^dagger] - trace
*/   
void qcd_gluonField(qcd_gaugeField *a, qcd_gaugeField *u)
{
   qcd_int_2 c1,c2,mu;
   qcd_int_4 x;
   qcd_real_8 tr;
   
   for(x=0; x<(u->geo->lV+u->geo->boundarySize); x++)
   for(mu=0; mu<4; mu++)
   {
      tr= 0.333333333333333333333*(u->D[x][mu][0][0].im + u->D[x][mu][1][1].im +u->D[x][mu][2][2].im);
      for(c1=0; c1<3; c1++)
      for(c2=0; c2<3; c2++)
      {
         a->D[x][mu][c1][c2] = (qcd_complex_16) {0.5*(u->D[x][mu][c2][c1].im + u->D[x][mu][c1][c2].im),
                                                 0.5*(u->D[x][mu][c2][c1].re - u->D[x][mu][c1][c2].re)};
      }
      for(c1=0; c1<3; c1++)
         a->D[x][mu][c1][c1].re -= tr;
   }
   return;
}//end qcd_gluonField



/* returns the lattice gluon field a_mu(x), calculated from
   SU(3) gauge-field u
   a = 1/(2i) [u-u^dagger] - trace
*/   
void qcd_localGluonField(qcd_complex_16 a[3][3], qcd_gaugeField *u,qcd_uint_8 x, qcd_uint_2 mu)
{
   qcd_uint_2 c1,c2;
   qcd_real_8 tr;
   
   tr= 0.333333333333333333333*(u->D[x][mu][0][0].im + u->D[x][mu][1][1].im +u->D[x][mu][2][2].im);
   for(c1=0; c1<3; c1++)
   for(c2=0; c2<3; c2++)
   {
      a[c1][c2] = (qcd_complex_16) {0.5*(u->D[x][mu][c2][c1].im + u->D[x][mu][c1][c2].im),
                                    0.5*(u->D[x][mu][c2][c1].re - u->D[x][mu][c1][c2].re)};
   }
   for(c1=0; c1<3; c1++)
      a[c1][c1].re -= tr;

   return;
}//end qcd_localGluonField


/* stochastic over-relaxation algorithm to fix u to landau-gauge.
   a local maximum of the functional
   sum_{mu,x} Re Tr [ U_mu (x) ]
   is returned
   The iteration is stopped when the precision qcd_LANDAU_GAUGE_FIXING_PRECISION is reached
   for the quality-measure
   1/V sum_x Tr(A_mu(x)-A_mu(x-mu)) ^ 2
   or the iteration number surpasses qcd_LANDAU_GAUGE_FIXING_MAX_ITER
*/   
int qcd_landauGauge(qcd_gaugeField *landauu, qcd_gaugeField *u, qcd_real_8 overparam)
{
   qcd_real_8 prec=1.0,tr,trold=0.0;
   qcd_int_4 iter=0,x,xeo;
   qcd_int_2 c1,c2,mu,eo;
   qcd_complex_16 C;
   qcd_uint_2 coord[4];
     
   qcd_complex_16 tmp[3][3],tmp2[3][3];
   qcd_gaugeTransformation g;
   qcd_complex_16 glu1[3][3],glu2[3][3];
   int myid = u->geo->myid;   
   
   /* put some tests for the input parameters here */
   

   qcd_copyGaugeField(landauu,u);   
   //if(myid==0) printf("gauge field copied\n");
   
   qcd_initGaugeTransformation(&g,u->geo);
      
   while(prec>qcd_LANDAU_GAUGE_FIXING_PRECISION && iter<qcd_LANDAU_GAUGE_FIXING_MAX_ITER)
   {  
      for(eo=0; eo<2; eo++)
      {    
         //start gauge field communication      
         qcd_communicateGaugeP(landauu);               
         qcd_waitall(landauu->geo);          //version 1 without latency hiding
         //if(myid==0) printf("gauge field communicated\n");
         
         /*
         //estimate precision
         if(eo==1)
         {
            tr=0;
            prec=0;
            for(x=0; x<u->geo->lV; x++)
            {
               qcd_zero3x3(tmp2);
               for(mu=0; mu<4; mu++)
               {
                  qcd_localGluonField(glu1,landauu,x,mu);
                  qcd_localGluonField(glu2,landauu,u->geo->minus[x][mu],mu);
                  qcd_sub3x3(tmp,glu1,glu2);
                  qcd_add3x3(tmp2,tmp2,tmp);
               }
               qcd_mul3x3(tmp,tmp2,tmp2);
               tr += tmp[0][0].re + tmp[1][1].re + tmp[2][2].re;
            }            
         }
         */
      
         for(xeo=0; xeo<u->geo->lV/2; xeo++)
         {
             x=u->geo->eo[xeo][eo];             

             /*/////////////////////////////////////// transformation matrix as in Dina's code ////////////////*/
             qcd_zero3x3(g.D[x]);
             for(mu=0; mu<4; mu++)
             {
                qcd_add3x3(g.D[x],g.D[x],landauu->D[x][mu]);
                qcd_addAdjoint3x3(g.D[x],g.D[x],landauu->D[u->geo->minus[x][mu]][mu]);
             }

             //project transformation matrix g back to SU(3)
             qcd_projectSU33x3(g.D[x]);

             if(drand48()<overparam)
             {
                qcd_MUL3x3(tmp,g.D[x],g.D[x]);
                qcd_copy3x3(g.D[x],tmp);
             }
             
             qcd_dagger3x3(g.D[x]);

         }//end volume loop
         
         //if(myid==0) printf("transformation calculated\n");
         
         
         //communicate transformation field
         qcd_communicateTransformationM(&g);
         qcd_waitall(g.geo);
         //if(myid==0) printf("transformation communicated\n");
                  
         for(xeo=0; xeo<u->geo->lV/2; xeo++)
         {
             // 1) apply transformation to links even (odd) links
             x=u->geo->eo[xeo][eo];               
             for(mu=0; mu<4; mu++)
             {
                qcd_MUL3x3(tmp,g.D[x],landauu->D[x][mu]);
                qcd_copy3x3(landauu->D[x][mu],tmp);
             }
             
             // 2) apply transformation to links odd (even) links
             x=u->geo->eo[xeo][(eo+1)%2];
             for(mu=0; mu<4; mu++)
             {
                qcd_MULADJOINT3x3(tmp,landauu->D[x][mu],g.D[u->geo->plus[x][mu]]);
                qcd_copy3x3(landauu->D[x][mu],tmp);
             }
         }//end volume loop         

      }//end even-odd loop                  
      
      //if(myid==0) printf("local gauge trafo on all sites done\n"); 
      
               
      //MPI_Allreduce(&tr, &prec, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      //prec/=4*4*u->geo->V;
      
      iter++;
      //if(myid==0) printf("Iteration %i, 1/4V tr( (d_mu A_mu)(d_mu A_mu)^+ ): %e\n",iter,prec);
            
      //calculate dina's stopping criterion
      prec=0;
      for(x=0; x<u->geo->lV; x++)
      for(mu=0; mu<4; mu++)
      {
         prec += qcd_trace3x3(landauu->D[x][mu]).re;
      }
      MPI_Allreduce(&prec, &tr, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      prec = abs(trold-tr)/tr;
      trold=tr;
      //if(myid==0) printf("real trace: %e\n", prec/(4*u->geo->V));      
   }//end iteration-loop 
   
   if(myid==0) printf("Gauge Fixed in %i iterations. Real trace: %e\n",iter, tr/(4*u->geo->V));
   if(prec>qcd_LANDAU_GAUGE_FIXING_PRECISION)
   {
      if(u->geo->myid==0)
         fprintf(stderr,"warning in qcd_landauGauge: required precision not reached after %i iterations\n",iter);
   }

   qcd_destroyGaugeTransformation(&g);
   return(iter);
}//end qcd_landauGauge
