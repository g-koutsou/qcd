#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <math.h>
#include <lime.h>
#include "qcd.h"
#include "qcd_init.c"
#include "qcd_communication.c"
#include "qcd_blas.c"
#include "qcd_io.c"
#include "qcd_wilson.c"
#include "qcd_smearing.c"


int main(int argc,char* argv[])
{
   qcd_geometry geo;
   qcd_uint_2 X[4];
   qcd_uint_2 P[4]={2, 2, 1, 2};
   qcd_uint_2 L[4]={48, 24, 24, 24};
   qcd_real_8 theta[4]={M_PI,0,0,0}; // boundary conditions
   qcd_vector vec,vec2,vec3;
   qcd_gaugeField u;
   qcd_propagator prop;
   int c1,c2,t,x,y,z,mu,nu;
   qcd_real_8 kappa=0.160856;
   qcd_real_8 residue,nv,nv2;
   double t1,t2;
   qcd_real_8 *tst1,*tst2,*uu;
   
   int myid,numprocs, namelen;    
   char processor_name[MPI_MAX_PROCESSOR_NAME];
   
   //set up MPI
   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD,&numprocs);         // num. of processes taking part in the calculation
   MPI_Comm_rank(MPI_COMM_WORLD,&myid);             // each process gets its ID
   MPI_Get_processor_name(processor_name,&namelen); // 
   
   qcd_initGeometry(&geo,L,P, theta, myid, numprocs);
   
   qcd_initVector(&vec,&geo);
   qcd_initVector(&vec2,&geo);
   qcd_initVector(&vec3,&geo);
   qcd_initGaugeField(&u,&geo);
   qcd_initPropagator(&prop,&geo);
   
      
   if(myid==0) printf("Local lattice: %i x %i x %i x %i\n",geo.lL[0],geo.lL[1],geo.lL[2],geo.lL[3]);
   
   t1=MPI_Wtime();   
   if(qcd_getGaugeField("/work/hch02h/conf.1000",qcd_GF_LIME,&u))
   {
      fprintf(stderr,"process %i: Error reading gauge field!\n",myid);
      exit(1);  
   }
   t2=MPI_Wtime();
   if(myid==0) printf("got gauge-field in %lf sek\n",t2-t1);
   
   t1=MPI_Wtime();
   if(qcd_getPropagator("sollist.1000_NoAPE",qcd_PROP_CMI,&prop))
   {
      fprintf(stderr,"process %i: Error reading propagator!\n",myid);
      exit(1);  
   }   
   t2=MPI_Wtime();
   if(myid==0) printf("got smeared sources in %lf sek\n",t2-t1);
   
   qcd_copyVectorPropagator(&vec,&prop,0,0);
   qcd_copyVector(&vec2,&vec);
   tst1= (qcd_real_8*) &(vec.D[123][0][0].re);
   tst2= (qcd_real_8*) &(vec2.D[123][0][0].re);
   uu  = (qcd_real_8*) &(u.D[123][2][0][0].re);
   qcd_APPLY_U(uu,tst1,tst2);
   qcd_APPLY_U_DAGGER(uu,tst2,tst1);
   qcd_copyVectorPropagator(&vec,&prop,0,0);
   qcd_subVector(&vec3,&vec2,&vec);
   nv2=qcd_normVector(&vec3);
   if(myid==0) printf("|Udagger U vec - vec| = %lf\n",nv2);
   
   for(nu=0; nu<4; nu++)
   for(c2=0; c2<3; c2++)
   {
      t1=MPI_Wtime();
      qcd_zeroVector(&vec2);
      nv2=qcd_normVector(&vec2);
      //if(myid==0) printf("|zero-vector| = %lf\n",nv2);
      if(myid==0) vec2.D[0][nu][c2] = (qcd_complex_16) {1,0};
      for(mu=0; mu<50; mu++)
      {
         if(qcd_gaussIteration3d(&vec2,&u,4.0,0))
         {
            fprintf(stderr,"Error while smearing!\n");
            exit(1);
         }
         //nv2 = qcd_normVector(&vec2);
         //if(myid==0) printf("iteration %i, |vec2| = %lf\n",mu,nv2);  
      }   
         
      t2=MPI_Wtime();
      if(myid==0) printf("Created smeared source in %lf sek\n",t2-t1);   
      qcd_copyVectorPropagator(&vec,&prop,nu,c2);
      
      qcd_subVector(&vec3,&vec2,&vec);
      nv=qcd_normVector(&vec);
      nv2=qcd_normVector(&vec2);
      residue = qcd_normVector(&vec3) / nv;
      
      if(myid==0)
      {
         printf("smeared vector nu=%i c=%i: |parallel - sequential| / |sequential|= %lf\n",nu,c2,residue);
         printf("|sequential| = %lf\n",nv);
         printf("|parallel|   = %lf\n",nv2);
      }         
   }      


   qcd_destroyVector(&vec);
   qcd_destroyVector(&vec2);
   qcd_destroyVector(&vec3);
   qcd_destroyPropagator(&prop);
   qcd_destroyGaugeField(&u);
   qcd_destroyGeometry(&geo);
   MPI_Finalize();         
}   
