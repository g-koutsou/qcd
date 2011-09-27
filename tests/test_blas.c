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


int main(int argc,char* argv[])
{
   qcd_geometry geo;
   qcd_uint_2 X[4];
   qcd_uint_2 P[4]={2, 2, 2, 2};
   qcd_uint_2 L[4]={8, 8, 8, 8};
   qcd_real_8 theta[4]={M_PI,0,0,0}; // boundary conditions
   qcd_vector vec,vec2,vec3;
   qcd_gaugeField u,u2,u3;
   int c1,c2,t,x,y,z,mu,nu;
   qcd_uint_8 i,j;
   qcd_real_8 rest;
   
   int myid,numprocs, namelen;    
   char processor_name[MPI_MAX_PROCESSOR_NAME];
   
   //set up MPI
   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD,&numprocs);         // num. of processes taking part in the calculation
   MPI_Comm_rank(MPI_COMM_WORLD,&myid);             // each process gets its ID
   MPI_Get_processor_name(processor_name,&namelen); // 
   
   srand48(myid*137);
   if(qcd_initGeometry(&geo,L,P, theta, myid, numprocs))
   {
      fprintf(stderr,"fuck!\n");
      exit(1);
   }
   
   qcd_initVector(&vec,&geo);
   qcd_initVector(&vec2,&geo);
   qcd_initVector(&vec3,&geo);
   qcd_initGaugeField(&u,&geo);
   qcd_initGaugeField(&u2,&geo);
   qcd_initGaugeField(&u3,&geo);      
      
   if(myid==0) printf("Local lattice: %i x %i x %i x %i\n",geo.lL[0],geo.lL[1],geo.lL[2],geo.lL[3]);
      
   for(i=0; i<geo.lV; i++)
   for(mu=0;mu<4;mu++)
   for(c1=0;c1<3;c1++)
   {
      vec.D[i][mu][c1] = (qcd_complex_16) {drand48(), drand48()};
      for(c2=0;c2<3;c2++)
         u.D[i][mu][c1][c2] = (qcd_complex_16) {drand48(), drand48()};
   }
   
   qcd_copyVector(&vec3,&vec);
   qcd_copyGaugeField(&u3,&u);
   

   for(t=0; t<geo.L[0];t++)
      qcd_scaleVector3d(&vec,2.0,t);
      
   qcd_scaleVector(&vec,0.5);      
   qcd_subVector(&vec2,&vec,&vec3);
   rest = qcd_normVector(&vec2);
   if(myid==0)
   {
      if(rest>1e-14)
         printf("something wrong in scale-vector, rest = %e\n",rest);
      else
         printf("scaleVector & scaleVector3d OK\n");               
   }
   
   
   qcd_addVector(&vec3,&vec2,&vec);
   qcd_subVector(&vec3,&vec3,&vec);
   qcd_subVector(&vec,&vec3,&vec2);
   rest = qcd_normVector(&vec);
   if(myid==0)
   {
      if(rest>1e-14)
         printf("something wrong in add/sub Vector, rest = %e\n",rest);
      else
         printf("add/sub vector OK\n");               
   }   
   
   qcd_copyVector(&vec3,&vec2);   
   qcd_addVector,(&vec,&vec2,&vec3);
   qcd_scaleVector(&vec,0.5);
   qcd_subVector(&vec2,&vec,&vec3);
   rest = qcd_normVector(&vec2);
   if(myid==0)
   {
      if(rest>1e-14)
         printf("something wrong in add-vector, rest = %e\n",rest);
      else
         printf("addVector OK\n");               
   }

   qcd_destroyVector(&vec);
   qcd_destroyVector(&vec2);
   qcd_destroyVector(&vec3);
   qcd_destroyGaugeField(&u);
   qcd_destroyGaugeField(&u2);
   qcd_destroyGaugeField(&u3);      
   qcd_destroyGeometry(&geo);
   MPI_Finalize();         
}   
