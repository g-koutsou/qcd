#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include "qcd.h"
#include "qcd_init.c"
#include "qcd_blas.c"
#include "qcd_io.c"

int main(int argc,char* argv[])
{
   qcd_geometry geo;
   qcd_uint_2 X[4];
   qcd_uint_2 P[4]={1, 2, 2, 2};
   qcd_uint_2 L[4]={48, 24, 24, 24};
   qcd_vector vec,vec2;
   qcd_gaugeField u;
   qcd_propagator prop;
   int c1,c2,t,x,y,z,mu,nu;
   double t1,t2;
   
   
   int myid,numprocs, namelen;    
   char processor_name[MPI_MAX_PROCESSOR_NAME];
   
   //set up MPI
   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD,&numprocs);         // num. of processes taking part in the calculation
   MPI_Comm_rank(MPI_COMM_WORLD,&myid);             // each process gets its ID
   MPI_Get_processor_name(processor_name,&namelen); // 
   
   qcd_initGeometry(&geo,L,P, myid, numprocs);
   
   if(myid==0)
      printf("sizeof(qcd_real_4) = %i\n",sizeof(qcd_real_4));
   t1=MPI_Wtime();
   qcd_initGaugeField(&u,&geo);
   t2=MPI_Wtime();
   if(myid==0)
      printf("time for initializing gauge field %lf\n",t2-t1);
   t1=MPI_Wtime();   
   qcd_initPropagator(&prop,&geo);
   t2=MPI_Wtime();
   if(myid==0)
      printf("time for initializing propagator %lf\n",t2-t1);   
   
   t1=MPI_Wtime();
   if(qcd_getGaugeField("/work/hch02h/qcd_test/conf.1000",qcd_GF_LIME,&u))
   {
      fprintf(stderr,"process %i: Error reading gauge field!\n",myid);
      exit(1);  
   }
   t2=MPI_Wtime();
   if(myid==0)
      printf("time for reading gauge field %lf\n",t2-t1);
   
   t1=MPI_Wtime();
   if(qcd_getPropagator("/work/hch02h/qcd_test/sollist.1000",qcd_PROP_CMI,&prop))
   {
      fprintf(stderr,"process %i: Error reading propagator!\n",myid);
      exit(1);  
   }
   t2=MPI_Wtime();
   if(myid==0)
      printf("time for reading propagator %lf\n",t2-t1);   
   
   
   qcd_destroyGaugeField(&u);
   qcd_destroyPropagator(&prop);   
   qcd_destroyGeometry(&geo);
   MPI_Finalize();         
}   
