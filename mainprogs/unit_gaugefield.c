#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <qcd.h>


int main(int argc,char* argv[])
{
   qcd_geometry geo;
   qcd_uint_2 P[4]={2, 2, 2, 2};
   qcd_uint_2 L[4]={16, 8, 8, 8};
   qcd_real_8 theta[4]={M_PI,0.,0.,0.}; // boundary conditions
   qcd_gaugeField u;
   qcd_uint_2 mu,c1;
   qcd_uint_4 i;
   qcd_real_8 p;
      
   int myid,numprocs, namelen;    
   char processor_name[MPI_MAX_PROCESSOR_NAME];
   
   //set up MPI
   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD,&numprocs);         // num. of processes taking part in the calculation
   MPI_Comm_rank(MPI_COMM_WORLD,&myid);             // each process gets its ID
   MPI_Get_processor_name(processor_name,&namelen); //
   
   qcd_initGeometry(&geo,L,P, theta, myid, numprocs);
   qcd_initGaugeField(&u,&geo);
   
   //create unit gauge-field
   qcd_zeroGaugeField(&u);
   for(i=0; i<geo.lV; i++)
   for(mu=0; mu<4; mu++)
   for(c1=0; c1<3; c1++)
      u.D[i][mu][c1][c1] = (qcd_complex_16) {1.0,0};
   
   p=qcd_calculatePlaquette(&u);
   if(myid==0) printf("plaquette of unit gauge-field: %e\n",p);
   
   qcd_writeGaugeField("unit_conf.0000",qcd_GF_LIME,&u,"plaquette = 1\ntrajectory nr = 0\nbeta = 0.0, kappa = 0.0, mu = 0.0, c2_rec = 0.0");
   
   qcd_destroyGaugeField(&u);
   qcd_destroyGeometry(&geo);
}   
