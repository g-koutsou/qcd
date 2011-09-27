#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <qcd.h>


int main(int argc,char* argv[])
{
   qcd_geometry geo;
   qcd_uint_2 P[4]={2, 2, 2, 2};
   qcd_uint_2 L[4]={48, 24, 24, 24};
   qcd_real_8 theta[4]={M_PI,0.,0.,0.}; // boundary conditions
   qcd_vector vec,vec2,vec3;
   qcd_gaugeField u;
   int t;
   qcd_uint_2 c2,nu;
   qcd_real_8 kappa=0.160856;
   qcd_real_8 residue;
   char srcname[255];
   
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
      
   if(myid==0) printf("Local lattice: %i x %i x %i x %i\n",geo.lL[0],geo.lL[1],geo.lL[2],geo.lL[3]);
      
   if(qcd_getGaugeField("conf.1000",qcd_GF_LIME,&u))
   {
      fprintf(stderr,"process %i: Error reading gauge field!\n",myid);
      exit(EXIT_FAILURE);  
   }
   if(myid==0) printf("got gauge-field\n");
   residue = qcd_calculatePlaquette(&u);
   if(myid==0) printf("plaquette = %e\n",residue);
   
   for(nu=0; nu<4; nu++)
   for(c2=0; c2<3; c2++)
   {
      if(qcd_getVector("prop.mass00.1000.430413_SM",qcd_PROP_HMC,nu,c2,&vec))
      {
         fprintf(stderr,"process %i: Error reading vector!\n",myid);
         exit(EXIT_FAILURE);
      }
      residue = qcd_normVector(&vec);
      if(myid==0) printf("|x|=%e\n", residue);
     
      sprintf(srcname,"src.%02i",nu*3+c2);
      if(qcd_getVector(srcname, qcd_PROP_CMI, nu,c2,&vec2))
      {
         fprintf(stderr,"process %i: Error reading vector!\n",myid);
         exit(EXIT_FAILURE);
      } 
      residue = qcd_normVector(&vec2);
      if(myid==0) printf("|b|=%e\n",residue);

      qcd_applyWilsonTMOp(&vec3,&vec,&u, 1.0/(2.0*kappa)-4.0 ,0.0064);

      residue = qcd_normVector(&vec3);
      if(myid==0) printf("|Dx|=%e\n",residue);      

      qcd_subVector(&vec,&vec3,&vec2);

      residue = qcd_normVector(&vec)  / qcd_normVector(&vec2);
      if(myid==0) printf("residue calculated\n");
      if(myid==0) printf("solution nu=%i c=%i: |Dx-b| / |b|= %e\n",nu,c2, (float) residue);
   }
   
   
   
   if(myid==0) printf("cleaning up:\n");          


   qcd_destroyVector(&vec);
   qcd_destroyVector(&vec2);
   qcd_destroyVector(&vec3);
   qcd_destroyGaugeField(&u);
   qcd_destroyGeometry(&geo);
   MPI_Finalize();         
   return(EXIT_SUCCESS);
}   
