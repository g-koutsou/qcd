#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <math.h>
#include <lime.h>
#include <qcd.h>


int main(int argc,char* argv[])
{
   qcd_geometry geo;
   qcd_uint_2 X[4],c1,c2,mu;
   qcd_uint_8 x;
   qcd_uint_2 P[4]={1, 1, 1, 1};
   qcd_uint_2 L[4]={8, 8, 8, 8};
   qcd_real_8 theta[4]={M_PI,0,0,0}; // boundary conditions
   qcd_real_8 plaquette;
   long t1,t2;   
   
   qcd_complex_16 tmp1[3][3], tmp2[3][3], tmp3[3][3];

   qcd_gaugeField u,lu;
   
   int myid,numprocs, namelen;    
   char processor_name[MPI_MAX_PROCESSOR_NAME];
   
   //set up MPI
   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD,&numprocs);         // num. of processes taking part in the calculation
   MPI_Comm_rank(MPI_COMM_WORLD,&myid);             // each process gets its ID
   MPI_Get_processor_name(processor_name,&namelen); // 
   
   qcd_initGeometry(&geo,L,P, theta, myid, numprocs);
   
   qcd_initGaugeField(&u,&geo);
   qcd_initGaugeField(&lu,&geo);

   if(myid==0) printf("Local lattice: %i x %i x %i x %i\n",geo.lL[0],geo.lL[1],geo.lL[2],geo.lL[3]);
   
   t1=MPI_Wtime();   
   if(qcd_getGaugeField("conf88.0000",qcd_GF_LIME,&u))
   {
      fprintf(stderr,"process %i: Error reading gauge field!\n",myid);
      exit(1);  
   }

   
   t2=MPI_Wtime();
   if(myid==0) printf("got gauge-field in %lf sek\n",t2-t1);
   
   //just to make sure
   qcd_projectSU3(&u);
/*
   for(x=0; x<geo.lV; x++)
   for(mu=0;mu<4;mu++)
   for(c1=0;c1<3;c1++)
   for(c2=0;c2<3;c2++)
   if(isnan(u.D[x][mu][c1][c2].re) || isnan(u.D[x][mu][c1][c2].im))
      printf("corrupt gauge field in x=%i mu=%i c1=%i c2=%i\n",x,mu,c1,c2);
   */
   
/*   //tests of 3x3 blas
   for(c1=0; c1<3; c1++)
   {
      for(c2=0; c2<3; c2++)
         printf("%f%+f ",u.D[12][2][c1][c2]);
      printf("\n");
   }      
   printf("\n");
   qcd_copy3x3(tmp1,u.D[12][2]);
   qcd_copy3x3(tmp2,tmp1);
   qcd_zero3x3(tmp3);
   qcd_sub3x3(tmp1,tmp1,tmp2);
   for(c1=0; c1<3; c1++)
   {
      for(c2=0; c2<3; c2++)
         printf("%f%+f ",tmp1[c1][c2]);
      printf("\n");
   }
   printf("\n");
   qcd_addAdjoint3x3(tmp1,tmp2,tmp2);
   for(c1=0; c1<3; c1++)
   {
      for(c2=0; c2<3; c2++)
         printf("%f%+f ",tmp1[c1][c2]);
      printf("\n");
   }
   printf("\n");  
   qcd_mulAdjoint3x3(tmp1,tmp2,tmp2); 
   for(c1=0; c1<3; c1++)
   {
      for(c2=0; c2<3; c2++)
         printf("%f%+f ",tmp1[c1][c2]);
      printf("\n");
   }
   printf("\n");  
   
   exit(0);
*/   
   
   
   
   
   
   //qcd_randomGaugeField(&u);
   
   
   t1=MPI_Wtime();
   qcd_landauGauge(&lu,&u);
   t2=MPI_Wtime();
   if(myid==0) printf("gauge-fixed in %lf sek\n",t2-t1);   
 
   plaquette = qcd_calculatePlaquette(&lu);
   if(myid==0) printf("plaquette %e\n",plaquette);  
   qcd_writeGaugeField("landauconf.0000",qcd_GF_LIME,&lu,"landau gauged fixed conf. yay!");
   if(qcd_getGaugeField("landauconf.0000",qcd_GF_LIME,&u))
   {
      fprintf(stderr,"process %i: Error reading gauge field!\n",myid);
      exit(1);  
   }
   plaquette = qcd_calculatePlaquette(&u);
   if(myid==0) printf("plaquette %e\n",plaquette);  
   
 
   qcd_destroyGaugeField(&u);
   qcd_destroyGaugeField(&lu);
   qcd_destroyGeometry(&geo);
   MPI_Finalize();         
}   
