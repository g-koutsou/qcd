/* 
 *
 * reads solution vectors
 * and creates a disconnected loop via
 * the one-end trick
 *
 * G.K 2011
 ****************************************/
 
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <mpi.h>
#include <qcd.h>
#include "projectors.h"


int main(int argc,char* argv[])
{
   qcd_uint_2 mu,nu,ku,lu,c1,c2,c3,c1p,c2p,c3p;// various loop variables
   qcd_uint_2 id1,id2,id3,cc1,cc2,al,be;
   qcd_uint_4 i,j,k, v,lx,ly,lz,ip1,im1,v3; 
   qcd_int_4 x,y,z;
   qcd_uint_2 ic1,ic2,ic3;                    //
   qcd_uint_4 t_sink, t_start, t_stop, t,lt;
   qcd_real_8 tmp;                            // general purpuse
   FILE *fp_momlist;
  
   FILE *fp_loop;                           // output file
   FILE *fp_sol1_list;
   FILE *fp_sol2_list;
  
   int params_len;                            // needed to read inputfiles
   char *params;                              // needed to read inputfiles

   char loop_name[qcd_MAX_STRING_LENGTH];
   
   char param_name[qcd_MAX_STRING_LENGTH];      // name of parameter file  
   char momlist_name[qcd_MAX_STRING_LENGTH];    // name of momenta-list file
   char sol1_list_name[qcd_MAX_STRING_LENGTH];
   char sol2_list_name[qcd_MAX_STRING_LENGTH]; 

   qcd_geometry geo;                            // geometry structure
   qcd_vector sol1, sol2;                              


   qcd_real_8 theta[4] = {M_PI,0.0,0.0,0.0};    // antiperiodic b.c. in time
   qcd_uint_2 L[4];
   qcd_uint_2 P[4];
   qcd_complex_16 phase_factor;         
   qcd_complex_16 z1, z2;                       // temp variables
   qcd_complex_16 C, C2;   
   qcd_complex_16 loop, loop2;
   qcd_real_8 plaq;
   qcd_int_4 ctr, ctr2;   
   qcd_complex_16 *block[16];                       // to store the block

   qcd_int_4 (*mom)[3];                         // momenta-list

   int myid,numprocs, namelen;    
   char processor_name[MPI_MAX_PROCESSOR_NAME];
             
   //set up MPI
   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD,&numprocs);         // num. of processes taking part in the calculation
   MPI_Comm_rank(MPI_COMM_WORLD,&myid);             // each process gets its ID
   MPI_Get_processor_name(processor_name,&namelen); // 
               
   
   //////////////////// READ INPUT FILE /////////////////////////////////////////////
      
   if(argc!=2)
   {
      if(myid==0) fprintf(stderr,"No input file specified\n");
      exit(EXIT_FAILURE);
   }

   strcpy(param_name,argv[1]);
   if(myid==0)
   {
      i=0;
      printf("opening input file %s\n",param_name);
      params=qcd_getParams(param_name,&params_len);
      if(params==NULL)
      {
         i=1;
      }
   }
   MPI_Bcast(&i,1,MPI_INT, 0, MPI_COMM_WORLD);
   if(i==1) exit(EXIT_FAILURE);
   MPI_Bcast(&params_len, 1, MPI_INT, 0, MPI_COMM_WORLD);
   if(myid!=0) params = (char*) malloc(params_len*sizeof(char));
   MPI_Bcast(params, params_len, MPI_CHAR, 0, MPI_COMM_WORLD);
   
   sscanf(qcd_getParam("<processors_txyz>",params,params_len),"%hd %hd %hd %hd",&P[0], &P[1], &P[2], &P[3]);
   sscanf(qcd_getParam("<lattice_txyz>",params,params_len),"%hd %hd %hd %hd",&L[0], &L[1], &L[2], &L[3]);
   if(P[0] != 1)
   {
      if(myid==0) fprintf(stderr,"Error! Number of processors in t-direction must be one.\n");
      exit(EXIT_FAILURE);
   }
   if(qcd_initGeometry(&geo,L,P, theta, myid, numprocs)) exit(EXIT_FAILURE);   
   
   if(myid==0) printf(" Local lattice: %i x %i x %i x %i\n",geo.lL[0],geo.lL[1],geo.lL[2],geo.lL[3]);
      
   sscanf(qcd_getParam("<t>",params,params_len),"%d %d",&t_start, &t_stop);
   if(myid==0) printf("Got loop time slices: %d ... %d\n",t_start,t_stop);
         
   int n_sol = -1;
   sscanf(qcd_getParam("<numb_sol>",params,params_len),"%d", &n_sol);
   if(myid==0) printf("Got number of solutions: %d\n", n_sol);
            	  
   strcpy(sol1_list_name,qcd_getParam("<solution_list_1>",params,params_len));
   if(myid==0) printf("Got first solution list: %s\n",sol1_list_name);

   strcpy(sol2_list_name,qcd_getParam("<solution_list_2>",params,params_len));
   if(myid==0) printf("Got first solution list: %s\n",sol2_list_name);
                  
   strcpy(loop_name,qcd_getParam("<loop_name>",params,params_len));
   if(myid==0) printf("Got output file name: %s\n",loop_name);
           
   strcpy(momlist_name,qcd_getParam("<momenta_list>",params,params_len));
   if(myid==0) printf("Got momenta-list file name: %s\n",momlist_name);
    
   free(params);

   
   qcd_initVector(&sol1, &geo);
   qcd_initVector(&sol2, &geo);       
   
   //################################################################################   
   // calculate the 'blocks'      
   
   k = 0;
   //open output file to write in
   if(myid==0)
   {
      fp_loop = fopen(loop_name,"w");   
      if(fp_loop==NULL)
      {
         printf("failed to open %s for writing\n",loop_name);
         k=1;
      }
      if(myid == 0)
	fprintf(fp_loop, "# i_sol p(x, y, z)  g   t            re im\n");

      fp_momlist = fopen(momlist_name,"r");
      if(fp_momlist==NULL)
      {
	printf("failed to open %s for reading\n",momlist_name);
         k=1;
      }
   }
   MPI_Bcast(&k,1,MPI_INT, 0, MPI_COMM_WORLD);
   if(k==1) exit(EXIT_FAILURE);
  
   //load momenta-list   
   int n_mom; 
   if(myid==0) fscanf(fp_momlist,"%i\n",&n_mom);
   MPI_Bcast(&n_mom,1,MPI_INT, 0, MPI_COMM_WORLD);
   if(myid==0) printf("will read %i momenta combinations\n",n_mom);

   mom = malloc(n_mom*3*sizeof(qcd_int_4));

   if(myid==0)
   {
      for(j=0; j<n_mom; j++)
	{
         fscanf(fp_momlist,"%i %i %i\n",&(mom[j][0]),&(mom[j][1]),&(mom[j][2]));
      }
      fclose(fp_momlist);   
   }
   MPI_Bcast(&(mom[0][0]),n_mom*3,MPI_INT,0, MPI_COMM_WORLD);
   if(myid==0) printf("momenta list read and broadcasted\n");   
     
   
   if((fp_sol2_list = fopen(sol2_list_name, "r")) == NULL)
     {
       if(myid == 0)
	 {
	   fprintf(stderr,"%s: failed to open file for reading\n", sol2_list_name);
	   exit(EXIT_FAILURE);
	 }
     }

   if((fp_sol1_list = fopen(sol1_list_name, "r")) == NULL)
     {
       if(myid == 0)
	 {
	   fprintf(stderr,"%s: failed to open file for reading\n", sol1_list_name);
	   exit(EXIT_FAILURE);
	 }
     }

   for(i=0; i<16; i++)
     block[i]= malloc(sizeof(qcd_complex_16)*geo.lV3);

   for(int isol=0; isol<n_sol; isol++)
     {
       char sol1_name[qcd_MAX_STRING_LENGTH];
       char sol2_name[qcd_MAX_STRING_LENGTH];

       fgets(sol1_name, qcd_MAX_STRING_LENGTH, fp_sol1_list);
       fgets(sol2_name, qcd_MAX_STRING_LENGTH, fp_sol2_list);

       *(index(sol1_name, '\n')) = '\0';
       *(index(sol2_name, '\n')) = '\0';

       qcd_getVectorLime(sol1_name, &sol1);
       qcd_getVectorLime(sol2_name, &sol2);


       for(t=t_start; t<=t_stop; t++)
	 {
	   lt = t % geo.lL[0];
	   int t_proc = t / geo.lL[0];
	   for(int ig=0; ig<16; ig++)
	     {
	       for(v3=0; v3<geo.lV3; v3++)   //set blocks to zero
		 block[ig][v3]= (qcd_complex_16) {0,0};
	       
	     }
	   if(geo.Pos[0] == t_proc)  //inside the local lattice, otherwise nothing to calculate
	     {
	       if(myid==0) 
		 printf("t=%i\n",t);

	       for(lx=0; lx<geo.lL[1]; lx++)
		 for(ly=0; ly<geo.lL[2]; ly++)
		   for(lz=0; lz<geo.lL[3]; lz++)
		     for(int igamma=0; igamma<16; igamma++)
		     {
		       int iv = qcd_LEXIC(lt, lx, ly, lz, geo.lL);
		       int is = qcd_LEXIC0(lx, ly, lz, geo.lL);
		       int icol = igamma % 4;
		       int irow = igamma / 4;
		       for(int c=0; c<3; c++)
			 for(int mu=0; mu<4; mu++)
			   {
			     block[igamma][is].re += qcd_CMULR(
							       qcd_CMUL(qcd_CONJ(sol1.D[iv][mu][c]), qcd_GAMMA[5][mu][irow]),
							       sol2.D[iv][icol][c]);

			     block[igamma][is].im += qcd_CMULI(
							       qcd_CMUL(qcd_CONJ(sol1.D[iv][mu][c]), qcd_GAMMA[5][mu][irow]),
							       sol2.D[iv][icol][c]);
			   }
		     }
	       //Fourier transform time-slice
	       for(int igamma=0; igamma<16; igamma++)
		 for(j=0; j<n_mom; j++)
		   {
		     /* if(myid==0) */
		     /*   { */
		     /*     printf("%i %+i %+i %+i\n",t,mom[j][0],mom[j][1],mom[j][2]); */
		     /*   } */
		     loop = (qcd_complex_16) {0,0};
		     
		     for(lx=0; lx<geo.lL[1]; lx++)
		       for(ly=0; ly<geo.lL[2]; ly++)
			 for(lz=0; lz<geo.lL[3]; lz++)
			   {
			     v3 = qcd_LEXIC0(lx,ly,lz,geo.lL);
			     x=lx+geo.Pos[1]*geo.lL[1];
			     y=ly+geo.Pos[2]*geo.lL[2];
			     z=lz+geo.Pos[3]*geo.lL[3];
			     tmp = (((double) mom[j][0]*x)/geo.L[1] + 
				    ((double) mom[j][1]*y)/geo.L[2] + 
				    ((double) mom[j][2]*z)/geo.L[3])*2*M_PI;
			     
			     C2=(qcd_complex_16) {cos(tmp), -sin(tmp)}; //TABULATE FOR LARGE SPEEDUP!!!
			     loop = qcd_CADD(loop, qcd_CMUL(block[igamma][v3],C2));
			   }
		     
		     MPI_Reduce(&(loop.re), &(loop2.re), 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		     if(myid == 0)
		       //                                         # i_sol p(x, y, z) g t    re im
		       fprintf(fp_loop, "  %5d  %+2d %+2d %+2d  %2d %3d %+e %+e\n", 
			       isol, mom[j][0], mom[j][1], mom[j][2], igamma, t, loop2.re, loop2.im);
		   }
	       
	     }
	 }
     }

   if(myid==0)
   {
      fclose(fp_loop);
   }   
   
                  
   
   //#####################################################################   
   // clean up
   if(myid==0) printf("cleaning up...\n");
   
   for(int igamma=0; igamma<16; igamma++)
     free(block[igamma]);
   free(mom);
   qcd_destroyVector(&sol1);
   qcd_destroyVector(&sol2);
   qcd_destroyGeometry(&geo);
   MPI_Finalize();
}//end main
