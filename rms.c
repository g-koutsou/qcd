/* source.c
 *
 * creates smeared sources
 * Gaussian smearing with 
 * APE-smeared gauge fields
 *
 * Tomasz Korzec 2009
 ****************************************/
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <math.h>
#include <qcd.h>
 
int main(int argc,char* argv[])
{
   FILE*  pfile;
   char*  params;
   char   gauge_name[qcd_MAX_STRING_LENGTH];
   char   param_name[qcd_MAX_STRING_LENGTH];
   char   out_name[qcd_MAX_STRING_LENGTH];
   qcd_int_4   x_src[4],lx_src[4],i,j,t,isource,nsmear,nsmearAPE;
   qcd_uint_2   mu,nu,col,c1,c2,s;
   qcd_real_8   alpha_i, alpha_f, d_alpha ,alphaAPE,plaq;
   int params_len;   

   qcd_geometry geo;
   qcd_gaugeField u, uAPE;
   qcd_gaugeField *u_ptr, *uAPE_ptr, *utmp_ptr;
   qcd_propagator source;
   qcd_vector vec;
   qcd_uint_2 P[4];
   qcd_uint_2 L[4];
   qcd_real_8 theta[4]={M_PI,0.,0.,0.}; // boundary conditions

   qcd_uint_4    ismear,nsources;
   qcd_uint_4    ape_ismear;
   
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
   if(qcd_initGeometry(&geo,L,P, theta, myid, numprocs)) exit(EXIT_FAILURE);
   
   if(myid==0) printf(" Local lattice: %i x %i x %i x %i\n",geo.lL[0],geo.lL[1],geo.lL[2],geo.lL[3]);  

   sscanf(qcd_getParam("<source_pos_txyz>",params,params_len),"%d %d %d %d",&x_src[0], &x_src[1], &x_src[2], &x_src[3]);
   if(myid==0) printf(" Got source coords: %d %d %d %d\n",x_src[0],x_src[1],x_src[2],x_src[3]);
   
   sscanf(qcd_getParam("<alpha_gauss>",params,params_len),"%lf %lf %lf",&alpha_i, &d_alpha, &alpha_f);
   if(myid==0) printf(" Got alpha_gauss: from %lf, to %lf, every %lf\n",alpha_i, alpha_f, d_alpha);
   sscanf(qcd_getParam("<nsmear_gauss>",params,params_len),"%d",&nsmear);
   if(myid==0) printf(" Got nsmear_gauss: %d\n",nsmear);
   sscanf(qcd_getParam("<alpha_APE>",params,params_len),"%lf",&alphaAPE);
   if(myid==0) printf(" Got alpha_APE: %lf\n",alphaAPE);
   sscanf(qcd_getParam("<nsmear_APE>",params,params_len),"%d",&nsmearAPE);
   if(myid==0) printf(" Got nsmear_APE: %d\n",nsmearAPE);   
   strcpy(gauge_name,qcd_getParam("<cfg_name>",params,params_len));
   if(myid==0) printf(" Got conf name: %s\n",gauge_name);
  

   free(params);      
   ///////////////////////////////////////////////////////////////////////////////////////////////////
   
   if(nsmear != 0)
     {
       qcd_initGaugeField(&u,&geo);
       qcd_initGaugeField(&uAPE,&geo);

       if(qcd_getGaugeField(gauge_name,qcd_GF_LIME,&u))
	 {
	   fprintf(stderr,"process %i: Error reading gauge field!\n",myid);
	   exit(EXIT_FAILURE);
	 }
       
       if(myid==0) printf("gauge-field loaded\n");
       plaq = qcd_calculatePlaquette(&u);
       if(myid==0) printf("plaquette = %e\n",plaq);
       
       u_ptr = &u;
       uAPE_ptr = &uAPE;   
       for(i=0; i<nsmearAPE; i++)
	 {
	   qcd_apeSmear3d(uAPE_ptr, u_ptr, alphaAPE);
	   utmp_ptr=u_ptr; u_ptr=uAPE_ptr; uAPE_ptr=utmp_ptr;   
	 }
       utmp_ptr=u_ptr; u_ptr=uAPE_ptr; uAPE_ptr=utmp_ptr; //reverse the last swap. Also needed when nsmearAPE=0
       qcd_destroyGaugeField(u_ptr);
       uAPE = *uAPE_ptr;
    
       if(myid==0) printf("gauge-field APE-smeared\n");
       plaq = qcd_calculatePlaquette(&uAPE);
       if(myid==0) printf("plaquette = %e\n",plaq); 
   
       //qcd_initPropagator(&source,&geo);
     }
   
   qcd_initVector(&vec,&geo);
  
   // which process has the source coords?
   for(i=0; i<4; i++)
      lx_src[i] = x_src[i]/geo.lL[i];
          
   for(double alpha=alpha_i; alpha <= alpha_f; alpha+=d_alpha)
     {       
       qcd_zeroVector(&vec); 
       if( (lx_src[0]==geo.Pos[0]) && 
	   (lx_src[1]==geo.Pos[1]) && 
	   (lx_src[2]==geo.Pos[2]) && 
	   (lx_src[3]==geo.Pos[3]) )
	 {
	   vec.D[qcd_LEXIC((x_src[0]%geo.lL[0]),
	   		   (x_src[1]%geo.lL[1]),
	   		   (x_src[2]%geo.lL[2]),
	   		   (x_src[3]%geo.lL[3]), geo.lL)][0][0].re=1.;
	 }
       for(i=0; i<nsmear; i++)
	 {
	   if(qcd_gaussIteration3d(&vec,&uAPE,alpha,x_src[0]))
	     {
	       fprintf(stderr,"process %i: Error while smearing!\n",geo.myid);
	       exit(EXIT_FAILURE);
	     }

	   double denom = 0, sum_r2 = 0;
	   if( lx_src[0]==geo.Pos[0] )
	     {
	       for(int z=0; z<geo.lL[3]; z++)
		 for(int y=0; y<geo.lL[2]; y++)
		   for(int x=0; x<geo.lL[1]; x++)
		     {
		       int gx = abs(geo.lL[1]*geo.Pos[1] + x - x_src[1]);
		       int gy = abs(geo.lL[2]*geo.Pos[2] + y - x_src[2]);
		       int gz = abs(geo.lL[3]*geo.Pos[3] + z - x_src[3]);

		       gx = gx>geo.L[1]/2 ? geo.L[1] - gx : gx;
		       gy = gy>geo.L[2]/2 ? geo.L[2] - gy : gy;
		       gz = gz>geo.L[3]/2 ? geo.L[3] - gz : gz;

		       int lv = qcd_LEXIC((x_src[0]%geo.lL[0]), x, y, z, geo.lL);
		       double loc_norm = 0;
		       for(int mu=0; mu<4; mu++)
			 for(int c=0; c<3; c++)
			   {
			     loc_norm += qcd_NORMSQUARED(vec.D[lv][mu][c]);
			   }
		       denom += loc_norm;
		       sum_r2 += loc_norm*(gx*gx + gy*gy + gz*gz);
		     }
	     }
	   MPI_Barrier(MPI_COMM_WORLD);
	   double glob_denom, glob_sum_r2;
	   MPI_Reduce(&denom, &glob_denom, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	   MPI_Reduce(&sum_r2, &glob_sum_r2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	   if(geo.myid == 0)
	     printf(" Nsmear = %3d, alpha = %10.5f, <r^2> = %+e\n", i+1, alpha, glob_sum_r2/glob_denom);
	 }
       if(geo.myid == 0)
	 printf("\n\n");
     }
   qcd_destroyVector(&vec);

   if(nsmear != 0)
     qcd_destroyGaugeField(&uAPE); 
   //qcd_writePropagator(out_name, qcd_PROP_CMI, &source);
   
   qcd_rng_finalize();
   
   ////////////////////////////////////// CLEAN UP AND EXIT ///////////////////////////////////////////
   //qcd_destroyPropagator(&source);
   qcd_destroyGeometry(&geo);
   MPI_Finalize();
   return(EXIT_SUCCESS);
}//end main
