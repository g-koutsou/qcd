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
   char*  params = NULL;
   char   gauge_name[qcd_MAX_STRING_LENGTH];
   char   param_name[qcd_MAX_STRING_LENGTH];
   char   out_name[qcd_MAX_STRING_LENGTH];
   qcd_int_4   x_src[4],lx_src[4],i,nsmear,nsmearAPE;
   qcd_real_8   alpha ,alphaAPE,plaq;
   int params_len;   

   qcd_geometry geo;
   qcd_gaugeField u, uAPE;
   qcd_gaugeField *u_ptr, *uAPE_ptr, *utmp_ptr;
   qcd_vector vec;
   qcd_uint_2 P[4];
   qcd_uint_2 L[4];
   qcd_real_8 theta[4]={M_PI,0.,0.,0.}; // boundary conditions

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
   if(P[0] != 1)
     {
       fprintf(stderr, " Must use only 1 process along t-direction, exiting...\n");
       exit(EXIT_FAILURE);
     }
   sscanf(qcd_getParam("<lattice_txyz>",params,params_len),"%hd %hd %hd %hd",&L[0], &L[1], &L[2], &L[3]);
   if(qcd_initGeometry(&geo,L,P, theta, myid, numprocs)) exit(EXIT_FAILURE);
   
   if(myid==0) printf(" Local lattice: %i x %i x %i x %i\n",geo.lL[0],geo.lL[1],geo.lL[2],geo.lL[3]);  

   sscanf(qcd_getParam("<source_pos_txyz>",params,params_len),"%d %d %d %d",&x_src[0], &x_src[1], &x_src[2], &x_src[3]);
   if(myid==0) printf(" Got source coords: %d %d %d %d\n",x_src[0],x_src[1],x_src[2],x_src[3]);
   
   sscanf(qcd_getParam("<alpha_gauss>",params,params_len),"%lf", &alpha);
   if(myid==0) printf(" Got alpha_gauss: %lf\n",alpha);
   sscanf(qcd_getParam("<nsmear_gauss>",params,params_len),"%d",&nsmear);
   if(myid==0) printf(" Got nsmear_gauss: %d\n",nsmear);
   sscanf(qcd_getParam("<alpha_APE>",params,params_len),"%lf",&alphaAPE);
   if(myid==0) printf(" Got alpha_APE: %lf\n",alphaAPE);
   sscanf(qcd_getParam("<nsmear_APE>",params,params_len),"%d",&nsmearAPE);
   if(myid==0) printf(" Got nsmear_APE: %d\n",nsmearAPE);   
   strcpy(gauge_name,qcd_getParam("<cfg_name>",params,params_len));
   if(myid==0) printf(" Got conf name: %s\n",gauge_name);
   strcpy(out_name,qcd_getParam("<src_block_name>",params,params_len));
   if(myid==0) printf(" Got out name: %s\n",out_name);
  

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
     }

   qcd_real_8 *src = malloc(sizeof(qcd_real_8)*geo.lL[1]*geo.lL[2]*geo.lL[3]);
   if(src == NULL)
     {
       fprintf(stderr, "process %i: malloc returned NULL!\n",geo.myid);
       exit(EXIT_FAILURE);
     }
   if( lx_src[0]==geo.Pos[0] )
     {
       for(int z=0; z<geo.lL[3]; z++)
	 for(int y=0; y<geo.lL[2]; y++)
	   for(int x=0; x<geo.lL[1]; x++)
	     {       
	       int lv = qcd_LEXIC((x_src[0]%geo.lL[0]), x, y, z, geo.lL);
	       double loc_norm = 0;
	       for(int mu=0; mu<4; mu++)
		 for(int c=0; c<3; c++)
		   {
		     loc_norm += qcd_NORMSQUARED(vec.D[lv][mu][c]);
		   }
	       src[x+geo.lL[1]*(y+z*geo.lL[2])] = loc_norm;
	     }
     }
   qcd_destroyVector(&vec);
   /*
    * This assumes only 1 process in time
    * If not, anything may happen
    */

   MPI_Datatype fileview;
   MPI_File fh;
   MPI_Status status;
   int globv3[] = {geo.L[3], geo.L[2], geo.L[1]};
   int locv3[] = {geo.lL[3], geo.lL[2], geo.lL[1]};
   int starts[] = {geo.Pos[3]*locv3[0], geo.Pos[2]*locv3[1], geo.Pos[1]*locv3[2]};

   MPI_Type_create_subarray(3, globv3, locv3, starts, MPI_ORDER_C, MPI_DOUBLE, &fileview);
   MPI_Type_commit(&fileview);
   MPI_File_open(MPI_COMM_WORLD, out_name, MPI_MODE_WRONLY|MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
   MPI_File_set_size(fh, 0);
   MPI_File_set_view(fh, 0, MPI_DOUBLE, fileview, "native", MPI_INFO_NULL);
   if(!qcd_isBigEndian())
     {
       qcd_swap_8(src, (geo.lL[1]*geo.lL[2]*geo.lL[3]));
     }
   MPI_File_write_all(fh, src, (geo.lL[1]*geo.lL[2]*geo.lL[3]), MPI_DOUBLE, &status);
   MPI_File_close(&fh);
   free(src);
   if(nsmear != 0)
     qcd_destroyGaugeField(&uAPE); 

   
   ////////////////////////////////////// CLEAN UP AND EXIT ///////////////////////////////////////////
   //qcd_destroyPropagator(&source);
   qcd_destroyGeometry(&geo);
   MPI_Finalize();
   return(EXIT_SUCCESS);
}//end main
