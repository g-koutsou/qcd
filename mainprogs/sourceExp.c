/* sourceExp.c
 *
 * creates sources with the phase Exp[I p x]
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
   char   param_name[qcd_MAX_STRING_LENGTH];
   char   out_name[qcd_MAX_STRING_LENGTH];
   qcd_uint_4   mom[4],i,j,t,isource;
   qcd_uint_2   mu,nu,col,c1,c2,s,lt,lx,ly,lz;
   qcd_uint_2   gt,gx,gy,gz;
   qcd_int_4   params_len;

   qcd_geometry geo;
   qcd_propagator prop;
   qcd_uint_2 P[4];
   qcd_uint_2 L[4];
   qcd_real_8 theta[4]={M_PI,0.,0.,0.}; // boundary conditions
   qcd_complex_16 expo;
   qcd_vector vec;                     

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
     
   sscanf(qcd_getParam("<momentum_txyz>",params,params_len),"%d %d %d %d",&mom[0], &mom[1], &mom[2], &mom[3]);
   if(myid==0) printf(" Got momentum coords: %d %d %d %d\n",mom[0],mom[1],mom[2],mom[3]);
  
   //   strcpy(gauge_name,qcd_getParam("<cfg_name>",params,params_len));
   //   if(myid==0) printf(" Got conf name: %s\n",gauge_name);
   strcpy(out_name,qcd_getParam("<source>",params,params_len));
   if(myid==0) printf(" Got source name %s\n",out_name);

   

   free(params);      
///////////////////////////////////////////////////////////////////////////////////////////////////



   qcd_initPropagator(&prop,&geo);
  
   //   for(i=0; i<4; i++)
   //      lx_src[i] = x_src[i]-geo.Pos[i]*geo.lL[i];  //source_pos in local lattice
        
      qcd_zeroPropagator(&prop); 

   for(mu=0;mu<4;mu++)
   for(col=0;col<3;col++)
   {
     for(lt=0; lt<geo.lL[0]; lt++)
     for(lx=0; lx<geo.lL[1]; lx++)
     for(ly=0; ly<geo.lL[2]; ly++)
     for(lz=0; lz<geo.lL[3]; lz++)
       {
	 gt = geo.Pos[0] * geo.lL[0] + lt ;
	 gx = geo.Pos[1] * geo.lL[1] + lx ;
	 gy = geo.Pos[2] * geo.lL[2] + ly ;
	 gz = geo.Pos[3] * geo.lL[3] + lz ;

         expo = (qcd_complex_16){cos(2*M_PI*mom[0]*gt/geo.L[0]+2*M_PI*mom[1]*gx/geo.L[1]+2*M_PI*mom[2]*gy/geo.L[2]+2*M_PI*mom[3]*gz/geo.L[3]),sin(2*M_PI*mom[0]*gt/geo.L[0]+2*M_PI*mom[1]*gx/geo.L[1]+2*M_PI*mom[2]*gy/geo.L[2]+2*M_PI*mom[3]*gz/geo.L[3])};

         prop.D[qcd_LEXIC(lt,lx,ly,lz,geo.lL)][mu][mu][col][col]=expo;

   }
      if(myid==0) printf("Propagator: col = %1d, spin = %1d\n",col,mu);  
   }

   //
   //   WRITING IN THE OLD FORMAT 
   //   qcd_writePropagator(out_name, qcd_PROP_CMI, &prop);
   //
 
   qcd_initVector(&vec,&geo);

   char tmp_name[qcd_MAX_STRING_LENGTH];
   for(mu=0;mu<4;mu++)
   for(c1=0;c1<3;c1++)
   {
    sprintf(tmp_name, "%s.%00005d", out_name, c1 + mu*3);
    qcd_copyVectorPropagator(&vec, &prop, mu, c1);
    //    qcd_writeVector(tmp_name, qcd_PROP_LIME, &vec);
    qcd_writeVector(tmp_name, qcd_PROP_LIME, mu, c1, &vec);
   }   
   
   ////////////////////////////////////// CLEAN UP AND EXIT ///////////////////////////////////////////
   qcd_destroyPropagator(&prop);
   qcd_destroyGeometry(&geo);
   qcd_destroyVector(&vec);
   MPI_Finalize();
   return(EXIT_SUCCESS);
}//end main
