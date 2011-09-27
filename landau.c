/* landau.c
 *
 * uses qcd-lib to parallely gauge fix configurations
 * to landau gauge
 *
 * Tomasz Korzec 2009
 ********************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <qcd.h>





int main(int argc,char* argv[])
{
   qcd_uint_4 i;
   int myid,numprocs, namelen;    
   char processor_name[MPI_MAX_PROCESSOR_NAME];
   int params_len;                               // needed to read inputfiles
   char *params;                                 // needed to read inputfiles				 
	char param_name[qcd_MAX_STRING_LENGTH];       // name of parameter file 
   qcd_real_8 pover; 		                      // stochastic over-relaxation probability	 
   char gauge_name[qcd_MAX_STRING_LENGTH];       // name of gauge-configuration file
   char landau_gauge_name[qcd_MAX_STRING_LENGTH];// name of gauge-configuration output file
   qcd_real_8 theta[4] = {M_PI,0.0,0.0,0.0};     // antiperiodic b.c. in time
   qcd_uint_2 L[4];
   qcd_uint_2 P[4];          
   qcd_geometry geo;
   qcd_gaugeField u, lu;
   qcd_real_8 plaquette;
             
   //set up MPI
   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD,&numprocs);         // num. of processes taking part in the calculation
   MPI_Comm_rank(MPI_COMM_WORLD,&myid);             // each process gets its ID
   MPI_Get_processor_name(processor_name,&namelen); //   
   char *message;
   char *message2;
           
   
   //////////////////// READ INPUT FILE /////////////////////////////////////////////
            
   if(argc!=2)
   {
      if(myid==0) fprintf(stderr,"No input file specified\n");
      MPI_Finalize();
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
   if(i==1)
   {
      MPI_Finalize();
      exit(EXIT_FAILURE);
   }
   MPI_Bcast(&params_len, 1, MPI_INT, 0, MPI_COMM_WORLD);
   if(myid!=0) params = (char*) malloc(params_len*sizeof(char));
   MPI_Bcast(params, params_len, MPI_CHAR, 0, MPI_COMM_WORLD);
   
   sscanf(qcd_getParam("<processors_txyz>",params,params_len),"%hd %hd %hd %hd",&P[0], &P[1], &P[2], &P[3]);
   sscanf(qcd_getParam("<lattice_txyz>",params,params_len),"%hd %hd %hd %hd",&L[0], &L[1], &L[2], &L[3]);
   if(qcd_initGeometry(&geo,L,P, theta, myid, numprocs))
   {
      MPI_Finalize();
      exit(EXIT_FAILURE);
   }   
   
   if(myid==0) printf(" Local lattice: %i x %i x %i x %i\n",geo.lL[0],geo.lL[1],geo.lL[2],geo.lL[3]);
      
   sscanf(qcd_getParam("<over_relaxation>",params,params_len),"%lf",&pover);
   if(myid==0) printf("Got over-relaxation parameter: %e\n",pover);
                     	    
   strcpy(gauge_name,qcd_getParam("<cfg_name>",params,params_len));
   if(myid==0) printf("Got conf name: %s\n",gauge_name);

   strcpy(landau_gauge_name,qcd_getParam("<cfg_output_name>",params,params_len));
   if(myid==0) printf("Got conf output name: %s\n",landau_gauge_name);
              
   free(params);
   
   
   ///////////////////////////////////////////////////////////////////////////
   
   qcd_initGaugeField(&u,&geo);
   qcd_initGaugeField(&lu,&geo);
   if(qcd_getGaugeField(gauge_name,qcd_GF_LIME,&u))
   {
      fprintf(stderr,"process %i: Error reading gauge field!\n",myid);
      MPI_Finalize();
      exit(EXIT_FAILURE);  
   }
   plaquette = qcd_calculatePlaquette(&u);
   if(myid==0) printf("Gauge field loaded, plaquette = %f\n",plaquette);
   
   qcd_landauGauge(&lu,&u,pover);
   plaquette = qcd_calculatePlaquette(&lu);
   if(myid==0) printf("Fixed to landau gauge, plaquette = %f\n",plaquette);
   
   if(qcd_getLimeMessage(gauge_name,&geo,&message))
   {
      fprintf(stderr,"process %i: Error reading gauge field!\n",myid);
      MPI_Finalize();
      exit(EXIT_FAILURE);     
   }
   message2=(char*) malloc(strlen(message)+strlen("\n gauge transformed to Landau gauge")+1);
   sprintf(message2,"%s\n gauge transformed to Landau gauge",message);
   if(myid==0) printf("message written: %s\n",message2);
   if(qcd_writeGaugeField(landau_gauge_name,qcd_GF_LIME,&lu,message2))
   {
      fprintf(stderr,"process %i: Error writing gauge field!\n",myid);
      MPI_Finalize();
      exit(EXIT_FAILURE);  
   }
   
   free(message);
   free(message2);
   qcd_destroyGaugeField(&u);
   qcd_destroyGaugeField(&lu);
   qcd_destroyGeometry(&geo);
   MPI_Finalize();
}
