/* smear_APE4D.c
 *
 * uses qcd-lib to parallely APE-smear configurations
 * in 4d
 *
 * Christos Kallidonis, Dec 2015
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
   char *params = NULL;                          // needed to read inputfiles				 
   char param_name[qcd_MAX_STRING_LENGTH];       // name of parameter file 
   char gauge_base_name[qcd_MAX_STRING_LENGTH];       // name of gauge-configuration file
   char gauge_name[qcd_MAX_STRING_LENGTH];       // name of gauge-configuration file
   char APE_gauge_name[qcd_MAX_STRING_LENGTH];// name of gauge-configuration output file
   char APE_base_name[qcd_MAX_STRING_LENGTH];// name of gauge-configuration output file
   qcd_real_8 theta[4] = {M_PI,0.0,0.0,0.0};     // antiperiodic b.c. in time
   qcd_uint_2 L[4];
   qcd_uint_2 P[4];          
   qcd_geometry geo;
   qcd_gaugeField u, uAPE;
   qcd_gaugeField *u_ptr, *uAPE_ptr, *utmp_ptr;
   qcd_real_8 plaq;
   qcd_real_8 alphaAPE,alphaAPE_orig;
   int nAPE;
   int trans_alpha;
   int gauge_traj;
             
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
                     	    
   strcpy(gauge_base_name,qcd_getParam("<cfg_name>",params,params_len));
   if(myid==0) printf("Got conf name: %s\n",gauge_base_name);
   sscanf(qcd_getParam("<cfg_traj>",params,params_len),"%d",&gauge_traj);
   if(myid==0) printf("Got conf trajectory: %d\n",gauge_traj);

   sprintf(gauge_name,"%s.%04d",gauge_base_name,gauge_traj);
   if(myid==0) printf("Will read conf from %s\n",gauge_name);

   strcpy(APE_base_name,qcd_getParam("<cfg_smeared_name>",params,params_len));
   if(myid==0) printf("Got smeared conf base name: %s\n",APE_base_name);
              
   sscanf(qcd_getParam("<alpha_APE>",params,params_len),"%lf",&alphaAPE);
   if(myid==0) printf(" Got alpha_APE: %lf\n",alphaAPE);
   alphaAPE_orig = alphaAPE;

   sscanf(qcd_getParam("<nsmear_APE>",params,params_len),"%d",&nAPE);
   if(myid==0) printf(" Got nsmear_APE: %d\n",nAPE);

   sscanf(qcd_getParam("<transfrom_alphaAPE>",params,params_len),"%d",&trans_alpha);
   if(trans_alpha){
     alphaAPE = alphaAPE/(6.0*(1.0-alphaAPE));
     if(myid==0) printf("* Will transform alphaAPE to qcd convention.\n");
   }
   else{
     if(myid==0) printf("* Will NOT transform alphaAPE to qcd convention\n");
   }

   if(myid==0) printf("alphaAPE in use is: %lf\n",alphaAPE);

   free(params);
      
   ///////////////////////////////////////////////////////////////////////////
   
   int j,k;
   j = 0;
   j += qcd_initGaugeField(&u, &geo);
   j += qcd_initGaugeField(&uAPE,&geo);
   MPI_Allreduce(&j, &k, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   if(k>0){
     if(myid==0) fprintf(stderr,"Error, not enough memory\n");
     exit(EXIT_FAILURE);
   }

   if(qcd_getGaugeField(gauge_name,qcd_GF_LIME,&u)){
     fprintf(stderr,"process %i: Error reading gauge field!\n",myid);
     MPI_Finalize();
     exit(EXIT_FAILURE);
   }
   if(myid==0) printf("gauge-field loaded\n");
   plaq = qcd_calculatePlaquette(&u);
   if(myid==0) printf("plaquette = %e\n",plaq);
   
   qcd_communicateGaugePM(&u);
   qcd_waitall(&geo);

   u_ptr = &u;
   uAPE_ptr = &uAPE;
   for(i=0; i<nAPE; i++){
     qcd_apeSmear4d(uAPE_ptr, u_ptr, alphaAPE);
     utmp_ptr=u_ptr; u_ptr=uAPE_ptr; uAPE_ptr=utmp_ptr;
   }
   utmp_ptr=u_ptr; u_ptr=uAPE_ptr; uAPE_ptr=utmp_ptr; //reverse the last swap. Also needed when nsmearAPE=0
   qcd_destroyGaugeField(u_ptr);
   uAPE = *uAPE_ptr;
   
   if(myid==0) printf("gauge-field 4D-APE smeared\n");
   plaq = qcd_calculatePlaquette(&uAPE);
   if(myid==0) printf("new plaquette = %e\n",plaq);
   
   if(qcd_getLimeMessage(gauge_name,&geo,&message)){
     fprintf(stderr,"process %i: Error reading gauge field!\n",myid);
     MPI_Finalize();
     exit(EXIT_FAILURE);     
   }
   

   sprintf(APE_gauge_name,"%s_nAPE%d_alphaAPE0p%02.0f.%04d",APE_base_name,nAPE,alphaAPE_orig*100,gauge_traj);
   if(myid==0) printf("Complete path to the smeared gauge field: %s\n",APE_gauge_name);

   message2=(char*) malloc(strlen(message)+strlen("\n gauge 4D-APE smeared. New plaquette: ")+16+1);
   sprintf(message2,"%s\n gauge 4D-APE smeared. New plaquette: %ef",message,plaq);
   if(myid==0) printf("message written: %s\n",message2);

   if(qcd_writeGaugeField(APE_gauge_name,qcd_GF_LIME,&uAPE,message2)){
     fprintf(stderr,"process %i: Error writing gauge field!\n",myid);
     MPI_Finalize();
     exit(EXIT_FAILURE);  
   }
   
   free(message);
   free(message2);

   //   qcd_destroyGaugeField(uAPE_ptr);
   //   qcd_destroyGaugeField(utmp_ptr);
   //   qcd_destroyGaugeField(&u);
   //   qcd_destroyGaugeField(&uAPE);

   qcd_destroyGeometry(&geo);
   MPI_Finalize();
}
