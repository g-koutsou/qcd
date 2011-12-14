/* b_minus_Dx.c
 *
 * loads source and solution
 * applies Dirac op on solution 
 * checks whether solution was created 
 * correctly from source
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
   char   src_name[qcd_MAX_STRING_LENGTH];
   char   srcname[qcd_MAX_STRING_LENGTH];
   char   sol_name[qcd_MAX_STRING_LENGTH];
   char   out_name[qcd_MAX_STRING_LENGTH];
   char   param_name[qcd_MAX_STRING_LENGTH];
   qcd_uint_4   i;
   qcd_uint_2   nu,c2;
   int   params_len;

   qcd_geometry geo;
   qcd_gaugeField u;

   qcd_vector vec,vec2,vec3,vec4;
   qcd_uint_2 P[4];
   qcd_uint_2 L[4];
   qcd_real_8 theta[4]={M_PI,0.,0.,0.},plaquette; // boundary conditions
   qcd_real_8 kappa, muT,residue, target_eps;
   
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
       
   strcpy(gauge_name,qcd_getParam("<cfg_name>",params,params_len));
   if(myid==0) printf(" Got conf name: %s\n",gauge_name);
   
   strcpy(src_name,qcd_getParam("<src_name>",params,params_len));
   if(myid==0) printf(" Got source name %s\n",src_name);
   
   strcpy(sol_name,qcd_getParam("<sol_name>",params,params_len));
   if(myid==0) printf(" Got solution name %s\n",sol_name);

   sscanf(qcd_getParam("<src_col>",params,params_len),"%hd",&c2);
   if(myid==0) printf(" Got source col=%u\n",c2); 

   sscanf(qcd_getParam("<src_spin>",params,params_len),"%hd",&nu);
   if(myid==0) printf(" Got source spin=%u\n",nu); 

   sscanf(qcd_getParam("<kappa>",params,params_len),"%lf",&kappa);
   if(myid==0) printf(" Got kappa=%f\n",kappa); 

   sscanf(qcd_getParam("<muT>",params,params_len),"%lf",&muT);
   if(myid==0) printf(" Got mu=%f\n",muT);

   sscanf(qcd_getParam("<eps>",params,params_len),"%lf",&target_eps);
   if(myid==0) printf(" Got eps=%e\n",target_eps);
   
   strcpy(out_name,qcd_getParam("<out_name>",params,params_len));
   if(myid==0) printf(" Got output name %s\n",out_name);

   free(params);      
///////////////////////////////////////////////////////////////////////////////////////////////////



   qcd_initVector(&vec,&geo);
   qcd_initVector(&vec2,&geo);
   qcd_initVector(&vec3,&geo);
   qcd_initVector(&vec4,&geo);
   qcd_initGaugeField(&u,&geo);
            
   if(qcd_getGaugeField(gauge_name,qcd_GF_LIME,&u))
   {
      fprintf(stderr,"process %i: Error reading gauge field!\n",myid);
      exit(EXIT_FAILURE);  
   }
   if(myid==0) printf("got gauge-field\n");
   plaquette = qcd_calculatePlaquette(&u);
   if(myid==0) printf("plaquette = %e\n",plaquette);
   if(plaquette> 1 || plaquette<0.3)
      if(myid==0) printf("POSSIBLE ERROR! SUSPICIOUS PLAQUETTE! PLEASE CHECK!\n");
   
   if(qcd_getVector(sol_name,qcd_PROP_HMC,nu,c2,&vec))
     {
       fprintf(stderr,"process %i: Error reading vector!\n",myid);
       exit(EXIT_FAILURE);
     }
   residue = qcd_normVector(&vec);
   if(myid==0) printf("|x|=%e\n", residue);
     
   sprintf(srcname,"%s",src_name);//, nu*3+c2);
   if(qcd_getVector(srcname, qcd_PROP_LIME, 0,0,&vec2))
     {
       fprintf(stderr,"process %i: Error reading vector!\n",myid);
       exit(EXIT_FAILURE);
     } 
   residue = qcd_normVector(&vec2);
   if(myid==0) printf("|b|=%e\n",residue);

   qcd_applyWilsonTMOp(&vec3,&vec,&u, 1.0/(2.0*kappa)-4.0 ,muT);

   residue = qcd_normVector(&vec3);
   if(myid==0) printf("|Dx|=%e\n",residue);  
   
   /* */
   //qcd_scaleVector (&vec3, );
   /* */
   qcd_subVector(&vec4,&vec3,&vec2);
		    
   residue = qcd_normVector(&vec4);
   if(myid==0) printf("|Dx-b|=%e\n",residue);  
   residue /= qcd_normVector(&vec2);

   if(myid==0) printf("residue calculated\n");
   if(myid==0) printf("solution nu=%i c=%i: |Dx-b| / |b|= %e\n",nu,c2, (float) residue);
   
   if(residue>target_eps){
     if(myid==0) printf("POSSIBLE ERROR! HIGH RESIDUE! PLEASE CHECK!\n");    
     return EXIT_FAILURE;
   }

   qcd_writeVector(out_name, qcd_PROP_LIME, 0, 0, &vec);
   
   
   if(myid==0) printf("cleaning up:\n");          


   qcd_destroyVector(&vec);
   qcd_destroyVector(&vec2);
   qcd_destroyVector(&vec3);
   qcd_destroyVector(&vec4);
   qcd_destroyGaugeField(&u);
   qcd_destroyGeometry(&geo);
   MPI_Finalize();         
   return(EXIT_SUCCESS);
}//end main
