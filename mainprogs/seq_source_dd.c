/* seq_source_dd.c
 *
 * creates sequential sources
 * for proton-operator-proton correlators with d-quarks in the current
 *
 * N sources are created, with the projectors defined in projectors.h 
 *
 *
 * sample input file for execution on 8 CPUs could contain:
 *
 * <processors_txyz>1 2 2 2</processors_txyz>
 * <lattice_txyz>8 8 8 8</lattice_txyz>
 * <t_sink>3</t_sink>
 * <nsources>2</nsources>
 * <projector_1>13</projector_1>
 * <seq_src_name_1>seq_src_Proj13_list</seq_src_name_1>
 * <projector_2>3</projector_2>
 * <seq_src_name_2>seq_src_Proj3_list</seq_src_name_2>
 * <alpha_gauss>4.0</alpha_gauss>
 * <nsmear_gauss>50</nsmear_gauss>
 * <alpha_APE>0.5</alpha_APE>
 * <nsmear_APE>20</nsmear_APE>
 * <cfg_name>conf88.0000</cfg_name>
 * <propagator_u>sollist_u.0000</propagator_u>
 *
 *
 *
 * Tomasz Korzec 2009
 ****************************************/
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <math.h>
#include <qcd.h>
#include "projectors.h"
 
 
 
 
int main(int argc,char* argv[])
{  
   qcd_uint_4 i,mu,nu,ku,lu,c1,c2,v,x,y,z;   // loop variables
   qcd_uint_4 c3,c1p,c2p,c3p,ctr,ctr2;
   qcd_uint_4 cc1,cc2,a,b;
   qcd_uint_4 isource;                 // ..
   qcd_uint_4 t_sink; 
   int lt_sink;    // sink/source time-slice

   qcd_uint_2 nsources;                // number of different sources

   qcd_uint_4 nsmear, nsmearAPE;       // gaussian and APE smearing: n
   qcd_real_8 alpha, alphaAPE;         // gaussian and APE smearing: alpha
   int params_len;                     // needed to read inputfiles
   char *params = NULL;                       // needed to read inputfiles
   char tmp_string[qcd_MAX_STRING_LENGTH]; // general purpuse
   char param_name[qcd_MAX_STRING_LENGTH];

   char gauge_name[qcd_MAX_STRING_LENGTH]; // name of gauge-config file
   char **source_name;                 // names of output files
   qcd_uint_2 *source_type;            // types of sources
      
   qcd_uint_2 L[4], P[4];              // lattice size and subdivision of lattice   
   qcd_geometry geo;                   // geometry structure   
   qcd_propagator prop_u;              // u-propagator
   qcd_propagator prop_tmp;            // needed when rotating etc.
   qcd_vector vec;                     // needed when smearing
   qcd_gaugeField u;                   // gauge field
   qcd_gaugeField uAPE;                // APE smeared gaugeField
   qcd_gaugeField *u_ptr, *uAPE_ptr, *utmp_ptr;
   
   char prop_u_name[qcd_MAX_STRING_LENGTH]; // file names of up and down quark propagators
   
   qcd_real_8 theta[4] = {M_PI, 0.0, 0.0, 0.0}; // antiperiodic b.c. in time
                
   // C gamma_5 x \bar C gamma_5
   // nonzero entries of [C gamma_5]_ab [\bar C gamma_5]_cd are in
   // a=cg5cg5_ind[..][0], b=cg5cg5_ind[..][1], c=cg5cg5_ind[..][2], d=cg5cg5_ind[..][3]
   qcd_uint_2     cg5cg5_ind[16][4];
   qcd_complex_16 cg5cg5_val[16];
   
   qcd_complex_16 C, factor;              
   qcd_real_8 plaq;
   
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
   
                       
   sscanf(qcd_getParam("<t_sink>",params,params_len),"%d",&t_sink);
   if(myid==0) printf("process %i: Got sink time slice: %d\n",myid, t_sink);
   
   if(t_sink >= L[0])
     {
       if(myid==0) fprintf(stderr, " Error: t_sink (=%d) >= L[0] (=%d),\n t_sink should be in [0, L[0]), did you forget to mod(t_sink, L[0]) ?\n", t_sink, L[0]);
       exit(EXIT_FAILURE);
     }

//   sscanf(qcd_getParam("<t_src>",params,params_len),"%d",&t_src);
//   if(myid==0) printf("process %i: Got source time slice: %d\n",myid, t_src);
      
   sscanf(qcd_getParam("<nsources>",params,params_len),"%hu",&nsources);
   if(myid==0) printf("process %i: Got number of sources: %d\n",myid, nsources);
   
   source_name= malloc(nsources*sizeof(*source_name));
   source_type= malloc(nsources*sizeof(qcd_uint_2));
   
   for(i=0; i<nsources; i++)
   {
      source_name[i]= malloc(qcd_MAX_STRING_LENGTH*sizeof(char));
      sprintf(tmp_string,"<projector_%d>",i+1);
      sscanf(qcd_getParam(tmp_string,params,params_len),"%hu",&source_type[i]);
      if(myid==0) printf("Sequential source #%d type: Proj%d\n",i, source_type[i]);
      sprintf(tmp_string,"<seq_src_name_%d>",i+1);
      strcpy(source_name[i],qcd_getParam(tmp_string,params,params_len));
      if(myid==0) printf("Sequential source #%d name: %s\n",i, source_name[i]);
   }
   
   sscanf(qcd_getParam("<alpha_gauss>",params,params_len),"%lf",&alpha);
   if(myid==0) printf(" Got alpha_gauss: %lf\n",alpha);
   sscanf(qcd_getParam("<nsmear_gauss>",params,params_len),"%d",&nsmear);
   if(myid==0) printf(" Got nsmear_gauss: %d\n",nsmear);
   sscanf(qcd_getParam("<alpha_APE>",params,params_len),"%lf",&alphaAPE);
   if(myid==0) printf(" Got alpha_APE: %lf\n",alphaAPE);
   sscanf(qcd_getParam("<nsmear_APE>",params,params_len),"%d",&nsmearAPE);
   if(myid==0) printf(" Got nsmear_APE: %d\n",nsmearAPE);   
   strcpy(gauge_name,qcd_getParam("<cfg_name>",params,params_len));
   if(myid==0) printf(" Got conf name: %s\n",gauge_name);
  
  
   strcpy(prop_u_name,qcd_getParam("<propagator_u>",params,params_len));
   if(myid==0) printf(" Got u-propagator name: %s\n",prop_u_name);

   free(params);


///////////////////////////////////////////////////////////////////////////////////////////////////

   /* load gauge field and APE-smear it */
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
  
   /* initialize propagators and vectors */ 
   qcd_initPropagator(&prop_u, &geo);
   qcd_initPropagator(&prop_tmp, &geo);
   qcd_initVector(&vec, &geo);

   
   /* load propagators */
   if(qcd_getPropagator(prop_u_name,qcd_PROP_LIME, &prop_u)) exit(EXIT_FAILURE);
   if(myid==0) printf("propagator loaded\n");    
      
          
   //################################################################################
   // smear the propagators on the sink side
   
   for(mu=0;mu<4;mu++)
   for(c1=0;c1<3;c1++)
   {
      qcd_copyVectorPropagator(&vec,&prop_u,mu,c1);
      for(i=0; i<nsmear; i++)
      {
         if(qcd_gaussIteration3d(&vec,&uAPE,alpha,t_sink))
         {
            fprintf(stderr,"process %i: Error while smearing!\n",geo.myid);
            exit(EXIT_FAILURE);
         }
      }
      qcd_copyPropagatorVector(&prop_u,&vec,mu,c1);
   }   
   if(myid==0) printf("propagator smeared\n");
   
   
   
   
   // tabulate non-vanishing entries of C*gamma5 x \bar C*gamma5
   ctr = 0;
   for(mu=0;mu<4;mu++) 
   for(nu=0;nu<4;nu++)
   for(ku=0;ku<4;ku++)
   for(lu=0;lu<4;lu++)
   {
      C = qcd_CMUL(qcd_CGAMMA[5][mu][nu],qcd_BAR_CGAMMA[5][ku][lu]);
      if(qcd_NORM(C)>1e-3)
      {
         cg5cg5_val[ctr] = C;
         cg5cg5_ind[ctr][0] = mu;
         cg5cg5_ind[ctr][1] = nu;
         cg5cg5_ind[ctr][2] = ku;
         cg5cg5_ind[ctr][3] = lu;                                        
         ctr++;
      }
   }
   
   for(isource=0; isource<nsources; isource++)
   {
      qcd_zeroPropagator(&prop_tmp);
   
               
      if(myid==0) printf("creating source with Proj%i.\n",source_type[isource]);
      lt_sink = t_sink - geo.Pos[0]*geo.lL[0]; // local t_sink
    
      // construct M^{c3 c3p}_{nu lu} with ubar u current and protons   
      if((lt_sink>=0) && (lt_sink<geo.lL[0]))      // otherwise this CPU has nothing to compute
      for(cc1=0;cc1<6;cc1++)
      {
         c1=qcd_EPS[cc1][0];
         c2=qcd_EPS[cc1][1];
         c3=qcd_EPS[cc1][2];
         for(cc2=0;cc2<6;cc2++)
         {     
            c1p=qcd_EPS[cc2][0];
            c2p=qcd_EPS[cc2][1];
            c3p=qcd_EPS[cc2][2];
            for(ctr2=0;ctr2<ctr;ctr2++)
            {
               mu = cg5cg5_ind[ctr2][0];
               nu = cg5cg5_ind[ctr2][1];
               ku = cg5cg5_ind[ctr2][2];
               lu = cg5cg5_ind[ctr2][3];
               for(a=0;a<4;a++)
               for(b=0;b<4;b++)
               if(qcd_NORM(PROJECTOR[source_type[isource]][a][b])>1e-3)
               { 
                    factor = qcd_CMUL(qcd_CSCALE(PROJECTOR[source_type[isource]][a][b],-qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),
                                    cg5cg5_val[ctr2]);
                  
                  for(z=0; z<geo.lL[3]; z++)
                  for(y=0; y<geo.lL[2]; y++)
                  for(x=0; x<geo.lL[1]; x++)
                  {
                     v=qcd_LEXIC(lt_sink,x,y,z,geo.lL);
                     prop_tmp.D[v][nu][lu][c3][c3p] = qcd_CADD(prop_tmp.D[v][nu][lu][c3][c3p],
                                                               qcd_CMUL(factor,
                                                                        qcd_CADD(qcd_CMUL(prop_u.D[v][mu][b][c1][c1p],
                                                                                          prop_u.D[v][a][ku][c2][c2p]),
                                                                                 qcd_CMUL(prop_u.D[v][mu][ku][c1][c1p],
                                                                                          prop_u.D[v][a][b][c2][c2p]))));
                  }//end volume loop
               }//end a/b loop
            }//end mu/gu/ku/ju loop
         }//end second color loop       
      }//end first color loop
             
   
      qcd_conjPropagator(&prop_tmp);
      qcd_gamma5Propagator(&prop_tmp);             
   
      for(mu=0;mu<4;mu++)
      for(c1=0;c1<3;c1++)
      {
         qcd_copyVectorPropagator(&vec,&prop_tmp,mu,c1);
         for(i=0; i<nsmear; i++)
         {
            if(qcd_gaussIteration3d(&vec,&uAPE,alpha,t_sink))
            {
               fprintf(stderr,"process %i: Error while smearing!\n",geo.myid);
               exit(EXIT_FAILURE);
            }
         }
         qcd_copyPropagatorVector(&prop_tmp,&vec,mu,c1);
	 char tmp_name[qcd_MAX_STRING_LENGTH];
	 sprintf(tmp_name, "%s.%00005d", source_name[isource], c1 + mu*3);
	 qcd_writeVector(tmp_name, qcd_PROP_LIME, mu, c1, &vec);
      }   
      if(myid==0) printf("seq. source smeared\n");
   
      //qcd_writePropagator(source_name[isource], qcd_PROP_CMI, &prop_tmp);
   }//end isource-loop
   
   
   
   

   
   //#################################################################################
   // clean up

   qcd_destroyPropagator(&prop_u);
   qcd_destroyPropagator(&prop_tmp);
   qcd_destroyGaugeField(&uAPE);
   qcd_destroyVector(&vec);
   for(i=0; i<nsources; i++)
      free(source_name[i]);
   free(source_name);
   free(source_type);
   qcd_destroyGeometry(&geo);
   
   MPI_Finalize();      
}//end main
