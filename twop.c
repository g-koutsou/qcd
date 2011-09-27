/* twop.c
 *
 * reads forward propagators
 * and creates nucleon two point functions
 *
 * Tomasz Korzec 2010
 ****************************************/
 
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
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
   qcd_uint_4 x_src[4];                       // source and sink coordinates
   qcd_uint_4 t_sink, t_start, t_stop, t,lt;
   qcd_real_8 tmp;                            // general purpuse
   FILE *fp_momlist;
  
   FILE *fp_corr_p;                           // output file
  
   int params_len;                            // needed to read inputfiles
   char *params;                              // needed to read inputfiles

   char gauge_name[qcd_MAX_STRING_LENGTH];      // name of gauge-configuration file
   char corr_p_name[qcd_MAX_STRING_LENGTH];     // name of output file proton 2pt function
   
   char param_name[qcd_MAX_STRING_LENGTH];      // name of parameter file  
   char momlist_name[qcd_MAX_STRING_LENGTH];    // name of momenta-list file
   char uprop_name[qcd_MAX_STRING_LENGTH];      // file names of up and down quark propagators
   char dprop_name[qcd_MAX_STRING_LENGTH];      

   qcd_geometry geo;                            // geometry structure
   qcd_propagator uprop;                        // propagator
   qcd_propagator dprop;                        // propagator
   qcd_vector vec;                              // needed when smearing
   qcd_gaugeField u;                            // gauge field 
   qcd_gaugeField uAPE;                         // APE smeared gaugeField
   qcd_gaugeField *u_ptr, *uAPE_ptr, *utmp_ptr;

   qcd_uint_4 nsmear, nsmearAPE;       // gaussian and APE smearing: n
   qcd_real_8 alpha, alphaAPE;         // gaussian and APE smearing: alpha

   qcd_real_8 theta[4] = {M_PI,0.0,0.0,0.0};    // antiperiodic b.c. in time
   qcd_uint_2 L[4];
   qcd_uint_2 P[4];
   qcd_complex_16 phase_factor;         
   qcd_complex_16 z1, z2;                       // temp variables
   qcd_complex_16 C, C2;   
   qcd_complex_16 corr, corr2;
   qcd_real_8 plaq;
   qcd_int_4 ctr, ctr2;
   qcd_int_2 cg5cg5_ind[16*16][4];
   qcd_complex_16 cg5cg5_val[16*16];
   
   qcd_complex_16 *block;                       // to store the block (2pt function before FT)

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
   if(myid==0) printf("Got sink time slices: %d ... %d\n",t_start,t_stop);
                     	  
   sscanf(qcd_getParam("<source_pos_txyz>",params,params_len),"%d %d %d %d",&x_src[0],&x_src[1],&x_src[2],&x_src[3]);
   if(myid==0) printf("Got source coords: %d %d %d %d\n",x_src[0],x_src[1],x_src[2],x_src[3]);
     
   strcpy(uprop_name,qcd_getParam("<propagator_u>",params,params_len));
   if(myid==0) printf("Got propagator file name: %s\n",uprop_name);
   strcpy(dprop_name,qcd_getParam("<propagator_d>",params,params_len));
   if(myid==0) printf("Got propagator file name: %s\n",dprop_name);
   
   strcpy(gauge_name,qcd_getParam("<cfg_name>",params,params_len));
   if(myid==0) printf("Got conf name: %s\n",gauge_name);
          
          
          
   strcpy(corr_p_name,qcd_getParam("<corr_name>",params,params_len));
   if(myid==0) printf("Got output file name: %s\n",corr_p_name);
           
   strcpy(momlist_name,qcd_getParam("<momenta_list>",params,params_len));
   if(myid==0) printf("Got momenta-list file name: %s\n",momlist_name);
    
   sscanf(qcd_getParam("<alpha_gauss>",params,params_len),"%lf",&alpha);
   if(myid==0) printf(" Got alpha_gauss: %lf\n",alpha);
   sscanf(qcd_getParam("<nsmear_gauss>",params,params_len),"%d",&nsmear);
   if(myid==0) printf(" Got nsmear_gauss: %d\n",nsmear);
   sscanf(qcd_getParam("<alpha_APE>",params,params_len),"%lf",&alphaAPE);
   if(myid==0) printf(" Got alpha_APE: %lf\n",alphaAPE);
   sscanf(qcd_getParam("<nsmear_APE>",params,params_len),"%d",&nsmearAPE);
   if(myid==0) printf(" Got nsmear_APE: %d\n",nsmearAPE); 
   free(params);



         
   //#####################################################################   
   // allocate memory
   // load gauge-field and APE-smear it
   j=0;
   j += qcd_initGaugeField(&u, &geo);
   j += qcd_initGaugeField(&uAPE,&geo);
   MPI_Allreduce(&j, &k, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   if(k>0)
   {
      if(myid==0) fprintf(stderr,"Error, not enough memory\n");
      exit(EXIT_FAILURE);
   }
      
   if(qcd_getGaugeField(gauge_name,qcd_GF_LIME,&u)) 
   {
      if(myid==0) fprintf(stderr,"Error reading gauge field\n");
      exit(EXIT_FAILURE);
   }
   if(myid==0) printf("gauge-field loaded\n");   
   plaq = qcd_calculatePlaquette(&u);
   if(myid==0) printf("plaquette = %e\n",plaq);
   qcd_communicateGaugePM(&u);
   qcd_waitall(&geo);
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
  
   j = 0;
   j += qcd_initPropagator(&uprop, &geo);
   j += qcd_initPropagator(&dprop, &geo);
   j += qcd_initVector(&vec, &geo);
   
   block = (qcd_complex_16*) malloc(geo.lV3*sizeof(qcd_complex_16));
            
   MPI_Allreduce(&j, &k, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   if(k>0)
   {
      if(myid==0) printf("not enough memory\n");
      exit(EXIT_FAILURE);
   }
   if(myid==0) printf("memory for propagators and gauge-field allocated\n");
         
   
   //##############################################################################
   
   
   
   // load propagators
   
   if(qcd_getPropagator(uprop_name,qcd_PROP_LIME, &uprop)) exit(EXIT_FAILURE);
   if(myid==0) printf("up propagator loaded\n");
   if(qcd_getPropagator(dprop_name,qcd_PROP_LIME, &dprop)) exit(EXIT_FAILURE);
   if(myid==0) printf("down propagator loaded\n");   

   //################################################################################
   // transform propagators to basis with theta-periodic boundaries in the temporal direction
   for(lt=0; lt<geo.lL[0]; lt++)
   {
      t = lt + geo.Pos[0] * geo.lL[0];
      phase_factor   = (qcd_complex_16) {cos(theta[0]*t/geo.L[0]),sin(theta[0]*t/geo.L[0])};
      qcd_mulPropagatorC3d(&uprop, phase_factor, (t+x_src[0]) % geo.L[0]);
      qcd_mulPropagatorC3d(&dprop, phase_factor, (t+x_src[0]) % geo.L[0]);
   }
   if(myid==0) printf("propagators transformed to basis with theta-periodic boundary conditions\n");
   
   
   //################################################################################
   // gaussian smearing of propagators (only the time-slices that will be used)
   for(mu=0;mu<4;mu++)
   for(c1=0;c1<3;c1++)
   {
      qcd_copyVectorPropagator(&vec,&uprop,mu,c1);
      for(i=0; i<nsmear; i++)
      {
         if(qcd_gaussIteration3dAll(&vec,&uAPE,alpha,i==0))
         {
            fprintf(stderr,"process %i: Error while smearing!\n",geo.myid);
            exit(EXIT_FAILURE);
         }
      }
      qcd_copyPropagatorVector(&uprop,&vec,mu,c1);
      qcd_copyVectorPropagator(&vec,&dprop,mu,c1);
      for(i=0; i<nsmear; i++)
      {
         if(qcd_gaussIteration3dAll(&vec,&uAPE,alpha,i==0))
         {
            fprintf(stderr,"process %i: Error while smearing!\n",geo.myid);
            exit(EXIT_FAILURE);
         }
      }
      qcd_copyPropagatorVector(&dprop,&vec,mu,c1);
   }
   qcd_destroyGaugeField(&uAPE);
   qcd_destroyVector(&vec);   
   if(myid==0) printf("propagators smeared\n");
         
   
   //################################################################################   
   // calculate the 'blocks'      
   
   ctr = 0;
   for(mu=0;mu<4;mu++)
   for(nu=0;nu<4;nu++)
   for(ku=0;ku<4;ku++)
   for(lu=0;lu<4;lu++)
   {
      C = qcd_CMUL(qcd_CGAMMA[5][mu][nu],qcd_BAR_CGAMMA[5][ku][lu]);
      if(qcd_NORM(C)>1e-3)
      {
         cg5cg5_val[ctr].re = C.re;
         cg5cg5_val[ctr].im = C.im;
         cg5cg5_ind[ctr][0] = mu;
         cg5cg5_ind[ctr][1] = nu;
         cg5cg5_ind[ctr][2] = ku;
         cg5cg5_ind[ctr][3] = lu;                                                            
         ctr++;
      }
   }
      
      
   //open output file to write in
   if(myid==0)
   {
      fp_corr_p = fopen(corr_p_name,"w");   
      if(fp_corr_p==NULL)
      {
         printf("failed to open %s for writing\n",corr_p_name);
         k=1;
      }
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
   if(myid==0) fscanf(fp_momlist,"%i\n",&i);
   MPI_Bcast(&i,1,MPI_INT, 0, MPI_COMM_WORLD);
   if(myid==0) printf("will read %i momenta combinations\n",i);

   mom = malloc(i*3*sizeof(qcd_int_4));

   if(myid==0)
   {
      for(j=0; j<i; j++)
      {
         fscanf(fp_momlist,"%i %i %i\n",&(mom[j][0]),&(mom[j][1]),&(mom[j][2]));
         //printf("got combination %i %i %i\n",mom[j][0],mom[j][1],mom[j][2]);  
      }
      fclose(fp_momlist);   
   }
   MPI_Bcast(&(mom[0][0]),i*3,MPI_INT,0, MPI_COMM_WORLD);
   if(myid==0) printf("momenta list read and broadcasted\n");   
     
             
   for(t=t_start; t<=t_stop; t++)
   {
      lt = ((t+x_src[0])%geo.L[0]) - geo.Pos[0]*geo.lL[0];
      if(lt>=0 && lt<geo.lL[0])  //inside the local lattice, otherwise nothing to calculate
      {
         if(myid==0) printf("t=%i\n",t);
         for(v3=0; v3<geo.lV3; v3++)   //set blocks to zero
            block[v3]= (qcd_complex_16) {0,0};

         for(al=0; al<4; al++)
         for(be=0; be<4; be++)
         if(qcd_NORM(PROJECTOR[13][al][be])>0.0001)
         {         
            /* alpha-beta-component*/
	   if(myid==0) printf("process %i: calculating %i-%i component\n",myid,al,be);
            for(ctr2=0; ctr2<ctr; ctr2++)
            {          
               mu=cg5cg5_ind[ctr2][0];
               nu=cg5cg5_ind[ctr2][1];
               ku=cg5cg5_ind[ctr2][2];
               lu=cg5cg5_ind[ctr2][3];                              
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
                     for(lx=0; lx<geo.lL[1]; lx++)
                     for(ly=0; ly<geo.lL[2]; ly++)
                     for(lz=0; lz<geo.lL[3]; lz++)
                     {
                        v3 = qcd_LEXIC0(lx,ly,lz,geo.lL);
                        v =  qcd_LEXIC(lt,lx,ly,lz,geo.lL);
                        block[v3] = qcd_CSUB(block[v3]
                                            ,qcd_CMUL(PROJECTOR[13][al][be] 
                                                     ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(cg5cg5_val[ctr2]
                                                                                           ,uprop.D[v][mu][ku][c1][c1p])
                                                                                  ,dprop.D[v][nu][lu][c2][c2p])
                                                                         ,uprop.D[v][be][al][c3][c3p])
                                                                ,qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2])));
                        block[v3] = qcd_CSUB(block[v3]
                                            ,qcd_CMUL(PROJECTOR[13][al][be]
                                                     ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(cg5cg5_val[ctr2]
                                                                                           ,uprop.D[v][mu][al][c1][c1p])
                                                                                  ,dprop.D[v][nu][lu][c2][c2p])
                                                                         ,uprop.D[v][be][ku][c3][c3p])
                                                    ,qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2])));

                     }//space loop
                  }//color2 loop    
               }//color1 loop
            }//nonvanishing cg5cg5 loop
         }//nonvanishing projector condition
         
         //Fourier transform time-slice
         
         for(j=0; j<i; j++)
         {
            if(myid==0) 
            {
               fprintf(fp_corr_p,"%i %+i %+i %+i ",t,mom[j][0],mom[j][1],mom[j][2]);
            }   
            corr = (qcd_complex_16) {0,0};
            
            for(lx=0; lx<geo.lL[1]; lx++)
            for(ly=0; ly<geo.lL[2]; ly++)
            for(lz=0; lz<geo.lL[3]; lz++)
            {
               v3 = qcd_LEXIC0(lx,ly,lz,geo.lL);
               x=lx+geo.Pos[1]*geo.lL[1] - x_src[1];
               y=ly+geo.Pos[2]*geo.lL[2] - x_src[2];
               z=lz+geo.Pos[3]*geo.lL[3] - x_src[3];
               tmp = (((double) mom[j][0]*x)/geo.L[1] + ((double) mom[j][1]*y)/geo.L[2] + ((double) mom[j][2]*z)/geo.L[3])*2*M_PI;
               C2=(qcd_complex_16) {cos(tmp), -sin(tmp)}; //TABULATE FOR LARGE SPEEDUP!!!
               corr=qcd_CADD(corr, qcd_CMUL(block[v3],C2));
            }
            //printf("process %i: corr = %f %+fi\n",myid,0.5*corr.re,0.5*corr.im);
            MPI_Reduce(&(corr.re), &(corr2.re), 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            if(myid==0) 
            {
               fprintf(fp_corr_p,"%+e %+e\n",corr2.re*0.5,corr2.im*0.5);
            }
         }
         
      }//end lt inside local block condition
      

   }//end t-loop   
         
   
   if(myid==0)
   {
      fclose(fp_corr_p);
   }   
   
                  
   
   //#####################################################################   
   // clean up
   if(myid==0) printf("cleaning up...\n");
   
   free(block);
   free(mom);
   qcd_destroyPropagator(&uprop);
   qcd_destroyPropagator(&dprop);
   qcd_destroyGeometry(&geo);
   MPI_Finalize();
}//end main
