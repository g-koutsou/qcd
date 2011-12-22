/* Proton nucleon three point function using all to all Propagator using TSM
 * This program takes the forward propagators                                                    
 * And the stochastic sources with the smeared stochastic solution vector
 * Then we calculate the three point function in momentum space
 * Kyriakos Hadjiyiannakou 2011                                                                 
 ****************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <qcd.h>
#include "projectors.h"
#include <qcd_gamma_up.h>

#define Nmu 4
#define Nic 3
#define MAX_NUM_PHI 1000
#define Noperators 16
int main(int argc, char* argv[]){

  qcd_uint_2 mu,nu,ku,lu,lambda,xita,beta,rho,phita, delta, sigma,zita, kappa; // dirac indices
  qcd_uint_2 c1,c2,c3,c4,c5,c6,tt;      // color indices
  qcd_uint_2 cc1,cc2;                   // two more color indices
  qcd_uint_2 i,j,k,l,m,it;                  // for temporal use only
  qcd_uint_4 num_momenta;                // store the number of total momenta
  qcd_int_4 kk,proj_index;               // this index takes the gamma matrix fo insertion operator
  FILE *fp_momlist;                     // file pointer to the list of momenta
  FILE *fp_threep;                      // file pointer to the output file where i will store the thee point functions
  FILE *fp_test;
  int params_len;                       // the string legth of the parameters
  char *params;                         // pointer to the params
  qcd_uint_2 aa,bb,cc,ff,gg,hh;         // also color indices
  qcd_int_4 imom;                       // index for momenta

  char gauge_name[qcd_MAX_STRING_LENGTH];       // store the path of gauge configuration
  char param_name[qcd_MAX_STRING_LENGTH];      // name of parameter file                        
  char momlist_name[qcd_MAX_STRING_LENGTH];    // name of momenta-list file                     
  char uprop_name[qcd_MAX_STRING_LENGTH];      // file with the filenames of up propagators   
  char dprop_name[qcd_MAX_STRING_LENGTH];      // file with the filenames of down propagators
  char output_name[qcd_MAX_STRING_LENGTH];     // name of output file

  char phi_name[qcd_MAX_STRING_LENGTH];        // file with the filenames of phi low precision solution vector
  char xi_1_name[qcd_MAX_STRING_LENGTH];       // file list of the noise vector that gives the low precision solution vector
 
  FILE *fp_phi_name,*fp_xi_1_name; // file pointer to the different input files

  char phi_files[MAX_NUM_PHI][qcd_MAX_STRING_LENGTH];        // for each number of phi solution vector we store a string for the path
  char xi_1_files[MAX_NUM_PHI][qcd_MAX_STRING_LENGTH];       // same as xi_1

  qcd_uint_4 num_phi, num_xi_1; // the num_phi == num_xi_1, num_chi == num_xi_2 check this must be equal 
  qcd_uint_4  t,lt;        // i will use only t for global time and lt for local time
  qcd_complex_16 phase_factor,sum;                     // complex number where we store the phase for fourier transform
  qcd_uint_4 v3,v;                                 // variables for v3 = lx*ly*lz , v = lx*ly*lz*lt
  qcd_uint_2 L[4];                                 // store global lattice dimensions
  qcd_uint_2 P[4];                                 // store the number of proccesors in each direction
  qcd_uint_4 x_src[4],t_noise;                       // x_src the coordinates of source and t_c time coordinate for insertion
  qcd_propagator uprop;                            // for up propagator
  qcd_propagator dprop;                           // for down propagator
  qcd_propagator uprop_unsmear;
  qcd_int_4 x,y,z;                                 // indices for spatial coordinates
  qcd_real_8 tmp;                                  // temporary complex number
  qcd_uint_4 lx,ly,lz;                             // local lattice variables
  qcd_vector phi,xi_1;                    // solution and noice vectors

  qcd_geometry geo;                                // a variable to start geometry

  qcd_int_4 (*mom)[3];                             // a list where we will store the momenta
  qcd_real_8 theta[4] = {M_PI,0.0,0.0,0.0};    // antiperiodic b.c. in time take back normal fermion fields                                                                                    
  qcd_real_8 plaq;                             // calculate the plaquette                    
  int myid,numprocs, namelen;                  // for mpi use
  char processor_name[MPI_MAX_PROCESSOR_NAME]; //mpi use

  // smearing variables
  qcd_vector vec;                             // temporary vector using for smearing
  qcd_gaugeField u, uAPE;                     // store gauge field before and after smearing
  qcd_gaugeField *u_ptr, *uAPE_ptr, *utmp_ptr; // some pointer to gauge fields
  qcd_uint_4 nsmear, nsmearAPE;               // smearing parameters
  qcd_real_8 alpha, alphaAPE;                  // smearing parameters

  qcd_int_4 ctr, ctr2, ctr11[Noperators], ctr22, ctr33;           // special spin indices trick for gammas multiplications  
  qcd_int_2 cg5cg5g5_ind[16*16*16][6], gammag5_ind[Noperators][16][2], cg5g5cg5_ind[16*16][4]; // store the non zero indices and gamma values
  qcd_complex_16 cg5cg5g5_val[16*16*16], gammag5_val[Noperators][16], cg5g5cg5_val[16*16];  // %
  qcd_complex_16 gammag5[Noperators][4][4], g5cg5[4][4];
  qcd_complex_16 z1, z2;                       // temp complex variables                                                                                                                                     
  qcd_complex_16 C, C2;                        // temp complex variables                                                                           
  qcd_complex_16 *W1[Nmu][Nic], *W2[Noperators][Nmu][Nic], *W3[Nmu][Nmu][Nmu][Nic], *W5[Nmu][Nmu][Nmu][Nic], *W7[Nmu][Nic]; //W8=W6=W4=W2         W terms for contractions    
  qcd_complex_16 *W1p[Nmu][Nic], *W2p[Noperators][Nmu][Nic], *W3p[Nmu][Nmu][Nmu][Nic], *W5p[Nmu][Nmu][Nmu][Nic], *W7p[Nmu][Nic]; //W8=W6=W4=W2    W terms after fourier
  qcd_complex_16 *W1p_r[Nmu][Nic], *W2p_r[Noperators][Nmu][Nic], *W3p_r[Nmu][Nmu][Nmu][Nic], *W5p_r[Nmu][Nmu][Nmu][Nic], *W7p_r[Nmu][Nic];        // W after reduce to root
  qcd_complex_16 *W12[Noperators][Nmu][Nmu], *W34[Noperators][Nmu][Nmu], *W56[Noperators][Nmu][Nmu], *W78[Noperators][Nmu][Nmu];                                                    // after do all the color-spin sums
  qcd_complex_16 *threep[Noperators][Nmu][Nmu], *threep_proj[Noperators],*threep_test[Noperators][Nmu][Nmu];                                                                                                 // store three point function


 //set up MPI                                                                                 
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);           
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);            
  MPI_Get_processor_name(processor_name,&namelen);  

  if(argc != 2)
    {
      if(myid == 0) fprintf(stderr,"No input file specified \n"); // check for input file
      exit(EXIT_FAILURE);
    }
  
  strcpy(param_name, argv[1]); // copy the name of input file to a character variable

  if(myid == 0 )
    {
      i=0; // boolean check
      printf("Opening input file for reading parameters %s \n", param_name);
      params = qcd_getParams( param_name, &params_len); // retrieve params in string and params_len only for root
      if(params == NULL) i=1; // if fail to take params return 1
    }

  MPI_Bcast(&i, 1, MPI_INT, 0, MPI_COMM_WORLD); // broadcast success or fail to all
  if(i == 1) exit(EXIT_FAILURE); 

  MPI_Bcast(&params_len, 1, MPI_INT, 0, MPI_COMM_WORLD); // broadcast params_len to all

  if(myid != 0) params = (char*)malloc(params_len*sizeof(char)); // allocate memory for params for all exect root beacause a previous routine do that for root

  MPI_Bcast(params, params_len, MPI_CHAR, 0, MPI_COMM_WORLD); // Broadcast the params to all processors

  /////////// reading from strings for all processors //////////
  
  sscanf(qcd_getParam("<processors_txyz>",params,params_len),"%hd %hd %hd %hd",&P[0], &P[1], &P[2], &P[3]); // read number of processors in each direction
  sscanf(qcd_getParam("<lattice_txyz>",params,params_len),"%hd %hd %hd %hd",&L[0], &L[1], &L[2], &L[3]);    // read the dimensions of global lattice

  if(P[0] != 1)
    {
      if( myid == 0) fprintf(stderr, " Error! , Number of processors in t-direction must be always 1"); // only use one processor in time direction
      exit(EXIT_FAILURE);
    }

  if(qcd_initGeometry(&geo, L, P, theta, myid, numprocs)) exit(EXIT_FAILURE);                           // start the geometry each processor takes a segment of lattice
  if(myid == 0) printf( " Local Lattice : %i x %i x %i x %i \n", geo.lL[0], geo.lL[1], geo.lL[2], geo.lL[3]); // print the local lattice

  strcpy(output_name,qcd_getParam("<output_name>",params,params_len)); //get  output filename                                                                
  if(myid==0) printf("Got output file name: %s\n",output_name);

  strcpy(uprop_name,qcd_getParam("<propagator_u>",params,params_len)); //get  u propagator file list
  if(myid==0) printf("Got propagator file name: %s\n",uprop_name);

  strcpy(dprop_name,qcd_getParam("<propagator_d>",params,params_len)); //get d propagator file list
  if(myid==0) printf("Got propagator file name: %s\n",dprop_name);

  strcpy(momlist_name,qcd_getParam("<momenta_list>",params,params_len)); //get the filename for the momenta list
  if(myid==0) printf("Got momenta-list file name: %s\n",momlist_name);

  strcpy(phi_name,qcd_getParam("<phi_list>",params,params_len)); //get the filename of phi list                                                       
  if(myid==0) printf("Got phi-list file name: %s\n",phi_name);
  sscanf(qcd_getParam("<num_phi>",params,params_len),"%d",&num_phi); // read number of phi vectors 
  if(myid == 0) printf("Number of phi vectors is: %d\n",num_phi);

  strcpy(xi_1_name,qcd_getParam("<xi_1_list>",params,params_len)); //get the filename of xi_1 list                                   
  if(myid==0) printf("Got xi_1-list file name: %s\n",xi_1_name);
  sscanf(qcd_getParam("<num_xi_1>",params,params_len),"%d",&num_xi_1); // read number of xi_1 vectors                                                                  
  if(myid == 0)printf("Number of xi_1 vectors is: %d\n",num_xi_1);

  if( num_phi != num_xi_1 )
    {
      if(myid == 0)fprintf(stderr,"Error, num_phi must be the same with num_xi_1");
      exit(EXIT_FAILURE);
    }
  

  strcpy(gauge_name,qcd_getParam("<cfg_name>",params,params_len));                // read gauge filename
  if(myid==0) printf("Got conf name: %s\n",gauge_name);
  
  sscanf(qcd_getParam("<alpha_gauss>",params,params_len),"%lf",&alpha);           // read alpha parameter for wuppertal smearing
  if(myid==0) printf(" Got alpha_gauss: %lf\n",alpha);
  
  sscanf(qcd_getParam("<nsmear_gauss>",params,params_len),"%d",&nsmear);          // number of iterations for wuppertal smearing
  if(myid==0) printf(" Got nsmear_gauss: %d\n",nsmear);
  
  sscanf(qcd_getParam("<alpha_APE>",params,params_len),"%lf",&alphaAPE);         // read alpha parameter for APE smearing
  if(myid==0) printf(" Got alpha_APE: %lf\n",alphaAPE);
  
  sscanf(qcd_getParam("<nsmear_APE>",params,params_len),"%d",&nsmearAPE);         // number of iterations for APE smearing
  if(myid==0) printf(" Got nsmear_APE: %d\n",nsmearAPE);

  sscanf(qcd_getParam("<source_pos_txyz>",params,params_len),"%d %d %d %d",&x_src[0],&x_src[1],&x_src[2],&x_src[3]);  // read the position where we put the source
  if(myid==0) printf("Got source coords: %d %d %d %d\n",x_src[0],x_src[1],x_src[2],x_src[3]);

  sscanf(qcd_getParam("<current_pos_t>",params,params_len),"%d",&t_noise);                                              // read the time position of insertion
  if(myid==0) printf("Got noise vector time-slice dilution: %d \n",t_noise);

  sscanf(qcd_getParam("<proj_index>",params,params_len),"%d",&proj_index);                                         // index to choose the projector                                   
  if(myid==0) printf("Got index for projector matrix : %d \n",proj_index);


  free(params); 

  //////////////////////// finish with parameters //////////////////////////////////////////////////
  
  //allocate memory for gauge-fields////////////////////////////////////////// no comments i just copied it from twop.c
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
  if(myid == 0 ) printf("start APE smearing\n");
  for(i=0; i<nsmearAPE; i++)
    {
      qcd_apeSmear3d(uAPE_ptr, u_ptr, alphaAPE);
      utmp_ptr=u_ptr; u_ptr=uAPE_ptr; uAPE_ptr=utmp_ptr;
    }
  utmp_ptr=u_ptr; u_ptr=uAPE_ptr; uAPE_ptr=utmp_ptr; //reverse the last swap. Also needed when nsmearAPE=0                                                                                           
  qcd_destroyGaugeField(u_ptr);
  uAPE = *uAPE_ptr;
  if(myid == 0) printf("printf finish APE smearing\n");
  //////////////////////////////////////////////////////////////////////////////////

  ///////////////////// start initialize the resources ///////////////////// only for propagator and vectors

  j = 0; // flag for ok allocation
  j += qcd_initPropagator(&uprop, &geo); // initialize memory for u propagator
  j += qcd_initPropagator(&dprop, &geo); // initialize memory for d propagator
  j += qcd_initPropagator(&uprop_unsmear, &geo);
  j += qcd_initVector(&vec, &geo);

  j += qcd_initVector(&phi, &geo);       // initialize phi vector
  j += qcd_initVector(&xi_1, &geo);      // initialize xi_1 vector

  MPI_Allreduce(&j, &k, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); // readuce k for all to see if the allocation is ok!

  if(k>0) // check for memory allocation flag
    {
      if(myid==0) printf("not enough memory\n");
      exit(EXIT_FAILURE);
    }
  if(myid==0) printf("memory for propagators and tsm vectors allocated\n");

  //############################################################///////////////////

  // load propagators

  if(qcd_getPropagator(uprop_name,qcd_PROP_LIME, &uprop)) exit(EXIT_FAILURE); // read u propagators from binary in ram
  if(myid==0) printf("up propagator loaded\n");
  if(qcd_getPropagator(dprop_name,qcd_PROP_LIME, &dprop)) exit(EXIT_FAILURE); // read d propagators from binary in ram
  if(myid==0) printf("down propagator loaded\n");

  //################################################################################                                                                                                                 
  // transform propagators to basis with theta-periodic boundaries in the temporal direction  !!! remember to do the same fo solution vectors
  for(int lt=0; lt<geo.lL[0]; lt++)
    {
      t = lt + geo.Pos[0] * geo.lL[0];
      phase_factor   = (qcd_complex_16) {cos(theta[0]*t/geo.L[0]),sin(theta[0]*t/geo.L[0])};
      qcd_mulPropagatorC3d(&uprop, phase_factor, (t+x_src[0]) % geo.L[0]);                         // !!!!!!!!!!! ask giannis why he use  (t+x_src[0]) % geo.L[0])
      qcd_mulPropagatorC3d(&dprop, phase_factor, (t+x_src[0]) % geo.L[0]);
    }
  if(myid==0) printf("propagators transformed to basis with theta-periodic boundary conditions\n");

  //################################################################################                                                                                                                 
  qcd_copyPropagatorPropagator(&uprop_unsmear, &uprop); // copy to one no smeared propagator before smearing

  int Tseperation = t_noise > x_src[0] ? t_noise - x_src[0]+1 : t_noise - x_src[0]+1 + geo.L[0] ;

  // gaussian smearing of propagators (only the time-slices that will be used)                                                                                                                       
    
  for(mu=0;mu<4;mu++)
    for(c1=0;c1<3;c1++)
      {
	if(myid == 0)printf("smearing mu = %d , c1 = %d \n",mu,c1);
			qcd_copyVectorPropagator(&vec,&uprop,mu,c1);             // copy each of 12 vectors for propagator to do smearing
		
	for(i=0; i<nsmear; i++)
	  {
	    
	    if(qcd_gaussIteration3d(&vec,&uAPE,alpha,t_noise))   // smearing for up propagator
	      {
		fprintf(stderr,"process %i: Error while smearing!\n",geo.myid);
		exit(EXIT_FAILURE);
	      }
	    
	  }
	qcd_copyPropagatorVector(&uprop,&vec,mu,c1);
	qcd_copyVectorPropagator(&vec,&dprop,mu,c1);
	for(i=0; i<nsmear; i++)
	  {
	   
	    if(qcd_gaussIteration3d(&vec,&uAPE,alpha,t_noise))        // smearing for down propagator
	      {
		fprintf(stderr,"process %i: Error while smearing!\n",geo.myid);
		exit(EXIT_FAILURE);
	      }
	    
	  }
	qcd_copyPropagatorVector(&dprop,&vec,mu,c1);
		
      }
  
  qcd_destroyVector(&vec);
                                               // destroy only the temporary vector i dont nead any more
  if(myid==0) printf("propagators smeared\n");
  /////////////////////////////////////////////////////////////////////////////////////////
      
  //calculate the dirac non zero value and indices
  
 

  //
  ctr = 0;
  for(beta = 0; beta < 4 ; beta++)                // this are only for C * Gamma5 * \bar{C * Gamma5}
    for(nu = 0; nu < 4 ; nu++)
      for(lambda = 0 ; lambda < 4 ; lambda++)
	for(xita = 0 ; xita < 4 ; xita++)
	  for(delta = 0 ; delta < 4 ; delta++)
	    for(kappa = 0 ; kappa < 4 ; kappa++)
	      {
		C = qcd_CMUL(qcd_GAMMA[5][delta][kappa],qcd_CMUL(qcd_CGAMMA[5][beta][nu],qcd_BAR_CGAMMA[5][lambda][xita]));
		if(qcd_NORM(C) > 1e-3)
		  {
		    cg5cg5g5_val[ctr].re = C.re;                      // store values
		    cg5cg5g5_val[ctr].im = C.im;
		    cg5cg5g5_ind[ctr][0] = beta;                      // store indices
		    cg5cg5g5_ind[ctr][1] = nu;
		    cg5cg5g5_ind[ctr][2] = lambda;
		    cg5cg5g5_ind[ctr][3] = xita;
		    cg5cg5g5_ind[ctr][4] = delta;
		    cg5cg5g5_ind[ctr][5] = kappa;
		    ctr++;                                         // important dont forget to increase this index
		  }
	      }
  
 
  for(kappa = 0 ; kappa < 4 ; kappa++)
    for(nu = 0 ; nu < 4 ; nu++){
      sum = (qcd_complex_16) {0,0};
      for(beta = 0 ; beta < 4 ; beta++){
	sum = qcd_CADD(sum,qcd_CMUL(qcd_GAMMA[5][beta][kappa],qcd_CGAMMA[5][beta][nu]));
      }
      g5cg5[kappa][nu] = sum ;
    }

  ctr33 =0 ;
  for(lambda = 0 ; lambda < 4 ; lambda++)
    for(xita = 0 ; xita <4 ; xita++)
      for(kappa =0; kappa < 4 ; kappa++)                               // now for gamma inserion matrix only                                       
	for(nu=0; nu < 4 ; nu++){
	  C = qcd_CMUL(qcd_BAR_CGAMMA[5][lambda][xita],g5cg5[kappa][nu]) ;
	  if(qcd_NORM(C) > 1e-3)
	    {
	      cg5g5cg5_val[ctr33].re = C.re;
	      cg5g5cg5_val[ctr33].im = C.im;
	      cg5g5cg5_ind[ctr33][0] = lambda;
	      cg5g5cg5_ind[ctr33][1] = xita;
	      cg5g5cg5_ind[ctr33][2] = kappa;
	      cg5g5cg5_ind[ctr33][3] = nu;
	      ctr33++;
	    }
	}

  for(int gamma_index = 0 ; gamma_index < Noperators ; gamma_index++){
    for(zita = 0 ; zita < 4 ; zita++)
      for(phita = 0 ; phita < 4 ; phita++){
	sum = (qcd_complex_16) {0,0};
	for( rho = 0 ; rho < 4 ; rho++){
	  sum = qcd_CADD(sum,qcd_CMUL(qcd_GAMMA[5][zita][rho],qcd_GAMMA_OPERATORS[gamma_index][rho][phita]));
	}
	gammag5[gamma_index][zita][phita] = sum;
      }
  }

  for(int gamma_index = 0 ; gamma_index < Noperators ; gamma_index++){
    ctr11[gamma_index] = 0 ;
    for(zita =0; zita < 4 ; zita++)                               // now for gamma inserion matrix only
      for(phita=0; phita < 4 ; phita++){
	C = gammag5[gamma_index][zita][phita] ;
	if(qcd_NORM(C) > 1e-3)
	  {
	    gammag5_val[gamma_index][ctr11[gamma_index]].re = C.re;
	    gammag5_val[gamma_index][ctr11[gamma_index]].im = C.im;
	    gammag5_ind[gamma_index][ctr11[gamma_index]][0] = zita;
	    gammag5_ind[gamma_index][ctr11[gamma_index]][1] = phita;
	    ctr11[gamma_index]++;
	  }
      }
  }

  if(myid == 0)printf("finish all spin collections\n");
  //open output file to write in and file for momenta                                                                                                            
  if(myid==0)
    {
      fp_threep = fopen(output_name,"w");
      if(fp_threep==NULL)
	{
	  printf("failed to open %s for writing\n",output_name);
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
  if(k==1) exit(EXIT_FAILURE); // check fo success
  //////////////////////////////

  // load the momenta list///////////////////////////////
  if(myid==0) fscanf(fp_momlist,"%d\n",&num_momenta);
  
  MPI_Bcast(&num_momenta,1,MPI_INT, 0, MPI_COMM_WORLD);
  if(myid==0) printf("will read %d momenta combinations\n",num_momenta);
  mom = malloc(num_momenta*3*sizeof(qcd_int_4));
  if(myid==0)
    {
      for(j=0; j<num_momenta; j++)
	{
	  fscanf(fp_momlist,"%d %d %d\n",&(mom[j][0]),&(mom[j][1]),&(mom[j][2]));
	  //printf("got combination %i %i %i\n",mom[j][0],mom[j][1],mom[j][2]);                                                                                                                        
	}
      fclose(fp_momlist);
    }
  MPI_Bcast(&(mom[0][0]),num_momenta*3,MPI_INT,0, MPI_COMM_WORLD);
  if(myid==0) printf("momenta list read and broadcasted\n");
  ////////////////////// succesfuly pass the momenta to all cpu
  
  //allocate memory for all the W terms ///////////////////////////////////////////////
  for(i=0; i < Nmu ; i++)
    for(j=0; j < Nic ; j++){
      W1[i][j] = (qcd_complex_16*)malloc(geo.lV3 * sizeof(qcd_complex_16));
      W7[i][j] = (qcd_complex_16*)malloc(geo.lV3 * sizeof(qcd_complex_16));
      W1p[i][j] = (qcd_complex_16*)malloc(num_momenta * sizeof(qcd_complex_16));
      W7p[i][j] =(qcd_complex_16*)malloc(num_momenta *sizeof(qcd_complex_16));
      W1p_r[i][j] =(qcd_complex_16*)malloc(num_momenta *sizeof(qcd_complex_16));
      W7p_r[i][j] =(qcd_complex_16*)malloc(num_momenta *sizeof(qcd_complex_16));
      for(int gamma_index = 0; gamma_index < Noperators ; gamma_index++){
	W2[gamma_index][i][j] = (qcd_complex_16*)malloc(geo.lV3 * sizeof(qcd_complex_16));
	W2p[gamma_index][i][j] =(qcd_complex_16*)malloc(num_momenta *sizeof(qcd_complex_16));
	W2p_r[gamma_index][i][j] =(qcd_complex_16*)malloc(num_momenta *sizeof(qcd_complex_16));
      }
    }

  for( i=0; i < Nmu ; i++)
    for( j=0; j<Nmu ; j++)
      for(k=0; k<Nmu ; k++)
	for(l=0; l<Nic ; l++){
	  W3[i][j][k][l] = (qcd_complex_16*)malloc(geo.lV3 * sizeof(qcd_complex_16));
	  W5[i][j][k][l] = (qcd_complex_16*)malloc(geo.lV3 * sizeof(qcd_complex_16));
	  W3p[i][j][k][l] = (qcd_complex_16*)malloc(num_momenta * sizeof(qcd_complex_16));
	  W5p[i][j][k][l] =(qcd_complex_16*)malloc(num_momenta * sizeof(qcd_complex_16));
	  W3p_r[i][j][k][l] = (qcd_complex_16*)malloc(num_momenta * sizeof(qcd_complex_16));
	  W5p_r[i][j][k][l] =(qcd_complex_16*)malloc(num_momenta * sizeof(qcd_complex_16));
	}

  for(int gamma_index = 0 ;gamma_index < Noperators ; gamma_index++){
    threep_proj[gamma_index] = (qcd_complex_16*)malloc(num_momenta*geo.L[0]*sizeof(qcd_complex_16));
    for( i = 0 ; i < Nmu ; i++)
      for(j = 0 ; j <Nmu ; j++)
	{
	  W12[gamma_index][i][j] = (qcd_complex_16*)malloc(num_momenta*geo.L[0]*sizeof(qcd_complex_16));
	  W34[gamma_index][i][j] = (qcd_complex_16*)malloc(num_momenta*geo.L[0]*sizeof(qcd_complex_16));
	  W56[gamma_index][i][j] = (qcd_complex_16*)malloc(num_momenta*geo.L[0]*sizeof(qcd_complex_16));
	  W78[gamma_index][i][j] = (qcd_complex_16*)malloc(num_momenta*geo.L[0]*sizeof(qcd_complex_16));
	  threep[gamma_index][i][j] = (qcd_complex_16*)malloc(num_momenta*geo.L[0]*sizeof(qcd_complex_16));
          threep_test[gamma_index][i][j] = (qcd_complex_16*)malloc(num_momenta*geo.L[0]*sizeof(qcd_complex_16));
	  //	if(W56[i][j] == NULL)printf("malaka\n");
	}
  }

  for(int gamma_index =0 ;gamma_index < Noperators ; gamma_index++){
    for(kk=0; kk < num_momenta*geo.L[0] ; kk++){
      threep_proj[gamma_index][kk] = (qcd_complex_16) {0,0};
      for( i = 0 ; i < Nmu ; i++)
	for(j =0 ; j < Nmu ; j++)
	  {
	    threep[gamma_index][i][j][kk] = (qcd_complex_16) {0,0}; // initialize to zero the three point function
	    threep_test[gamma_index][i][j][kk] = (qcd_complex_16) {0,0}; // initialize to zero the three point function                                  
	  }
    }
  }
  
 

  ///////////////////////////////// finish with the allocation for W terms
  
  // read the names of files from a list////////////////////////////////////////////////////////////
  fp_phi_name = fopen(phi_name,"r");
  fp_xi_1_name = fopen(xi_1_name,"r");
  
  if(fp_phi_name == NULL){
    if(myid == 0)fprintf(stderr,"Error opening list with phi vectors\n");
    exit(EXIT_FAILURE);
  }
  if(fp_xi_1_name == NULL){
    if(myid == 0)fprintf(stderr,"Error opening list with xi_1 vectors\n");
    exit(EXIT_FAILURE);
  }
  
  for(int i=0;i<num_phi;i++)
    {
      fscanf(fp_phi_name,"%s",&(phi_files[i][0]));
    }  
    
fclose(fp_phi_name);


  for(int i=0;i<num_xi_1;i++)
    {
      imom++;
      fscanf(fp_xi_1_name,"%s",&(xi_1_files[i][0]));
    }
  
  fclose(fp_xi_1_name);
  
  //###################################### MAIN CALCULATIONS ####################################################
  for(int r = 0 ; r < num_phi ; r++){        // a summation over the ensemble of noise vectors

    for(int gamma_index = 0 ;gamma_index < Noperators ; gamma_index++){
      for( i = 0 ; i < Nmu ; i++)
	for(j = 0 ; j <Nmu ; j++)
	  for(kk = 0 ; kk < num_momenta*geo.L[0] ; kk++)
	    {
	      W12[gamma_index][i][j][kk] = (qcd_complex_16) {0,0};         // for each noise vector we have to set the values of them to zero
	      W34[gamma_index][i][j][kk] = (qcd_complex_16) {0,0};
	      W56[gamma_index][i][j][kk] = (qcd_complex_16) {0,0};
	      W78[gamma_index][i][j][kk] = (qcd_complex_16) {0,0};
	    }
    }

    //load the current phi and xi_1
    qcd_getVector(&(phi_files[r][0]),qcd_PROP_LIME,0,0,&phi);
    qcd_getVector(&(xi_1_files[r][0]),qcd_PROP_LIME,0,0,&xi_1);

    //we have to multiply phi by the exponential with theta
    for(int llt=0; llt<geo.lL[0]; llt++)
      {
	t = llt + geo.Pos[0] * geo.lL[0];
	phase_factor   = (qcd_complex_16) {cos(theta[0]*t/geo.L[0]),sin(theta[0]*t/geo.L[0])};
	qcd_mulVectorC3d(&phi, phase_factor, (t+t_noise) % geo.L[0]);                               
      }
    
    //we have to do smearing on the xi
    
    for(i=0; i<nsmear; i++)
      {
    	if(qcd_gaussIteration3d(&xi_1,&uAPE,alpha,t_noise))
    	  {
    	    fprintf(stderr,"process %i: Error while smearing!\n",geo.myid);
    	    exit(EXIT_FAILURE);
    	  }
      }
    	
    // set to zero the W // this variables must set to zero for each time slice
    for(i=0; i < Nmu ; i++)
      for(j=0; j < Nic ; j++)
	for(v=0; v<geo.lV3; v++){
	  W1[i][j][v] = (qcd_complex_16) {0,0};
	  W7[i][j][v] = (qcd_complex_16) {0,0};
	}
    
    for( i=0; i < Nmu ; i++)
      for( j=0; j<Nmu ; j++)
	for(k=0; k<Nmu ; k++)
	  for(l=0; l<Nic ; l++)
	    for(v=0; v<geo.lV3 ; v++){
	      W3[i][j][k][l][v] = (qcd_complex_16) {0,0};
	      W5[i][j][k][l][v] = (qcd_complex_16) {0,0};
	    }

	    ///////////////////////////////////////////////////////////////////////
    if(myid == 0)printf("start W1 W3 W5 W7 for r = %d \n",r);
   
    t = t_noise;                                         // very important this happens only for the sink
    for(lx = 0 ; lx < geo.lL[1] ; lx++)                  // each processors do his work for its local lattice
      for(ly = 0 ; ly< geo.lL[2] ; ly++)
	for(lz = 0 ; lz < geo.lL[3] ; lz++){
	  v3 = qcd_LEXIC0(lx,ly,lz,geo.lL);
	  v =  qcd_LEXIC(t,lx,ly,lz,geo.lL);
	  	      
	  for(cc1 = 0 ; cc1 < 6 ; cc1++){
	    aa = qcd_EPS[cc1][0];
	    bb = qcd_EPS[cc1][1];
	    cc = qcd_EPS[cc1][2];
	    
	    for(cc2 = 0 ; cc2 < 6 ; cc2++){
	      ff = qcd_EPS[cc2][0];
	      gg = qcd_EPS[cc2][1];
	      hh = qcd_EPS[cc2][2];
	      
	      for(ctr2 = 0; ctr2 < ctr ; ctr2++){
		beta = cg5cg5g5_ind[ctr2][0];
		nu = cg5cg5g5_ind[ctr2][1];
		lambda = cg5cg5g5_ind[ctr2][2];
		xita = cg5cg5g5_ind[ctr2][3];
		delta = cg5cg5g5_ind[ctr2][4];
		kappa = cg5cg5g5_ind[ctr2][5];       
		
		W1[delta][hh][v3] = qcd_CADD(W1[delta][hh][v3],qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( cg5cg5g5_val[ctr2] , 
				    qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),xi_1.D[v][kappa][cc]),uprop.D[v][beta][lambda][aa][ff]),dprop.D[v][nu][xita][bb][gg]));
		
		for(sigma =0 ; sigma < 4 ; sigma++)
		  W3[delta][sigma][lambda][ff][v3] = qcd_CADD(W3[delta][sigma][lambda][ff][v3],
							      qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( cg5cg5g5_val[ctr2], qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]), xi_1.D[v][kappa][cc]),uprop.D[v][beta][sigma][aa][hh]),
								       dprop.D[v][nu][xita][bb][gg]));
	
	      } // close dirac loop for terms W1 & W3

	      for(delta = 0 ; delta < 4 ; delta++)
		for(ctr2 = 0 ; ctr2 < ctr33 ; ctr2++){
		  lambda = cg5g5cg5_ind[ctr2][0];
		  xita = cg5g5cg5_ind[ctr2][1];
		  kappa = cg5g5cg5_ind[ctr2][2];
		  nu = cg5g5cg5_ind[ctr2][3];
		  W7[delta][hh][v3] = qcd_CADD(W7[delta][hh][v3],qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( cg5g5cg5_val[ctr2] ,
				      qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),xi_1.D[v][kappa][aa]),uprop.D[v][delta][lambda][cc][ff]),dprop.D[v][nu][xita][bb][gg]));
		  
		  for(sigma =0 ; sigma < 4 ; sigma++)
		    W5[delta][sigma][lambda][ff][v3] = qcd_CADD(W5[delta][sigma][lambda][ff][v3],qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_CSCALE( cg5g5cg5_val[ctr2] ,
						       qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]),xi_1.D[v][kappa][aa]),uprop.D[v][delta][sigma][cc][hh]), 
						       dprop.D[v][nu][xita][bb][gg]));
		} // close dirac loop for terms W5 & W7
	      
	    } // close color
	  } //close color
	} //close spatials

    if(myid ==0) printf("finish W1 W3 W5 W7 for r = %d \n",r);
    
    for( it = 0 ; it < Tseperation ; it++)                                         
      {
	if(myid == 0)printf("start W2 r = %d , t = %d\n",r,it);
	t = ((it+x_src[0])%geo.L[0]);
	for(int gamma_index = 0 ;gamma_index < Noperators ; gamma_index++){
	  for(i=0; i < Nmu ; i++)
	    for(j=0; j < Nic ; j++)
	      for(v=0; v<geo.lV3; v++)
		W2[gamma_index][i][j][v] = (qcd_complex_16) {0,0};          // this term must be zero and for noise loop & for time loop
	}
	
	for(int gamma_index = 0 ;gamma_index < Noperators ; gamma_index++){
	  for(lx = 0 ; lx < geo.lL[1] ; lx++)                  // each processors do his work for its local lattice
	    for(ly = 0 ; ly< geo.lL[2] ; ly++)
	      for(lz = 0 ; lz < geo.lL[3] ; lz++){
		v3 = qcd_LEXIC0(lx,ly,lz,geo.lL);
		v =  qcd_LEXIC(t,lx,ly,lz,geo.lL);

		// for W2 term
		for(tt=0 ; tt< Nic ;tt++)
		  for( ctr2 =0 ; ctr2 < ctr11[gamma_index] ; ctr2++){
		    zita = gammag5_ind[gamma_index][ctr2][0];
		    phita = gammag5_ind[gamma_index][ctr2][1];
		    
		    for(hh=0; hh<Nic ;hh++)
		      for(sigma =0 ;sigma<Nmu ; sigma++){
			W2[gamma_index][sigma][hh][v3] = qcd_CADD(W2[gamma_index][sigma][hh][v3], 
								  qcd_CMUL(gammag5_val[gamma_index][ctr2],qcd_CMUL(qcd_CONJ(phi.D[v][zita][tt]),uprop_unsmear.D[v][phita][sigma][tt][hh])));
		    
		      }
		  } //close color and dirac
	      } // close spatial for W2
	} // close for all operators
	    // Fourier Transform //////////////////////////////////

	for(i=0; i < Nmu ; i++) // set to zero
	  for(j=0; j < Nic ; j++)
	    for(v=0; v< num_momenta; v++){
	      W1p[i][j][v] = (qcd_complex_16) {0,0};
	      W7p[i][j][v] = (qcd_complex_16) {0,0};
	      W1p_r[i][j][v] = (qcd_complex_16) {0,0};
	      W7p_r[i][j][v] = (qcd_complex_16) {0,0};
	      for(int gamma_index = 0 ;gamma_index < Noperators ; gamma_index++){
		W2p[gamma_index][i][j][v] = (qcd_complex_16) {0,0};
		W2p_r[gamma_index][i][j][v] = (qcd_complex_16) {0,0};
	      }
	    }

	for( i=0; i < Nmu ; i++) // set to zero
	  for( j=0; j<Nmu ; j++)
	    for(k=0; k<Nmu ; k++)
	      for(l=0; l<Nic ; l++)
		for(v=0; v<num_momenta ; v++){
		  W3p[i][j][k][l][v] = (qcd_complex_16) {0,0};
		  W5p[i][j][k][l][v] = (qcd_complex_16) {0,0};
		  W3p_r[i][j][k][l][v] = (qcd_complex_16) {0,0};
		  W5p_r[i][j][k][l][v] = (qcd_complex_16) {0,0};
		}
	if(myid == 0)printf("start fourier for r = %d, t = %d\n",r,t);

	
	for( j=0 ; j< num_momenta ; j++){ // begin fourier
	    
	  for(lx=0; lx<geo.lL[1]; lx++)
	    for(ly=0; ly<geo.lL[2]; ly++)
	      for(lz=0; lz<geo.lL[3]; lz++)
		{
		  
		  v3 = qcd_LEXIC0(lx,ly,lz,geo.lL);
		  x=lx+geo.Pos[1]*geo.lL[1] - x_src[1];                   /// !!!!!!!!!! ask giannis why substract x_src
		  y=ly+geo.Pos[2]*geo.lL[2] - x_src[2];
		  z=lz+geo.Pos[3]*geo.lL[3] - x_src[3];
		  tmp = (((double) mom[j][0]*x)/geo.L[1] + ((double) mom[j][1]*y)/geo.L[2] + ((double) mom[j][2]*z)/geo.L[3])*2*M_PI;
		  C2=(qcd_complex_16) {cos(tmp),sin(tmp)}; //TABULATE FOR LARGE SPEEDUP!!!                                                                                                   
		  
		  for( hh = 0 ; hh < Nic ; hh++)
		    for( delta = 0 ; delta < Nmu ; delta++){
		      W1p[delta][hh][j]=qcd_CADD(W1p[delta][hh][j], qcd_CMUL(W1[delta][hh][v3],C2));
		      for(int gamma_index = 0 ;gamma_index < Noperators ; gamma_index++)
			W2p[gamma_index][delta][hh][j]=qcd_CADD(W2p[gamma_index][delta][hh][j], qcd_CMUL(W2[gamma_index][delta][hh][v3],C2));
		      W7p[delta][hh][j]=qcd_CADD(W7p[delta][hh][j], qcd_CMUL(W7[delta][hh][v3],C2));

		      for( mu = 0 ; mu< Nmu ; mu++)
			for( nu = 0 ; nu <Nmu ; nu++){
			  W3p[delta][mu][nu][hh][j] = qcd_CADD(W3p[delta][mu][nu][hh][j], qcd_CMUL(W3[delta][mu][nu][hh][v3],C2));
			  W5p[delta][mu][nu][hh][j] = qcd_CADD(W5p[delta][mu][nu][hh][j], qcd_CMUL(W5[delta][mu][nu][hh][v3],C2));
			}
		    }
		}
	    
	  for( hh = 0 ; hh < Nic ; hh++) // reduce the values to root
	    for(delta = 0 ;delta <Nmu ; delta++){
	      MPI_Reduce(&(W1p[delta][hh][j].re), &(W1p_r[delta][hh][j].re), 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	      for(int gamma_index = 0 ;gamma_index < Noperators ; gamma_index++)
		MPI_Reduce(&(W2p[gamma_index][delta][hh][j].re), &(W2p_r[gamma_index][delta][hh][j].re), 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	      MPI_Reduce(&(W7p[delta][hh][j].re), &(W7p_r[delta][hh][j].re), 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	      for( mu = 0 ; mu < Nmu ; mu++)
		for(nu = 0 ; nu < Nmu ; nu++){
		  MPI_Reduce(&(W3p[delta][mu][nu][hh][j].re), &(W3p_r[delta][mu][nu][hh][j].re), 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		  MPI_Reduce(&(W5p[delta][mu][nu][hh][j].re), &(W5p_r[delta][mu][nu][hh][j].re), 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		}
	    }
	  

	  
	} // close loop for momentum for fourier

	


	    ////////////////////////////////////////////////////////////////////////
	  
	for(int gamma_index = 0 ;gamma_index < Noperators ; gamma_index++){
	  for( int q = 0 ; q < num_momenta; q++)
	    for( delta = 0 ; delta < Nmu ; delta++)
	      for( sigma = 0 ; sigma < Nmu ; sigma++){
		
		for( hh = 0 ; hh < Nic ; hh++){
		  W12[gamma_index][delta][sigma][q*geo.L[0]+t] = qcd_CADD(W12[gamma_index][delta][sigma][q*geo.L[0]+t],qcd_CMUL(W1p_r[delta][hh][0] , W2p_r[gamma_index][sigma][hh][q]));
		  W78[gamma_index][delta][sigma][q*geo.L[0]+t] = qcd_CADD(W78[gamma_index][delta][sigma][q*geo.L[0]+t],qcd_CMUL(W7p_r[delta][hh][0] , W2p_r[gamma_index][sigma][hh][q]));
		  
		  for(int lambda = 0 ; lambda < Nmu ; lambda++){
		     W34[gamma_index][delta][sigma][q*geo.L[0]+t] = qcd_CADD(W34[gamma_index][delta][sigma][q*geo.L[0]+t],qcd_CMUL(W3p_r[delta][sigma][lambda][hh][0],W2p_r[gamma_index][lambda][hh][q]));
		     W56[gamma_index][delta][sigma][q*geo.L[0]+t] = qcd_CADD(W56[gamma_index][delta][sigma][q*geo.L[0]+t],qcd_CMUL(W5p_r[delta][sigma][lambda][hh][0],W2p_r[gamma_index][lambda][hh][q]));
		  }
		  
		}
	      }

	  for( int q = 0 ; q < num_momenta; q++)
	    for( delta = 0 ; delta < Nmu ; delta++)
	      for( sigma = 0 ; sigma < Nmu ; sigma++){
		threep[gamma_index][delta][sigma][q*geo.L[0]+t] = qcd_CADD(threep[gamma_index][delta][sigma][q*geo.L[0]+t],qcd_CADD(qcd_CSUB(W12[gamma_index][delta][sigma][q*geo.L[0]+t],
																	     W34[gamma_index][delta][sigma][q*geo.L[0]+t]),qcd_CSUB(W56[gamma_index][delta][sigma][q*geo.L[0]+t], W78[gamma_index][delta][sigma][q*geo.L[0]+t])));
		threep_test[gamma_index][delta][sigma][q*geo.L[0]+t] = threep[gamma_index][delta][sigma][q*geo.L[0]+t];
	      }
	}

	
      } //close time

    
    if(myid == 0){
      for(int gamma_index = 0 ;gamma_index < Noperators ; gamma_index++){
	for( int q = 0 ; q < num_momenta; q++)
	  for(int it = 0 ; it < Tseperation ; it++)
	    {

	      t = ((it+x_src[0])%geo.L[0]);
	      for(int delta = 0 ; delta < Nmu ; delta++)
		for(int sigma = 0 ; sigma < Nmu ; sigma++){
		  threep_test[gamma_index][delta][sigma][q*geo.L[0]+t].re = -(threep_test[gamma_index][delta][sigma][q*geo.L[0]+t].re/(r+1));
		  threep_test[gamma_index][delta][sigma][q*geo.L[0]+t].im = -(threep_test[gamma_index][delta][sigma][q*geo.L[0]+t].im/(r+1));
		}

	      for(int delta = 0 ; delta < Nmu ; delta++)
		for(int sigma = 0 ; sigma < Nmu ; sigma++){
		  threep_proj[gamma_index][q*geo.L[0]+t] = qcd_CADD(threep_proj[gamma_index][q*geo.L[0]+t],qcd_CMUL(PROJECTOR[proj_index][delta][sigma],threep_test[gamma_index][sigma][delta][q*geo.L[0]+t]));
		}

	      fprintf(fp_threep,"%6d %2d %3d %4d %+e %+e\n",r,gamma_index,it,q,threep_proj[gamma_index][q*geo.L[0]+t].re,threep_proj[gamma_index][q*geo.L[0]+t].im);

	    }
      }
    }
    
  } // close r
 



  if(myid == 0)printf("finish program \n");
  //finish MPI library
  qcd_waitall(&geo);
  MPI_Finalize();


  return 0;
}
