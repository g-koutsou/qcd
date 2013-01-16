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


#define MAX_STRING 1024

int 
main(int argc,char* argv[])
{
  int numprocs, myid, namelen;
  char processor_name[MAX_STRING];
  //set up MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);         // num. of processes taking part in the calculation
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);             // each process gets its ID
  MPI_Get_processor_name(processor_name,&namelen); // 
  
  
  //////////////////// READ INPUT FILE /////////////////////////////////////////////      
  if(argc!=2) {
    if(myid==0) 
      fprintf(stderr,"No input file specified\n");
     
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
   
  char param_name[MAX_STRING];
  char *params = NULL;
  int params_len;
  strcpy(param_name,argv[1]);
  if(myid==0) {
    printf("opening input file %s\n",param_name);
    params = qcd_getParams(param_name, &params_len);
    if(params==NULL) {
      if(myid == 0) 
	fprintf(stderr, " error loading parameter file\n");

      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
  }

  MPI_Bcast(&params_len, 1, MPI_INT, 0, MPI_COMM_WORLD);   
  if(myid!=0) 
    params = (char*) malloc(params_len*sizeof(char));
   
  MPI_Bcast(params, params_len, MPI_CHAR, 0, MPI_COMM_WORLD);
   
  unsigned short int P[4], L[4];
  sscanf(qcd_getParam("<processors_txyz>",params,params_len),"%hu %hu %hu %hu",&P[0], &P[1], &P[2], &P[3]);
  sscanf(qcd_getParam("<lattice_txyz>",params,params_len),"%hu %hu %hu %hu",&L[0], &L[1], &L[2], &L[3]);

  if(P[0] != 1) {
    if(myid==0) 
      fprintf(stderr,"Error! Number of processors in t-direction must be one.\n");
     
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
   
  qcd_real_8 theta[4] = {M_PI,0.0,0.0,0.0};   
  qcd_geometry geo;
  if(qcd_initGeometry(&geo,L,P, theta, myid, numprocs))
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
   
  if(myid==0) 
    printf(" Local lattice: %i x %i x %i x %i\n",geo.lL[0],geo.lL[1],geo.lL[2],geo.lL[3]);
      
  int t_start, t_stop;
  sscanf(qcd_getParam("<t>", params, params_len),"%d %d", &t_start, &t_stop);
  if(myid==0) 
    printf(" Got sink time slices: %d ... %d\n", t_start, t_stop);
   
  int nt_sink = t_stop - t_start + 1;
   
  int x_src[4];
  sscanf(qcd_getParam("<source_pos_txyz>", params, params_len),"%d %d %d %d",&x_src[0],&x_src[1],&x_src[2],&x_src[3]);
  if(myid==0) 
    printf("Got source coords: %d %d %d %d\n", x_src[0], x_src[1], x_src[2], x_src[3]);
      
  char uprop_name[MAX_STRING];
  strcpy(uprop_name,qcd_getParam("<propagator_u>", params, params_len));
  if(myid==0) 
    printf("Got propagator file name: %s\n", uprop_name);
   
  char dprop_name[MAX_STRING];
  strcpy(dprop_name,qcd_getParam("<propagator_d>",params,params_len));
  if(myid==0) 
    printf("Got propagator file name: %s\n", dprop_name);
   
  char gauge_name[MAX_STRING];
  strcpy(gauge_name,qcd_getParam("<cfg_name>",params,params_len));
  if(myid==0) 
    printf("Got conf name: %s\n",gauge_name);
          
  char corr_p_name[MAX_STRING];
  strcpy(corr_p_name,qcd_getParam("<corr_name>",params,params_len));
  if(myid==0) 
    printf("Got output file name: %s\n",corr_p_name);
           
  char momlist_name[MAX_STRING];
  strcpy(momlist_name,qcd_getParam("<momenta_list>",params,params_len));
  if(myid==0) 
    printf("Got momenta-list file name: %s\n",momlist_name);

  double alpha;
  sscanf(qcd_getParam("<alpha_gauss>",params,params_len),"%lf",&alpha);
  if(myid==0) 
    printf(" Got alpha_gauss: %lf\n",alpha);

  int nsmear;
  sscanf(qcd_getParam("<nsmear_gauss>",params,params_len),"%d",&nsmear);
  if(myid==0) 
    printf(" Got nsmear_gauss: %d\n",nsmear);

  double alphaAPE;
  sscanf(qcd_getParam("<alpha_APE>",params,params_len),"%lf",&alphaAPE);
  if(myid==0) 
    printf(" Got alpha_APE: %lf\n",alphaAPE);

  int nsmearAPE;
  sscanf(qcd_getParam("<nsmear_APE>",params,params_len),"%d",&nsmearAPE);
  if(myid==0) 
    printf(" Got nsmear_APE: %d\n",nsmearAPE); 

  free(params);
         
  //#####################################################################   
  // allocate memory
  // load gauge-field and APE-smear it
  int ierr;
  qcd_gaugeField u, uAPE;
  ierr = qcd_initGaugeField(&u, &geo);
  if(ierr) {
    if(myid == 0)
      fprintf(stderr, " error allocating gauge field\n");
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  ierr = qcd_initGaugeField(&uAPE,&geo);
  if(ierr) {
    if(myid == 0)
      fprintf(stderr, " error allocating APE gauge field\n");
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
      
  ierr = qcd_getGaugeField(gauge_name,qcd_GF_LIME,&u);
  if(ierr) {
    if(myid==0) 
      fprintf(stderr,"Error reading gauge field\n");

    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  if(myid==0) 
    printf("gauge-field loaded\n");   

  double plaq = qcd_calculatePlaquette(&u);
  if(myid==0) 
    printf("plaquette = %e\n",plaq);

  qcd_communicateGaugePM(&u);
  qcd_waitall(&geo);
  qcd_gaugeField *u_ptr[] = {&u, &uAPE};
  for(int i=0; i<nsmearAPE; i++) {
    qcd_apeSmear3d(u_ptr[(i+1)%2], u_ptr[i%2], alphaAPE);
  }
  if(nsmearAPE % 2 == 0)
    qcd_copyGaugeField(&uAPE, &u);
   
  qcd_communicateGaugePM(&uAPE);
  qcd_waitall(&geo);

  plaq = qcd_calculatePlaquette(&uAPE);
  if(myid==0) 
    printf("plaquette = %e\n",plaq);

  qcd_destroyGaugeField(&u);
   
  qcd_propagator uprop, dprop;
  qcd_vector vec;
  ierr = qcd_initPropagator(&uprop, &geo);
  if(ierr) {
    if(myid == 0)
      fprintf(stderr, " error allocating up-prop memory\n");

    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  ierr = qcd_initPropagator(&dprop, &geo);
  if(ierr) {
    if(myid == 0)
      fprintf(stderr, " error allocating down-prop memory\n");

    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  ierr = qcd_initVector(&vec, &geo);
  if(ierr) {
    if(myid == 0)
      fprintf(stderr, " error allocating temporary vector memory\n");

    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  //##############################################################################
  // load propagators
  ierr = qcd_getPropagator(uprop_name,qcd_PROP_LIME, &uprop);
  if(ierr) {
    if(myid == 0)
      fprintf(stderr, " %s: error reading file\n", uprop_name);

    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
  if(myid==0) 
    printf("up propagator loaded\n");

  ierr = qcd_getPropagator(dprop_name,qcd_PROP_LIME, &dprop);
  if(ierr) {
    if(myid == 0)
      fprintf(stderr, " %s: error reading file\n", dprop_name);

    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
  if(myid==0) 
    printf("down propagator loaded\n");   

  //################################################################################
  // transform propagators to basis with theta-periodic boundaries in the temporal direction
  for(int lt=0; lt<geo.lL[0]; lt++) {
    int t = lt + geo.Pos[0] * geo.lL[0];     
    qcd_complex_16 phase_factor = (qcd_complex_16) {
      cos(theta[0]*t/geo.L[0]),
      sin(theta[0]*t/geo.L[0])
    };     
    qcd_mulPropagatorC3d(&uprop, phase_factor, (t+x_src[0]) % geo.L[0]);
    qcd_mulPropagatorC3d(&dprop, phase_factor, (t+x_src[0]) % geo.L[0]);
  }
  if(myid==0) 
    printf("propagators transformed to basis with theta-periodic boundary conditions\n");
      
  //################################################################################
  // gaussian smearing of propagators (only the time-slices that will be used)
  for(int mu=0;mu<4;mu++)
    for(int c1=0;c1<3;c1++) {
      qcd_copyVectorPropagator(&vec,&uprop,mu,c1);
      for(int i=0; i<nsmear; i++) {
	qcd_gaussIteration3dAll(&vec,&uAPE,alpha,i==0);
      }
      qcd_copyPropagatorVector(&uprop,&vec,mu,c1);

      qcd_copyVectorPropagator(&vec,&dprop,mu,c1);
      for(int i=0; i<nsmear; i++) {
	qcd_gaussIteration3dAll(&vec,&uAPE,alpha,i==0);
      }
      qcd_copyPropagatorVector(&dprop,&vec,mu,c1);
    }
  qcd_destroyVector(&vec);
  qcd_destroyGaugeField(&uAPE);
   
  if(myid==0) 
    printf("propagators smeared\n");         
   
  //################################################################################   
  // calculate the 'blocks'         
  int ctr = 0;
  unsigned int cg5cg5_ind[16*16][4];
  qcd_complex_16 cg5cg5_val[16*16];
  for(int mu=0;mu<4;mu++)
    for(int nu=0;nu<4;nu++)
      for(int ku=0;ku<4;ku++)
	for(int lu=0;lu<4;lu++) {
	  qcd_complex_16 C = qcd_CMUL(qcd_CGAMMA[5][mu][nu],
				      qcd_BAR_CGAMMA[5][ku][lu]);
	  if(qcd_NORM(C)>1e-3) {
	    cg5cg5_val[ctr].re = C.re;
	    cg5cg5_val[ctr].im = C.im;
	    cg5cg5_ind[ctr][0] = mu;
	    cg5cg5_ind[ctr][1] = nu;
	    cg5cg5_ind[ctr][2] = ku;
	    cg5cg5_ind[ctr][3] = lu;                                                            
	    ctr++;
	  }
	}   
   
  //load momenta-list   
  int n_mom;
  int (*mom)[3] = NULL;
  if(myid == 0) {
    FILE *fp_momlist = fopen(momlist_name,"r");
    if(fp_momlist==NULL) {
      fprintf(stderr, "%s: failed to open for reading\n", momlist_name);
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    fscanf(fp_momlist,"%i\n",&n_mom);
    mom = malloc(n_mom*3*sizeof(int));
    for(int m=0; m<n_mom; m++)
      fscanf(fp_momlist, "%d %d %d\n", &(mom[m][0]), &(mom[m][1]), &(mom[m][2]));
     
    fclose(fp_momlist);
  }
   
  MPI_Bcast(&n_mom, 1 ,MPI_INT, 0, MPI_COMM_WORLD);
  if(myid != 0)
    mom = malloc(n_mom*3*sizeof(int));
   
  MPI_Bcast(&(mom[0][0]), 3*n_mom, MPI_INT, 0, MPI_COMM_WORLD);

  qcd_complex_16 (*block)[2][16] = malloc(geo.lV3*nt_sink*2*16*sizeof(qcd_complex_16));
  if(block == NULL) {
    if(myid == 0)                                                    
      fprintf(stderr, " error allocating block memory\n");
     
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
  for(int v=0; v<geo.lV3*nt_sink; v++)
    for(int ch=0; ch<2; ch++)
      for(int mu=0; mu<16; mu++)
	block[v][ch][mu] = (qcd_complex_16) {0,0};

  for(int t=t_start; t<=t_stop; t++) {
    int lt = ((t+x_src[0])%geo.L[0]) - geo.Pos[0]*geo.lL[0];
    int it = t-t_start;

    qcd_propagator *q1[2] = {&uprop, &dprop};
    qcd_propagator *q2[2] = {&dprop, &uprop};

    if(myid==0) 
      printf("t = %3d\n",t);
     
    for(int al=0; al<4; al++)
      for(int be=0; be<4; be++)
	for(int ctr2=0; ctr2<ctr; ctr2++) {          
	  int mu = cg5cg5_ind[ctr2][0];
	  int nu = cg5cg5_ind[ctr2][1];
	  int ku = cg5cg5_ind[ctr2][2];
	  int lu = cg5cg5_ind[ctr2][3];                              
	  for(int cc1=0;cc1<6;cc1++) {
	    int c1 = qcd_EPS[cc1][0];
	    int c2 = qcd_EPS[cc1][1];
	    int c3 = qcd_EPS[cc1][2];
	    for(int cc2=0;cc2<6;cc2++) {          
	      int c1p = qcd_EPS[cc2][0];
	      int c2p = qcd_EPS[cc2][1];
	      int c3p = qcd_EPS[cc2][2];
	      for(int ch=0; ch<2; ch++) {

#pragma omp parallel for
		for(int v3=0; v3<geo.lV3; v3++) {
		  int v = lt + v3*geo.lL[0];
		  int iv = v3 + it*geo.lV3;
		  block[iv][ch][al + be*4] = qcd_CSUB(block[iv][ch][al + be*4]
						      ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(cg5cg5_val[ctr2]
											     ,q1[ch]->D[v][mu][ku][c1][c1p])
										    ,q2[ch]->D[v][nu][lu][c2][c2p])
									   ,q1[ch]->D[v][be][al][c3][c3p])
								  ,qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]));

		  block[iv][ch][al + be*4] = qcd_CSUB(block[iv][ch][al + be*4]
						      ,qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(cg5cg5_val[ctr2]
											     ,q1[ch]->D[v][mu][al][c1][c1p])
										    ,q2[ch]->D[v][nu][lu][c2][c2p])
									   ,q1[ch]->D[v][be][ku][c3][c3p])
								  ,qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]));
		}//space loop
	      }
	    }//color2 loop    
	  }//color1 loop
	}//nonvanishing cg5cg5 loop
  }
   
  qcd_complex_16 (*corr_p)[2][16] = NULL, (*corr)[2][16] = NULL;
  corr_p = malloc(sizeof(qcd_complex_16)*2*n_mom*16);
  if(myid == 0)
    corr = malloc(sizeof(qcd_complex_16)*2*n_mom*nt_sink*16);
   
  for(int it=0; it<nt_sink; it++) {
#pragma omp parallel for
    for(int m=0; m<n_mom; m++) {
      for(int mu=0; mu<16; mu++)
	for(int ch=0; ch<2; ch++)
	  corr_p[m][ch][mu] = (qcd_complex_16) {0,0};

      for(int lx=0; lx<geo.lL[1]; lx++)
	for(int ly=0; ly<geo.lL[2]; ly++)
	  for(int lz=0; lz<geo.lL[3]; lz++) {
	    int v3 = qcd_LEXIC0(lx,ly,lz,geo.lL);
	    int x = lx + geo.Pos[1]*geo.lL[1] - x_src[1];
	    int y = ly + geo.Pos[2]*geo.lL[2] - x_src[2];
	    int z = lz + geo.Pos[3]*geo.lL[3] - x_src[3];
	    double angle = 
	      ((double) mom[m][0]*x/(double)geo.L[1] + 
	       (double) mom[m][1]*y/(double)geo.L[2] + 
	       (double) mom[m][2]*z/(double)geo.L[3])*2*M_PI;
       
	    qcd_complex_16 phase = (qcd_complex_16) {cos(angle), -sin(angle)}; //TABULATE FOR LARGE SPEEDUP!!!
	    for(int ch=0; ch<2; ch++)
	      for(int mu=0; mu<16; mu++)
		corr_p[m][ch][mu] = qcd_CADD(corr_p[m][ch][mu], qcd_CMUL(block[v3+it*geo.lV3][ch][mu], phase));
	  }
    }
    MPI_Reduce(corr_p, ((double *)corr+2*it*n_mom*2*16), 2*n_mom*2*16, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);     
  }//end t-loop   
            
  if(myid==0) {
    FILE *fp = fopen(corr_p_name, "w");
    if(fp == NULL) {
      fprintf(stderr, " %s: error opening file for writing\n", corr_p_name);
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    for(int it=0; it<nt_sink; it++)
      for(int m=0; m<n_mom; m++) 
	for(int ch=0; ch<2; ch++) 
	  for(int mu=0; mu<4; mu++) {
	    fprintf(fp, " %3d %+d %+d %+d %+e %+e  %+e %+e  %+e %+e  %+e %+e %3s\n", 
		    it % geo.L[0], 
		    mom[m][0], mom[m][1], mom[m][2],
		    corr[m + it*n_mom][ch][mu*4 + 0].re*0.5, corr[m + it*n_mom][ch][mu*4 + 0].im*0.5,
		    corr[m + it*n_mom][ch][mu*4 + 1].re*0.5, corr[m + it*n_mom][ch][mu*4 + 1].im*0.5,
		    corr[m + it*n_mom][ch][mu*4 + 2].re*0.5, corr[m + it*n_mom][ch][mu*4 + 2].im*0.5,
		    corr[m + it*n_mom][ch][mu*4 + 3].re*0.5, corr[m + it*n_mom][ch][mu*4 + 3].im*0.5,
		    ch == 0 ? "ppn" : "pnn");
	  }

    fclose(fp);
  }
  //#####################################################################   
  // clean up
  if(myid == 0)
    free(corr);

  free(corr_p);
  free(block);
  free(mom);
  qcd_destroyPropagator(&uprop);
  qcd_destroyPropagator(&dprop);
  qcd_destroyGeometry(&geo);
  MPI_Finalize();
}//end main
