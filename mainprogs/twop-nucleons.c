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

#define NC 3
#define NS 4

#include <qcd_prop_contract.h>
#include <qcd_prop_ops.h>
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

  /* { */
  /*   qcd_complex_16 (*y)[NS][NC][NC] = uprop.D[0]; */
  /*   qcd_complex_16 (*x)[NS][NC][NC] = dprop.D[0]; */
  /*   qcd_complex_16 q[4][4][3][3], p[4][4][3][3]; */
    
  /*   prop_contract_02(q,y,x); */

  /*   memset(p, '\0', sizeof(qcd_complex_16)*3*3*4*4); */

  /*   for(int eps0=0; eps0<6; eps0++) { */
  /*     int a0 = qcd_EPS[eps0][0]; */
  /*     int b0 = qcd_EPS[eps0][1]; */
  /*     int c0 = qcd_EPS[eps0][2];     */
  /*     int sign0 = eps0 < 3 ? +1 : -1; */
  /*     for(int eps1=0; eps1<6; eps1++) { */
  /* 	int a1 = qcd_EPS[eps1][0]; */
  /* 	int b1 = qcd_EPS[eps1][1]; */
  /* 	int c1 = qcd_EPS[eps1][2]; */
  /* 	int sign1 = eps1 < 3 ? +1 : -1; */
      
  /* 	for(int mu=0; mu<4; mu++) */
  /* 	  for(int nu=0; nu<4; nu++) {	     */
  /* 	    for(int ku=0; ku<4; ku++) { */
  /* 	      p[mu][nu][a0][a1].re += sign1*sign0*qcd_CMULR(uprop.D[0][ku][mu][b0][b1], dprop.D[0][ku][nu][c0][c1]); */
  /* 	      p[mu][nu][a0][a1].im += sign1*sign0*qcd_CMULI(uprop.D[0][ku][mu][b0][b1], dprop.D[0][ku][nu][c0][c1]); */
  /* 	    } */
  /* 	  } */
  /*     } */
  /*   } */

  /*   int c0=2,c1=0; */
  /*   for(int mu=0; mu<4; mu++) { */
  /*     for(int nu=0; nu<4; nu++) */
  /* 	printf(" %+6.2e%+6.2e*J ", p[mu][nu][c0][c1].re, p[mu][nu][c0][c1].im); */
      
  /*     printf("\n"); */
  /*   } */
  /*   printf("\n"); */
  /*   for(int mu=0; mu<4; mu++) { */
  /*     for(int nu=0; nu<4; nu++) */
  /* 	printf(" %+6.2e%+6.2e*J ", q[mu][nu][c0][c1].re, q[mu][nu][c0][c1].im); */
      
  /*     printf("\n"); */
  /*   } */
  /* } */

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

  qcd_complex_16 (*block)[4][2][16] = malloc(geo.lV3*nt_sink*4*2*16*sizeof(qcd_complex_16));
  if(block == NULL) {
    if(myid == 0)
      fprintf(stderr, " error allocating block memory\n");
    
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
  for(int v=0; v<geo.lV3*nt_sink; v++)
    for(int chan=0; chan<4; chan++)
      for(int isosp=0; isosp<2; isosp++)
	for(int mu=0; mu<16; mu++)
	  block[v][chan][isosp][mu] = (qcd_complex_16) {0,0};

  for(int t=t_start; t<=t_stop; t++) {
    int lt = ((t+x_src[0])%geo.L[0]) - geo.Pos[0]*geo.lL[0];
    int it = t-t_start;

    qcd_propagator *q1[2] = {&uprop, &dprop};
    qcd_propagator *q2[2] = {&dprop, &uprop};

    if(myid==0)
      printf("t = %3d\n",t);

    for(int isosp=0; isosp<2; isosp++) {
#pragma omp parallel for
      for(int v3=0; v3<geo.lV3; v3++) {
	int v = lt + v3*geo.lL[0];
	int iv = v3 + it*geo.lV3;
	qcd_complex_16 (*y)[NS][NC][NC] = q1[isosp]->D[v];
	qcd_complex_16 (*x)[NS][NC][NC] = q2[isosp]->D[v];
	qcd_complex_16 z[NS][NS][NC][NC];
	qcd_complex_16 aux[NS][NS][NC][NC];
	
	qcd_complex_16 A[NS][NS], B[NS][NS];

	/*
	 * \chi_1 - to - \chi_1
	 */
	qcd_complex_16 Cg5xCg5[NS][NS][NC][NC];
	prop_Cg5_G(aux, x);
	prop_G_Cg5(Cg5xCg5, aux);
	prop_contract_02(z, Cg5xCg5, y);

	for(int s0=0; s0<NS; s0++)
	  for(int s1=0; s1<NS; s1++){
	    A[s0][s1] = (qcd_complex_16){0.,0.};
	    B[s0][s1] = (qcd_complex_16){0.,0.};
	  }
	
	for(int c1=0; c1<NC; c1++)
	  for(int c0=0; c0<NC; c0++)
	    for(int s0=0; s0<NS; s0++)
	      for(int s1=0; s1<NS; s1++)
		for(int s2=0; s2<NS; s2++) {
		  A[s0][s1].re += qcd_CMULR(y[s0][s1][c0][c1], z[s2][s2][c0][c1]);
		  A[s0][s1].im += qcd_CMULI(y[s0][s1][c0][c1], z[s2][s2][c0][c1]);

		  B[s0][s1].re += qcd_CMULR(y[s0][s2][c0][c1], z[s2][s1][c0][c1]);
		  B[s0][s1].im += qcd_CMULI(y[s0][s2][c0][c1], z[s2][s1][c0][c1]);
		}

	for(int s0=0; s0<NS; s0++)
	  for(int s1=0; s1<NS; s1++) {
	    block[iv][0][isosp][s1 + s0*NS].re = A[s0][s1].re + B[s0][s1].re;
	    block[iv][0][isosp][s1 + s0*NS].im = A[s0][s1].im + B[s0][s1].im;
	  }

	/*
	 * \chi_1 - to - \chi_2
	 */
	qcd_complex_16 yg5[NS][NS][NC][NC];
	prop_G_gamma_5(yg5, y);

	qcd_complex_16 Cg5xC[NS][NS][NC][NC];
	prop_Cg5_G(aux, x);
	prop_G_C(Cg5xC, aux);
	prop_contract_13(z, Cg5xC, y);

	for(int s0=0; s0<NS; s0++)
	  for(int s1=0; s1<NS; s1++){
	    A[s0][s1] = (qcd_complex_16){0.,0.};
	    B[s0][s1] = (qcd_complex_16){0.,0.};
	  }
	
	for(int c1=0; c1<NC; c1++)
	  for(int c0=0; c0<NC; c0++)
	    for(int s0=0; s0<NS; s0++)
	      for(int s1=0; s1<NS; s1++)
		for(int s2=0; s2<NS; s2++) {
		  A[s0][s1].re += qcd_CMULR(yg5[s0][s1][c0][c1], z[s2][s2][c0][c1]);
		  A[s0][s1].im += qcd_CMULI(yg5[s0][s1][c0][c1], z[s2][s2][c0][c1]);
					       
		  B[s0][s1].re += qcd_CMULR(yg5[s2][s1][c0][c1], z[s2][s0][c0][c1]);
		  B[s0][s1].im += qcd_CMULI(yg5[s2][s1][c0][c1], z[s2][s0][c0][c1]);
		}

	for(int s0=0; s0<NS; s0++)
	  for(int s1=0; s1<NS; s1++) {
	    block[iv][1][isosp][s1 + s0*NS].re = A[s0][s1].re + B[s0][s1].re;
	    block[iv][1][isosp][s1 + s0*NS].im = A[s0][s1].im + B[s0][s1].im;
	  }

	/*
	 * \chi_2 - to - \chi_1
	 */
	qcd_complex_16 g5y[NS][NS][NC][NC];
	prop_gamma_5_G(g5y, y);

	qcd_complex_16 CxCg5[NS][NS][NC][NC];
	prop_C_G(aux, x);
	prop_G_Cg5(CxCg5, aux);
	prop_contract_02(z, CxCg5, y);

	for(int s0=0; s0<NS; s0++)
	  for(int s1=0; s1<NS; s1++){
	    A[s0][s1] = (qcd_complex_16){0.,0.};
	    B[s0][s1] = (qcd_complex_16){0.,0.};
	  }
	
	for(int c1=0; c1<NC; c1++)
	  for(int c0=0; c0<NC; c0++)
	    for(int s0=0; s0<NS; s0++)
	      for(int s1=0; s1<NS; s1++)
		for(int s2=0; s2<NS; s2++) {
		  A[s0][s1].re += qcd_CMULR(g5y[s0][s1][c0][c1], z[s2][s2][c0][c1]);
		  A[s0][s1].im += qcd_CMULI(g5y[s0][s1][c0][c1], z[s2][s2][c0][c1]);

		  B[s0][s1].re += qcd_CMULR(g5y[s0][s2][c0][c1], z[s2][s1][c0][c1]);
		  B[s0][s1].im += qcd_CMULI(g5y[s0][s2][c0][c1], z[s2][s1][c0][c1]);
		}

	for(int s0=0; s0<NS; s0++)
	  for(int s1=0; s1<NS; s1++) {
	    block[iv][2][isosp][s1 + s0*NS].re = A[s0][s1].re + B[s0][s1].re;
	    block[iv][2][isosp][s1 + s0*NS].im = A[s0][s1].im + B[s0][s1].im;
	  }

	/*
	 * \chi_2 - to - \chi_2
	 */
	qcd_complex_16 g5yg5[NS][NS][NC][NC];
	prop_gamma_5_G(g5yg5, yg5);

	qcd_complex_16 CxC[NS][NS][NC][NC];
	qcd_complex_16 w[NS][NS][NC][NC];
	prop_C_G(aux, x);
	prop_G_C(CxC, aux);
	prop_contract_02(z, CxC, y);
	prop_contract_02(w, CxC, yg5);

	for(int s0=0; s0<NS; s0++)
	  for(int s1=0; s1<NS; s1++){
	    A[s0][s1] = (qcd_complex_16){0.,0.};
	    B[s0][s1] = (qcd_complex_16){0.,0.};
	  }
	
	for(int c1=0; c1<NC; c1++)
	  for(int c0=0; c0<NC; c0++)
	    for(int s0=0; s0<NS; s0++)
	      for(int s1=0; s1<NS; s1++)
		for(int s2=0; s2<NS; s2++) {
		  A[s0][s1].re += qcd_CMULR(g5yg5[s0][s1][c0][c1], z[s2][s2][c0][c1]);
		  A[s0][s1].im += qcd_CMULI(g5yg5[s0][s1][c0][c1], z[s2][s2][c0][c1]);

		  B[s0][s1].re += qcd_CMULR(g5y[s0][s2][c0][c1], w[s2][s1][c0][c1]);
		  B[s0][s1].im += qcd_CMULI(g5y[s0][s2][c0][c1], w[s2][s1][c0][c1]);
		}

	for(int s0=0; s0<NS; s0++)
	  for(int s1=0; s1<NS; s1++) {
	    block[iv][3][isosp][s1 + s0*NS].re = A[s0][s1].re + B[s0][s1].re;
	    block[iv][3][isosp][s1 + s0*NS].im = A[s0][s1].im + B[s0][s1].im;
	  }
      
	
      }
    }
  }
  qcd_complex_16 (*corr_p)[4][2][16] = NULL, (*corr)[4][2][16] = NULL;
  corr_p = malloc(sizeof(qcd_complex_16)*2*n_mom*16*4);
  if(myid == 0)
    corr = malloc(sizeof(qcd_complex_16)*2*n_mom*nt_sink*16*4);
   
  for(int it=0; it<nt_sink; it++) {
#pragma omp parallel for
    for(int m=0; m<n_mom; m++) {
      for(int mu=0; mu<16; mu++)
	for(int chan=0; chan<4; chan++)
	  for(int isosp=0; isosp<2; isosp++)
	    corr_p[m][chan][isosp][mu] = (qcd_complex_16) {0,0};

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
	    for(int chan=0; chan<4; chan++)
	      for(int isosp=0; isosp<2; isosp++)
		for(int mu=0; mu<16; mu++)
		  corr_p[m][chan][isosp][mu] = qcd_CADD(corr_p[m][chan][isosp][mu], 
							qcd_CMUL(block[v3+it*geo.lV3][chan][isosp][mu], phase));
	  }
    }
    MPI_Reduce(corr_p, ((double *)corr+2*it*n_mom*2*4*16), 2*n_mom*2*4*16, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  }//end t-loop
            
  if(myid==0) {
    FILE *fp = fopen(corr_p_name, "w");
    if(fp == NULL) {
      fprintf(stderr, " %s: error opening file for writing\n", corr_p_name);
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    char *chan_tag[] = {"1-1",
			"1-2",
			"2-1",
			"2-2"};
    for(int isosp=0; isosp<2; isosp++)
      for(int it=0; it<nt_sink; it++)	
	for(int m=0; m<n_mom; m++)
	  for(int chan=0; chan<4; chan++)
	    for(int mu=0; mu<4; mu++) {
	      fprintf(fp, " %3d %+d %+d %+d %+e %+e  %+e %+e  %+e %+e  %+e %+e %3s %3s\n",
		      it % geo.L[0],
		      mom[m][0], mom[m][1], mom[m][2],
		      corr[m + it*n_mom][chan][isosp][mu*4 + 0].re, corr[m + it*n_mom][chan][isosp][mu*4 + 0].im,
		      corr[m + it*n_mom][chan][isosp][mu*4 + 1].re, corr[m + it*n_mom][chan][isosp][mu*4 + 1].im,
		      corr[m + it*n_mom][chan][isosp][mu*4 + 2].re, corr[m + it*n_mom][chan][isosp][mu*4 + 2].im,
		      corr[m + it*n_mom][chan][isosp][mu*4 + 3].re, corr[m + it*n_mom][chan][isosp][mu*4 + 3].im,
		      isosp == 0 ? "ppm" : "pmm",
		      chan_tag[chan]);
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
