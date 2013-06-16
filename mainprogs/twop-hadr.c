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

#define N_MESONS 10      // 1, \gamma_5, \gamma_\mu, \gamma_5 \gamma_\mu
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
      
  int t_start, t_stop;
  sscanf(qcd_getParam("<t>", params, params_len),"%d %d", &t_start, &t_stop);
  if(myid==0) 
    printf(" Got sink time slices: %d ... %d\n", t_start, t_stop);

  char gauge_name[MAX_STRING];
  strcpy(gauge_name,qcd_getParam("<cfg_name>",params,params_len));
  if(myid==0) 
    printf("Got conf name: %s\n",gauge_name);
          
   
  int nt_sink = t_stop - t_start + 1;

  int n_two_points;
  sscanf(qcd_getParam("<numb_twop>", params, params_len),"%d",&n_two_points);

  
  int x_src[n_two_points][4];
  char uprop_name[n_two_points][MAX_STRING];  
  char dprop_name[n_two_points][MAX_STRING];
  char corr_nucleons_name[n_two_points][MAX_STRING];
  char corr_mesons_name[n_two_points][MAX_STRING];
  for(int i=0; i<n_two_points; i++) {
    char tag[MAX_STRING];
    sprintf(tag,"<source_pos_%d_txyz>", i);
    sscanf(qcd_getParam(tag, params, params_len),"%d %d %d %d",&x_src[i][0],&x_src[i][1],&x_src[i][2],&x_src[i][3]);
    if(myid==0) 
      printf("Got source coords %d: %d %d %d %d\n", i, x_src[i][0], x_src[i][1], x_src[i][2], x_src[i][3]);

    sprintf(tag,"<propagator_u_%d>",i);
    strcpy(uprop_name[i],qcd_getParam(tag, params, params_len));
    if(myid==0) 
      printf("Got propagator file name %d: %s\n", i, uprop_name[i]);
   
    sprintf(tag,"<propagator_d_%d>",i);
    strcpy(dprop_name[i],qcd_getParam(tag,params,params_len));
    if(myid==0) 
      printf("Got propagator file name %d: %s\n", i, dprop_name[i]);

    sprintf(tag,"<nucleon_corr_name_%d>",i);
    strcpy(corr_nucleons_name[i],qcd_getParam(tag,params,params_len));
    if(myid==0) 
      printf("Got output file name %d: %s\n",i,corr_nucleons_name[i]);

    sprintf(tag,"<mesons_corr_name_%d>",i);
    strcpy(corr_mesons_name[i],qcd_getParam(tag,params,params_len));
    if(myid==0) 
      printf("Got output file name %d: %s\n",i,corr_mesons_name[i]);
  }
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

  double t0 = MPI_Wtime();
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
  t0 = MPI_Wtime() - t0;
  if(myid==0)
    printf("done APE3d in %g sec\n",t0);  

  plaq = qcd_calculatePlaquette(&uAPE);
  if(myid==0) 
    printf("plaquette = %e\n",plaq);

  qcd_destroyGaugeField(&u);

  qcd_propagator uprop, dprop;
  qcd_propagator uprop_sm, dprop_sm;
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

  ierr = qcd_initPropagator(&uprop_sm, &geo);
  if(ierr) {
    if(myid == 0)
      fprintf(stderr, " error allocating smeared up-prop memory\n");
      
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  ierr = qcd_initPropagator(&dprop_sm, &geo);
  if(ierr) {
    if(myid == 0)
      fprintf(stderr, " error allocating smeared down-prop memory\n");
      
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  ierr = qcd_initVector(&vec, &geo);
  if(ierr) {
    if(myid == 0)
      fprintf(stderr, " error allocating temporary vector memory\n");
      
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
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
  
  qcd_complex_16 (*block_nucl)[4][2][16] = malloc(geo.lV3*nt_sink*4*2*16*sizeof(qcd_complex_16));
  if(block_nucl == NULL) {
    if(myid == 0)
      fprintf(stderr, " error allocating block_nucl memory\n");
    
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  qcd_complex_16 (*block_meson)[2][2][N_MESONS] = malloc(geo.lV3*nt_sink*N_MESONS*2*2*sizeof(qcd_complex_16));
  if(block_meson == NULL) {
    if(myid == 0)
      fprintf(stderr, " error allocating block_meson memory\n");
    
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
  
  for(int itwop=0; itwop<n_two_points; itwop++) {
    //##############################################################################
    // load propagators
    t0 = MPI_Wtime();
    ierr = qcd_getPropagator(uprop_name[itwop],qcd_PROP_LIME, &uprop);
    if(ierr) {
      if(myid == 0)
	fprintf(stderr, " %s: error reading file\n", uprop_name[itwop]);
      
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    t0 = MPI_Wtime() - t0;
    if(myid==0) 
      printf("up propagator loaded in %g sec\n", t0);
    
    t0 = MPI_Wtime();
    ierr = qcd_getPropagator(dprop_name[itwop],qcd_PROP_LIME, &dprop);
    if(ierr) {
      if(myid == 0)
	fprintf(stderr, " %s: error reading file\n", dprop_name[itwop]);
      
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    t0 = MPI_Wtime() - t0;
    if(myid==0) 
      printf("down propagator loaded in %g sec\n", t0);
    
    t0 = MPI_Wtime(); 
    //################################################################################
    // transform propagators to basis with theta-periodic boundaries in the temporal direction
    
    for(int lt=0; lt<geo.lL[0]; lt++) {
      int t = lt + geo.Pos[0] * geo.lL[0];     
      qcd_complex_16 phase_factor = (qcd_complex_16) {
	cos(theta[0]*t/geo.L[0]),
	sin(theta[0]*t/geo.L[0])
      };     
      qcd_mulPropagatorC3d(&uprop, phase_factor, (t+x_src[itwop][0]) % geo.L[0]);
      qcd_mulPropagatorC3d(&dprop, phase_factor, (t+x_src[itwop][0]) % geo.L[0]);
    }
    t0 = MPI_Wtime() - t0;
    if(myid==0)
      printf("propagators transformed to basis with theta-periodic boundary conditions in %g sec\n",t0);
    
    //################################################################################
    // gaussian smearing of propagators (only the time-slices that will be used)
    t0 = MPI_Wtime();
    for(qcd_uint_2 mu=0;mu<4;mu++)
      for(qcd_uint_2 c1=0;c1<3;c1++) {
	qcd_copyVectorPropagator(&vec,&uprop,mu,c1);
	for(int i=0; i<nsmear; i++) {
	  qcd_gaussIteration3dAll(&vec,&uAPE,alpha,(qcd_uint_2)(i==0));
	}
	qcd_copyPropagatorVector(&uprop_sm,&vec,mu,c1);
	
	qcd_copyVectorPropagator(&vec,&dprop,mu,c1);
	for(int i=0; i<nsmear; i++) {
	  qcd_gaussIteration3dAll(&vec,&uAPE,alpha,(qcd_uint_2)(i==0));
	}
	qcd_copyPropagatorVector(&dprop_sm,&vec,mu,c1);
      }

    t0 = MPI_Wtime()-t0;
    if(myid==0)
      printf("propagators smeared in %g sec\n", t0);
    
    //################################################################################
    // calculate the 'blocks'
    
    for(int v=0; v<geo.lV3*nt_sink; v++)
      for(int chan=0; chan<4; chan++)
	for(int isosp=0; isosp<2; isosp++)
	  for(int mu=0; mu<16; mu++)
	    block_nucl[v][chan][isosp][mu] = (qcd_complex_16) {0,0};

    for(int v=0; v<geo.lV3*nt_sink; v++)
      for(int ud=0; ud<2; ud++)
	for(int ism=0; ism<2; ism++)
	  for(int imes=0; imes<N_MESONS; imes++)
	    block_meson[v][ism][ud][imes] = (qcd_complex_16) {0,0};

    t0 = MPI_Wtime();    
    for(int t=t_start; t<=t_stop; t++) {
      int lt = ((t+x_src[itwop][0])%geo.L[0]) - geo.Pos[0]*geo.lL[0];
      int it = t-t_start;
      
      qcd_propagator *q1[4] = {&uprop_sm, &dprop_sm, &uprop, &dprop};
      qcd_propagator *q2[2] = {&dprop_sm, &uprop_sm};
      
      for(int ism=0; ism<2; ism++)
	for(int ud=0; ud<2; ud++) {
#pragma omp parallel for
	  for(int v3=0; v3<geo.lV3; v3++) {
	    int v = lt + v3*(int)geo.lL[0];
	    int iv = v3 + it*(int)geo.lV3;
	    qcd_complex_16 C[NS][NS][NC][NC];
	    qcd_complex_16 (*x)[NS][NC][NC] = q1[ud + ism*2]->D[v];
	    qcd_complex_16 aux[NS][NS][NC][NC];
	    qcd_complex_16 y[NS][NS][NC][NC];
	    qcd_complex_16 y0[NS][NS][NC][NC];
	    qcd_complex_16 y1[NS][NS][NC][NC];

	    for(int imes=0; imes<N_MESONS; imes++) {
	      prop_G_dag(y, x);
	      switch(imes) {
	      case 0:  /* 1 */
		memcpy(y0, x, sizeof(qcd_complex_16)*NS*NS*NC*NC);
		memcpy(y1, y, sizeof(qcd_complex_16)*NS*NS*NC*NC);
		break;
	      case 1:  /* g5 */
		prop_gamma_5_G(y0, x);
		prop_gamma_5_G(y1, y);
		break;
	      case 2:  /* gx */
		prop_gamma_x_G(y0, x);
		prop_gamma_x_G(y1, y);
		break;
	      case 3:  /* gy */
		prop_gamma_y_G(y0, x);
		prop_gamma_y_G(y1, y);
		break;
	      case 4:  /* gz */
		prop_gamma_z_G(y0, x);
		prop_gamma_z_G(y1, y);
		break;
	      case 5:  /* gt */
		prop_gamma_t_G(y0, x);
		prop_gamma_t_G(y1, y);
		break;
	      case 6:  /* g5gx */
		prop_gamma_x_G(aux, x);
		prop_gamma_5_G(y0, aux);
		prop_gamma_x_G(aux, y);
		prop_gamma_5_G(y1, aux);
		break;
	      case 7:  /* g5gy */
		prop_gamma_y_G(aux, x);
		prop_gamma_5_G(y0, aux);
		prop_gamma_y_G(aux, y);
		prop_gamma_5_G(y1, aux);
		break;
	      case 8:  /* g5gz */
		prop_gamma_z_G(aux, x);
		prop_gamma_5_G(y0, aux);
		prop_gamma_z_G(aux, y);
		prop_gamma_5_G(y1, aux);
		break;
	      case 9:  /* g5gt */
		prop_gamma_t_G(aux, x);
		prop_gamma_5_G(y0, aux);
		prop_gamma_t_G(aux, y);
		prop_gamma_5_G(y1, aux);
		break;
	      }
	
	      prop_G_G(C, y1 ,y0);
	      for(int sp=0; sp<NS; sp++)
		for(int col=0; col<NC; col++) {
		  block_meson[iv][ism][ud][imes].re += C[sp][sp][col][col].re;
		  block_meson[iv][ism][ud][imes].im += C[sp][sp][col][col].im;
		}	      
	    }
	  }
	}
	  
    
      for(int isosp=0; isosp<2; isosp++) {
#pragma omp parallel for
	for(int v3=0; v3<geo.lV3; v3++) {
	  int v = lt + v3*(int)geo.lL[0];
	  int iv = v3 + it*(int)geo.lV3;
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
	      block_nucl[iv][0][isosp][s1 + s0*NS].re = A[s0][s1].re + B[s0][s1].re;
	      block_nucl[iv][0][isosp][s1 + s0*NS].im = A[s0][s1].im + B[s0][s1].im;
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
	      block_nucl[iv][1][isosp][s1 + s0*NS].re = A[s0][s1].re + B[s0][s1].re;
	      block_nucl[iv][1][isosp][s1 + s0*NS].im = A[s0][s1].im + B[s0][s1].im;
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
	      block_nucl[iv][2][isosp][s1 + s0*NS].re = A[s0][s1].re + B[s0][s1].re;
	      block_nucl[iv][2][isosp][s1 + s0*NS].im = A[s0][s1].im + B[s0][s1].im;
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
	      block_nucl[iv][3][isosp][s1 + s0*NS].re = A[s0][s1].re + B[s0][s1].re;
	      block_nucl[iv][3][isosp][s1 + s0*NS].im = A[s0][s1].im + B[s0][s1].im;
	    }
	  
	  
	}
      }
    }
    t0 = MPI_Wtime() - t0;
    if(myid==0)
      printf(" Done twop-%d contractions in %g sec\n", itwop, t0);



    qcd_complex_16 (*corr_nucl_p)[4][2][16] = NULL, (*corr_nucl)[4][2][16] = NULL;
    corr_nucl_p = malloc(sizeof(qcd_complex_16)*2*n_mom*16*4);
    if(myid == 0)
      corr_nucl = malloc(sizeof(qcd_complex_16)*2*n_mom*nt_sink*16*4);
    
    for(int it=0; it<nt_sink; it++) {
#pragma omp parallel for
      for(int m=0; m<n_mom; m++) {
	for(int mu=0; mu<16; mu++)
	  for(int chan=0; chan<4; chan++)
	    for(int isosp=0; isosp<2; isosp++)
	      corr_nucl_p[m][chan][isosp][mu] = (qcd_complex_16) {0,0};
	
	for(int lx=0; lx<geo.lL[1]; lx++)
	  for(int ly=0; ly<geo.lL[2]; ly++)
	    for(int lz=0; lz<geo.lL[3]; lz++) {
	      unsigned long int v3 = qcd_LEXIC0(lx,ly,lz,geo.lL);
	      int x = lx + geo.Pos[1]*geo.lL[1] - x_src[itwop][1];
	      int y = ly + geo.Pos[2]*geo.lL[2] - x_src[itwop][2];
	      int z = lz + geo.Pos[3]*geo.lL[3] - x_src[itwop][3];
	      double angle =
		((double) mom[m][0]*x/(double)geo.L[1] +
		 (double) mom[m][1]*y/(double)geo.L[2] +
		 (double) mom[m][2]*z/(double)geo.L[3])*2*M_PI;
	      
	      qcd_complex_16 phase = (qcd_complex_16) {cos(angle), -sin(angle)}; //TABULATE FOR LARGE SPEEDUP!!!
	      for(int chan=0; chan<4; chan++)
		for(int isosp=0; isosp<2; isosp++)
		  for(int mu=0; mu<16; mu++)
		    corr_nucl_p[m][chan][isosp][mu] = qcd_CADD(corr_nucl_p[m][chan][isosp][mu], 
							       qcd_CMUL(block_nucl[v3+it*geo.lV3][chan][isosp][mu], phase));
	    }
      }
      MPI_Reduce(corr_nucl_p, &corr_nucl[n_mom*it][0][0][0].re, 2*n_mom*2*4*16, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }//end t-loop
    
    if(myid==0) {
      FILE *fp = fopen(corr_nucleons_name[itwop], "w");
      if(fp == NULL) {
	fprintf(stderr, " %s: error opening file for writing\n", corr_nucleons_name[itwop]);
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
			corr_nucl[m + it*n_mom][chan][isosp][mu*4 + 0].re, corr_nucl[m + it*n_mom][chan][isosp][mu*4 + 0].im,
			corr_nucl[m + it*n_mom][chan][isosp][mu*4 + 1].re, corr_nucl[m + it*n_mom][chan][isosp][mu*4 + 1].im,
			corr_nucl[m + it*n_mom][chan][isosp][mu*4 + 2].re, corr_nucl[m + it*n_mom][chan][isosp][mu*4 + 2].im,
			corr_nucl[m + it*n_mom][chan][isosp][mu*4 + 3].re, corr_nucl[m + it*n_mom][chan][isosp][mu*4 + 3].im,
			isosp == 0 ? "ppm" : "pmm",
			chan_tag[chan]);
	      }
      fclose(fp);
    }
    //#####################################################################
    // clean up
    if(myid == 0)
      free(corr_nucl);
    
    free(corr_nucl_p);
    
  
    qcd_complex_16 (*corr_meson_p)[2][2][N_MESONS] = NULL, (*corr_meson)[2][2][N_MESONS] = NULL;
    corr_meson_p = malloc(sizeof(qcd_complex_16)*2*n_mom*2*N_MESONS);
    if(myid == 0)
      corr_meson = malloc(sizeof(qcd_complex_16)*2*n_mom*nt_sink*2*N_MESONS);
    
    for(int it=0; it<nt_sink; it++) {
#pragma omp parallel for
      for(int m=0; m<n_mom; m++) {
	for(int ism=0; ism<2; ism++)
	  for(int ud=0; ud<2; ud++)
	    for(int imes=0; imes<N_MESONS; imes++)
	      corr_meson_p[m][ism][ud][imes] = (qcd_complex_16) {0,0};
	
	for(int lx=0; lx<geo.lL[1]; lx++)
	  for(int ly=0; ly<geo.lL[2]; ly++)
	    for(int lz=0; lz<geo.lL[3]; lz++) {
	      unsigned long int v3 = qcd_LEXIC0(lx,ly,lz,geo.lL);
	      int x = lx + geo.Pos[1]*geo.lL[1] - x_src[itwop][1];
	      int y = ly + geo.Pos[2]*geo.lL[2] - x_src[itwop][2];
	      int z = lz + geo.Pos[3]*geo.lL[3] - x_src[itwop][3];
	      double angle =
		((double) mom[m][0]*x/(double)geo.L[1] +
		 (double) mom[m][1]*y/(double)geo.L[2] +
		 (double) mom[m][2]*z/(double)geo.L[3])*2*M_PI;
	      
	      qcd_complex_16 phase = (qcd_complex_16) {cos(angle), -sin(angle)}; //TABULATE FOR LARGE SPEEDUP!!!
	      for(int ism=0; ism<2; ism++)
		for(int ud=0; ud<2; ud++)
		  for(int imes=0; imes<N_MESONS; imes++)
		    corr_meson_p[m][ism][ud][imes] = qcd_CADD(corr_meson_p[m][ism][ud][imes], 
							      qcd_CMUL(block_meson[v3+it*geo.lV3][ism][ud][imes], phase));
	    }
      }
      MPI_Reduce(corr_meson_p, &corr_meson[n_mom*it][0][0][0].re, 2*n_mom*2*N_MESONS*2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }//end t-loop
    
    if(myid==0) {
      FILE *fp = fopen(corr_mesons_name[itwop], "w");
      if(fp == NULL) {
	fprintf(stderr, " %s: error opening file for writing\n", corr_mesons_name[itwop]);
	MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
      }
      char *gamma_tag[] = {"=1=",
			   "=g5=",
			   "=gx=",
			   "=gy=",
			   "=gz=",
			   "=gt=",
			   "=g5gx=",
			   "=g5gy=",
			   "=g5gz=",
			   "=g5gt="};
      char *smear_tag[] = {"SS", "LS"};
      char *ud_tag[] = {"up", "dn"};
      for(int it=0; it<nt_sink; it++)	
	for(int m=0; m<n_mom; m++)
	  for(int imes=0; imes<N_MESONS; imes++)
	    for(int ism=0; ism<2; ism++) 
	      for(int ud=0; ud<2; ud++) {
		fprintf(fp, " %3d %+d %+d %+d %+e %+e %7s %2s %2s\n",
			it % geo.L[0],
			mom[m][0], mom[m][1], mom[m][2],
			corr_meson[m + it*n_mom][ism][ud][imes].re, corr_meson[m + it*n_mom][ism][ud][imes].im,
			gamma_tag[imes],
			smear_tag[ism],
			ud_tag[ud]);
	      }
      fclose(fp);
    }
    //#####################################################################
    // clean up
    if(myid == 0)
      free(corr_meson);
    
    free(corr_meson_p);
  }
  
  free(block_nucl);
  free(block_meson);
  free(mom);
  qcd_destroyVector(&vec);
  qcd_destroyGaugeField(&uAPE);
  qcd_destroyPropagator(&uprop);
  qcd_destroyPropagator(&dprop);
  qcd_destroyPropagator(&uprop_sm);
  qcd_destroyPropagator(&dprop_sm);
  qcd_destroyGeometry(&geo);
  MPI_Finalize();
  return 0;
}//end main
