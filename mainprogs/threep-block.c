/* threep_idris.c
 *
 * reads forward and backward propagators
 * and creates three point functions
 *
 * 0 and 1 Derivative operators
 * vector, axial, tensor
 *
 * Tomasz Korzec 2009
 ****************************************/
 
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <qcd.h>

void
writeBlock(char fname[], qcd_complex_16 **block, int x_src[4], int nt_ins, int t_ins[nt_ins], int site_size, char *site_tags[], qcd_geometry geo)
{
  int lx = geo.lL[1], ly = geo.lL[2], lz = geo.lL[3];
  int Lx = geo.L[1], Ly = geo.L[2], Lz = geo.L[3], Lt = geo.L[0];
  int myid = geo.myid;
  unsigned long int offset;
  if(myid == 0) {
    FILE *fp = fopen(fname, "w");
    if(fp == NULL) {
      fprintf(stderr, " %s: error opening file for writing\n", fname);
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    char header[1024], line[1024];
    strcpy(header, "begin-header;\n");

    /* lattice dimensions */
    sprintf(line, "lattice dimensions (x,y,z,t) = %3d,%3d,%3d,%3d;\n", Lx, Ly, Lz, Lt);
    strcat(header, line);

    /* source coordinates */
    sprintf(line, "source position (x,y,z,t) = %3d,%3d,%3d,%3d;\n", x_src[1], x_src[2], x_src[3], x_src[0]);
    strcat(header, line);

    /* t-insertions lines of header */
    sprintf(line, "number of insertions = %d;\n", nt_ins);
    strcat(header, line);
    sprintf(line, "t-insertions =");
    for(int t=0; t<nt_ins; t++) {
      char st[1024];
      sprintf(st, " %d,", t_ins[t]);
      strcat(line, st);
    }
    line[strlen(line)-1] = '\0'; /* deletes last comma */
    strcat(line, ";\n");
    strcat(header, line);
    
    /* block-sites lines of header */
    sprintf(line, "lattice-site size = %d;\n", site_size);
    strcat(header, line);
    sprintf(line, "lattice-site tags =");
    for(int s=0; s<site_size; s++) {
      char st[1024];
      sprintf(st, " %s,", site_tags[s]);
      strcat(line, st);
    }
    line[strlen(line)-1] = '\0'; /* deletes last comma */
    strcat(line, ";\n");
    strcat(header, line);
    
    strcat(header, "end-header;\n");
    fprintf(fp, header);

    offset = ftell(fp);
  }
  MPI_Bcast(&offset, sizeof(unsigned long int), MPI_BYTE, 0, MPI_COMM_WORLD);

  //header written. Now use MPI-I/O to write binary data
  int sizes[5], lsizes[5], starts[5];
  sizes[0]=nt_ins;
  sizes[1]=Lz;
  sizes[2]=Ly;
  sizes[3]=Lx;
  sizes[4]=site_size*2;
  lsizes[0]=nt_ins;
  lsizes[1]=lz;
  lsizes[2]=ly;
  lsizes[3]=lx;
  lsizes[4]=sizes[4];
  starts[0]=0;
  starts[1]=geo.Pos[3]*lsizes[1];
  starts[2]=geo.Pos[2]*lsizes[2];
  starts[3]=geo.Pos[1]*lsizes[3];                  
  starts[4]=0;
  
  MPI_Datatype subblock;
  MPI_Type_create_subarray(5,sizes,lsizes,starts,MPI_ORDER_C,MPI_DOUBLE,&subblock);
  MPI_Type_commit(&subblock);

  MPI_File fh;
  MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
  MPI_File_set_view(fh, offset, MPI_DOUBLE, subblock, "native", MPI_INFO_NULL);

  int chunksize = site_size*sizeof(qcd_complex_16);
  double *buffer = malloc(chunksize*geo.lV3*nt_ins);
  if(buffer==NULL) {
    fprintf(stderr," malloc() returned NULL in %s, out of memory?\n",__func__);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  for(int t=0; t<nt_ins; t++)
    for(int z=0; z<lz; z++)
      for(int y=0; y<ly; y++)
	for(int x=0; x<lx; x++)	  
	  for(int s=0; s<site_size; s++) {
	    int ix = 2*(s + site_size*(x + lx*(y + ly*(z + lz*t))));
	    int iv = x + lx*(y + ly*(z + lz*t));
	    buffer[ix+0] = block[s][iv].re;
	    buffer[ix+1] = block[s][iv].im;
	  }

  if(!qcd_isBigEndian())
    qcd_swap_8((double *) buffer, site_size*2*lx*ly*lz*nt_ins);
   
  MPI_File_write_all(fh, buffer, site_size*2*lx*ly*lz*nt_ins, MPI_DOUBLE, MPI_STATUS_IGNORE);
  free(buffer);
  MPI_File_close(&fh);
  MPI_Type_free(&subblock);
  return;
}



int 
main(int argc,char* argv[])
{
  int numprocs, myid;             
  //set up MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);         // num. of processes taking part in the calculation
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);             // each process gets its ID
  
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int namelen;
  MPI_Get_processor_name(processor_name,&namelen); // 
  
  qcd_complex_16 g5sigmu0[5][4][4];            // gamma_5 * [gamma_mu, gamma_0] *1/2   
  for(int i=0; i<5; i++)
    for(int mu=0; mu<4; mu++)
      for(int nu=0; nu<4; nu++)	{           
	g5sigmu0[i][mu][nu]= (qcd_complex_16){0,0};
	for(int ku=0; ku<4; ku++)
	  for(int lu=0; lu<4; lu++) {
	    g5sigmu0[i][mu][nu] = qcd_CADD(g5sigmu0[i][mu][nu],
					   qcd_CMUL(qcd_CMUL(qcd_GAMMA[5][mu][ku],
							     qcd_GAMMA[i][ku][lu]),
						    qcd_GAMMA[0][lu][nu]));
	    g5sigmu0[i][mu][nu] = qcd_CSUB(g5sigmu0[i][mu][nu],
					   qcd_CMUL(qcd_CMUL(qcd_GAMMA[5][mu][ku],
							     qcd_GAMMA[0][ku][lu]),
						    qcd_GAMMA[i][lu][nu]));
	  }
	g5sigmu0[i][mu][nu] = qcd_CSCALE(g5sigmu0[i][mu][nu],0.5);
      }

  //////////////////// READ INPUT FILE /////////////////////////////////////////////
      
  if(argc!=2) {
    if(myid==0) fprintf(stderr,"No input file specified\n");
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  char param_name[qcd_MAX_STRING_LENGTH];
  int params_len;
  char *params = NULL;
  strcpy(param_name,argv[1]);
  if(myid==0) {
    printf("opening input file %s\n",param_name);
    params = qcd_getParams(param_name, &params_len);
    if(params == NULL) {
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
  }

  MPI_Bcast(&params_len, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if(myid!=0) 
    params = (char*) malloc(params_len*sizeof(char));

  MPI_Bcast(params, params_len, MPI_CHAR, 0, MPI_COMM_WORLD);

  qcd_uint_2 L[4];
  qcd_uint_2 P[4];   
  sscanf(qcd_getParam("<processors_txyz>",params,params_len),"%hd %hd %hd %hd",&P[0], &P[1], &P[2], &P[3]);
  sscanf(qcd_getParam("<lattice_txyz>",params,params_len),"%hd %hd %hd %hd",&L[0], &L[1], &L[2], &L[3]);

  qcd_real_8 theta[4] = {M_PI,0.0,0.0,0.0};
  qcd_geometry geo;
  if(qcd_initGeometry(&geo, L, P, theta, myid, numprocs)) 
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
   
  if(P[0] != 1) {
    if(myid==0)
      fprintf(stderr, " Number of processes along T should be 1 (got %d)\n", P[0]);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  if(myid==0) 
    printf(" Local lattice: %i x %i x %i x %i\n",geo.lL[0],geo.lL[1],geo.lL[2],geo.lL[3]);
   
  int t_start, t_stop;
  sscanf(qcd_getParam("<t>",params,params_len),"%d %d",&t_start, &t_stop);
  if(myid==0) 
    printf("Got insertion time slices: %d ... %d\n",t_start,t_stop);
                     	  
  int x_src[4];
  sscanf(qcd_getParam("<source_pos_txyz>",params,params_len),"%d %d %d %d",&x_src[0],&x_src[1],&x_src[2],&x_src[3]);
  if(myid==0) 
    printf("Got source coords: %d %d %d %d\n",x_src[0],x_src[1],x_src[2],x_src[3]);
   
  int t_sink;
  sscanf(qcd_getParam("<t_sink>",params,params_len),"%d",&t_sink);
  if(myid==0) 
    printf("Got sink time slice: %d\n",t_sink);
  
  if(t_sink >= L[0]) {
    if(myid==0) 
      fprintf(stderr, 
	      " Error: t_sink (=%d) >= L[0] (=%d),\n t_sink should be in [0, L[0]), did you forget to mod(t_sink, L[0]) ?\n", 
	      t_sink, L[0]);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
   
  char prop_name[qcd_MAX_STRING_LENGTH];
  strcpy(prop_name,qcd_getParam("<propagator>",params,params_len));
  if(myid==0) 
    printf("Got propagator file name: %s\n",prop_name);

  char bprop_name[qcd_MAX_STRING_LENGTH];
  strcpy(bprop_name,qcd_getParam("<seq_prop>",params,params_len));
  if(myid==0) 
    printf("Got sequential propagator file name: %s\n",bprop_name);

  char gauge_name[qcd_MAX_STRING_LENGTH];
  strcpy(gauge_name,qcd_getParam("<cfg_name>",params,params_len));
  if(myid==0) 
    printf("Got conf name: %s\n",gauge_name);
          
  char blockloc_p_name[qcd_MAX_STRING_LENGTH];  
  strcpy(blockloc_p_name,qcd_getParam("<block_name_p_local>",params,params_len));
  if(myid==0)
    printf("Got output file name: %s\n",blockloc_p_name);

  char blockloc_s_name[qcd_MAX_STRING_LENGTH];  
  strcpy(blockloc_s_name,qcd_getParam("<block_name_s_local>",params,params_len));
  if(myid==0)
    printf("Got output file name: %s\n",blockloc_s_name);

  char blockloc_v_name[qcd_MAX_STRING_LENGTH];  
  strcpy(blockloc_v_name,qcd_getParam("<block_name_v_local>",params,params_len));
  if(myid==0)
    printf("Got output file name: %s\n",blockloc_v_name);

  char blocknoe_v_name[qcd_MAX_STRING_LENGTH];  
  strcpy(blocknoe_v_name,qcd_getParam("<block_name_v_noether>",params,params_len));
  if(myid==0)
    printf("Got output file name: %s\n",blocknoe_v_name);

  char blockloc_a_name[qcd_MAX_STRING_LENGTH];  
  strcpy(blockloc_a_name,qcd_getParam("<block_name_a_local>",params,params_len));
  if(myid==0)
    printf("Got output file name: %s\n",blockloc_a_name);

  char block_vD_name[qcd_MAX_STRING_LENGTH];    
  strcpy(block_vD_name,qcd_getParam("<block_name_vD>",params,params_len));
  if(myid==0)
    printf("Got output file name: %s\n",block_vD_name);

  char block_aD_name[qcd_MAX_STRING_LENGTH];     
  strcpy(block_aD_name,qcd_getParam("<block_name_aD>",params,params_len));
  if(myid==0)
    printf("Got output file name: %s\n",block_aD_name);

  char block_d1_name[qcd_MAX_STRING_LENGTH];        
  strcpy(block_d1_name,qcd_getParam("<block_name_d1>",params,params_len));
  if(myid==0)
    printf("Got output file name: %s\n",block_d1_name);

  free(params);         
  //#####################################################################   
  int ierr = 0;

  qcd_propagator prop;
  ierr = qcd_initPropagator(&prop, &geo);
  if(ierr != 0) {
    fprintf(stderr, " qcd_initPropagator() returned in error, out of memory?\n");
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  qcd_propagator backprop;
  ierr = qcd_initPropagator(&backprop, &geo);
  if(ierr != 0) {
    fprintf(stderr, " qcd_initPropagator() returned in error, out of memory?\n");
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  qcd_gaugeField u;
  ierr = qcd_initGaugeField(&u, &geo);
  if(ierr != 0) {
    fprintf(stderr, " qcd_initGaugeField() returned in error, out of memory?\n");
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
  int nt_ins = t_stop - t_start + 1;
  qcd_complex_16 **block_s = (qcd_complex_16 **) malloc(1*sizeof(qcd_complex_16 *));
  qcd_complex_16 **block_p = (qcd_complex_16 **) malloc(1*sizeof(qcd_complex_16 *));
  qcd_complex_16 **block_n = (qcd_complex_16 **) malloc(4*sizeof(*block_n));
  qcd_complex_16 **block_l = (qcd_complex_16 **) malloc(4*sizeof(*block_l));
  qcd_complex_16 **block_a = (qcd_complex_16 **) malloc(4*sizeof(*block_a));
  qcd_complex_16 **block_vD = (qcd_complex_16 **) malloc(16*sizeof(*block_vD));
  qcd_complex_16 **block_aD = (qcd_complex_16 **) malloc(16*sizeof(*block_aD));
  qcd_complex_16 **block_d1 = (qcd_complex_16 **) malloc(16*sizeof(*block_d1));

  block_s[0] = (qcd_complex_16 *) malloc(geo.lV3*nt_ins*sizeof(qcd_complex_16));
  if(block_s[0] == NULL) {
    fprintf(stderr," malloc() returned NULL, out of memmory?\n");
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }  

  block_p[0] = (qcd_complex_16 *) malloc(geo.lV3*nt_ins*sizeof(qcd_complex_16));
  if(block_p[0] == NULL) {
    fprintf(stderr," malloc() returned NULL, out of memmory?\n");
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }  

  for(int i=0; i<4; i++) {
    block_n[i] = (qcd_complex_16 *) malloc(geo.lV3*nt_ins*sizeof(qcd_complex_16));
    if(block_n[i] == NULL) {
      fprintf(stderr," malloc() returned NULL, out of memmory?\n");
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    block_l[i] = (qcd_complex_16 *) malloc(geo.lV3*nt_ins*sizeof(qcd_complex_16));
    if(block_l[i]==NULL) {
      fprintf(stderr," malloc() returned NULL, out of memmory?\n");
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    block_a[i] = (qcd_complex_16 *) malloc(geo.lV3*nt_ins*sizeof(qcd_complex_16));
    if(block_a[i]==NULL) {
      fprintf(stderr," malloc() returned NULL, out of memmory?\n");
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
  }
   
  for(int i=0; i<16; i++) {
    block_vD[i] = (qcd_complex_16*) malloc(geo.lV3*nt_ins*sizeof(qcd_complex_16));
    if(block_vD[i]==NULL) {
      fprintf(stderr," malloc() returned NULL, out of memmory?\n");
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
     
    block_aD[i] = (qcd_complex_16*) malloc(geo.lV3*nt_ins*sizeof(qcd_complex_16));
    if(block_aD[i]==NULL) {
      fprintf(stderr," malloc() returned NULL, out of memmory?\n");
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
     
    block_d1[i] = (qcd_complex_16*) malloc(geo.lV3*nt_ins*sizeof(qcd_complex_16));
    if(block_d1[i]==NULL) {
      fprintf(stderr," malloc() returned NULL, out of memmory?\n");
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
  }

  //##############################################################################
  // load gauge-field
  if(myid == 0)
    printf("reading gauge-field: %s\n", gauge_name);

  ierr = qcd_getGaugeField(gauge_name,qcd_GF_LIME,&u);
  if(ierr != 0) {
    if(myid == 0)
      fprintf(stderr, " %s: error reading gauge-filed\n", gauge_name);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  qcd_real_8 plaq = qcd_calculatePlaquette(&u);
  if(myid==0)
    printf("plaquette = %e\n",plaq);
  
  qcd_communicateGaugePM(&u);
  qcd_waitall(&geo);

  if(myid == 0)
    printf("reading propagator: %s\n", prop_name);

  ierr = qcd_getPropagator(prop_name,qcd_PROP_LIME, &prop);
  if(ierr != 0) {
    if(myid == 0)
      fprintf(stderr, " %s: error reading propagator\n", prop_name);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
  
  if(myid == 0)
    printf("reading propagator: %s\n", bprop_name);

  ierr = qcd_getPropagator(bprop_name,qcd_PROP_LIME, &backprop);
  if(ierr != 0) {
    if(myid == 0)
      fprintf(stderr, " %s: error reading backward propagator\n", bprop_name);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
  
  //################################################################################
  // transform propagators to basis with theta-periodic boundaries in the temporal direction
  for(int lt=0; lt<geo.lL[0]; lt++) {
    int t = lt + geo.Pos[0] * geo.lL[0];
    qcd_complex_16 phase_factor = (qcd_complex_16) {cos(theta[0]*t/geo.L[0]),sin(theta[0]*t/geo.L[0])};
    // the backward-props get the complex-conjugated phases, so after the conjugation
    // in the next step they will be correct
    qcd_complex_16 phase_factor_b = (qcd_complex_16) {
      cos(theta[0]*(t_sink-(t+x_src[0])+2*(t_sink-x_src[0]))/geo.L[0]),
      -sin(theta[0]*(t_sink-(t+x_src[0])+2*(t_sink-x_src[0]))/geo.L[0])
    };
    qcd_mulPropagatorC3d(&prop, phase_factor, (t+x_src[0]) % geo.L[0]);
    qcd_mulPropagatorC3d(&backprop, phase_factor_b, (t+x_src[0]) % geo.L[0]);
  }
  
  qcd_waitall(&geo);
  qcd_communicatePropagatorPM(&prop);
  
  //#####################################################################   
  // complex-conjugate the backward propagators, multiply by gamma_5. Works only in chiral gamma basis
  qcd_conjPropagator(&backprop);
  qcd_gamma5Propagator(&backprop);   
  
  qcd_waitall(&geo);
  qcd_communicatePropagatorPM(&backprop);
  qcd_waitall(&geo);

  for(int i=0; i<geo.lV3*nt_ins; i++) {   //set blocks to zero  
    block_s[0][i] = (qcd_complex_16) {0,0};
    block_p[0][i] = (qcd_complex_16) {0,0};
    for(int mu=0; mu<4; mu++) {
      block_n[mu][i]= (qcd_complex_16) {0,0};
      block_l[mu][i]= (qcd_complex_16) {0,0};
      block_a[mu][i]= (qcd_complex_16) {0,0};
    }
    for(int mu=0; mu<16; mu++) {      
      block_vD[mu][i]= (qcd_complex_16) {0,0};
      block_aD[mu][i]= (qcd_complex_16) {0,0};
      block_d1[mu][i]= (qcd_complex_16) {0,0};
    }
  }


  for(int t=t_start; t<=t_stop; t++) {
    int it = t-t_start;
    if(myid==0) 
      printf("t = %2d\n", t);
        
    /* 
       Dirty fix to flip sign in case the
       three-point function crosses the
       end of the lattice.
       
       A better fix would be to ask for "t_src" and "source-sink separation"
       in the input parameters, and work it out from there.
    */
    double sign = +1;
    if(t_sink < x_src[0]) {
      sign = -1;
    }
    
    /* 
       This should just evaluate to lt = (t+x_src[0])%geo.L[0], since we've guaranteed
       that there is only 1 process along t-axis (P[0] = 1) 
    */
    int lt = ((t + x_src[0]) % geo.L[0]) - geo.Pos[0]*geo.lL[0];
#pragma omp parallel for
    for(int v3=0; v3<geo.lV3; v3++) {   
      qcd_complex_16 backfor;                      // backward-prop x forward-prop partially traced
      qcd_complex_16 bdfmu[4][4][4];               // stores backward-prop D_mu forward-prop
      
      int v = lt + v3*geo.lL[0];         
      int iv = v3 + it*geo.lV3;
      for(int id1=0; id1<4; id1++)
	for(int id3=0; id3<4; id3++) { 
	  for(int mu=0; mu<4; mu++) {
	    int vp1 = geo.plus[v][mu];              
	    int vm1 = geo.minus[v][mu];

	    bdfmu[mu][id1][id3] = (qcd_complex_16) {0,0};
	    for(int id2=0; id2<4; id2++)
	      for(int ic1=0; ic1<3; ic1++)
		for(int ic2=0; ic2<3; ic2++)
                  for(int ic3=0; ic3<3; ic3++) {
		    // x x x+mu
		    bdfmu[mu][id1][id3] = qcd_CADD(bdfmu[mu][id1][id3],
						   qcd_CMUL(backprop.D[v][id1][id2][ic1][ic2],
							    qcd_CMUL(u.D[v][mu][ic1][ic3],
								     prop.D[vp1][id3][id2][ic3][ic2])));
		    // x  x-mu  x-mu
		    bdfmu[mu][id1][id3] = qcd_CSUB(bdfmu[mu][id1][id3],
						   qcd_CMUL(backprop.D[v][id1][id2][ic1][ic2],
							    qcd_CMUL(qcd_CONJ(u.D[vm1][mu][ic3][ic1]),
								     prop.D[vm1][id3][id2][ic3][ic2])));
		    // x+mu  x  x
		    bdfmu[mu][id1][id3] = qcd_CSUB(bdfmu[mu][id1][id3],
						   qcd_CMUL(backprop.D[vp1][id1][id2][ic1][ic2],
							    qcd_CMUL(qcd_CONJ(u.D[v][mu][ic3][ic1]),
								     prop.D[v][id3][id2][ic3][ic2])));
		    // x-mu  x-mu  x
		    bdfmu[mu][id1][id3] = qcd_CADD(bdfmu[mu][id1][id3],
						   qcd_CMUL(backprop.D[vm1][id1][id2][ic1][ic2],
							    qcd_CMUL(u.D[vm1][mu][ic1][ic3],
								     prop.D[v][id3][id2][ic3][ic2])));

                  }

	    if(qcd_NORM(qcd_ONE_PLUS_GAMMA[mu][id1][id3])>1e-4)
	      for(int id2=0; id2<4; id2++)
		for(int ic1=0; ic1<3; ic1++)
                  for(int ic2=0; ic2<3; ic2++)
		    for(int ic3=0; ic3<3; ic3++) {
                      // x+mu  x  x
                      block_n[mu][iv] = qcd_CADD(block_n[mu][iv],
                                                qcd_CMUL(backprop.D[vp1][id1][id2][ic1][ic2],
                                                         qcd_CMUL(qcd_ONE_PLUS_GAMMA[mu][id1][id3],
                                                                  qcd_CMUL(qcd_CONJ(u.D[v][mu][ic3][ic1]),
                                                                           prop.D[v][id3][id2][ic3][ic2]))));
                      // x  x-mu  x-mu
                      block_n[mu][iv] = qcd_CADD(block_n[mu][iv],
                                                qcd_CMUL(backprop.D[v][id1][id2][ic1][ic2],
                                                         qcd_CMUL(qcd_ONE_PLUS_GAMMA[mu][id1][id3],
                                                                  qcd_CMUL(qcd_CONJ(u.D[vm1][mu][ic3][ic1]),
                                                                           prop.D[vm1][id3][id2][ic3][ic2]))));
		    }
	  
	    if(qcd_NORM(qcd_ONE_MINUS_GAMMA[mu][id1][id3])>1e-4)
	      for(int id2=0; id2<4; id2++)
		for(int ic1=0; ic1<3; ic1++)
                  for(int ic2=0; ic2<3; ic2++)
		    for(int ic3=0; ic3<3; ic3++) {
                      // x  x  x+mu
                      block_n[mu][iv] = qcd_CSUB(block_n[mu][iv],
                                                qcd_CMUL(backprop.D[v][id1][id2][ic1][ic2],
                                                         qcd_CMUL(qcd_ONE_MINUS_GAMMA[mu][id1][id3],
                                                                  qcd_CMUL(u.D[v][mu][ic1][ic3],
                                                                           prop.D[vp1][id3][id2][ic3][ic2]))));
		      
                      // x-mu  x-mu  x
                      block_n[mu][iv] = qcd_CSUB(block_n[mu][iv],
                                                qcd_CMUL(backprop.D[vm1][id1][id2][ic1][ic2],
                                                         qcd_CMUL(qcd_ONE_MINUS_GAMMA[mu][id1][id3],
                                                                  qcd_CMUL(u.D[vm1][mu][ic1][ic3],
                                                                           prop.D[v][id3][id2][ic3][ic2]))));
		    }         
	  }//end mu loop

	  //now local operators
	  //pre-calculate a partial trace of Backward x Forward prop:
	  backfor = (qcd_complex_16) {0,0};
	  for(int id2=0; id2<4; id2++)
	    for(int ic1=0; ic1<3; ic1++)
	      for(int ic2=0; ic2<3; ic2++) {
		backfor = qcd_CADD(backfor, qcd_CMUL(backprop.D[v][id1][id2][ic1][ic2],
						     prop.D[v][id3][id2][ic1][ic2]));
	      }
	  
	  
	  /*********** pseudo-scalar current ***********/
	  if(qcd_NORM(qcd_GAMMA[5][id1][id3])>1e-4) {
	    block_p[0][iv] = qcd_CADD(block_p[0][iv], qcd_CMUL(backfor,
							       qcd_GAMMA[5][id1][id3]));
	  }

	  /*********** scalar current ***********/
	  if(id1 == id3) {
	    block_s[0][iv] = qcd_CADD(block_s[0][iv], backfor);
	  }
	  
	  /*********** local vector current ***********/
	  for(int mu=0; mu<4; mu++)
	    if(qcd_NORM(qcd_GAMMA[mu][id1][id3])>1e-4) {
	      block_l[mu][iv] = qcd_CADD(block_l[mu][iv], qcd_CMUL(backfor,
								   qcd_GAMMA[mu][id1][id3]));                             
	    }

	  /*********** local axial current ***********/     
	  for(int mu=0; mu<4; mu++)
	    if(qcd_NORM(qcd_G5GAMMA[mu][id1][id3])>1e-4) {
	      block_a[mu][iv] = qcd_CADD(block_a[mu][iv], qcd_CMUL(backfor,                                                      
								   qcd_G5GAMMA[mu][id1][id3]));
	    }
	}//end id1,id3 loop

      for(int id1=0; id1<4; id1++)
	for(int id3=0; id3<4; id3++) {
	  for(int mu=0; mu<4; mu++)
	    for(int nu=0; nu<=mu; nu++) {                       
	      /*********** one derivative vector operator ***********/
	      if(qcd_NORM(qcd_GAMMA[mu][id1][id3])>1e-4) {
		block_vD[mu*4+nu][iv] = qcd_CADD(block_vD[mu*4+nu][iv],
						 qcd_CMUL(qcd_GAMMA[mu][id1][id3],
							  bdfmu[nu][id1][id3]));                                          
	      }
	      if(qcd_NORM(qcd_GAMMA[nu][id1][id3])>1e-4) {  
		block_vD[mu*4+nu][iv] = qcd_CADD(block_vD[mu*4+nu][iv],
						 qcd_CMUL(qcd_GAMMA[nu][id1][id3],
							  bdfmu[mu][id1][id3]));
	      }

	      /*********** one derivative axial and axial-antisymmetric operators ***********/
	      if(qcd_NORM(qcd_G5GAMMA[mu][id1][id3])>1e-4) {
		block_aD[mu*4+nu][iv] = qcd_CADD(block_aD[mu*4+nu][iv],
						qcd_CMUL(qcd_G5GAMMA[mu][id1][id3],
							 bdfmu[nu][id1][id3]));
		
		block_d1[mu*4+nu][iv] = qcd_CADD(block_d1[mu*4+nu][iv],
						qcd_CMUL(qcd_G5GAMMA[mu][id1][id3],
							 bdfmu[nu][id1][id3]));
	      }
	      
	      if(qcd_NORM(qcd_G5GAMMA[nu][id1][id3])>1e-4) {
		block_aD[mu*4+nu][iv] = qcd_CADD(block_aD[mu*4+nu][iv],
						qcd_CMUL(qcd_G5GAMMA[nu][id1][id3],
							 bdfmu[mu][id1][id3]));
		
		block_d1[mu*4+nu][iv] = qcd_CSUB(block_d1[mu*4+nu][iv],
						qcd_CMUL(qcd_G5GAMMA[nu][id1][id3],
							 bdfmu[mu][id1][id3]));
	      }
	    }//end mu nu loop
	}//end id1 id3 loop   
    }//end volume loop    
    
#pragma omp parallel for
    for(int v3=0; v3<geo.lV3; v3++) {         
      int iv = v3 + it*geo.lV3;
      block_s[0][iv].re *= sign;
      block_s[0][iv].im *= sign;

      block_p[0][iv].re *= sign;
      block_p[0][iv].im *= sign;

      for(int mu=0; mu<4; mu++) {
	block_l[mu][iv].re *= sign;
	block_l[mu][iv].im *= sign;

	block_n[mu][iv].re *= sign;
	block_n[mu][iv].im *= sign;

	block_a[mu][iv].re *= sign;
	block_a[mu][iv].im *= sign;
      }

      for(int mu=0; mu<16; mu++) {
	block_vD[mu][iv].re *= sign;
	block_vD[mu][iv].im *= sign;

	block_aD[mu][iv].re *= sign;
	block_aD[mu][iv].im *= sign;

	block_d1[mu][iv].re *= sign;
	block_d1[mu][iv].im *= sign;	
      }
    }


  }// end t-insertion loop

  /*
   * Scalar
   */  
  {
    int site_size = 1;
    char *site_tags[] = {"(scalar: 0)"};
    int t_ins[nt_ins];
    for(int t=t_start; t<=t_stop; t++){
      t_ins[t-t_start] = t;
    }  
    writeBlock(blockloc_s_name, block_s, x_src, nt_ins, t_ins, site_size, site_tags, geo);
  }

  /*
   * Pseudo-scalar
   */  
  {
    int site_size = 1;
    char *site_tags[] = {"(pseudoscalar: 0)"};
    int t_ins[nt_ins];
    for(int t=t_start; t<=t_stop; t++){
      t_ins[t-t_start] = t;
    }  
    writeBlock(blockloc_p_name, block_p, x_src, nt_ins, t_ins, site_size, site_tags, geo);
  }

  /*
   * Axial-vector
   */  
  {
    int site_size = 4;
    char *site_tags[] = {"(axial-vector-t: 0)",
			 "(axial-vector-x: 1)",
			 "(axial-vector-y: 2)",
			 "(axial-vector-z: 3)"};
    int t_ins[nt_ins];
    for(int t=t_start; t<=t_stop; t++){
      t_ins[t-t_start] = t;
    }  
    writeBlock(blockloc_a_name, block_a, x_src, nt_ins, t_ins, site_size, site_tags, geo);
  }

  /*
   * Vector-local
   */  
  {
    int site_size = 4;
    char *site_tags[] = {"(vector-local-t: 0)",
			 "(vector-local-x: 1)",
			 "(vector-local-y: 2)",
			 "(vector-local-z: 3)"};
    int t_ins[nt_ins];
    for(int t=t_start; t<=t_stop; t++){
      t_ins[t-t_start] = t;
    }  
    writeBlock(blockloc_v_name, block_l, x_src, nt_ins, t_ins, site_size, site_tags, geo);
  }

  /*
   * Vector-noether
   */  
  {
    int site_size = 4;
    char *site_tags[] = {"(vector-noether-t: 0)",
			 "(vector-noether-x: 1)",
			 "(vector-noether-y: 2)",
			 "(vector-noether-z: 3)"};
    int t_ins[nt_ins];
    for(int t=t_start; t<=t_stop; t++){
      t_ins[t-t_start] = t;
    }  
    for(int mu=0; mu<4; mu++)
      for(int v=0; v<geo.lV3*nt_ins; v++) {
	block_n[mu][v].re /= 4.0;
	block_n[mu][v].im /= 4.0;
      }                                    
    writeBlock(blocknoe_v_name, block_n, x_src, nt_ins, t_ins, site_size, site_tags, geo);
  }

  /*
   * Vector derivative
   */
  {
    int site_size = 4*(4+1)/2;
    char *site_tags[] = {"(vector-deriv-tt: 0 0)",
			 "(vector-deriv-xt: 1 0)",
			 "(vector-deriv-xx: 1 1)",
			 "(vector-deriv-yt: 2 0)",
			 "(vector-deriv-yx: 2 1)",
			 "(vector-deriv-yy: 2 2)",
			 "(vector-deriv-zt: 3 0)",
			 "(vector-deriv-zx: 3 1)",
			 "(vector-deriv-zy: 3 2)",
			 "(vector-deriv-zz: 3 3)"};
    int t_ins[nt_ins];
    for(int t=t_start; t<=t_stop; t++){
      t_ins[t-t_start] = t;
    }  
    int x = 0;
    qcd_complex_16 *block_vD_swap[site_size];
    for(int mu=0; mu<4; mu++)
      for(int nu=0; nu<=mu; nu++) {
	block_vD_swap[x] = block_vD[nu + mu*4];
	for(int v=0; v<geo.lV3*nt_ins; v++) {
	  block_vD_swap[x][v].re /= 8.0;
	  block_vD_swap[x][v].im /= 8.0;
	}
	x++;
      }
    writeBlock(block_vD_name, block_vD_swap, x_src, nt_ins, t_ins, site_size, site_tags, geo);
  }

  /*
   * Axial-vector derivative
   */
  {
    int site_size = 4*(4+1)/2;
    char *site_tags[] = {"(axial-vector-deriv-tt: 0 0)",
			 "(axial-vector-deriv-xt: 1 0)",
			 "(axial-vector-deriv-xx: 1 1)",
			 "(axial-vector-deriv-yt: 2 0)",
			 "(axial-vector-deriv-yx: 2 1)",
			 "(axial-vector-deriv-yy: 2 2)",
			 "(axial-vector-deriv-zt: 3 0)",
			 "(axial-vector-deriv-zx: 3 1)",
			 "(axial-vector-deriv-zy: 3 2)",
			 "(axial-vector-deriv-zz: 3 3)"};
    int t_ins[nt_ins];
    for(int t=t_start; t<=t_stop; t++){
      t_ins[t-t_start] = t;
    }  
    int x = 0;
    qcd_complex_16 *block_aD_swap[site_size];
    for(int mu=0; mu<4; mu++)
      for(int nu=0; nu<=mu; nu++) {
	block_aD_swap[x] = block_aD[nu + mu*4];
	for(int v=0; v<geo.lV3*nt_ins; v++) {
	  block_aD_swap[x][v].re /= 8.0;
	  block_aD_swap[x][v].im /= 8.0;
	}
	x++;
      }
    writeBlock(block_aD_name, block_aD_swap, x_src, nt_ins, t_ins, site_size, site_tags, geo);
  }

  /*
   * Axial-vector derivative antisymmetric
   */
  {
    int site_size = (4-1)*4/2;
    char *site_tags[] = {"(axial-vector-deriv-antisymmetric-xt: 1 0)",
			 "(axial-vector-deriv-antisymmetric-yt: 2 0)",
			 "(axial-vector-deriv-antisymmetric-yx: 2 1)",
			 "(axial-vector-deriv-antisymmetric-zt: 3 0)",
			 "(axial-vector-deriv-antisymmetric-zx: 3 1)",
			 "(axial-vector-deriv-antisymmetric-zy: 3 2)"};
    
    int t_ins[nt_ins];
    for(int t=t_start; t<=t_stop; t++){
      t_ins[t-t_start] = t;
    }  
    int x = 0;
    qcd_complex_16 *block_d1_swap[site_size];
    for(int mu=0; mu<4; mu++)
      for(int nu=0; nu<mu; nu++) {
	block_d1_swap[x] = block_d1[nu + mu*4];
	for(int v=0; v<geo.lV3*nt_ins; v++) {
	  block_d1_swap[x][v].re /= 8.0;
	  block_d1_swap[x][v].im /= 8.0;
	}
	x++;
      }
    writeBlock(block_d1_name, block_d1_swap, x_src, nt_ins, t_ins, site_size, site_tags, geo);
  }
  
  for(int i=0; i<4; i++) {
    free(block_n[i]);
    free(block_l[i]);
    free(block_a[i]);
  }
  for(int i=0; i<16; i++) {
    free(block_vD[i]);
    free(block_aD[i]);
    free(block_d1[i]);
  }
  free(block_s[0]);
  free(block_s);
  free(block_p[0]);
  free(block_p);

  free(block_n);
  free(block_l);
  free(block_a);

  free(block_vD);
  free(block_aD);
  free(block_d1);

  qcd_destroyPropagator(&prop);
  qcd_destroyPropagator(&backprop);
  qcd_destroyGaugeField(&u);
  qcd_destroyGeometry(&geo);
  MPI_Finalize();
  return 0;
}
