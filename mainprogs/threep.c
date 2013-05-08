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





int main(int argc,char* argv[])
{
   qcd_uint_2 mu,nu,ku,lu,c1,c2,id1,id2,id3;  // various loop variables
   qcd_uint_4 i,k,v,lx,ly,lz,ip1,im1; 
   qcd_int_4 x,y,z;
   qcd_uint_4 numOfMom;                       // number of momenta
   qcd_uint_2 ic1,ic2,ic3;                    //
   qcd_uint_4 x_src[4];                       // source and sink coordinates
   qcd_uint_4 t_sink, t_start, t_stop, t;
   int ierr;
   int lt;
   qcd_real_8 tmp;                            // general purpuse
   FILE *fp_momlist;
  
   FILE *fp_corrnoe_v;                      // output files
   FILE *fp_corrloc_s;      
   FILE *fp_corrloc_p;      
   FILE *fp_corrloc_v;      
   FILE *fp_corrloc_a;      
   FILE *fp_corrloc_t;   
   FILE *fp_corr_vD;
   FILE *fp_corr_aD;
//   FILE *fp_corr_tD;
   FILE *fp_corr_d1;
  
   int params_len;               // needed to read inputfiles
   char *params;                 // needed to read inputfiles

   char gauge_name[qcd_MAX_STRING_LENGTH];      // name of gauge-configuration file
   char corrloc_s_name[qcd_MAX_STRING_LENGTH];  // name of output file name scalar current
   char corrloc_p_name[qcd_MAX_STRING_LENGTH];  // name of output file name pseudo-scalar current
   char corrloc_v_name[qcd_MAX_STRING_LENGTH];  // name of output file name local vector current
   char corrnoe_v_name[qcd_MAX_STRING_LENGTH];  // name of output file name noether vector current
   char corrloc_a_name[qcd_MAX_STRING_LENGTH];  // name of output file name local axial current
   char corrloc_t_name[qcd_MAX_STRING_LENGTH];// name of output file name local tensor current
   char corr_vD_name[qcd_MAX_STRING_LENGTH];    // name of output file name one derivative vector
   char corr_aD_name[qcd_MAX_STRING_LENGTH];    // name of output file name one derivative axial
//   char corr_tD_name[qcd_MAX_STRING_LENGTH];  // name of output file name one derivative tensor
   char corr_d1_name[qcd_MAX_STRING_LENGTH];    // name of output file name one derivative axial antisymmetric
   
   char param_name[qcd_MAX_STRING_LENGTH];      // name of parameter file  
   char momlist_name[qcd_MAX_STRING_LENGTH];    // name of momenta-list file
   char prop_name[qcd_MAX_STRING_LENGTH];       // file names of up and down quark propagators
   char bprop_name[qcd_MAX_STRING_LENGTH]; 

   qcd_geometry geo;                            // geometry structure
   qcd_propagator prop;                         // propagator
   qcd_propagator backprop;                     // backward prop. 
   qcd_gaugeField u;                            // gauge field 

   qcd_real_8 theta[4] = {M_PI,0.0,0.0,0.0};    // antiperiodic b.c. in time
   qcd_uint_2 L[4];
   qcd_uint_2 P[4];
   qcd_complex_16 phase_factor, phase_factor_b;         
   qcd_complex_16 z1, z2;                       // temp variables
  
   qcd_complex_16 *block_p;			// to store the blocks (pseudo-scalar)
   qcd_complex_16 *block_s;			// to store the blocks (scalar)
   qcd_complex_16 **block_n, **block_l;         // to store the blocks (local & noether vector)
   qcd_complex_16 **block_a,  **block_t;        // (local axial & local tensor)
   qcd_complex_16 **block_vD, **block_aD;       // (one derivative vector & axial)
   qcd_complex_16 /* **block_tD ,*/ **block_d1; // (one derivative tensor & antisym. vector)
   
   qcd_complex_16 g5sigmu0[5][4][4];            // gamma_5 * [gamma_mu, gamma_0] *1/2
   qcd_complex_16 g5sigmu1[5][4][4];            // gamma_5 * [gamma_mu, gamma_1] *1/2
   qcd_complex_16 g5sigmu2[5][4][4];            // gamma_5 * [gamma_mu, gamma_2] *1/2
   qcd_complex_16 g5sigmu3[5][4][4];            // gamma_5 * [gamma_mu, gamma_3] *1/2
   
   qcd_int_4 (*mom)[3];                         // momenta-list

   int myid,numprocs, namelen;    
   char processor_name[MPI_MAX_PROCESSOR_NAME];
   				 
				 
             
             
             
   //set up MPI
   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD,&numprocs);         // num. of processes taking part in the calculation
   MPI_Comm_rank(MPI_COMM_WORLD,&myid);             // each process gets its ID
   MPI_Get_processor_name(processor_name,&namelen); // 
   
   
   for(i=0; i<5; i++)
   for(mu=0; mu<4; mu++)
   for(nu=0; nu<4; nu++)
   {           
      g5sigmu0[i][mu][nu]= (qcd_complex_16){0,0};
      g5sigmu1[i][mu][nu]= (qcd_complex_16){0,0};
      g5sigmu2[i][mu][nu]= (qcd_complex_16){0,0};
      g5sigmu3[i][mu][nu]= (qcd_complex_16){0,0};
      for(ku=0; ku<4; ku++)
      for(lu=0; lu<4; lu++)
      {
	/* 0 */
	g5sigmu0[i][mu][nu] = qcd_CADD(g5sigmu0[i][mu][nu],qcd_CMUL(qcd_CMUL(qcd_GAMMA[5][mu][ku],
									     qcd_GAMMA[i][ku][lu]),
								    qcd_GAMMA[0][lu][nu]));
	g5sigmu0[i][mu][nu] = qcd_CSUB(g5sigmu0[i][mu][nu],qcd_CMUL(qcd_CMUL(qcd_GAMMA[5][mu][ku],
									     qcd_GAMMA[0][ku][lu]),
								    qcd_GAMMA[i][lu][nu]));
	/* 1 */
	g5sigmu1[i][mu][nu] = qcd_CADD(g5sigmu1[i][mu][nu],qcd_CMUL(qcd_CMUL(qcd_GAMMA[5][mu][ku],
									     qcd_GAMMA[i][ku][lu]),
								    qcd_GAMMA[1][lu][nu]));
	g5sigmu1[i][mu][nu] = qcd_CSUB(g5sigmu1[i][mu][nu],qcd_CMUL(qcd_CMUL(qcd_GAMMA[5][mu][ku],
									     qcd_GAMMA[1][ku][lu]),
								    qcd_GAMMA[i][lu][nu]));
	/* 2 */
	g5sigmu2[i][mu][nu] = qcd_CADD(g5sigmu2[i][mu][nu],qcd_CMUL(qcd_CMUL(qcd_GAMMA[5][mu][ku],
									     qcd_GAMMA[i][ku][lu]),
								    qcd_GAMMA[2][lu][nu]));
	g5sigmu2[i][mu][nu] = qcd_CSUB(g5sigmu2[i][mu][nu],qcd_CMUL(qcd_CMUL(qcd_GAMMA[5][mu][ku],
									     qcd_GAMMA[2][ku][lu]),
								    qcd_GAMMA[i][lu][nu]));
	/* 3 */
	g5sigmu3[i][mu][nu] = qcd_CADD(g5sigmu3[i][mu][nu],qcd_CMUL(qcd_CMUL(qcd_GAMMA[5][mu][ku],
									     qcd_GAMMA[i][ku][lu]),
								    qcd_GAMMA[3][lu][nu]));
	g5sigmu3[i][mu][nu] = qcd_CSUB(g5sigmu3[i][mu][nu],qcd_CMUL(qcd_CMUL(qcd_GAMMA[5][mu][ku],
									     qcd_GAMMA[3][ku][lu]),
								    qcd_GAMMA[i][lu][nu]));
      }
      g5sigmu0[i][mu][nu] = qcd_CSCALE(g5sigmu0[i][mu][nu],0.5);
      g5sigmu1[i][mu][nu] = qcd_CSCALE(g5sigmu1[i][mu][nu],0.5);
      g5sigmu2[i][mu][nu] = qcd_CSCALE(g5sigmu2[i][mu][nu],0.5);
      g5sigmu3[i][mu][nu] = qcd_CSCALE(g5sigmu3[i][mu][nu],0.5);
   }
   
   
   
   
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
      
   sscanf(qcd_getParam("<t>",params,params_len),"%d %d",&t_start, &t_stop);
   if(myid==0) printf("Got insertion time slices: %d ... %d\n",t_start,t_stop);
                     	  
   sscanf(qcd_getParam("<source_pos_txyz>",params,params_len),"%d %d %d %d",&x_src[0],&x_src[1],&x_src[2],&x_src[3]);
   if(myid==0) printf("Got source coords: %d %d %d %d\n",x_src[0],x_src[1],x_src[2],x_src[3]);
   
   sscanf(qcd_getParam("<t_sink>",params,params_len),"%d",&t_sink);
   if(myid==0) printf("Got sink time slice: %d\n",t_sink);
  
   if(t_sink >= L[0])
     {
       if(myid==0) fprintf(stderr, " Error: t_sink (=%d) >= L[0] (=%d),\n t_sink should be in [0, L[0]), did you forget to mod(t_sink, L[0]) ?\n", t_sink, L[0]);
       exit(EXIT_FAILURE);
     }

   strcpy(prop_name,qcd_getParam("<propagator>",params,params_len));
   if(myid==0) printf("Got propagator file name: %s\n",prop_name);
   strcpy(bprop_name,qcd_getParam("<seq_prop>",params,params_len));
   if(myid==0) printf("Got sequential propagator file name: %s\n",bprop_name);

   strcpy(gauge_name,qcd_getParam("<cfg_name>",params,params_len));
   if(myid==0) printf("Got conf name: %s\n",gauge_name);
          
   strcpy(corrloc_p_name,qcd_getParam("<corr_name_p_local>",params,params_len));
   if(myid==0) printf("Got output file name: %s\n",corrloc_p_name);
   
   strcpy(corrloc_s_name,qcd_getParam("<corr_name_s_local>",params,params_len));
   if(myid==0) printf("Got output file name: %s\n",corrloc_s_name);
   
   strcpy(corrloc_v_name,qcd_getParam("<corr_name_v_local>",params,params_len));
   if(myid==0) printf("Got output file name: %s\n",corrloc_v_name);
   
   strcpy(corrnoe_v_name,qcd_getParam("<corr_name_v_noether>",params,params_len));
   if(myid==0) printf("Got output file name: %s\n",corrnoe_v_name);
   
   strcpy(corrloc_a_name,qcd_getParam("<corr_name_a_local>",params,params_len));
   if(myid==0) printf("Got output file name: %s\n",corrloc_a_name);
   
  strcpy(corrloc_t_name,qcd_getParam("<corr_name_t_local>",params,params_len));
   if(myid==0) printf("Got output file name: %s\n",corrloc_t_name);
   
   strcpy(corr_vD_name,qcd_getParam("<corr_name_vD>",params,params_len));
   if(myid==0) printf("Got output file name: %s\n",corr_vD_name);
      
   strcpy(corr_aD_name,qcd_getParam("<corr_name_aD>",params,params_len));
   if(myid==0) printf("Got output file name: %s\n",corr_aD_name);
   
/*   strcpy(corr_tD_name,qcd_getParam("<corr_name_tD>",params,params_len));
   if(myid==0) printf("Got output file name: %s\n",corr_tD_name); */
      
   strcpy(corr_d1_name,qcd_getParam("<corr_name_d1>",params,params_len));
   if(myid==0) printf("Got output file name: %s\n",corr_d1_name);
      
   strcpy(momlist_name,qcd_getParam("<momenta_list>",params,params_len));
   if(myid==0) printf("Got momenta-list file name: %s\n",momlist_name);
    
   free(params);



         
   //#####################################################################   
   // allocate memory
  
   ierr = 0;
   ierr += qcd_initPropagator(&prop, &geo);
   ierr += qcd_initPropagator(&backprop, &geo);
   ierr += qcd_initGaugeField(&u, &geo);
   
   block_n = (qcd_complex_16**)malloc(4*sizeof(*block_n));
   block_l = (qcd_complex_16**)malloc(4*sizeof(*block_l));
   block_a = (qcd_complex_16**)malloc(4*sizeof(*block_a));
   block_t = (qcd_complex_16**)malloc(16*sizeof(*block_t));
   block_vD = (qcd_complex_16**)malloc(16*sizeof(*block_vD));
   block_aD = (qcd_complex_16**)malloc(16*sizeof(*block_aD));
//   block_tD = (qcd_complex_16**)malloc(16*sizeof(*block_tD));
   block_d1 = (qcd_complex_16**)malloc(16*sizeof(*block_d1));

   block_s = (qcd_complex_16 *)malloc(geo.lV3*sizeof(qcd_complex_16));
   block_p = (qcd_complex_16 *)malloc(geo.lV3*sizeof(qcd_complex_16));
   for(i=0; i<4; i++)
   {
      block_n[i] = (qcd_complex_16*) malloc(geo.lV3*sizeof(qcd_complex_16));
      if(block_n[i]==NULL)
      {
         fprintf(stderr,"process %i: out of memmory (for noether isovector block)\n",myid);
         ierr++;
      }
      block_l[i] = (qcd_complex_16*) malloc(geo.lV3*sizeof(qcd_complex_16));
      if(block_l[i]==NULL)
      {
         fprintf(stderr,"process %i: out of memmory (for local isovector block)\n",myid);
         ierr++;
      }
      block_a[i] = (qcd_complex_16*) malloc(geo.lV3*sizeof(qcd_complex_16));
      if(block_a[i]==NULL)
      {
         fprintf(stderr,"process %i: out of memmory (for local axial block)\n",myid);
         ierr++;
      }        
   }
   for(i=0; i<16; i++)
   {
      block_vD[i] = (qcd_complex_16*) malloc(geo.lV3*sizeof(qcd_complex_16));
      if(block_vD[i]==NULL)
      {
         fprintf(stderr,"process %i: out of memmory (for vD block)\n",myid);
         ierr++;
      }
      block_aD[i] = (qcd_complex_16*) malloc(geo.lV3*sizeof(qcd_complex_16));
      if(block_aD[i]==NULL)
      {
         fprintf(stderr,"process %i: out of memmory (for aD block)\n",myid);
         ierr++;
      }

      block_t[i] = (qcd_complex_16*) malloc(geo.lV3*sizeof(qcd_complex_16));
      if(block_t[i]==NULL)
      {
         fprintf(stderr,"process %i: out of memmory (for local tensor block)\n",myid);
         ierr++;
      }

/*      
      block_tD[i] = (qcd_complex_16*) malloc(geo.lV3*sizeof(qcd_complex_16));
      if(block_tD[i]==NULL)
      {
         fprintf(stderr,"process %i: out of memmory (for tD block)\n",myid);
         ierr++;
      }   
*/
      block_d1[i] = (qcd_complex_16*) malloc(geo.lV3*sizeof(qcd_complex_16));
      if(block_d1[i]==NULL)
      {
         fprintf(stderr,"process %i: out of memmory (for d1 block)\n",myid);
         ierr++;
      }         
   }
         
   MPI_Allreduce(&ierr, &k, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   if(k>0)
   {
      if(myid==0) printf("not enough memory\n");
      exit(EXIT_FAILURE);
   }
   if(myid==0) printf("memory for propagators and gauge-field allocated\n");
         
   
   //##############################################################################
   // load gauge-field
   if(qcd_getGaugeField(gauge_name,qcd_GF_LIME,&u)) exit(EXIT_FAILURE);
   if(myid==0) printf("gauge-field loaded\n");   
   
   qcd_communicateGaugePM(&u);
   
   // load propagators
   
   if(qcd_getPropagator(prop_name,qcd_PROP_LIME, &prop)) exit(EXIT_FAILURE);
   if(myid==0) printf("propagator loaded\n");
   if(qcd_getPropagator(bprop_name,qcd_PROP_LIME, &backprop)) exit(EXIT_FAILURE);
   if(myid==0) printf("sequential propagator loaded\n");   
   
   
   //################################################################################
   // transform propagators to basis with theta-periodic boundaries in the temporal direction
   for(lt=0; lt<geo.lL[0]; lt++)
   {
      t = lt + geo.Pos[0] * geo.lL[0];
      phase_factor   = (qcd_complex_16) {cos(theta[0]*t/geo.L[0]),sin(theta[0]*t/geo.L[0])};
      // the backward-props get the complex-conjugated phases, so after the conjugation
      // in the next step they will be correct
      phase_factor_b = (qcd_complex_16) {cos(theta[0]*(t_sink-(t+x_src[0])+2*(t_sink-x_src[0]))/geo.L[0]),-sin(theta[0]*(t_sink-(t+x_src[0])+2*(t_sink-x_src[0]))/geo.L[0])};
      qcd_mulPropagatorC3d(&prop, phase_factor, (t+x_src[0]) % geo.L[0]);
      qcd_mulPropagatorC3d(&backprop, phase_factor_b, (t+x_src[0]) % geo.L[0]);
   }
   if(myid==0) printf("propagators transformed to basis with theta-periodic boundary conditions\n");
   
   qcd_waitall(&geo);
   qcd_communicatePropagatorPM(&prop);
         
   //#####################################################################   
   // complex-conjugate the backward propagators, multiply by gamma_5. Works only in chiral gamma basis
   qcd_conjPropagator(&backprop);
   qcd_gamma5Propagator(&backprop);   
   if(myid==0) printf("backward propagators complex-conjugated and multiplied by gamma_5\n");

   qcd_waitall(&geo);
   qcd_communicatePropagatorPM(&backprop);
   qcd_waitall(&geo);
   if(myid==0) printf("communication done.\n");
   //#####################################################################   
   // calculate the 'blocks'
   
   
   
   qcd_complex_16 *corr[4], *corr_master[4];
   for(t=t_start; t<=t_stop; t++)
   {
      if(myid==0) printf("t=%i\n",t);
      
      for(mu=0; mu<4; mu++)
      for(i=0; i<geo.lV3; i++)   //set blocks to zero
      {
         block_n[mu][i]= (qcd_complex_16) {0,0};
         block_l[mu][i]= (qcd_complex_16) {0,0};
         block_a[mu][i]= (qcd_complex_16) {0,0};
	 block_s[i] = (qcd_complex_16) {0,0};
	 block_p[i] = (qcd_complex_16) {0,0};
      }
      for(mu=0; mu<16; mu++)
      for(i=0; i<geo.lV3; i++)   //set blocks to zero
      {
         block_vD[mu][i]= (qcd_complex_16) {0,0};
         block_aD[mu][i]= (qcd_complex_16) {0,0};
//         block_tD[mu][i]= (qcd_complex_16) {0,0};
         block_d1[mu][i]= (qcd_complex_16) {0,0};
         block_t[mu][i]= (qcd_complex_16) {0,0};
      }

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

      lt = ((t+x_src[0])%geo.L[0]) - geo.Pos[0]*geo.lL[0];
      if(lt>=0 && lt<geo.lL[0])  //inside the local lattice, otherwise nothing to calculate
      {
	int j;
#pragma omp parallel for private(id1,id2,id3,ip1,im1,ic1,ic2,ic3,i,mu,nu)
         for(j=0; j<geo.lV3; j++)
	   {   
	     qcd_complex_16 backfor;                      // backward-prop x forward-prop partially traced
	     qcd_complex_16 bdfmu[4][4][4];               // stores backward-prop D_mu forward-prop
	     
            i=lt + j*geo.lL[0];         
            for(id1=0; id1<4; id1++)
            for(id3=0; id3<4; id3++)
            { 
               for(mu=0; mu<4; mu++) 
               {
                  ip1 = geo.plus[i][mu];              
                  im1 = geo.minus[i][mu];

                  bdfmu[mu][id1][id3] = (qcd_complex_16) {0,0};
                  for(id2=0; id2<4; id2++)
                  for(ic1=0; ic1<3; ic1++)
                  for(ic2=0; ic2<3; ic2++)
                  for(ic3=0; ic3<3; ic3++)
                  {
                     // x x x+mu
                     bdfmu[mu][id1][id3] = qcd_CADD(bdfmu[mu][id1][id3],
                                                    qcd_CMUL(backprop.D[i][id1][id2][ic1][ic2],
                                                             qcd_CMUL(u.D[i][mu][ic1][ic3],
                                                                      prop.D[ip1][id3][id2][ic3][ic2])));
                     // x  x-mu  x-mu
                     bdfmu[mu][id1][id3] = qcd_CSUB(bdfmu[mu][id1][id3],
                                                    qcd_CMUL(backprop.D[i][id1][id2][ic1][ic2],
                                                             qcd_CMUL(qcd_CONJ(u.D[im1][mu][ic3][ic1]),
                                                                      prop.D[im1][id3][id2][ic3][ic2])));
                     // x+mu  x  x
                     bdfmu[mu][id1][id3] = qcd_CSUB(bdfmu[mu][id1][id3],
                                                    qcd_CMUL(backprop.D[ip1][id1][id2][ic1][ic2],
                                                             qcd_CMUL(qcd_CONJ(u.D[i][mu][ic3][ic1]),
                                                                      prop.D[i][id3][id2][ic3][ic2])));
                     // x-mu  x-mu  x
                     bdfmu[mu][id1][id3] = qcd_CADD(bdfmu[mu][id1][id3],
                                                    qcd_CMUL(backprop.D[im1][id1][id2][ic1][ic2],
                                                             qcd_CMUL(u.D[im1][mu][ic1][ic3],
                                                                      prop.D[i][id3][id2][ic3][ic2])));
                  }

                  




                  if(qcd_NORM(qcd_ONE_PLUS_GAMMA[mu][id1][id3])>1e-4)
                  for(id2=0; id2<4; id2++)
                  for(ic1=0; ic1<3; ic1++)
                  for(ic2=0; ic2<3; ic2++)
                  for(ic3=0; ic3<3; ic3++)
                  {
                      // x+mu  x  x
                      block_n[mu][j] = qcd_CADD(block_n[mu][j],
                                                qcd_CMUL(backprop.D[ip1][id1][id2][ic1][ic2],
                                                         qcd_CMUL(qcd_ONE_PLUS_GAMMA[mu][id1][id3],
                                                                  qcd_CMUL(qcd_CONJ(u.D[i][mu][ic3][ic1]),
                                                                           prop.D[i][id3][id2][ic3][ic2]))));
                      // x  x-mu  x-mu
                      block_n[mu][j] = qcd_CADD(block_n[mu][j],
                                                qcd_CMUL(backprop.D[i][id1][id2][ic1][ic2],
                                                         qcd_CMUL(qcd_ONE_PLUS_GAMMA[mu][id1][id3],
                                                                  qcd_CMUL(qcd_CONJ(u.D[im1][mu][ic3][ic1]),
                                                                           prop.D[im1][id3][id2][ic3][ic2]))));
                  }
                  if(qcd_NORM(qcd_ONE_MINUS_GAMMA[mu][id1][id3])>1e-4)
                  for(id2=0; id2<4; id2++)
                  for(ic1=0; ic1<3; ic1++)
                  for(ic2=0; ic2<3; ic2++)
                  for(ic3=0; ic3<3; ic3++)
                  {
                      // x  x  x+mu
                      block_n[mu][j] = qcd_CSUB(block_n[mu][j],
                                                qcd_CMUL(backprop.D[i][id1][id2][ic1][ic2],
                                                         qcd_CMUL(qcd_ONE_MINUS_GAMMA[mu][id1][id3],
                                                                  qcd_CMUL(u.D[i][mu][ic1][ic3],
                                                                           prop.D[ip1][id3][id2][ic3][ic2]))));

                      // x-mu  x-mu  x
                      block_n[mu][j] = qcd_CSUB(block_n[mu][j],
                                                qcd_CMUL(backprop.D[im1][id1][id2][ic1][ic2],
                                                         qcd_CMUL(qcd_ONE_MINUS_GAMMA[mu][id1][id3],
                                                                  qcd_CMUL(u.D[im1][mu][ic1][ic3],
                                                                           prop.D[i][id3][id2][ic3][ic2]))));
                  }         
               }//end mu loop

               //now local operators
               //pre-calculate a partial trace of Backward x Forward prop:
               backfor=(qcd_complex_16) {0,0};
               for(id2=0; id2<4; id2++)
               for(ic1=0; ic1<3; ic1++)
               for(ic2=0; ic2<3; ic2++)
               {
                  backfor = qcd_CADD(backfor, qcd_CMUL(backprop.D[i][id1][id2][ic1][ic2],
                                                       prop.D[i][id3][id2][ic1][ic2]));
               }


               /*********** pseudo-scalar current ***********/
               if(qcd_NORM(qcd_GAMMA[5][id1][id3])>1e-4)
               {
		 block_p[j] = qcd_CADD(block_p[j], qcd_CMUL(backfor,
							    qcd_GAMMA[5][id1][id3]));
               }

               /*********** scalar current ***********/
               if(id1 == id3)
               {
		 block_s[j] = qcd_CADD(block_s[j], backfor);
               }

               /*********** local vector current ***********/
               for(mu=0; mu<4; mu++)
               if(qcd_NORM(qcd_GAMMA[mu][id1][id3])>1e-4)
               {
                  block_l[mu][j] = qcd_CADD(block_l[mu][j], qcd_CMUL(backfor,
                                                                     qcd_GAMMA[mu][id1][id3]));                             
               }

               /*********** local axial current ***********/     
               for(mu=0; mu<4; mu++)
               if(qcd_NORM(qcd_G5GAMMA[mu][id1][id3])>1e-4)
               {
                  block_a[mu][j] = qcd_CADD(block_a[mu][j], qcd_CMUL(backfor,                                                      
                                                                     qcd_G5GAMMA[mu][id1][id3]));
               }
               /*********** local tensor current ***********/     
               
               for(mu=0; mu<4; mu++)
		 {
		   if(qcd_NORM(g5sigmu0[mu][id1][id3])>1e-4)
		     {
		       block_t[0*4 + mu][j] = qcd_CADD(block_t[0*4 + mu][j], qcd_CMUL(backfor,
										      g5sigmu0[mu][id1][id3]));
		     }                        
		   if(qcd_NORM(g5sigmu1[mu][id1][id3])>1e-4)
		     {
		       block_t[1*4 + mu][j] = qcd_CADD(block_t[1*4 + mu][j], qcd_CMUL(backfor,
										      g5sigmu1[mu][id1][id3]));
		     }                        
		   if(qcd_NORM(g5sigmu2[mu][id1][id3])>1e-4)
		     {
		       block_t[2*4 + mu][j] = qcd_CADD(block_t[2*4 + mu][j], qcd_CMUL(backfor,
										      g5sigmu2[mu][id1][id3]));
		     }                        
		   if(qcd_NORM(g5sigmu3[mu][id1][id3])>1e-4)
		     {
		       block_t[3*4 + mu][j] = qcd_CADD(block_t[3*4 + mu][j], qcd_CMUL(backfor,
										      g5sigmu3[mu][id1][id3]));
		     }                        
		 }
            }//end id1,id3 loop
            for(id1=0; id1<4; id1++)
            for(id3=0; id3<4; id3++)
            {
               for(mu=0; mu<4; mu++)
               for(nu=0; nu<=mu; nu++)
               {                       
                  /*********** one derivative vector operator ***********/
                  if(qcd_NORM(qcd_GAMMA[mu][id1][id3])>1e-4)
                  {
                      block_vD[mu*4+nu][j] = qcd_CADD(block_vD[mu*4+nu][j],
                                                      qcd_CMUL(qcd_GAMMA[mu][id1][id3],
                                                               bdfmu[nu][id1][id3]));                                          
                  }
                  if(qcd_NORM(qcd_GAMMA[nu][id1][id3])>1e-4)
                  {  
                      block_vD[mu*4+nu][j] = qcd_CADD(block_vD[mu*4+nu][j],
                                                      qcd_CMUL(qcd_GAMMA[nu][id1][id3],
                                                               bdfmu[mu][id1][id3]));
                  }

                  /*********** one derivative axial and axial-antisymmetric operators ***********/
                  if(qcd_NORM(qcd_G5GAMMA[mu][id1][id3])>1e-4)
                  {
                     block_aD[mu*4+nu][j] = qcd_CADD(block_aD[mu*4+nu][j],
                                                     qcd_CMUL(qcd_G5GAMMA[mu][id1][id3],
                                                              bdfmu[nu][id1][id3]));

                     block_d1[mu*4+nu][j] = qcd_CADD(block_d1[mu*4+nu][j],
                                                     qcd_CMUL(qcd_G5GAMMA[mu][id1][id3],
                                                              bdfmu[nu][id1][id3]));
                  }

                  if(qcd_NORM(qcd_G5GAMMA[nu][id1][id3])>1e-4)
                  {
		    block_aD[mu*4+nu][j] = qcd_CADD(block_aD[mu*4+nu][j],
						    qcd_CMUL(qcd_G5GAMMA[nu][id1][id3],
							     bdfmu[mu][id1][id3]));

                     block_d1[mu*4+nu][j] = qcd_CSUB(block_d1[mu*4+nu][j],
                                                     qcd_CMUL(qcd_G5GAMMA[nu][id1][id3],
                                                              bdfmu[mu][id1][id3]));
                  }
               }//end mu nu loop
            }//end id1 id3 loop   
            
	   }//end j loop   
      }//end lt inside local block condition
      //#####################################################################   
      // Fourier transform the blocks and write the correlators

      if(myid==0) printf("time-slice-blocks calculated\n");

      if(t==t_start)         
      {
         k=0;
         if(myid==0)
         {
            printf("t_start = %i, t = %i: opening output files\n",t_start, t); fflush(stdout);
            printf("opening %s\n",corrloc_s_name); fflush(stdout);
            fp_corrloc_s = fopen(corrloc_s_name,"w");   
            if(fp_corrloc_s==NULL)
            {
               printf("failed to open %s for writing\n",corrloc_s_name); fflush(stdout);
               k=1;
            }
            printf("opening %s\n",corrloc_p_name); fflush(stdout);
            fp_corrloc_p = fopen(corrloc_p_name,"w");   
            if(fp_corrloc_p==NULL)
            {
               printf("failed to open %s for writing\n",corrloc_p_name); fflush(stdout);
               k=1;
            }
            printf("opening %s\n",corrnoe_v_name); fflush(stdout);
            fp_corrnoe_v = fopen(corrnoe_v_name,"w");   
            if(fp_corrnoe_v==NULL)
            {
               printf("failed to open %s for writing\n",corrnoe_v_name); fflush(stdout);
               k=1;
            }
            printf("%s opened\n",corrnoe_v_name); fflush(stdout);
            printf("opening %s\n",corrloc_v_name); fflush(stdout);
            fp_corrloc_v = fopen(corrloc_v_name,"w"); fflush(stdout);
            if(fp_corrloc_v==NULL)
            {
               printf("failed to open %s for writing\n",corrloc_v_name); fflush(stdout);
               k=1;
            }            
            printf("%s opened\n",corrloc_v_name); fflush(stdout);
            printf("opening %s\n",corrloc_a_name); fflush(stdout);
            fp_corrloc_a = fopen(corrloc_a_name,"w");
            if(fp_corrloc_a==NULL)
            {
               printf("failed to open %s for writing\n",corrloc_a_name); fflush(stdout);
               k=1;
            }   
            printf("%s opened\n",corrloc_a_name); fflush(stdout);
           
            fp_corrloc_t = fopen(corrloc_t_name,"w");
            if(fp_corrloc_t==NULL)
            {
               printf("failed to open %s for writing\n",corrloc_t_name); fflush(stdout);
               k=1;
            }

            printf("opening %s\n",corr_vD_name); fflush(stdout);
            fp_corr_vD = fopen(corr_vD_name,"w");   
            if(fp_corr_vD==NULL)
            {
               printf("failed to open %s for writing\n",corr_vD_name); fflush(stdout);
               k=1;
            }
            printf("%s opened\n",corr_vD_name); fflush(stdout);
            printf("opening %s\n",corr_aD_name); fflush(stdout);
            fp_corr_aD = fopen(corr_aD_name,"w");
            if(fp_corr_aD==NULL)
            {
               printf("failed to open %s for writing\n",corr_aD_name); fflush(stdout);
               k=1;
            }
            printf("%s opened\n",corr_aD_name); fflush(stdout);
/*
            fp_corr_tD = fopen(corr_tD_name,"w");
            if(fp_corr_tD==NULL)
            {
               printf("failed to open %s for writing\n",corr_tD_name); fflush(stdout);
               k=1;
            }
*/
            printf("opening %s\n",corr_d1_name); fflush(stdout);
            fp_corr_d1 = fopen(corr_d1_name,"w");
            if(fp_corr_d1==NULL)
            {
               printf("failed to open %s for writing\n",corr_d1_name); fflush(stdout);
               k=1;
            }
            printf("%s opened\n",corr_d1_name); fflush(stdout);
            fp_momlist = fopen(momlist_name,"r");
            if(fp_momlist==NULL)
            {
               printf("error opening %s for reading\n",momlist_name); fflush(stdout);
               k=1;
            }            
            printf("output files opened\n");fflush(stdout);
         }
         
         MPI_Bcast(&k,1,MPI_INT, 0, MPI_COMM_WORLD);
         if(k==1) exit(EXIT_FAILURE);
         
         if(myid==0) fscanf(fp_momlist,"%u\n",&numOfMom);
         MPI_Bcast(&numOfMom,1,MPI_INT, 0, MPI_COMM_WORLD);
         if(myid==0) printf("will read %i momenta combinations\n",numOfMom);
         
         mom = malloc(numOfMom*3*sizeof(qcd_int_4));
         
         if(myid==0)
         {
	   for(int j=0; j<numOfMom; j++)
            {
               fscanf(fp_momlist,"%i %i %i\n",&(mom[j][0]),&(mom[j][1]),&(mom[j][2]));
               //printf("got combination %i %i %i\n",mom[j][0],mom[j][1],mom[j][2]);  
            }   
            fclose(fp_momlist);   
         }
         MPI_Bcast(&(mom[0][0]),numOfMom*3,MPI_INT,0, MPI_COMM_WORLD);
         if(myid==0) printf("momenta list read and broadcasted\n");   

	 for(int j=0; j<4; j++)
	   {
	     corr[j] = malloc(numOfMom*sizeof(qcd_complex_16));
	     if(myid == 0)
	       corr_master[j] = malloc(numOfMom*sizeof(qcd_complex_16));
	   }
      }
      
      int j;
#pragma omp parallel for private(lx,ly,lz,v,x,y,z,tmp)
      for(j=0; j<numOfMom; j++)
	{
	  corr[0][j] = (qcd_complex_16) {0,0};
	  corr[1][j] = (qcd_complex_16) {0,0};
          
	  if(lt>=0 && lt<geo.lL[0])  //inside the local lattice, otherwise nothing to calculate
            for(lx=0; lx<geo.lL[1]; lx++)
	      for(ly=0; ly<geo.lL[2]; ly++)
		for(lz=0; lz<geo.lL[3]; lz++)
		  {
		    v = qcd_LEXIC0(lx,ly,lz,geo.lL);
		    x=lx+geo.Pos[1]*geo.lL[1] - x_src[1];
		    y=ly+geo.Pos[2]*geo.lL[2] - x_src[2];
		    z=lz+geo.Pos[3]*geo.lL[3] - x_src[3];
		    tmp = (((double)mom[j][0]*x)/geo.L[1] + ((double)mom[j][1]*y)/geo.L[2] + ((double)mom[j][2]*z)/geo.L[3])*2*M_PI;
		    qcd_complex_16 C2 = (qcd_complex_16) {cos(tmp), sin(tmp)};
		    corr[0][j]=qcd_CADD(corr[0][j], qcd_CMUL(block_s[v],C2));
		    corr[1][j]=qcd_CADD(corr[1][j], qcd_CMUL(block_p[v],C2));
		  }
	}
      MPI_Reduce(&(corr[0][0].re), &(corr_master[0][0].re), numOfMom*2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(&(corr[1][0].re), &(corr_master[1][0].re), numOfMom*2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

      if(myid==0) 
	for(int j=0; j<numOfMom; j++)
	  {
	    fprintf(fp_corrloc_s,"%i %+i %+i %+i %+e %+e\n",t,mom[j][0],mom[j][1],mom[j][2],sign*corr_master[0][j].re,sign*corr_master[0][j].im);
	    fprintf(fp_corrloc_p,"%i %+i %+i %+i %+e %+e\n",t,mom[j][0],mom[j][1],mom[j][2],sign*corr_master[1][j].re,sign*corr_master[1][j].im);
	  }
      
      for(mu=0; mu<4; mu++)
      {
#pragma omp parallel for private(lx,ly,lz,v,x,y,z,tmp)
         for(j=0; j<numOfMom; j++)
         {
	    corr[0][j] = (qcd_complex_16) {0,0};
            corr[1][j] = (qcd_complex_16) {0,0};
            corr[2][j] = (qcd_complex_16) {0,0};
            
            if(lt>=0 && lt<geo.lL[0])  //inside the local lattice, otherwise nothing to calculate
            for(lx=0; lx<geo.lL[1]; lx++)
            for(ly=0; ly<geo.lL[2]; ly++)
            for(lz=0; lz<geo.lL[3]; lz++)
            {
               v = qcd_LEXIC0(lx,ly,lz,geo.lL);
               x=lx+geo.Pos[1]*geo.lL[1] - x_src[1];
               y=ly+geo.Pos[2]*geo.lL[2] - x_src[2];
               z=lz+geo.Pos[3]*geo.lL[3] - x_src[3];
               tmp = (((double)mom[j][0]*x)/geo.L[1] + ((double)mom[j][1]*y)/geo.L[2] + ((double)mom[j][2]*z)/geo.L[3])*2*M_PI;
               qcd_complex_16 C2=(qcd_complex_16) {cos(tmp), sin(tmp)};
               corr[0][j]=qcd_CADD(corr[0][j], qcd_CMUL(block_n[mu][v],C2));
               corr[1][j]=qcd_CADD(corr[1][j], qcd_CMUL(block_l[mu][v],C2));
               corr[2][j]=qcd_CADD(corr[2][j], qcd_CMUL(block_a[mu][v],C2));
            }
         }

	 for(int j=0; j<3; j++)
	   MPI_Reduce(&(corr[j][0].re), &(corr_master[j][0].re), numOfMom*2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	 if(myid==0) 
	   for(int j=0; j<numOfMom; j++)
	     {
	       fprintf(fp_corrnoe_v,"%i %+i %+i %+i %+e %+e %i\n",t,mom[j][0],mom[j][1],mom[j][2],
		       sign*corr_master[0][j].re*0.25,sign*corr_master[0][j].im*0.25,mu);
	       fprintf(fp_corrloc_v,"%i %+i %+i %+i %+e %+e %i\n",t,mom[j][0],mom[j][1],mom[j][2],
		       sign*corr_master[1][j].re,sign*corr_master[1][j].im,mu);
	       fprintf(fp_corrloc_a,"%i %+i %+i %+i %+e %+e %i\n",t,mom[j][0],mom[j][1],mom[j][2],
		       sign*corr_master[2][j].re,sign*corr_master[2][j].im,mu);
	     }
      }//end mu loop


      for(mu=0; mu<4; mu++)
      for(nu=mu+1; nu<4; nu++)
      {
#pragma omp parallel for private(lx,ly,lz,v,x,y,z,tmp)
         for(j=0; j<numOfMom; j++)
         {
	    corr[0][j] = (qcd_complex_16) {0,0};
            
            if(lt>=0 && lt<geo.lL[0])  //inside the local lattice, otherwise nothing to calculate
            for(lx=0; lx<geo.lL[1]; lx++)
            for(ly=0; ly<geo.lL[2]; ly++)
            for(lz=0; lz<geo.lL[3]; lz++)
            {
               v = qcd_LEXIC0(lx,ly,lz,geo.lL);
               x=lx+geo.Pos[1]*geo.lL[1] - x_src[1];
               y=ly+geo.Pos[2]*geo.lL[2] - x_src[2];
               z=lz+geo.Pos[3]*geo.lL[3] - x_src[3];
               tmp = (((double) mom[j][0]*x)/geo.L[1] + ((double) mom[j][1]*y)/geo.L[2] + ((double) mom[j][2]*z)/geo.L[3])*2*M_PI;
               qcd_complex_16 C2=(qcd_complex_16) {cos(tmp), sin(tmp)};
               corr[0][j]=qcd_CADD(corr[0][j], qcd_CMUL(block_t[mu*4+nu][v],C2));
            }
         }
	 
	 for(int j=0; j<1; j++)   
	   MPI_Reduce(&(corr[j][0].re), &(corr_master[j][0].re), numOfMom*2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	 
	 if(myid==0) 
	   for(int j=0; j<numOfMom; j++)
	     {
               fprintf(fp_corrloc_t,"%i %+i %+i %+i %+e %+e %i %i\n",t,mom[j][0],mom[j][1],mom[j][2],
		       sign*corr_master[0][j].re,sign*corr_master[0][j].im,mu,nu);
	     }   
      }

      for(mu=0; mu<4; mu++)
      for(nu=0; nu<=mu; nu++)
      {
#pragma omp parallel for private(lx,ly,lz,v,x,y,z,tmp)
         for(j=0; j<numOfMom; j++)
         {
	    corr[0][j] = (qcd_complex_16) {0,0};
            corr[1][j] = (qcd_complex_16) {0,0};
            corr[2][j] = (qcd_complex_16) {0,0};
            //corr[3][j] = (qcd_complex_16) {0,0};
            
            if(lt>=0 && lt<geo.lL[0])  //inside the local lattice, otherwise nothing to calculate
            for(lx=0; lx<geo.lL[1]; lx++)
            for(ly=0; ly<geo.lL[2]; ly++)
            for(lz=0; lz<geo.lL[3]; lz++)
            {
               v = qcd_LEXIC0(lx,ly,lz,geo.lL);
               x=lx+geo.Pos[1]*geo.lL[1] - x_src[1];
               y=ly+geo.Pos[2]*geo.lL[2] - x_src[2];
               z=lz+geo.Pos[3]*geo.lL[3] - x_src[3];
               tmp = (((double) mom[j][0]*x)/geo.L[1] + ((double) mom[j][1]*y)/geo.L[2] + ((double) mom[j][2]*z)/geo.L[3])*2*M_PI;
               qcd_complex_16 C2=(qcd_complex_16) {cos(tmp), sin(tmp)};
               corr[0][j]=qcd_CADD(corr[0][j], qcd_CMUL(block_vD[mu*4+nu][v],C2));
               corr[1][j]=qcd_CADD(corr[1][j], qcd_CMUL(block_aD[mu*4+nu][v],C2));
               corr[2][j]=qcd_CADD(corr[2][j], qcd_CMUL(block_d1[mu*4+nu][v],C2));
               //corr[3][j]=qcd_CADD(corr[3][j], qcd_CMUL(block_t[mu][v],C2));
            }
         }
	 
	 for(int j=0; j<3; j++)   
	   MPI_Reduce(&(corr[j][0].re), &(corr_master[j][0].re), numOfMom*2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	 
	 if(myid==0) 
	   for(int j=0; j<numOfMom; j++)
	     {
               fprintf(fp_corr_vD,"%i %+i %+i %+i %+e %+e %i %i\n",t,mom[j][0],mom[j][1],mom[j][2],
		       sign*corr_master[0][j].re*0.125,sign*corr_master[0][j].im*0.125,mu,nu);
               fprintf(fp_corr_aD,"%i %+i %+i %+i %+e %+e %i %i\n",t,mom[j][0],mom[j][1],mom[j][2],
		       sign*corr_master[1][j].re*0.125,sign*corr_master[1][j].im*0.125,mu,nu);
            }   

	 if(myid==0)
	   if(mu != nu)
	     for(int j=0; j<numOfMom; j++)
	       fprintf(fp_corr_d1,"%i %+i %+i %+i %+e %+e %i %i\n",t,mom[j][0],mom[j][1],mom[j][2],
		       sign*corr_master[2][j].re*0.125,sign*corr_master[2][j].im*0.125,mu,nu);
	 
      }//end mu/nu loop    

      if(t==t_stop && myid==0)
      {
         fclose(fp_corrnoe_v);
         fclose(fp_corrloc_v);
         fclose(fp_corrloc_s);
         fclose(fp_corrloc_p);
         fclose(fp_corrloc_a);   
         fclose(fp_corr_vD);
         fclose(fp_corr_aD);
         fclose(fp_corr_d1);
         fclose(fp_corrloc_t);

	 for(int j=0; j<4; j++)
	   {
	     free(corr[j]);
	     if(myid == 0)
	       free(corr_master[j]);
	   }

      }

   }//end t-loop   
         
   //#####################################################################   
   // clean up
   if(myid==0) printf("cleaning up...\n");
   for(i=0;i<4;i++)
   {
      free(block_n[i]);
      free(block_l[i]);
      free(block_a[i]);
      free(block_t[i]);
   }
   free(block_n);
   free(block_l);
   free(block_a);
   free(block_s);
   free(block_p);
   for(i=0; i<16; i++)
   {
      free(block_vD[i]);
      free(block_aD[i]);
      free(block_d1[i]);
   }
   free(block_vD);
   free(block_aD);
   free(block_d1);
   free(mom);
   qcd_destroyPropagator(&prop);
   qcd_destroyPropagator(&backprop);
   qcd_destroyGaugeField(&u);
   qcd_destroyGeometry(&geo);
   MPI_Finalize();
}//end main
