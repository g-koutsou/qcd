/* zfac.c
 *
 * reads half fourier transformed propagator and
 * writes out:
 *
 * 1) the propagator with momentum p (12x12 matrix)
 *    S(p) = sum_xy exp(-ip(x-y)) < d(x) \bar d(y) >
 *
 * 2) the bare vertex function with momentum p (12x12 matrix)
 *    for a set of local and 1-derivative bilinear operators
 *    G(p) = sum_xyz exp(-ip(x-y)) < u(x) \bar u(z) J d(z') \bar d(y) >
 *    where
 *    \bar u(z) J d(z') is the bilinear operator
 *
 *
 * Tomasz Korzec 2009
 ************************************************************************/
 
 
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <qcd.h>


int main(int argc,char* argv[])
{
   qcd_uint_2 mu,nu,rho,ku,lu,c1,c2;          // various loop variables
   qcd_uint_4 i,j,k,lt,x,y,z,t,ip1,im1; 
   qcd_uint_2 id1,id2,id3,id4;
   qcd_uint_2 ic1,ic2,ic3,ic4;                    
   qcd_real_8 tmp;                            // general purpuse
   qcd_uint_2 xx[4];
   int faciip1, faciim1;                      // factors related to b.c.
  
   FILE *fp_vfun_v = NULL, *fp_vfun_a = NULL;   // output files
   FILE *fp_vfun_s = NULL, *fp_vfun_p = NULL;
   FILE *fp_vfun_t = NULL, *fp_vfun_vD = NULL, *fp_vfun_aD = NULL;
   FILE *fp_vfun_tD = NULL, *fp_vfun_d1 = NULL;
   FILE *fp_pprop = NULL;      
  
   int params_len;               // needed to read inputfiles
   char *params = NULL;                 // needed to read inputfiles

   char gauge_name[qcd_MAX_STRING_LENGTH];      // name of gauge-configuration file
   char vfun_s_name[qcd_MAX_STRING_LENGTH];     // output file name, local scalar density vertex function
   char vfun_p_name[qcd_MAX_STRING_LENGTH];     // output file name, local pseudoscalar density vertex function
   char vfun_v_name[qcd_MAX_STRING_LENGTH];     // output file name, local vector current vertex function
   char vfun_a_name[qcd_MAX_STRING_LENGTH];     // output file name, local axial current vertex function
   char vfun_t_name[qcd_MAX_STRING_LENGTH];     // output file name, local tensor current vertex function
   char vfun_vD_name[qcd_MAX_STRING_LENGTH];    // output file name, 1 derivative vector operator vertex function
   char vfun_aD_name[qcd_MAX_STRING_LENGTH];    // output file name, 1 derivative axial operator vertex function
   char vfun_tD_name[qcd_MAX_STRING_LENGTH];    // output file name, 1 derivative tensor operator vertex function
   char vfun_d1_name[qcd_MAX_STRING_LENGTH];    // output file name, 1 derivative d1 operator vertex function
   char pprop_name[qcd_MAX_STRING_LENGTH];      // name of output file, momentum propagator
   
   char param_name[qcd_MAX_STRING_LENGTH];      // name of parameter file

   qcd_geometry geo;                            // geometry structure
   qcd_propagator prop, lprop, rprop;           // propagator
   qcd_gaugeField u;                            // gauge field 
   char prop_name[qcd_MAX_STRING_LENGTH];       // name of inverted Fourier source

   qcd_real_8 theta[4] = {M_PI,0.0,0.0,0.0};    // antiperiodic b.c. in time
   qcd_uint_2 L[4];
   qcd_uint_2 P[4];
   qcd_complex_16 phase_factor, C, Cphase_factor;
   
   qcd_complex_16 Sup[4][4][3][3];              // Fourier transformed up-quark & down-quark
   qcd_complex_16 Sdown[4][4][3][3];            // Propagators

   qcd_complex_16 vfun_s[4][4][3][3];           // vertex function noether current
   qcd_complex_16 vfun_p[4][4][3][3];           // vertex function noether current
   qcd_complex_16 vfun_v[4][4][4][3][3];        // vertex function local vector current
   qcd_complex_16 vfun_a[4][4][4][3][3];        // vertex function local axial current
   qcd_complex_16 vfun_t[16][4][4][3][3];        // vertex function local tensor current
   qcd_complex_16 vfun_vD[16][4][4][3][3];      // vertex function 1 derivative vector operator
   qcd_complex_16 vfun_aD[16][4][4][3][3];      // vertex function 1 derivative axial operator
   qcd_complex_16 vfun_tD[64][4][4][3][3];      // vertex function 1 derivative tensor operator
   qcd_complex_16 vfun_d1[16][4][4][3][3];      // vertex function 1 derivative d1 operator
   qcd_complex_16 vfun_tmp[64][4][4][3][3];     // temp variable
   
   qcd_complex_16 lxr, lDmur[4];
   
   qcd_complex_16 g5sig[5][5][4][4];            // gamma_5 * [gamma_mu, gamma_nu] *1/2
   
   qcd_int_2  pn[4];                            //integer momentum
   qcd_real_8 p[4];                             //= n*2*pi/L


   int myid,numprocs, namelen;    
   char processor_name[MPI_MAX_PROCESSOR_NAME];     
   				 
				 
             
             
   //////////////////////////////////////////////////////////////////////////////////////          
   //set up MPI
   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD,&numprocs);         // num. of processes taking part in the calculation
   MPI_Comm_rank(MPI_COMM_WORLD,&myid);             // each process gets its ID
   MPI_Get_processor_name(processor_name,&namelen); // 
   
   for(i=0; i<5; i++)
   for(mu=0; mu<4; mu++)
   for(nu=0; nu<4; nu++)
   {  
      for(j=0; j<5; j++)
      {         
         g5sig[i][j][mu][nu]= (qcd_complex_16){0,0};

         for(ku=0; ku<4; ku++)
         for(lu=0; lu<4; lu++)
         {
            g5sig[i][j][mu][nu] = qcd_CADD(g5sig[i][j][mu][nu],qcd_CMUL(qcd_CMUL(qcd_GAMMA[5][mu][ku],
                                                                                 qcd_GAMMA[i][ku][lu]),
                                                                        qcd_GAMMA[j][lu][nu]));
            g5sig[i][j][mu][nu] = qcd_CSUB(g5sig[i][j][mu][nu],qcd_CMUL(qcd_CMUL(qcd_GAMMA[5][mu][ku],
                                                                                 qcd_GAMMA[j][ku][lu]),
                                                                        qcd_GAMMA[i][lu][nu]));
         }
         g5sig[i][j][mu][nu] = qcd_CSCALE(g5sig[i][j][mu][nu],0.5);
      }
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
        
   strcpy(prop_name,qcd_getParam("<propagator>",params,params_len));
   if(myid==0) printf("Got propagator file name: %s\n",prop_name);

   strcpy(gauge_name,qcd_getParam("<cfg_name>",params,params_len));
   if(myid==0) printf("Got conf name: %s\n",gauge_name);
                                 
   strcpy(vfun_s_name,qcd_getParam("<vertex_function_scalar_name>",params,params_len));
   if(myid==0) printf("Got output file name: %s\n",vfun_s_name);
   
   strcpy(vfun_p_name,qcd_getParam("<vertex_function_pseudoscalar_name>",params,params_len));
   if(myid==0) printf("Got output file name: %s\n",vfun_p_name);

   strcpy(vfun_v_name,qcd_getParam("<vertex_function_vector_name>",params,params_len));
   if(myid==0) printf("Got output file name: %s\n",vfun_v_name);

   strcpy(vfun_a_name,qcd_getParam("<vertex_function_axial_name>",params,params_len));
   if(myid==0) printf("Got output file name: %s\n",vfun_a_name);

   strcpy(vfun_t_name,qcd_getParam("<vertex_function_tensor_name>",params,params_len));
   if(myid==0) printf("Got output file name: %s\n",vfun_t_name);

   strcpy(vfun_vD_name,qcd_getParam("<vertex_function_vectorD_name>",params,params_len));
   if(myid==0) printf("Got output file name: %s\n",vfun_vD_name);

   strcpy(vfun_aD_name,qcd_getParam("<vertex_function_axialD_name>",params,params_len));
   if(myid==0) printf("Got output file name: %s\n",vfun_aD_name);
   
   strcpy(vfun_tD_name,qcd_getParam("<vertex_function_tensorD_name>",params,params_len));
   if(myid==0) printf("Got output file name: %s\n",vfun_tD_name);
   
   strcpy(vfun_d1_name,qcd_getParam("<vertex_function_d1_name>",params,params_len));
   if(myid==0) printf("Got output file name: %s\n",vfun_d1_name);

   strcpy(pprop_name,qcd_getParam("<pprop_name>",params,params_len));
   if(myid==0) printf("Got output file name: %s\n",pprop_name);
         
   sscanf(qcd_getParam("<momentum>",params,params_len),"%hd %hd %hd %hd",&pn[0], &pn[1], &pn[2], &pn[3]);
   if(myid==0) printf("Got momentum: [(%i+0.5)*2*pi/%i, %i*2*pi/%i, %i*2*pi/%i, %i*2*pi/%i] \n",pn[0],L[0],pn[1],L[1],pn[2],L[2],pn[3],L[3]);
   p[0]=(pn[0]*2*M_PI+theta[0])/L[0];
   p[1]=(pn[1]*2*M_PI+theta[1])/L[1];
   p[2]=(pn[2]*2*M_PI+theta[2])/L[2];
   p[3]=(pn[3]*2*M_PI+theta[3])/L[3];
    
   free(params);
   
   
   //#####################################################################   
   // allocate memory
  
   j = 0;
   j += qcd_initPropagator(&prop, &geo);
   j += qcd_initPropagator(&lprop, &geo);
   j += qcd_initPropagator(&rprop, &geo);
   j += qcd_initGaugeField(&u, &geo);
   
   MPI_Allreduce(&j, &k, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
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
   
   // load propagator
   if(qcd_getPropagator(prop_name,qcd_PROP_CMI, &prop)) exit(EXIT_FAILURE);
   if(myid==0) printf("propagator loaded\n");   
   
   //################################################################################
   // transform propagators to basis with theta-periodic boundaries in the temporal direction
   for(lt=0; lt<geo.lL[0]; lt++)
   {
      t = lt + geo.Pos[0] * geo.lL[0];
      phase_factor   = (qcd_complex_16) {cos(theta[0]*t/geo.L[0]),sin(theta[0]*t/geo.L[0])};
      qcd_mulPropagatorC3d(&prop, phase_factor, lt);
   }
   if(myid==0) printf("propagators transformed to basis with theta-periodic boundary conditions\n");
   
   //transform propagators to physical basis
   for(i=0; i<geo.lV; i++)
   for(c1=0; c1<3; c1++)
   for(c2=0; c2<3; c2++)
   {
      //this works only with g5=diag(1 1 -1 -1)
      lprop.D[i][0][0][c1][c2] = (qcd_complex_16){-prop.D[i][0][0][c2][c1].im, -prop.D[i][0][0][c2][c1].re};
      rprop.D[i][0][0][c1][c2] = (qcd_complex_16){ prop.D[i][0][0][c1][c2].im, -prop.D[i][0][0][c1][c2].re};
      lprop.D[i][0][1][c1][c2] = (qcd_complex_16){-prop.D[i][1][0][c2][c1].im, -prop.D[i][1][0][c2][c1].re};
      rprop.D[i][0][1][c1][c2] = (qcd_complex_16){ prop.D[i][0][1][c1][c2].im, -prop.D[i][0][1][c1][c2].re};
      lprop.D[i][0][2][c1][c2] = (qcd_complex_16){ prop.D[i][2][0][c2][c1].re, -prop.D[i][2][0][c2][c1].im};
      rprop.D[i][0][2][c1][c2] = (qcd_complex_16){ prop.D[i][0][2][c1][c2].re,  prop.D[i][0][2][c1][c2].im};
      lprop.D[i][0][3][c1][c2] = (qcd_complex_16){ prop.D[i][3][0][c2][c1].re, -prop.D[i][3][0][c2][c1].im};
      rprop.D[i][0][3][c1][c2] = (qcd_complex_16){ prop.D[i][0][3][c1][c2].re,  prop.D[i][0][3][c1][c2].im};
      lprop.D[i][1][0][c1][c2] = (qcd_complex_16){-prop.D[i][0][1][c2][c1].im, -prop.D[i][0][1][c2][c1].re};
      rprop.D[i][1][0][c1][c2] = (qcd_complex_16){ prop.D[i][1][0][c1][c2].im, -prop.D[i][1][0][c1][c2].re};
      lprop.D[i][1][1][c1][c2] = (qcd_complex_16){-prop.D[i][1][1][c2][c1].im, -prop.D[i][1][1][c2][c1].re};
      rprop.D[i][1][1][c1][c2] = (qcd_complex_16){ prop.D[i][1][1][c1][c2].im, -prop.D[i][1][1][c1][c2].re};
      lprop.D[i][1][2][c1][c2] = (qcd_complex_16){ prop.D[i][2][1][c2][c1].re, -prop.D[i][2][1][c2][c1].im};
      rprop.D[i][1][2][c1][c2] = (qcd_complex_16){ prop.D[i][1][2][c1][c2].re,  prop.D[i][1][2][c1][c2].im};
      lprop.D[i][1][3][c1][c2] = (qcd_complex_16){ prop.D[i][3][1][c2][c1].re, -prop.D[i][3][1][c2][c1].im};
      rprop.D[i][1][3][c1][c2] = (qcd_complex_16){ prop.D[i][1][3][c1][c2].re,  prop.D[i][1][3][c1][c2].im};
      lprop.D[i][2][0][c1][c2] = (qcd_complex_16){ prop.D[i][0][2][c2][c1].re, -prop.D[i][0][2][c2][c1].im};
      rprop.D[i][2][0][c1][c2] = (qcd_complex_16){ prop.D[i][2][0][c1][c2].re,  prop.D[i][2][0][c1][c2].im};
      lprop.D[i][2][1][c1][c2] = (qcd_complex_16){ prop.D[i][1][2][c2][c1].re, -prop.D[i][1][2][c2][c1].im};
      rprop.D[i][2][1][c1][c2] = (qcd_complex_16){ prop.D[i][2][1][c1][c2].re,  prop.D[i][2][1][c1][c2].im};
      lprop.D[i][2][2][c1][c2] = (qcd_complex_16){ prop.D[i][2][2][c2][c1].im,  prop.D[i][2][2][c2][c1].re};
      rprop.D[i][2][2][c1][c2] = (qcd_complex_16){-prop.D[i][2][2][c1][c2].im,  prop.D[i][2][2][c1][c2].re};
      lprop.D[i][2][3][c1][c2] = (qcd_complex_16){ prop.D[i][3][2][c2][c1].im,  prop.D[i][3][2][c2][c1].re};
      rprop.D[i][2][3][c1][c2] = (qcd_complex_16){-prop.D[i][2][3][c1][c2].im,  prop.D[i][2][3][c1][c2].re};
      lprop.D[i][3][0][c1][c2] = (qcd_complex_16){ prop.D[i][0][3][c2][c1].re, -prop.D[i][0][3][c2][c1].im};
      rprop.D[i][3][0][c1][c2] = (qcd_complex_16){ prop.D[i][3][0][c1][c2].re,  prop.D[i][3][0][c1][c2].im};
      lprop.D[i][3][1][c1][c2] = (qcd_complex_16){ prop.D[i][1][3][c2][c1].re, -prop.D[i][1][3][c2][c1].im};
      rprop.D[i][3][1][c1][c2] = (qcd_complex_16){ prop.D[i][3][1][c1][c2].re,  prop.D[i][3][1][c1][c2].im};
      lprop.D[i][3][2][c1][c2] = (qcd_complex_16){ prop.D[i][2][3][c2][c1].im,  prop.D[i][2][3][c2][c1].re};
      rprop.D[i][3][2][c1][c2] = (qcd_complex_16){-prop.D[i][3][2][c1][c2].im,  prop.D[i][3][2][c1][c2].re};
      lprop.D[i][3][3][c1][c2] = (qcd_complex_16){ prop.D[i][3][3][c2][c1].im,  prop.D[i][3][3][c2][c1].re};
      rprop.D[i][3][3][c1][c2] = (qcd_complex_16){-prop.D[i][3][3][c1][c2].im,  prop.D[i][3][3][c1][c2].re};
   }
   if(myid==0) printf("propagators transformed to physical basis\n");
   qcd_waitall(&geo);

   
   ////////////////////////////////////////////////////////////////////////////////////////////////////////
   //create and print out the fourier transformed propagators
   memset(&(Sdown[0][0][0][0].re),0,12*12*sizeof(qcd_complex_16));
   memset(&(Sup[0][0][0][0].re),0,12*12*sizeof(qcd_complex_16));
   for(t=0; t<geo.lL[0]; t++)
   for(x=0; x<geo.lL[1]; x++)
   for(y=0; y<geo.lL[2]; y++)
   for(z=0; z<geo.lL[3]; z++)
   {
      tmp = (t+geo.lL[0]*geo.Pos[0])*p[0];
      tmp+= (x+geo.lL[1]*geo.Pos[1])*p[1];
      tmp+= (y+geo.lL[2]*geo.Pos[2])*p[2];
      tmp+= (z+geo.lL[3]*geo.Pos[3])*p[3];
      phase_factor = (qcd_complex_16) {cos(tmp), sin(tmp)};
      Cphase_factor= qcd_CONJ(phase_factor);
      i = qcd_LEXIC(t,x,y,z,geo.lL);
      for(mu=0; mu<4; mu++)
      for(nu=0; nu<4; nu++)
      for(c1=0; c1<3; c1++)
      for(c2=0; c2<3; c2++)
      {
         Sdown[mu][nu][c1][c2] = qcd_CADD(Sdown[mu][nu][c1][c2], qcd_CMUL(Cphase_factor,rprop.D[i][mu][nu][c1][c2]));
         Sup[mu][nu][c1][c2]   = qcd_CSUB(Sup[mu][nu][c1][c2],   qcd_CMUL( phase_factor,lprop.D[i][mu][nu][c1][c2]));
      }
   }   
   
   if(myid==0) printf("local momentum propagators calculated\n");

   if(myid==0)
   {
      j=0;
      if( (fp_vfun_s=fopen(vfun_s_name,"w"))==NULL) j++;
      if( (fp_vfun_p=fopen(vfun_p_name,"w"))==NULL) j++;
      if( (fp_vfun_v=fopen(vfun_v_name,"w"))==NULL) j++;
      if( (fp_vfun_a=fopen(vfun_a_name,"w"))==NULL) j++;
      if( (fp_vfun_t=fopen(vfun_t_name,"w"))==NULL) j++;
      if( (fp_vfun_vD=fopen(vfun_vD_name,"w"))==NULL) j++;
      if( (fp_vfun_aD=fopen(vfun_aD_name,"w"))==NULL) j++;
      if( (fp_vfun_tD=fopen(vfun_tD_name,"w"))==NULL) j++;
      if( (fp_vfun_d1=fopen(vfun_d1_name,"w"))==NULL) j++;
      if( (fp_pprop=fopen(pprop_name,"w"))==NULL) j++;
      if(j>0) fprintf(stderr,"Error while opening output files for writing!\n");
   }
   MPI_Bcast(&j,1,MPI_INT, 0, MPI_COMM_WORLD);
   if(j>0) exit(EXIT_FAILURE);
   
   if(myid==0) printf("all output files ready for writing\n");   
      
   for(mu=0; mu<4; mu++)
   for(nu=0; nu<4; nu++)
   for(c1=0; c1<3; c1++)
   for(c2=0; c2<3; c2++)
   {
      MPI_Reduce(&(Sdown[mu][nu][c1][c2].re), &(C.re), 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(&(Sdown[mu][nu][c1][c2].im), &(C.im), 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      if(myid==0)
         fprintf(fp_pprop,"d %i %i %i %i %e %e\n",mu, nu, c1, c2, C.re/geo.V, C.im/geo.V);
      MPI_Reduce(&(Sup[mu][nu][c1][c2].re), &(C.re), 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(&(Sup[mu][nu][c1][c2].im), &(C.im), 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      if(myid==0)
         fprintf(fp_pprop,"u %i %i %i %i %e %e\n",mu, nu, c1, c2, C.re/geo.V, C.im/geo.V);
   }

   if(myid==0) printf("global momentum propagators calculated and written\n");

   ////////////////////////////////////////////////////////////////////////////////////////////////////////
   //create and print out the vertex functions
   qcd_communicatePropagatorPM(&lprop);
   qcd_waitall(&geo);
   qcd_communicatePropagatorPM(&rprop);
   qcd_waitall(&geo);
     
   memset(&(vfun_s[0][0][0][0].re),0,4*4*3*3*sizeof(qcd_complex_16));
   memset(&(vfun_p[0][0][0][0].re),0,4*4*3*3*sizeof(qcd_complex_16));
   memset(&(vfun_v[0][0][0][0][0].re),0,4*4*4*3*3*sizeof(qcd_complex_16));
   memset(&(vfun_a[0][0][0][0][0].re),0,4*4*4*3*3*sizeof(qcd_complex_16));
   memset(&(vfun_t[0][0][0][0][0].re),0,16*4*4*3*3*sizeof(qcd_complex_16));
   memset(&(vfun_vD[0][0][0][0][0].re),0,16*4*4*3*3*sizeof(qcd_complex_16));
   memset(&(vfun_aD[0][0][0][0][0].re),0,16*4*4*3*3*sizeof(qcd_complex_16));
   memset(&(vfun_tD[0][0][0][0][0].re),0,64*4*4*3*3*sizeof(qcd_complex_16));
   memset(&(vfun_d1[0][0][0][0][0].re),0,16*4*4*3*3*sizeof(qcd_complex_16));
   
   for(i=0; i<geo.lV; i++)
   {   
      qcd_antilexic(xx, i, geo.lL);
      
      for(id1=0; id1<4; id1++)
      for(id2=0; id2<4; id2++)
      for(id3=0; id3<4; id3++)
      for(id4=0; id4<4; id4++)
      for(ic1=0; ic1<3; ic1++)
      for(ic4=0; ic4<3; ic4++)
      { 
         lxr   = (qcd_complex_16) {0,0};
         for(ic2=0; ic2<3; ic2++)
         {
            lxr = qcd_CADD(lxr, qcd_CMUL(lprop.D[i][id1][id2][ic1][ic2],
                                         rprop.D[i][id3][id4][ic2][ic4]));
         }
         
         for(mu=0; mu<4; mu++)
         {
            ip1 = geo.plus[i][mu];              
            im1 = geo.minus[i][mu];
            
            // handle cases where fields wrap around the T direction
            faciip1=1;
            faciim1=1;
            if(mu==0)
            {
               if((xx[0]+geo.Pos[0]*geo.lL[0]) == 0)
                  faciim1 = -1;
               if((xx[0]+geo.Pos[0]*geo.lL[0]) == geo.L[0]-1)
                  faciip1 = -1;   
            }
            
            lDmur[mu] = (qcd_complex_16){0,0};
            
            if((faciim1 == -1) || (faciip1 == -1))
            for(ic2=0; ic2<3; ic2++)   
            for(ic3=0; ic3<3; ic3++)
            {
               //calculate: propagator(i) D_mu propagator(i)

               // x x x+mu
               lDmur[mu] = qcd_CADD(lDmur[mu],
                                    qcd_CMUL(lprop.D[i][id1][id2][ic1][ic2],
                                             qcd_CMUL(u.D[i][mu][ic2][ic3],
                                                      qcd_CSCALE(rprop.D[ip1][id3][id4][ic3][ic4],faciip1))));
               // x  x-mu  x-mu
               lDmur[mu] = qcd_CSUB(lDmur[mu],
                                    qcd_CMUL(lprop.D[i][id1][id2][ic1][ic2],
                                             qcd_CMUL(qcd_CONJ(u.D[im1][mu][ic3][ic2]),
                                                      qcd_CSCALE(rprop.D[im1][id3][id4][ic3][ic4],faciim1))));
               // x+mu  x  x
               lDmur[mu] = qcd_CSUB(lDmur[mu],
                                    qcd_CMUL(qcd_CSCALE(lprop.D[ip1][id1][id2][ic1][ic2],faciip1),
                                             qcd_CMUL(qcd_CONJ(u.D[i][mu][ic3][ic2]),
                                                      rprop.D[i][id3][id4][ic3][ic4])));
               // x-mu  x-mu  x
               lDmur[mu] = qcd_CADD(lDmur[mu],
                                    qcd_CMUL(qcd_CSCALE(lprop.D[im1][id1][id2][ic1][ic2],faciim1),
                                             qcd_CMUL(u.D[im1][mu][ic2][ic3],
                                                      rprop.D[i][id3][id4][ic3][ic4])));
            }//end ic2, ic3 loops
            else
            for(ic2=0; ic2<3; ic2++)   
            for(ic3=0; ic3<3; ic3++)
            {
               //calculate: propagator(i) D_mu propagator(i)

               // x x x+mu
               lDmur[mu] = qcd_CADD(lDmur[mu],
                                    qcd_CMUL(lprop.D[i][id1][id2][ic1][ic2],
                                             qcd_CMUL(u.D[i][mu][ic2][ic3],
                                                      rprop.D[ip1][id3][id4][ic3][ic4])));
               // x  x-mu  x-mu
               lDmur[mu] = qcd_CSUB(lDmur[mu],
                                    qcd_CMUL(lprop.D[i][id1][id2][ic1][ic2],
                                             qcd_CMUL(qcd_CONJ(u.D[im1][mu][ic3][ic2]),
                                                      rprop.D[im1][id3][id4][ic3][ic4])));
               // x+mu  x  x
               lDmur[mu] = qcd_CSUB(lDmur[mu],
                                    qcd_CMUL(lprop.D[ip1][id1][id2][ic1][ic2],
                                             qcd_CMUL(qcd_CONJ(u.D[i][mu][ic3][ic2]),
                                                      rprop.D[i][id3][id4][ic3][ic4])));
               // x-mu  x-mu  x
               lDmur[mu] = qcd_CADD(lDmur[mu],
                                    qcd_CMUL(lprop.D[im1][id1][id2][ic1][ic2],
                                             qcd_CMUL(u.D[im1][mu][ic2][ic3],
                                                      rprop.D[i][id3][id4][ic3][ic4])));
            }//end ic2, ic3 loops
            
         }//end mu loop




         //local operators
         /*********** local scalar density ***********/
         vfun_s[id1][id4][ic1][ic4] = qcd_CADD(vfun_s[id1][id4][ic1][ic4], qcd_CMUL(lxr,
                                                                                    qcd_ONE[id2][id3]));         
         
         /*********** local pseudoscalar density ***********/
         vfun_p[id1][id4][ic1][ic4] = qcd_CADD(vfun_p[id1][id4][ic1][ic4], qcd_CMUL(lxr,
                                                                                    qcd_GAMMA[5][id2][id3]));                  
         /*********** local vector current ***********/
         for(mu=0; mu<4; mu++)
         if(qcd_NORM(qcd_GAMMA[mu][id2][id3])>1e-4)
         {
            vfun_v[mu][id1][id4][ic1][ic4] = qcd_CADD(vfun_v[mu][id1][id4][ic1][ic4], qcd_CMUL(lxr,
                                                                                               qcd_GAMMA[mu][id2][id3]));
         }
         /*********** local axial current ***********/     
         for(mu=0; mu<4; mu++)
         if(qcd_NORM(qcd_G5GAMMA[mu][id2][id3])>1e-4)
         {
            vfun_a[mu][id1][id4][ic1][ic4] = qcd_CADD(vfun_a[mu][id1][id4][ic1][ic4], qcd_CMUL(lxr,
                                                                                               qcd_G5GAMMA[mu][id2][id3]));
         }
         /*********** local tensor current ***********/     
         for(mu=0; mu<4; mu++)
         for(nu=0; nu<4; nu++)
         if(nu != mu)
         if(qcd_NORM(g5sig[mu][nu][id2][id3])>1e-4)
         {
            vfun_t[mu*4+nu][id1][id4][ic1][ic4] = qcd_CADD(vfun_t[mu*4+nu][id1][id4][ic1][ic4], qcd_CMUL(lxr,
                                                                                                g5sig[mu][nu][id2][id3]));
         }


         for(mu=0; mu<4; mu++)
         for(nu=0; nu<=mu; nu++)
         {                       
            /*********** one derivative vector operator ***********/
            if(qcd_NORM(qcd_GAMMA[mu][id2][id3])>1e-4)
            {
                vfun_vD[mu*4+nu][id1][id4][ic1][ic4] = qcd_CADD(vfun_vD[mu*4+nu][id1][id4][ic1][ic4],
                                                                qcd_CMUL(qcd_GAMMA[mu][id2][id3],
                                                                         lDmur[nu]));                                          
            }
            if(qcd_NORM(qcd_GAMMA[nu][id2][id3])>1e-4)
            {  
                vfun_vD[mu*4+nu][id1][id4][ic1][ic4] = qcd_CADD(vfun_vD[mu*4+nu][id1][id4][ic1][ic4],
                                                                qcd_CMUL(qcd_GAMMA[nu][id2][id3],
                                                                         lDmur[mu]));
            }

            /*********** one derivative axial and axial-antisymmetric operators ***********/
            if(qcd_NORM(qcd_G5GAMMA[mu][id2][id3])>1e-4)
            {
               vfun_aD[mu*4+nu][id1][id4][ic1][ic4] = qcd_CADD(vfun_aD[mu*4+nu][id1][id4][ic1][ic4],
                                                               qcd_CMUL(qcd_G5GAMMA[mu][id2][id3],
                                                                        lDmur[nu]));

               vfun_d1[mu*4+nu][id1][id4][ic1][ic4] = qcd_CADD(vfun_d1[mu*4+nu][id1][id4][ic1][ic4],
                                                               qcd_CMUL(qcd_G5GAMMA[mu][id2][id3],
                                                                        lDmur[nu]));
            }

            if(qcd_NORM(qcd_G5GAMMA[nu][id2][id3])>1e-4)
            {
               vfun_aD[mu*4+nu][id1][id4][ic1][ic4] = qcd_CADD(vfun_aD[mu*4+nu][id1][id4][ic1][ic4],
                                                               qcd_CMUL(qcd_G5GAMMA[nu][id2][id3],
                                                                        lDmur[mu]));

               vfun_d1[mu*4+nu][id1][id4][ic1][ic4] = qcd_CSUB(vfun_d1[mu*4+nu][id1][id4][ic1][ic4],
                                                               qcd_CMUL(qcd_G5GAMMA[nu][id2][id3],
                                                                        lDmur[mu]));
            }
         }//end mu nu loop

         for(mu=0; mu<4; mu++)
         for(nu=0; nu<4; nu++)
         for(rho=0; rho<4; rho++)
         {                       
            /*********** one derivative tensor operator ***********/
            if(qcd_NORM(g5sig[mu][nu][id2][id3])>1e-4)
            {        
               vfun_tD[mu*16+nu*4+rho][id1][id4][ic1][ic4] = qcd_CADD(vfun_tD[mu*16+nu*4+rho][id1][id4][ic1][ic4],
                                                                      qcd_CMUL(g5sig[mu][nu][id2][id3],
                                                                               lDmur[rho]));
            }

            if(qcd_NORM(g5sig[mu][rho][id2][id3])>1e-4)
            {
               vfun_tD[mu*16+nu*4+rho][id1][id4][ic1][ic4] = qcd_CADD(vfun_tD[mu*16+nu*4+rho][id1][id4][ic1][ic4],
                                                                      qcd_CMUL(g5sig[mu][rho][id2][id3],
                                                                               lDmur[nu]));
            }
         }//end mu/nu/rho loop
      
      
      }//end id1,id2,id3,id4,ic1,ic4 loops   
   }//end i loop

   if(myid==0) printf("local vertex functions calculated\n");


   //global sums and output into files   
   
   MPI_Reduce(&(vfun_s[0][0][0][0].re), &(vfun_tmp[0][0][0][0][0].re), 4*4*3*3*2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   if(myid==0)   
   for(id1=0; id1<4; id1++)
   for(id4=0; id4<4; id4++)
   for(ic1=0; ic1<3; ic1++)
   for(ic4=0; ic4<3; ic4++)
      fprintf(fp_vfun_s,"%i %i %i %i %+e %+e\n",id1,id4,ic1,ic4,vfun_tmp[0][id1][id4][ic1][ic4].re/geo.V,vfun_tmp[0][id1][id4][ic1][ic4].im/geo.V);

   MPI_Reduce(&(vfun_p[0][0][0][0].re), &(vfun_tmp[0][0][0][0][0].re), 4*4*3*3*2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   if(myid==0)   
   for(id1=0; id1<4; id1++)
   for(id4=0; id4<4; id4++)
   for(ic1=0; ic1<3; ic1++)
   for(ic4=0; ic4<3; ic4++)
      fprintf(fp_vfun_p,"%i %i %i %i %+e %+e\n",id1,id4,ic1,ic4,vfun_tmp[0][id1][id4][ic1][ic4].re/geo.V,vfun_tmp[0][id1][id4][ic1][ic4].im/geo.V);
      
   MPI_Reduce(&(vfun_v[0][0][0][0][0].re), &(vfun_tmp[0][0][0][0][0].re), 4*4*4*3*3*2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   if(myid==0)   
   for(mu=0; mu<4; mu++)
   for(id1=0; id1<4; id1++)
   for(id4=0; id4<4; id4++)
   for(ic1=0; ic1<3; ic1++)
   for(ic4=0; ic4<3; ic4++)
      fprintf(fp_vfun_v,"%i %i %i %i %i %+e %+e\n",mu,id1,id4,ic1,ic4,vfun_tmp[mu][id1][id4][ic1][ic4].re/geo.V,vfun_tmp[mu][id1][id4][ic1][ic4].im/geo.V);
      
   MPI_Reduce(&(vfun_a[0][0][0][0][0].re), &(vfun_tmp[0][0][0][0][0].re), 4*4*4*3*3*2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   if(myid==0)   
   for(mu=0; mu<4; mu++)
   for(id1=0; id1<4; id1++)
   for(id4=0; id4<4; id4++)
   for(ic1=0; ic1<3; ic1++)
   for(ic4=0; ic4<3; ic4++)
      fprintf(fp_vfun_a,"%i %i %i %i %i %+e %+e\n",mu,id1,id4,ic1,ic4,vfun_tmp[mu][id1][id4][ic1][ic4].re/geo.V,vfun_tmp[mu][id1][id4][ic1][ic4].im/geo.V);
      
   MPI_Reduce(&(vfun_t[0][0][0][0][0].re), &(vfun_tmp[0][0][0][0][0].re), 16*4*4*3*3*2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   if(myid==0)   
   for(mu=0; mu<4; mu++)
   for(nu=0; nu<4; nu++)
   if(mu != nu)
   for(id1=0; id1<4; id1++)
   for(id4=0; id4<4; id4++)
   for(ic1=0; ic1<3; ic1++)
   for(ic4=0; ic4<3; ic4++)
      fprintf(fp_vfun_t,"%i %i %i %i %i %i %+e %+e\n",mu,nu,id1,id4,ic1,ic4,vfun_tmp[mu*4+nu][id1][id4][ic1][ic4].re/geo.V,vfun_tmp[mu*4+nu][id1][id4][ic1][ic4].im/geo.V);
      
   if(myid==0) printf("global 0-derivative vertex functions calculated and written\n");
   
   MPI_Reduce(&(vfun_vD[0][0][0][0][0].re), &(vfun_tmp[0][0][0][0][0].re), 16*4*4*3*3*2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   if(myid==0)
   for(mu=0; mu<4; mu++)
   for(nu=0; nu<=mu; nu++)
   for(id1=0; id1<4; id1++)
   for(id4=0; id4<4; id4++)
   for(ic1=0; ic1<3; ic1++)
   for(ic4=0; ic4<3; ic4++)
      fprintf(fp_vfun_vD,"%i %i %i %i %i %i %+e %+e\n",mu,nu,id1,id4,ic1,ic4,vfun_tmp[mu*4+nu][id1][id4][ic1][ic4].re/(8*geo.V),vfun_tmp[mu*4+nu][id1][id4][ic1][ic4].im/(8*geo.V));
            
   MPI_Reduce(&(vfun_aD[0][0][0][0][0].re), &(vfun_tmp[0][0][0][0][0].re), 16*4*4*3*3*2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   if(myid==0)
   for(mu=0; mu<4; mu++)
   for(nu=0; nu<=mu; nu++)
   for(id1=0; id1<4; id1++)
   for(id4=0; id4<4; id4++)
   for(ic1=0; ic1<3; ic1++)
   for(ic4=0; ic4<3; ic4++)
      fprintf(fp_vfun_aD,"%i %i %i %i %i %i %+e %+e\n",mu,nu,id1,id4,ic1,ic4,vfun_tmp[mu*4+nu][id1][id4][ic1][ic4].re/(8*geo.V),vfun_tmp[mu*4+nu][id1][id4][ic1][ic4].im/(8*geo.V));
      
   MPI_Reduce(&(vfun_d1[0][0][0][0][0].re), &(vfun_tmp[0][0][0][0][0].re), 16*4*4*3*3*2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   if(myid==0)
   for(mu=0; mu<4; mu++)
   for(nu=0; nu<=mu; nu++)
   for(id1=0; id1<4; id1++)
   for(id4=0; id4<4; id4++)
   for(ic1=0; ic1<3; ic1++)
   for(ic4=0; ic4<3; ic4++)
      fprintf(fp_vfun_d1,"%i %i %i %i %i %i %+e %+e\n",mu,nu,id1,id4,ic1,ic4,vfun_tmp[mu*4+nu][id1][id4][ic1][ic4].re/(8*geo.V),vfun_tmp[mu*4+nu][id1][id4][ic1][ic4].im/(8*geo.V));

   MPI_Reduce(&(vfun_tD[0][0][0][0][0].re), &(vfun_tmp[0][0][0][0][0].re), 64*4*4*3*3*2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   if(myid==0)
   for(mu=0; mu<4; mu++)
   for(nu=0; nu<4; nu++)
   for(rho=0; rho<4; rho++)
   for(id1=0; id1<4; id1++)
   for(id4=0; id4<4; id4++)
   for(ic1=0; ic1<3; ic1++)
   for(ic4=0; ic4<3; ic4++)
      fprintf(fp_vfun_tD,"%i %i %i %i %i %i %i %+e %+e\n",mu,nu,rho,id1,id4,ic1,ic4,vfun_tmp[mu*16+nu*4+rho][id1][id4][ic1][ic4].re/(8*geo.V),vfun_tmp[mu*16+nu*4+rho][id1][id4][ic1][ic4].im/(8*geo.V));

   if(myid==0) printf("global 1-derivative vertex functions calculated and written\n");

   if(myid==0)
   {
      fclose(fp_vfun_s);
      fclose(fp_vfun_p);
      fclose(fp_vfun_v);
      fclose(fp_vfun_a);
      fclose(fp_vfun_t);
      fclose(fp_vfun_vD);
      fclose(fp_vfun_aD);
      fclose(fp_vfun_tD);
      fclose(fp_vfun_d1);
      fclose(fp_pprop);
   }
   
   
   qcd_destroyPropagator(&lprop);
   qcd_destroyPropagator(&rprop);
   qcd_destroyPropagator(&prop);
   qcd_destroyGaugeField(&u);
   qcd_destroyGeometry(&geo);
   MPI_Finalize();
}//end main 
