/* twop_all_hadr.c
 *
 * reads forward propagators
 * and creates two point functions for the 40 baryons
 *
 * Christos Kallidonis
 *
 * April 2012
 *
 * Mods: June 2015, Added contractions for pseudoscalar, scalar, vector and axial vector mesons
 ****************************************/
 
/* SAMPLE INPUT FILE

   echo "<processors_txyz>1 8 8 8</processors_txyz>" 		  > twop_0100.ini
   echo "<lattice_txyz>64 32 32 32</lattice_txyz>"          >> twop_0100.ini
   echo "<t>0 63</t>"                                       >> twop_0100.ini
   echo "<source_pos_txyz>0 0 0 0</source_pos_txyz>"        >> twop_0100.ini
   echo "<propagator_u>propu_list_0100.txt</propagator_u>"  >> twop_0100.ini
   echo "<propagator_d>propd_list_0100.txt</propagator_d>"  >> twop_0100.ini
   echo "<propagator_s>props_list_0100.txt</propagator_s>"  >> twop_0100.ini
   echo "<propagator_c>propc_list_0100.txt</propagator_c>"  >> twop_0100.ini
   echo "<propagators_udsc>1 1 0 0</propagators_udsc>"      >> twop_0100.ini
   echo "<unitary_sc>0</unitary_sc>"                        >> twop_0100.ini
   echo "<mesons>1</mesons>"                                >> twop_0100.ini
   echo "<cfg_name>conf.0100</cfg_name>"                    >> twop_0100.ini
   echo "<corr_name>twopt_0100</corr_name>"                 >> twop_0100.ini
   echo "<momenta_list>momentalist</momenta_list>"          >> twop_0100.ini
   echo "<alpha_gauss>4</alpha_gauss>"                      >> twop_0100.ini
   echo "<nsmear_gauss>110</nsmear_gauss>"                  >> twop_0100.ini
   echo "<alpha_APE>0.5</alpha_APE>"                        >> twop_0100.ini
   echo "<nsmear_APE>50</nsmear_APE>"                       >> twop_0100.ini
   echo "<particle_list>particle.list</particle_list>"      >> twop_0100.ini

*/ 
 
 
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <qcd.h>
#include <omp.h>

int main(int argc,char* argv[])
{

  qcd_uint_2 mu,nu,ku,lu,c1,c2,c3,c1p,c2p,c3p;// various loop variables
  qcd_uint_2 id1,id2,id3,cc1,cc2,al,be;
  qcd_uint_4 i,j,k, v,lx,ly,lz,ip1,im1,v3,lv3,tslices,iarr; 
  qcd_int_4 x,y,z;
  qcd_uint_2 ic1,ic2,ic3;                    //
  qcd_uint_4 x_src[4];                       // source and sink coordinates
  qcd_uint_4 t_sink, t_start, t_stop, t,lt;
  qcd_real_8 tmp;                            // general purpuse
  qcd_int_4 pc,p_id,p_num,*p_arr,*p_arr12,*p_arr32,m_id,nmom,p12,p32;
  qcd_uint_4 uyes,dyes,syes,cyes,mesyes;
  qcd_uint_4 unit_sc;
   
  FILE *fp_momlist,*fp_parlist;
   
  int params_len;                            // needed to read inputfiles
  char *params;                              // needed to read inputfiles

  char gauge_name[qcd_MAX_STRING_LENGTH];      // name of gauge-configuration file
  char corr_p_name[qcd_MAX_STRING_LENGTH];
  char param_name[qcd_MAX_STRING_LENGTH];      // name of parameter file  
  char momlist_name[qcd_MAX_STRING_LENGTH];    // name of momenta-list file
  char particle_list[qcd_MAX_STRING_LENGTH];
  char uprop_name[qcd_MAX_STRING_LENGTH];      // file names of up and down quark propagators
  char dprop_name[qcd_MAX_STRING_LENGTH];      
  char sprop_name[qcd_MAX_STRING_LENGTH];      // file names of up and down quark propagators
  char cprop_name[qcd_MAX_STRING_LENGTH];
   
   
  qcd_geometry geo;                            // geometry structure
  qcd_propagator uprop,sprop,uprop_pb,dprop_pb;                        // propagator
  qcd_propagator dprop,cprop,sprop_pb,cprop_pb;                        // propagator
  qcd_vector vec,vec_u,vec_d,vec_s,vec_c;                              // needed when smearing
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
  qcd_complex_16 ***bcorr12[16],***bcorr32[16][3],bcorrsum12,bcorrsum32[3];
  qcd_complex_16 **mcorrp[19],**mcorra[19],mcorrsum[38];
  qcd_real_8 plaq;
  qcd_int_4 ctr, ctr2;
  qcd_int_2 cg5cg5_ind[16*16][4];
  qcd_complex_16 cg5cg5_val[16*16];
  qcd_complex_16 one_plus_ig5[4],one_minus_ig5[4],g5[4]; //-for transformation purposes
   
  qcd_complex_16 *block12[4][4],*block32[4][4],*blocknp[4][4];                       // to store the block (2pt function before FT)
  qcd_complex_16 *udblock,*dublock,*sdblock,*sublock,*cdblock,*cublock;          // for mesons
  qcd_complex_16 *a0ublock,*rho1ublock,*rho2ublock,*rho3ublock,*a11ublock,*a12ublock,*a13ublock; // for mesons
  qcd_complex_16 *a0dblock,*rho1dblock,*rho2dblock,*rho3dblock,*a11dblock,*a12dblock,*a13dblock; // for mesons
  qcd_complex_16 *Ks1ublock,*Ks2ublock,*Ks3ublock,*Ds1ublock,*Ds2ublock,*Ds3ublock,*K11ublock,*K12ublock,*K13ublock; // for mesons
  qcd_complex_16 *Ks1dblock,*Ks2dblock,*Ks3dblock,*Ds1dblock,*Ds2dblock,*Ds3dblock,*K11dblock,*K12dblock,*K13dblock; // for mesons
  qcd_complex_16 scalarGamma,vectorGamma[3],axVectorGamma[3]; // for mesons


  qcd_int_4 (*mom)[3];                         // momenta-list

  int myid,numprocs, namelen;    
  char processor_name[MPI_MAX_PROCESSOR_NAME];
   				 
				 
             
  qcd_complex_16 i_im ; // Imaginary i         
             
  //set up MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);         // num. of processes taking part in the calculation
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);             // each process gets its ID
  MPI_Get_processor_name(processor_name,&namelen); // 

  int nthreads;

#pragma omp parallel
  {

    nthreads=omp_get_num_threads();
  }

  if(myid==0) printf("Running OpenMP with num_threads=%d\n",nthreads);

               
   
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
  tslices = t_stop-t_start+1;
  if(myid==0)  printf("Got sink time slices: %d ... %d\n Will run for a total of %d time slices\n",t_start,t_stop,tslices);
   
   
   
                     	  
  sscanf(qcd_getParam("<source_pos_txyz>",params,params_len),"%d %d %d %d",&x_src[0],&x_src[1],&x_src[2],&x_src[3]);
  if(myid==0) printf("Got source coords: %d %d %d %d\n",x_src[0],x_src[1],x_src[2],x_src[3]);
   
  sscanf(qcd_getParam("<propagators_udsc>",params,params_len),"%d %d %d %d",&uyes,&dyes,&syes,&cyes);
  if(myid==0){
    printf("Will use propagators:\n");  
    if(uyes) printf("up\n");
    if(dyes) printf("down\n");
    if(syes) printf("strange\n");
    if(cyes) printf("charm\n");
  }
     
  if(syes || cyes){
    sscanf(qcd_getParam("<unitary_sc>",params,params_len),"%d",&unit_sc);
    if(myid==0){
      if(unit_sc) printf("Using unitary setup strange and charm propagators\n");
      else printf("Using mixed action setup strange and charm propagators\n");
    }
  }

  sscanf(qcd_getParam("<mesons>",params,params_len),"%d",&mesyes);
  if(myid==0){
    if(mesyes) printf("Will do contractions for mesons as well\n");
    else printf("Will NOT do contractions for mesons\n");
  }

     
  strcpy(uprop_name,qcd_getParam("<propagator_u>",params,params_len));
  if(myid==0) printf("Got propagator file name: %s\n",uprop_name);
  strcpy(dprop_name,qcd_getParam("<propagator_d>",params,params_len));
  if(myid==0) printf("Got propagator file name: %s\n",dprop_name);
  strcpy(sprop_name,qcd_getParam("<propagator_s>",params,params_len));
  if(myid==0) printf("Got propagator file name: %s\n",sprop_name);   
  strcpy(cprop_name,qcd_getParam("<propagator_c>",params,params_len));
  if(myid==0) printf("Got propagator file name: %s\n",cprop_name);   
   
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
   
  strcpy(particle_list,qcd_getParam("<particle_list>",params,params_len));
  if(myid==0) printf("Got particle-list file name: %s\n",particle_list);
   
  free(params);

  lv3 = geo.lL[1]*geo.lL[2]*geo.lL[3];

         
  //#####################################################################   
  // allocate memory
  // load gauge-field and APE-smear it
   
  if(myid==0) printf("Memory allocations...\n");
   
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
   


  if(myid==0) printf("gauge-field APE smeared\n");
  plaq = qcd_calculatePlaquette(&uAPE);
  if(myid==0) printf("new plaquette = %e\n",plaq);  

  j = 0;
  if (uyes) j += qcd_initPropagator(&uprop, &geo);
  if (dyes) j += qcd_initPropagator(&dprop, &geo);
  if (syes) j += qcd_initPropagator(&sprop, &geo);
  if (cyes) j += qcd_initPropagator(&cprop, &geo);
   
  if (uyes) j += qcd_initPropagator(&uprop_pb, &geo);
  if (dyes) j += qcd_initPropagator(&dprop_pb, &geo);
  if (syes) j += qcd_initPropagator(&sprop_pb, &geo);
  if (cyes) j += qcd_initPropagator(&cprop_pb, &geo);
   
  for(mu=0;mu<4;mu++){
    for(nu=0;nu<4;nu++){
      block12[mu][nu] = (qcd_complex_16*) malloc(geo.lV3*sizeof(qcd_complex_16));
      block32[mu][nu] = (qcd_complex_16*) malloc(geo.lV3*sizeof(qcd_complex_16));
      blocknp[mu][nu] = (qcd_complex_16*) malloc(geo.lV3*sizeof(qcd_complex_16));
			
      if(block12[mu][nu]==NULL){
	if(myid==0) printf("Block12 %d %d not properly initialized\n",mu,nu);
	exit(EXIT_FAILURE);
      }
      if(block32[mu][nu]==NULL){
	if(myid==0) printf("Block32 %d %d not properly initialized\n",mu,nu);
	exit(EXIT_FAILURE);
      }
      if(blocknp[mu][nu]==NULL){
	if(myid==0) printf("Blocknp %d %d not properly initialized\n",mu,nu);
	exit(EXIT_FAILURE);
      }
    }
  }
	
  udblock = (qcd_complex_16*) malloc(geo.lV3*sizeof(qcd_complex_16));
  dublock = (qcd_complex_16*) malloc(geo.lV3*sizeof(qcd_complex_16));
  sublock = (qcd_complex_16*) malloc(geo.lV3*sizeof(qcd_complex_16));
  sdblock = (qcd_complex_16*) malloc(geo.lV3*sizeof(qcd_complex_16));
  cublock = (qcd_complex_16*) malloc(geo.lV3*sizeof(qcd_complex_16));	
  cdblock = (qcd_complex_16*) malloc(geo.lV3*sizeof(qcd_complex_16));

  a0ublock   = (qcd_complex_16*) malloc(geo.lV3*sizeof(qcd_complex_16));
  a0dblock   = (qcd_complex_16*) malloc(geo.lV3*sizeof(qcd_complex_16));
  rho1ublock = (qcd_complex_16*) malloc(geo.lV3*sizeof(qcd_complex_16));
  rho1dblock = (qcd_complex_16*) malloc(geo.lV3*sizeof(qcd_complex_16));
  rho2ublock = (qcd_complex_16*) malloc(geo.lV3*sizeof(qcd_complex_16)); 
  rho2dblock = (qcd_complex_16*) malloc(geo.lV3*sizeof(qcd_complex_16));
  rho3ublock = (qcd_complex_16*) malloc(geo.lV3*sizeof(qcd_complex_16));
  rho3dblock = (qcd_complex_16*) malloc(geo.lV3*sizeof(qcd_complex_16));
  a11ublock  = (qcd_complex_16*) malloc(geo.lV3*sizeof(qcd_complex_16));
  a11dblock  = (qcd_complex_16*) malloc(geo.lV3*sizeof(qcd_complex_16));
  a12ublock  = (qcd_complex_16*) malloc(geo.lV3*sizeof(qcd_complex_16));
  a12dblock  = (qcd_complex_16*) malloc(geo.lV3*sizeof(qcd_complex_16));
  a13ublock  = (qcd_complex_16*) malloc(geo.lV3*sizeof(qcd_complex_16));
  a13dblock  = (qcd_complex_16*) malloc(geo.lV3*sizeof(qcd_complex_16));
  Ks1ublock  = (qcd_complex_16*) malloc(geo.lV3*sizeof(qcd_complex_16));
  Ks1dblock  = (qcd_complex_16*) malloc(geo.lV3*sizeof(qcd_complex_16));
  Ks2ublock  = (qcd_complex_16*) malloc(geo.lV3*sizeof(qcd_complex_16));
  Ks2dblock  = (qcd_complex_16*) malloc(geo.lV3*sizeof(qcd_complex_16));
  Ks3ublock  = (qcd_complex_16*) malloc(geo.lV3*sizeof(qcd_complex_16));
  Ks3dblock  = (qcd_complex_16*) malloc(geo.lV3*sizeof(qcd_complex_16));
  Ds1ublock  = (qcd_complex_16*) malloc(geo.lV3*sizeof(qcd_complex_16));
  Ds1dblock  = (qcd_complex_16*) malloc(geo.lV3*sizeof(qcd_complex_16));
  Ds2ublock  = (qcd_complex_16*) malloc(geo.lV3*sizeof(qcd_complex_16));
  Ds2dblock  = (qcd_complex_16*) malloc(geo.lV3*sizeof(qcd_complex_16));
  Ds3ublock  = (qcd_complex_16*) malloc(geo.lV3*sizeof(qcd_complex_16));
  Ds3dblock  = (qcd_complex_16*) malloc(geo.lV3*sizeof(qcd_complex_16));
  K11ublock  = (qcd_complex_16*) malloc(geo.lV3*sizeof(qcd_complex_16));
  K11dblock  = (qcd_complex_16*) malloc(geo.lV3*sizeof(qcd_complex_16));
  K12ublock  = (qcd_complex_16*) malloc(geo.lV3*sizeof(qcd_complex_16));
  K12dblock  = (qcd_complex_16*) malloc(geo.lV3*sizeof(qcd_complex_16));
  K13ublock  = (qcd_complex_16*) malloc(geo.lV3*sizeof(qcd_complex_16));
  K13dblock  = (qcd_complex_16*) malloc(geo.lV3*sizeof(qcd_complex_16));

  MPI_Allreduce(&j, &k, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if(k>0)
    {
      if(myid==0) printf("not enough memory\n");
      exit(EXIT_FAILURE);
    }
  if(myid==0) printf("memory for propagators, blocks and gauge-field allocated\n");

  //-load momenta-list 	
  if(myid==0){     
    fp_momlist = fopen(momlist_name,"r");
    if(fp_momlist==NULL){
      printf("failed to open %s for reading\n",momlist_name);
      k=1;
    }
  }
  MPI_Bcast(&k,1,MPI_INT, 0, MPI_COMM_WORLD);
  if(k==1) exit(EXIT_FAILURE);
  
  
  if(myid==0) fscanf(fp_momlist,"%i\n",&nmom);
  MPI_Bcast(&nmom,1,MPI_INT, 0, MPI_COMM_WORLD);
  if(myid==0) printf("will read %i momenta combinations\n",nmom);

  mom = malloc(nmom*3*sizeof(qcd_int_4));

  if(myid==0)
    {
      for(j=0; j<nmom; j++)
	{
	  fscanf(fp_momlist,"%i %i %i\n",&(mom[j][0]),&(mom[j][1]),&(mom[j][2]));
	  //printf("got combination %i %i %i\n",mom[j][0],mom[j][1],mom[j][2]);  
	}
      fclose(fp_momlist);   
    }
  MPI_Bcast(&(mom[0][0]),nmom*3,MPI_INT,0, MPI_COMM_WORLD);
  if(myid==0) printf("momenta list read and broadcasted\n"); 
  //------------------------------------------------------------------

  //-load particle list
  if(myid==0){     
    fp_parlist = fopen(particle_list,"r");
    if(fp_parlist==NULL){
      printf("failed to open %s for reading\n",particle_list);
      k=1;
    }
  }
  MPI_Bcast(&k,1,MPI_INT, 0, MPI_COMM_WORLD);
  if(k==1) exit(EXIT_FAILURE);
  
  
  if(myid==0) fscanf(fp_parlist,"%i\n",&p_num);
  MPI_Bcast(&p_num,1,MPI_INT, 0, MPI_COMM_WORLD);
  if(myid==0){
    if(p_num==1) printf("will run only for particle ");
    else printf("will run for the following %d particles\n",p_num);
  }

  p_arr = malloc(p_num*sizeof(qcd_int_4));	
	
  if(myid==0){
    p32 = 0;
    p12 = 0;
    for(pc=0;pc<p_num;pc++){
      fscanf(fp_parlist,"%i\n",&p_arr[pc]);
      printf("%s\n",particle_names[p_arr[pc]]);
			
      if(particles32[p_arr[pc]]) p32++;
      else p12++;			
    }	
    fclose(fp_parlist);
    if(p_num>1) printf("\n%d of those are spin-1/2 and %d are spin-3/2\n",p12,p32);			
  }
	
  MPI_Bcast(&p_arr[0],p_num,MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&p12,1,MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&p32,1,MPI_INT, 0, MPI_COMM_WORLD);
	
  if(p12) p_arr12 = malloc(p12*sizeof(qcd_int_4));
  if(p32) p_arr32 = malloc(p32*sizeof(qcd_int_4));
	
  i=0; j=0;
  for(pc=0;pc<p_num;pc++){
    p_id = p_arr[pc];
		
    if(particles32[p_id]){
      p_arr32[i++]=p_id;
    }
    else{
      p_arr12[j++]=p_id;
    }
  }
	
  //------------------------------------------------------------------


  //-allocate memory for correlators

  if(mesyes){
    for(i=0;i<19;i++){	
      mcorrp[i]   = malloc(tslices*sizeof(qcd_complex_16*));
      mcorra[i]   = malloc(tslices*sizeof(qcd_complex_16*));
			
      if(mcorrp[i]==NULL || mcorra[i]==NULL){
	if(myid==0) printf("mcorrp or mcorra %d not properly initialized\n",i);
	exit(EXIT_FAILURE);
      }
			
						
      for(j=0;j<tslices;j++){
	mcorrp[i][j]   = malloc(nmom*sizeof(qcd_complex_16));
	mcorra[i][j]   = malloc(nmom*sizeof(qcd_complex_16));
				
	if(mcorrp[i][j]==NULL || mcorra[i][j]==NULL){
	  if(myid==0) printf("mcorrp or mcorra %d %d not properly initialized\n",i,j);
	  exit(EXIT_FAILURE);
	}
				
      }
    }
  }
		
  if(p12){
    for(i=0;i<16;i++){
      bcorr12[i]  = malloc(p12*sizeof(qcd_complex_16**));			
      if(bcorr12[i]==NULL){
	if(myid==0) printf("bcorr12 %d not properly initialized\n",i);
	exit(EXIT_FAILURE);
      }
			
      for(k=0;k<p12;k++){
	bcorr12[i][k] = malloc(tslices*sizeof(qcd_complex_16*));
	if(bcorr12[i][k]==NULL){
	  if(myid==0) printf("bcorr12 %d %d not properly initialized\n",i,k);
	  exit(EXIT_FAILURE);
	}				
				
	for(t=0;t<tslices;t++){
	  bcorr12[i][k][t] = malloc(nmom*sizeof(qcd_complex_16));
				
	  if(bcorr12[i][k][t]==NULL){
	    if(myid==0) printf("bcorr12 %d %d %d not properly initialized\n",i,k,t);
	    exit(EXIT_FAILURE);
	  }

	}
      }
    }
  }
	
  if(p32){
    for(i=0;i<16;i++){
      for(j=0;j<3;j++){
	bcorr32[i][j]  = malloc(p32*sizeof(qcd_complex_16**));
	if(bcorr32[i][j]==NULL){
	  if(myid==0) printf("bcorr32 %d %d not properly initialized\n",i,j);
	  exit(EXIT_FAILURE);
	}
			
	for(k=0;k<p32;k++){
	  bcorr32[i][j][k] = malloc(tslices*sizeof(qcd_complex_16*));
	  if(bcorr32[i][j][k]==NULL){
	    if(myid==0) printf("bcorr32 %d %d %d not properly initialized\n",i,j,k);
	    exit(EXIT_FAILURE);
	  }
	  for(t=0;t<tslices;t++){
	    bcorr32[i][j][k][t] = malloc(nmom*sizeof(qcd_complex_16));
	    if(bcorr32[i][j][k][t]==NULL){
	      if(myid==0) printf("bcorr32 %d %d %d %d not properly initialized\n",i,j,k,t);
	      exit(EXIT_FAILURE);
	    }
	  }
	}
      }
    }
  }

  if(myid==0) printf("\nmemory for all correlators allocated\n\n");
   
  //##############################################################################
  // load propagators

  if(myid==0) printf("Propagators loading...\n");  
	
  if (uyes){
    if(qcd_getPropagator(uprop_name,qcd_PROP_LIME, &uprop)) exit(EXIT_FAILURE);
    if(myid==0) printf("up propagator loaded\n");
  }
  if (dyes){
    if(qcd_getPropagator(dprop_name,qcd_PROP_LIME, &dprop)) exit(EXIT_FAILURE);
    if(myid==0) printf("down propagator loaded\n");  
  }
  if (syes){
    if(qcd_getPropagator(sprop_name,qcd_PROP_LIME, &sprop)) exit(EXIT_FAILURE);
    if(myid==0) printf("strange propagator loaded\n");
  }
  if (cyes){
    if(qcd_getPropagator(cprop_name,qcd_PROP_LIME, &cprop)) exit(EXIT_FAILURE);
    if(myid==0) printf("charm propagator loaded\n");   
  }   

  //################################################################################
   
  //transform propagators to the physical basis
   
  i_im.re = 0;
  i_im.im = 1;
 
  for(mu=0;mu<4;mu++){   
    one_plus_ig5[mu]  = qcd_CADD( qcd_ONE[mu][mu],qcd_CMUL(i_im,qcd_GAMMA[5][mu][mu]) );	   
    one_minus_ig5[mu] = qcd_CSUB( qcd_ONE[mu][mu],qcd_CMUL(i_im,qcd_GAMMA[5][mu][mu]) );
    g5[mu] = qcd_CSCALE(qcd_GAMMA[5][mu][mu],1.0);
  }

  if(!unit_sc){
    if(myid==0) printf("Transforming propagators according to mixed action setup\n");
    for(t=t_start; t<=t_stop; t++){
      lt = ((t+x_src[0])%geo.L[0]) - geo.Pos[0]*geo.lL[0];       
      for(lx=0; lx<geo.lL[1]; lx++)
	for(ly=0; ly<geo.lL[2]; ly++)
	  for(lz=0; lz<geo.lL[3]; lz++){
	    v = qcd_LEXIC(lt,lx,ly,lz,geo.lL);
        
	    for(c1=0;c1<3;c1++)
	      for(c2=0;c2<3;c2++)
		for(mu=0;mu<4;mu++)
		  for(nu=0;nu<4;nu++){
	     
		    if (uyes) uprop_pb.D[v][mu][nu][c1][c2] = qcd_CSCALE(qcd_CMUL(uprop.D[v][mu][nu][c1][c2],
										  qcd_CMUL(one_plus_ig5[mu],one_plus_ig5[nu])
										  ), 0.5);

		    if (dyes) dprop_pb.D[v][mu][nu][c1][c2] = qcd_CSCALE(qcd_CMUL(dprop.D[v][mu][nu][c1][c2],
										  qcd_CMUL(one_minus_ig5[mu],one_minus_ig5[nu])
										  ), 0.5);
	     

		    if (syes) sprop_pb.D[v][mu][nu][c1][c2] = qcd_CSCALE(qcd_CMUL(sprop.D[v][mu][nu][c1][c2],
										  qcd_CMUL(one_plus_ig5[mu],one_plus_ig5[nu])
										  ), 0.5);
	     
		    if (cyes) cprop_pb.D[v][mu][nu][c1][c2] = qcd_CSCALE(qcd_CMUL(cprop.D[v][mu][nu][c1][c2],
										  qcd_CMUL(one_plus_ig5[mu],one_plus_ig5[nu])
										  ), 0.5);	     	     
		  }//mu,nu,c2,c1	  
	  }//-space
    }//-time

  }
  else{
    if(myid==0) printf("Transforming propagators according to unitary setup\n");
    for(t=t_start; t<=t_stop; t++){
      lt = ((t+x_src[0])%geo.L[0]) - geo.Pos[0]*geo.lL[0];       
      for(lx=0; lx<geo.lL[1]; lx++)
	for(ly=0; ly<geo.lL[2]; ly++)
	  for(lz=0; lz<geo.lL[3]; lz++){
	    v = qcd_LEXIC(lt,lx,ly,lz,geo.lL);
        
	    for(c1=0;c1<3;c1++)
	      for(c2=0;c2<3;c2++)
		for(mu=0;mu<4;mu++)
		  for(nu=0;nu<4;nu++){
	     
		    if (uyes)uprop_pb.D[v][mu][nu][c1][c2] = qcd_CSCALE(qcd_CMUL(uprop.D[v][mu][nu][c1][c2],
										 qcd_CMUL(one_plus_ig5[mu],one_plus_ig5[nu])
										 ), 0.5);

		    if (dyes)dprop_pb.D[v][mu][nu][c1][c2] = qcd_CSCALE(qcd_CMUL(dprop.D[v][mu][nu][c1][c2],
										 qcd_CMUL(one_minus_ig5[mu],one_minus_ig5[nu])
										 ), 0.5);
	     

		    if (syes)sprop_pb.D[v][mu][nu][c1][c2] = qcd_CSCALE(qcd_CSUB(sprop.D[v][mu][nu][c1][c2],
										 qcd_CMUL(cprop.D[v][mu][nu][c1][c2],qcd_CMUL(qcd_GAMMA[5][mu][mu],qcd_GAMMA[5][nu][nu]))
										 ), 0.5);
	     
		    if (cyes)cprop_pb.D[v][mu][nu][c1][c2] = qcd_CSCALE(qcd_CSUB(cprop.D[v][mu][nu][c1][c2],
										 qcd_CMUL(sprop.D[v][mu][nu][c1][c2],qcd_CMUL(qcd_GAMMA[5][mu][mu],qcd_GAMMA[5][nu][nu]))
										 ), 0.5);
		  }//mu,nu,c2,c1	  
	  }//-space
    }//-time		
		
		
  }

  if (uyes) qcd_destroyPropagator(&uprop);
  if (dyes) qcd_destroyPropagator(&dprop);
  if (syes) qcd_destroyPropagator(&sprop);
  if (cyes) qcd_destroyPropagator(&cprop); 
   
  if(myid==0) printf("propagators transformed to physcial basis\n");  
   
  //################################################################################ 
  // transform propagators to basis with theta-periodic boundaries in the temporal direction
  for(lt=0; lt<geo.lL[0]; lt++)
    {
      t = lt + geo.Pos[0] * geo.lL[0];
      phase_factor   = (qcd_complex_16) {cos(theta[0]*t/geo.L[0]),sin(theta[0]*t/geo.L[0])};
      if (uyes) qcd_mulPropagatorC3d(&uprop_pb, phase_factor, (t+x_src[0]) % geo.L[0]);
      if (dyes) qcd_mulPropagatorC3d(&dprop_pb, phase_factor, (t+x_src[0]) % geo.L[0]);
      if (syes) qcd_mulPropagatorC3d(&sprop_pb, phase_factor, (t+x_src[0]) % geo.L[0]);
      if (cyes) qcd_mulPropagatorC3d(&cprop_pb, phase_factor, (t+x_src[0]) % geo.L[0]);    
    }
  if(myid==0) printf("propagators transformed to basis with theta-periodic boundary conditions\n");
   
   
  //################################################################################

  // gaussian smearing of propagators

  if(uyes && dyes && syes && cyes){
    if(myid==0) printf("Smearing the propagators...\n Will use optimized smearing for all propagators\n");
                
    if( qcd_gaussIteration3d_opt(&uprop_pb,&dprop_pb,&sprop_pb,&cprop_pb,&geo,&uAPE,nsmear,alpha,t_start,t_stop,x_src[0]) )
      {
	fprintf(stderr,"process %i: Error while smearing!\n",geo.myid);
	exit(EXIT_FAILURE);
      }

    qcd_destroyGaugeField(&uAPE);

    if(myid==0) printf("propagators smeared\n");
  }
  else{
    if(myid==0) printf("Smearing the propagators...\n Will use standard smearing function\n");
    j=0;
    j += qcd_initVector(&vec, &geo);
   
    MPI_Allreduce(&j, &k, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    if(k>0){
      if(myid==0) printf("not enough memory\n");
      exit(EXIT_FAILURE);
    }

    for(mu=0;mu<4;mu++)
      for(c1=0;c1<3;c1++){
		
	if(uyes){	
	  qcd_copyVectorPropagator(&vec,&uprop_pb,mu,c1);
	  for(i=0; i<nsmear; i++){
	    if(qcd_gaussIteration3dAll(&vec,&uAPE,alpha,i==0)){
	      fprintf(stderr,"process %i: Error while smearing!\n",geo.myid);
	      exit(EXIT_FAILURE);
	    }
	  }
	  qcd_copyPropagatorVector(&uprop_pb,&vec,mu,c1);
	}

	if(dyes){
	  qcd_copyVectorPropagator(&vec,&dprop_pb,mu,c1);
	  for(i=0; i<nsmear; i++){
	    if(qcd_gaussIteration3dAll(&vec,&uAPE,alpha,i==0)){
	      fprintf(stderr,"process %i: Error while smearing!\n",geo.myid);
	      exit(EXIT_FAILURE);
	    }
	  }
	  qcd_copyPropagatorVector(&dprop_pb,&vec,mu,c1);
	}
      
	if(syes){
	  qcd_copyVectorPropagator(&vec,&sprop_pb,mu,c1);
	  for(i=0; i<nsmear; i++){
	    if(qcd_gaussIteration3dAll(&vec,&uAPE,alpha,i==0)){
	      fprintf(stderr,"process %i: Error while smearing!\n",geo.myid);
	      exit(EXIT_FAILURE);
	    }
	  }
	  qcd_copyPropagatorVector(&sprop_pb,&vec,mu,c1);
	}
      
	if(cyes){
	  qcd_copyVectorPropagator(&vec,&cprop_pb,mu,c1);
	  for(i=0; i<nsmear; i++){
	    if(qcd_gaussIteration3dAll(&vec,&uAPE,alpha,i==0)){
	      fprintf(stderr,"process %i: Error while smearing!\n",geo.myid);
	      exit(EXIT_FAILURE);
	    }
	  }
	  qcd_copyPropagatorVector(&cprop_pb,&vec,mu,c1);  
	}
 	
	if(myid==0) printf("Indices mu=%d,c1=%d smeared\n",mu,c1);     
      }
	
    qcd_destroyGaugeField(&uAPE);
    qcd_destroyVector(&vec);   
    if(myid==0) printf("propagators smeared\n");
  }//-else
 

  //################################################################################   
  
  //-open output files to write in
	
  char corr_f_name[qcd_MAX_STRING_LENGTH],corr_m_name[qcd_MAX_STRING_LENGTH];
  FILE *fp_corr,*fp_corr_mes;
	
  k=0;
  if(myid==0){
    sprintf(corr_f_name,"%s_baryons.out",corr_p_name);
    fp_corr = fopen(corr_f_name,"w");		
    if(fp_corr==NULL){
      printf("failed to open %s for writing\n",corr_f_name);
      k=1;
    }
			
    if(mesyes){
      sprintf(corr_m_name,"%s_mesons.out",corr_p_name);
      fp_corr_mes = fopen(corr_m_name,"w");   
      if(fp_corr_mes==NULL){
	printf("failed to open %s for writing\n",corr_m_name);
	k=1;
      }
    }			
  }
	
  MPI_Bcast(&k,1,MPI_INT, 0, MPI_COMM_WORLD);
  if(k>=1) exit(EXIT_FAILURE);	
	
  if(myid==0) printf("Files opened properly\n");
	
  //################################################################################ 
  //################################################################################ 
  //################################################################################ 	

  //--------------------------- C O N T R A C T I O N S	

  //------------------------------------------------------------MESONS------------------------------
  if(mesyes){

    for(t=t_start; t<=t_stop; t++){
      lt = ((t+x_src[0])%geo.L[0]) - geo.Pos[0]*geo.lL[0];
      if(lt>=0 && lt<geo.lL[0]){  //inside the local lattice, otherwise nothing to calculate
      
	if(myid==0) printf("mesons t=%i\n",t);
	for(v3=0; v3<geo.lV3; v3++){
	  // set blocks to zero
	  udblock[v3] = (qcd_complex_16) {0,0};
	  dublock[v3] = (qcd_complex_16) {0,0};
	  sublock[v3] = (qcd_complex_16) {0,0};
	  sdblock[v3] = (qcd_complex_16) {0,0};
	  cublock[v3] = (qcd_complex_16) {0,0};
	  cdblock[v3] = (qcd_complex_16) {0,0};

	  a0ublock[v3]   = (qcd_complex_16) {0,0};
	  a0dblock[v3]   = (qcd_complex_16) {0,0};
	  rho1ublock[v3] = (qcd_complex_16) {0,0};
	  rho1dblock[v3] = (qcd_complex_16) {0,0};
	  rho2ublock[v3] = (qcd_complex_16) {0,0}; 
	  rho2dblock[v3] = (qcd_complex_16) {0,0};
	  rho3ublock[v3] = (qcd_complex_16) {0,0};
	  rho3dblock[v3] = (qcd_complex_16) {0,0};
	  a11ublock[v3]  = (qcd_complex_16) {0,0};
	  a11dblock[v3]  = (qcd_complex_16) {0,0};
	  a12ublock[v3]  = (qcd_complex_16) {0,0};
	  a12dblock[v3]  = (qcd_complex_16) {0,0};
	  a13ublock[v3]  = (qcd_complex_16) {0,0};
	  a13dblock[v3]  = (qcd_complex_16) {0,0};
	  Ks1ublock[v3]  = (qcd_complex_16) {0,0};
	  Ks1dblock[v3]  = (qcd_complex_16) {0,0};
	  Ks2ublock[v3]  = (qcd_complex_16) {0,0};
	  Ks2dblock[v3]  = (qcd_complex_16) {0,0};
	  Ks3ublock[v3]  = (qcd_complex_16) {0,0};
	  Ks3dblock[v3]  = (qcd_complex_16) {0,0};
	  Ds1ublock[v3]  = (qcd_complex_16) {0,0};
	  Ds1dblock[v3]  = (qcd_complex_16) {0,0};
	  Ds2ublock[v3]  = (qcd_complex_16) {0,0};
	  Ds2dblock[v3]  = (qcd_complex_16) {0,0};
	  Ds3ublock[v3]  = (qcd_complex_16) {0,0};
	  Ds3dblock[v3]  = (qcd_complex_16) {0,0};
	  K11ublock[v3]  = (qcd_complex_16) {0,0};
	  K11dblock[v3]  = (qcd_complex_16) {0,0};
	  K12ublock[v3]  = (qcd_complex_16) {0,0};
	  K12dblock[v3]  = (qcd_complex_16) {0,0};
	  K13ublock[v3]  = (qcd_complex_16) {0,0};
	  K13dblock[v3]  = (qcd_complex_16) {0,0};
	}

#pragma omp parallel for private(lz,ly,lx,v,mu,nu,c1,c2,ku,lu,c3p)
	for(v3=0; v3<lv3; v3++) {
	  lz = v3 % geo.lL[3];
	  ly = ((v3 - lz)/geo.lL[3]) % geo.lL[2];
	  lx = ((v3 - lz)/geo.lL[3] - ly) / geo.lL[2];
	  v =  qcd_LEXIC(lt,lx,ly,lz,geo.lL);
	
	  for(mu=0;mu<4;mu++){
	    for(nu=0;nu<4;nu++){
	      for(c1=0;c1<3;c1++){
		for(c2=0;c2<3;c2++){  
		  //pion
		  udblock[v3] = qcd_CSUB( udblock[v3],qcd_CMUL(uprop_pb.D[v][mu][nu][c1][c2],qcd_CONJ(uprop_pb.D[v][mu][nu][c1][c2])) );  // pi+ , u-dbar
		  dublock[v3] = qcd_CSUB( dublock[v3],qcd_CMUL(dprop_pb.D[v][mu][nu][c1][c2],qcd_CONJ(dprop_pb.D[v][mu][nu][c1][c2])) );  // pi- , d-ubar

		  //-kaon
		  sublock[v3] = qcd_CSUB( sublock[v3],qcd_CMUL(sprop_pb.D[v][mu][nu][c1][c2],qcd_CONJ(uprop_pb.D[v][mu][nu][c1][c2])) ); // K0 , s-dbar
		  sdblock[v3] = qcd_CSUB( sdblock[v3],qcd_CMUL(sprop_pb.D[v][mu][nu][c1][c2],qcd_CONJ(dprop_pb.D[v][mu][nu][c1][c2])) ); // K- , s-ubar

		  //-Dmeson
		  cublock[v3] = qcd_CSUB( cublock[v3],qcd_CMUL(cprop_pb.D[v][mu][nu][c1][c2],qcd_CONJ(uprop_pb.D[v][mu][nu][c1][c2])) ); // D+ , c-dbar
		  cdblock[v3] = qcd_CSUB( cdblock[v3],qcd_CMUL(cprop_pb.D[v][mu][nu][c1][c2],qcd_CONJ(dprop_pb.D[v][mu][nu][c1][c2])) ); // D0 , c-ubar

		  for(ku=0;ku<4;ku++){
		    for(lu=0;lu<4;lu++){
		      scalarGamma = qcd_CMUL(qcd_GAMMA[5][ku][mu],qcd_GAMMA[5][nu][lu]);
		      for(c3p=1;c3p<=3;c3p++){
			vectorGamma[c3p-1] = qcd_CMUL(qcd_G5GAMMA[c3p][ku][mu],qcd_G5GAMMA[c3p][nu][lu]);
			axVectorGamma[c3p-1]  = qcd_CMUL(qcd_GAMMA[c3p][ku][mu],qcd_GAMMA[c3p][nu][lu]);
		      }

		      //-a0 scalar meson
		      a0ublock[v3] = qcd_CSUB(a0ublock[v3], qcd_CMUL(scalarGamma,qcd_CMUL(uprop_pb.D[v][mu][nu][c1][c2],qcd_CONJ(uprop_pb.D[v][ku][lu][c1][c2]))));
		      a0dblock[v3] = qcd_CSUB(a0dblock[v3], qcd_CMUL(scalarGamma,qcd_CMUL(dprop_pb.D[v][mu][nu][c1][c2],qcd_CONJ(dprop_pb.D[v][ku][lu][c1][c2]))));

		      //-rho vector meson
		      rho1ublock[v3] = qcd_CADD(rho1ublock[v3], qcd_CMUL(vectorGamma[0],qcd_CMUL(uprop_pb.D[v][mu][nu][c1][c2],qcd_CONJ(uprop_pb.D[v][ku][lu][c1][c2]))));
		      rho1dblock[v3] = qcd_CADD(rho1dblock[v3], qcd_CMUL(vectorGamma[0],qcd_CMUL(dprop_pb.D[v][mu][nu][c1][c2],qcd_CONJ(dprop_pb.D[v][ku][lu][c1][c2]))));
		      rho2ublock[v3] = qcd_CADD(rho2ublock[v3], qcd_CMUL(vectorGamma[1],qcd_CMUL(uprop_pb.D[v][mu][nu][c1][c2],qcd_CONJ(uprop_pb.D[v][ku][lu][c1][c2]))));
		      rho2dblock[v3] = qcd_CADD(rho2dblock[v3], qcd_CMUL(vectorGamma[1],qcd_CMUL(dprop_pb.D[v][mu][nu][c1][c2],qcd_CONJ(dprop_pb.D[v][ku][lu][c1][c2]))));
		      rho3ublock[v3] = qcd_CADD(rho3ublock[v3], qcd_CMUL(vectorGamma[2],qcd_CMUL(uprop_pb.D[v][mu][nu][c1][c2],qcd_CONJ(uprop_pb.D[v][ku][lu][c1][c2]))));
		      rho3dblock[v3] = qcd_CADD(rho3dblock[v3], qcd_CMUL(vectorGamma[2],qcd_CMUL(dprop_pb.D[v][mu][nu][c1][c2],qcd_CONJ(dprop_pb.D[v][ku][lu][c1][c2]))));

		      //-a1 axial vector meson
		      a11ublock[v3] = qcd_CSUB(a11ublock[v3], qcd_CMUL(axVectorGamma[0],qcd_CMUL(uprop_pb.D[v][mu][nu][c1][c2],qcd_CONJ(uprop_pb.D[v][ku][lu][c1][c2]))));
		      a11dblock[v3] = qcd_CSUB(a11dblock[v3], qcd_CMUL(axVectorGamma[0],qcd_CMUL(dprop_pb.D[v][mu][nu][c1][c2],qcd_CONJ(dprop_pb.D[v][ku][lu][c1][c2]))));
		      a12ublock[v3] = qcd_CSUB(a12ublock[v3], qcd_CMUL(axVectorGamma[1],qcd_CMUL(uprop_pb.D[v][mu][nu][c1][c2],qcd_CONJ(uprop_pb.D[v][ku][lu][c1][c2]))));
		      a12dblock[v3] = qcd_CSUB(a12dblock[v3], qcd_CMUL(axVectorGamma[1],qcd_CMUL(dprop_pb.D[v][mu][nu][c1][c2],qcd_CONJ(dprop_pb.D[v][ku][lu][c1][c2]))));
		      a13ublock[v3] = qcd_CSUB(a13ublock[v3], qcd_CMUL(axVectorGamma[2],qcd_CMUL(uprop_pb.D[v][mu][nu][c1][c2],qcd_CONJ(uprop_pb.D[v][ku][lu][c1][c2]))));
		      a13dblock[v3] = qcd_CSUB(a13dblock[v3], qcd_CMUL(axVectorGamma[2],qcd_CMUL(dprop_pb.D[v][mu][nu][c1][c2],qcd_CONJ(dprop_pb.D[v][ku][lu][c1][c2]))));

		      //-Kstar vector meson
		      Ks1ublock[v3] = qcd_CADD(Ks1ublock[v3], qcd_CMUL(vectorGamma[0],qcd_CMUL(sprop_pb.D[v][ku][lu][c1][c2],qcd_CONJ(uprop_pb.D[v][mu][nu][c1][c2]))));
		      Ks1dblock[v3] = qcd_CADD(Ks1dblock[v3], qcd_CMUL(vectorGamma[0],qcd_CMUL(sprop_pb.D[v][ku][lu][c1][c2],qcd_CONJ(dprop_pb.D[v][mu][nu][c1][c2]))));
		      Ks2ublock[v3] = qcd_CADD(Ks2ublock[v3], qcd_CMUL(vectorGamma[1],qcd_CMUL(sprop_pb.D[v][ku][lu][c1][c2],qcd_CONJ(uprop_pb.D[v][mu][nu][c1][c2]))));
		      Ks2dblock[v3] = qcd_CADD(Ks2dblock[v3], qcd_CMUL(vectorGamma[1],qcd_CMUL(sprop_pb.D[v][ku][lu][c1][c2],qcd_CONJ(dprop_pb.D[v][mu][nu][c1][c2]))));
		      Ks3ublock[v3] = qcd_CADD(Ks3ublock[v3], qcd_CMUL(vectorGamma[2],qcd_CMUL(sprop_pb.D[v][ku][lu][c1][c2],qcd_CONJ(uprop_pb.D[v][mu][nu][c1][c2]))));
		      Ks3dblock[v3] = qcd_CADD(Ks3dblock[v3], qcd_CMUL(vectorGamma[2],qcd_CMUL(sprop_pb.D[v][ku][lu][c1][c2],qcd_CONJ(dprop_pb.D[v][mu][nu][c1][c2]))));

		      //-Dstar vector meson
		      Ds1ublock[v3] = qcd_CADD(Ds1ublock[v3], qcd_CMUL(vectorGamma[0],qcd_CMUL(cprop_pb.D[v][ku][lu][c1][c2],qcd_CONJ(uprop_pb.D[v][mu][nu][c1][c2]))));
		      Ds1dblock[v3] = qcd_CADD(Ds1dblock[v3], qcd_CMUL(vectorGamma[0],qcd_CMUL(cprop_pb.D[v][ku][lu][c1][c2],qcd_CONJ(dprop_pb.D[v][mu][nu][c1][c2]))));
		      Ds2ublock[v3] = qcd_CADD(Ds2ublock[v3], qcd_CMUL(vectorGamma[1],qcd_CMUL(cprop_pb.D[v][ku][lu][c1][c2],qcd_CONJ(uprop_pb.D[v][mu][nu][c1][c2]))));
		      Ds2dblock[v3] = qcd_CADD(Ds2dblock[v3], qcd_CMUL(vectorGamma[1],qcd_CMUL(cprop_pb.D[v][ku][lu][c1][c2],qcd_CONJ(dprop_pb.D[v][mu][nu][c1][c2]))));
		      Ds3ublock[v3] = qcd_CADD(Ds3ublock[v3], qcd_CMUL(vectorGamma[2],qcd_CMUL(cprop_pb.D[v][ku][lu][c1][c2],qcd_CONJ(uprop_pb.D[v][mu][nu][c1][c2]))));
		      Ds3dblock[v3] = qcd_CADD(Ds3dblock[v3], qcd_CMUL(vectorGamma[2],qcd_CMUL(cprop_pb.D[v][ku][lu][c1][c2],qcd_CONJ(dprop_pb.D[v][mu][nu][c1][c2]))));

		      //-K1 axial vector meson
		      K11ublock[v3] = qcd_CSUB(K11ublock[v3], qcd_CMUL(axVectorGamma[0],qcd_CMUL(sprop_pb.D[v][ku][lu][c1][c2],qcd_CONJ(uprop_pb.D[v][mu][nu][c1][c2]))));
		      K11dblock[v3] = qcd_CSUB(K11dblock[v3], qcd_CMUL(axVectorGamma[0],qcd_CMUL(sprop_pb.D[v][ku][lu][c1][c2],qcd_CONJ(dprop_pb.D[v][mu][nu][c1][c2]))));
		      K12ublock[v3] = qcd_CSUB(K12ublock[v3], qcd_CMUL(axVectorGamma[1],qcd_CMUL(sprop_pb.D[v][ku][lu][c1][c2],qcd_CONJ(uprop_pb.D[v][mu][nu][c1][c2]))));
		      K12dblock[v3] = qcd_CSUB(K12dblock[v3], qcd_CMUL(axVectorGamma[1],qcd_CMUL(sprop_pb.D[v][ku][lu][c1][c2],qcd_CONJ(dprop_pb.D[v][mu][nu][c1][c2]))));
		      K13ublock[v3] = qcd_CSUB(K13ublock[v3], qcd_CMUL(axVectorGamma[2],qcd_CMUL(sprop_pb.D[v][ku][lu][c1][c2],qcd_CONJ(uprop_pb.D[v][mu][nu][c1][c2]))));
		      K13dblock[v3] = qcd_CSUB(K13dblock[v3], qcd_CMUL(axVectorGamma[2],qcd_CMUL(sprop_pb.D[v][ku][lu][c1][c2],qcd_CONJ(dprop_pb.D[v][mu][nu][c1][c2]))));
		    }//-lu
		  }//-ku

		}}//-color c1,c2
	    }}//-dirac mu,nu
	}//-space       
 	 
	//Fourier transform time-slice
         
        for(j=0; j<nmom; j++){
         			 
	  for(i=0;i<38;i++) mcorrsum[i] = (qcd_complex_16) {0,0};
			           
	  for(lx=0; lx<geo.lL[1]; lx++)
            for(ly=0; ly<geo.lL[2]; ly++)
	      for(lz=0; lz<geo.lL[3]; lz++){  
		v3 = qcd_LEXIC0(lx,ly,lz,geo.lL);
		x=lx+geo.Pos[1]*geo.lL[1] - x_src[1];
		y=ly+geo.Pos[2]*geo.lL[2] - x_src[2];
		z=lz+geo.Pos[3]*geo.lL[3] - x_src[3];
		tmp = (((double) mom[j][0]*x)/geo.L[1] + ((double) mom[j][1]*y)/geo.L[2] + ((double) mom[j][2]*z)/geo.L[3])*2*M_PI;
		C2=(qcd_complex_16) {cos(tmp), -sin(tmp)}; //TABULATE FOR LARGE SPEEDUP!!!
		mcorrsum[0]  = qcd_CADD(mcorrsum[0],  qcd_CMUL(udblock[v3],C2)); // pi+
		mcorrsum[1]  = qcd_CADD(mcorrsum[1],  qcd_CMUL(dublock[v3],C2)); // pi-
		mcorrsum[2]  = qcd_CADD(mcorrsum[2],  qcd_CMUL(sdblock[v3],C2)); // K-
		mcorrsum[3]  = qcd_CADD(mcorrsum[3],  qcd_CMUL(sublock[v3],C2)); // K0
		mcorrsum[4]  = qcd_CADD(mcorrsum[4],  qcd_CMUL(cublock[v3],C2)); // D+
		mcorrsum[5]  = qcd_CADD(mcorrsum[5],  qcd_CMUL(cdblock[v3],C2)); // D0                                             
		mcorrsum[6]  = qcd_CADD(mcorrsum[6],  qcd_CMUL(a0ublock[v3],C2)); // a0+
		mcorrsum[7]  = qcd_CADD(mcorrsum[7],  qcd_CMUL(a0dblock[v3],C2)); // a0-
		mcorrsum[8]  = qcd_CADD(mcorrsum[8],  qcd_CMUL(rho1ublock[v3],C2)); // rho1+  
		mcorrsum[9]  = qcd_CADD(mcorrsum[9],  qcd_CMUL(rho1dblock[v3],C2)); // rho1- 
		mcorrsum[10] = qcd_CADD(mcorrsum[10], qcd_CMUL(rho2ublock[v3],C2)); // rho2+ 
		mcorrsum[11] = qcd_CADD(mcorrsum[11], qcd_CMUL(rho2dblock[v3],C2)); // rho2- 
		mcorrsum[12] = qcd_CADD(mcorrsum[12], qcd_CMUL(rho3ublock[v3],C2)); // rho3+ 
		mcorrsum[13] = qcd_CADD(mcorrsum[13], qcd_CMUL(rho3dblock[v3],C2)); // rho3- 
		mcorrsum[14] = qcd_CADD(mcorrsum[14], qcd_CMUL(a11ublock[v3],C2)); // a11+	 
		mcorrsum[15] = qcd_CADD(mcorrsum[15], qcd_CMUL(a11dblock[v3],C2)); // a11-	 
		mcorrsum[16] = qcd_CADD(mcorrsum[16], qcd_CMUL(a12ublock[v3],C2)); // a12+	 
		mcorrsum[17] = qcd_CADD(mcorrsum[17], qcd_CMUL(a12dblock[v3],C2)); // a12-	 
		mcorrsum[18] = qcd_CADD(mcorrsum[18], qcd_CMUL(a13ublock[v3],C2)); // a13+	 
		mcorrsum[19] = qcd_CADD(mcorrsum[19], qcd_CMUL(a13dblock[v3],C2)); // a13-   
		mcorrsum[20] = qcd_CADD(mcorrsum[20], qcd_CMUL(Ks1ublock[v3],C2)); // Kstar0_1
		mcorrsum[21] = qcd_CADD(mcorrsum[21], qcd_CMUL(Ks1dblock[v3],C2)); // Kstar-_1
		mcorrsum[22] = qcd_CADD(mcorrsum[22], qcd_CMUL(Ks2ublock[v3],C2)); // Kstar0_2
		mcorrsum[23] = qcd_CADD(mcorrsum[23], qcd_CMUL(Ks2dblock[v3],C2)); // Kstar-_2
		mcorrsum[24] = qcd_CADD(mcorrsum[24], qcd_CMUL(Ks3ublock[v3],C2)); // Kstar0_3
		mcorrsum[25] = qcd_CADD(mcorrsum[25], qcd_CMUL(Ks3dblock[v3],C2)); // Kstar-_3
		mcorrsum[26] = qcd_CADD(mcorrsum[26], qcd_CMUL(Ds1ublock[v3],C2)); // Dstar+_1
		mcorrsum[27] = qcd_CADD(mcorrsum[27], qcd_CMUL(Ds1dblock[v3],C2)); // Dstar0_1
		mcorrsum[28] = qcd_CADD(mcorrsum[28], qcd_CMUL(Ds2ublock[v3],C2)); // Dstar+_2
		mcorrsum[29] = qcd_CADD(mcorrsum[29], qcd_CMUL(Ds2dblock[v3],C2)); // Dstar0_2
		mcorrsum[30] = qcd_CADD(mcorrsum[30], qcd_CMUL(Ds3ublock[v3],C2)); // Dstar+_3
		mcorrsum[31] = qcd_CADD(mcorrsum[31], qcd_CMUL(Ds3dblock[v3],C2)); // Dstar0_3
		mcorrsum[32] = qcd_CADD(mcorrsum[32], qcd_CMUL(K11ublock[v3],C2)); // K1star+_1
		mcorrsum[33] = qcd_CADD(mcorrsum[33], qcd_CMUL(K11dblock[v3],C2)); // K1star0_1
		mcorrsum[34] = qcd_CADD(mcorrsum[34], qcd_CMUL(K12ublock[v3],C2)); // K1star+_2
		mcorrsum[35] = qcd_CADD(mcorrsum[35], qcd_CMUL(K12dblock[v3],C2)); // K1star0_2
		mcorrsum[36] = qcd_CADD(mcorrsum[36], qcd_CMUL(K13ublock[v3],C2)); // K1star+_3
		mcorrsum[37] = qcd_CADD(mcorrsum[37], qcd_CMUL(K13dblock[v3],C2)); // K1star0_3		      
	      }

	  MPI_Reduce(&(mcorrsum[0].re),  &(mcorrp[0][t][j].re),  2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	  MPI_Reduce(&(mcorrsum[1].re),  &(mcorra[0][t][j].re),  2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	  MPI_Reduce(&(mcorrsum[2].re),  &(mcorrp[1][t][j].re),  2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	  MPI_Reduce(&(mcorrsum[3].re),  &(mcorra[1][t][j].re),  2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	  MPI_Reduce(&(mcorrsum[4].re),  &(mcorrp[2][t][j].re),  2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	  MPI_Reduce(&(mcorrsum[5].re),  &(mcorra[2][t][j].re),  2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	  MPI_Reduce(&(mcorrsum[6].re),  &(mcorrp[3][t][j].re),  2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	  MPI_Reduce(&(mcorrsum[7].re),  &(mcorra[3][t][j].re),  2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	  MPI_Reduce(&(mcorrsum[8].re),  &(mcorrp[4][t][j].re),  2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	  MPI_Reduce(&(mcorrsum[9].re),  &(mcorra[4][t][j].re),  2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	  MPI_Reduce(&(mcorrsum[10].re), &(mcorrp[5][t][j].re),  2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	  MPI_Reduce(&(mcorrsum[11].re), &(mcorra[5][t][j].re),  2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	  MPI_Reduce(&(mcorrsum[12].re), &(mcorrp[6][t][j].re),  2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	  MPI_Reduce(&(mcorrsum[13].re), &(mcorra[6][t][j].re),  2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	  MPI_Reduce(&(mcorrsum[14].re), &(mcorrp[7][t][j].re),  2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	  MPI_Reduce(&(mcorrsum[15].re), &(mcorra[7][t][j].re),  2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	  MPI_Reduce(&(mcorrsum[16].re), &(mcorrp[8][t][j].re),  2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	  MPI_Reduce(&(mcorrsum[17].re), &(mcorra[8][t][j].re),  2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	  MPI_Reduce(&(mcorrsum[18].re), &(mcorrp[9][t][j].re),  2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	  MPI_Reduce(&(mcorrsum[19].re), &(mcorra[9][t][j].re),  2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	  MPI_Reduce(&(mcorrsum[20].re), &(mcorrp[10][t][j].re), 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	  MPI_Reduce(&(mcorrsum[21].re), &(mcorra[10][t][j].re), 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	  MPI_Reduce(&(mcorrsum[22].re), &(mcorrp[11][t][j].re), 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	  MPI_Reduce(&(mcorrsum[23].re), &(mcorra[11][t][j].re), 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	  MPI_Reduce(&(mcorrsum[24].re), &(mcorrp[12][t][j].re), 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	  MPI_Reduce(&(mcorrsum[25].re), &(mcorra[12][t][j].re), 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	  MPI_Reduce(&(mcorrsum[26].re), &(mcorrp[13][t][j].re), 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	  MPI_Reduce(&(mcorrsum[27].re), &(mcorra[13][t][j].re), 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	  MPI_Reduce(&(mcorrsum[28].re), &(mcorrp[14][t][j].re), 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	  MPI_Reduce(&(mcorrsum[29].re), &(mcorra[14][t][j].re), 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	  MPI_Reduce(&(mcorrsum[30].re), &(mcorrp[15][t][j].re), 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	  MPI_Reduce(&(mcorrsum[31].re), &(mcorra[15][t][j].re), 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	  MPI_Reduce(&(mcorrsum[32].re), &(mcorrp[16][t][j].re), 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	  MPI_Reduce(&(mcorrsum[33].re), &(mcorra[16][t][j].re), 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	  MPI_Reduce(&(mcorrsum[34].re), &(mcorrp[17][t][j].re), 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	  MPI_Reduce(&(mcorrsum[35].re), &(mcorra[17][t][j].re), 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	  MPI_Reduce(&(mcorrsum[36].re), &(mcorrp[18][t][j].re), 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	  MPI_Reduce(&(mcorrsum[37].re), &(mcorra[18][t][j].re), 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        }//-j
         
      }//end lt inside local block condition  
    }//end t-loop
  }//-mesons-if
	
  //################################################################################ 
  //################################################################################ 

  //-------------------------------------------------------------------- SPIN 1/2 BARYONS
    
  for(pc=0;pc<p12;pc++){
    p_id = p_arr12[pc]; 
                 
    for(t=t_start; t<=t_stop; t++){
      lt = ((t+x_src[0])%geo.L[0]) - geo.Pos[0]*geo.lL[0];
      if(lt>=0 && lt<geo.lL[0]){  //inside the local lattice, otherwise nothing to calculate
      
        for(v3=0; v3<geo.lV3; v3++)   //set blocks to zero
	  for(mu=0;mu<4;mu++)
	    for(nu=0;nu<4;nu++){
	      block12[mu][nu][v3]= (qcd_complex_16) {0,0};
	    }

	qcd_contractions2pt(p_id, block12, &uprop_pb, &dprop_pb, &sprop_pb, &cprop_pb, &geo, lt);
         
	//Fourier transform time-slice
        for(j=0; j<nmom; j++){
	  iarr = 0;
	  for(mu=0;mu<4;mu++){
	    for(nu=0;nu<4;nu++){		
	      bcorrsum12 = (qcd_complex_16) {0,0};

	      for(lx=0; lx<geo.lL[1]; lx++)
		for(ly=0; ly<geo.lL[2]; ly++)
		  for(lz=0; lz<geo.lL[3]; lz++){           
		    v3 = qcd_LEXIC0(lx,ly,lz,geo.lL);
		    x=lx+geo.Pos[1]*geo.lL[1] - x_src[1];
		    y=ly+geo.Pos[2]*geo.lL[2] - x_src[2];
		    z=lz+geo.Pos[3]*geo.lL[3] - x_src[3];
		    tmp = (((double) mom[j][0]*x)/geo.L[1] + ((double) mom[j][1]*y)/geo.L[2] + ((double) mom[j][2]*z)/geo.L[3])*2*M_PI;
		    C2=(qcd_complex_16) {cos(tmp), -sin(tmp)}; //TABULATE FOR LARGE SPEEDUP!!!
		    bcorrsum12 = qcd_CADD(bcorrsum12, qcd_CMUL(block12[mu][nu][v3],C2));
		  }
	      MPI_Reduce(&(bcorrsum12.re), &(bcorr12[iarr++][pc][t][j].re), 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	    }//-nu
	  }//-mu
	}//-j
             
      }//end lt inside local block condition
    }//end t-loop      
  }//-for pc
  
  //-------------------------------------------------------------------- SPIN 3/2 BARYONS   
  for(pc=0;pc<p32;pc++){
    p_id = p_arr32[pc]; 
    
             
    for(t=t_start; t<=t_stop; t++){
      lt = ((t+x_src[0])%geo.L[0]) - geo.Pos[0]*geo.lL[0];
      if(lt>=0 && lt<geo.lL[0]){  //inside the local lattice, otherwise nothing to calculate
		
	for(v3=0; v3<geo.lV3; v3++){   //set blocks to zero (again...)
	  for(mu=0;mu<4;mu++){
	    for(nu=0;nu<4;nu++){    
	      block12[mu][nu][v3]= (qcd_complex_16) {0,0};
	      block32[mu][nu][v3]= (qcd_complex_16) {0,0};
	      blocknp[mu][nu][v3]= (qcd_complex_16) {0,0};           
	    }}}

	qcd_contractions2pt_pr(p_id, block12, block32, blocknp, &uprop_pb, &dprop_pb, &sprop_pb, &cprop_pb, &geo, lt);
     
	//Fourier transform time-slice
        for(j=0; j<nmom; j++){
	  iarr = 0;
	  for(mu=0;mu<4;mu++){
	    for(nu=0;nu<4;nu++){			 
			 
	      for(i=0;i<3;i++) bcorrsum32[i] = (qcd_complex_16) {0,0};
                        
	      for(lx=0; lx<geo.lL[1]; lx++)
		for(ly=0; ly<geo.lL[2]; ly++)
		  for(lz=0; lz<geo.lL[3]; lz++){
		    v3 = qcd_LEXIC0(lx,ly,lz,geo.lL);
		    x=lx+geo.Pos[1]*geo.lL[1] - x_src[1];
		    y=ly+geo.Pos[2]*geo.lL[2] - x_src[2];
		    z=lz+geo.Pos[3]*geo.lL[3] - x_src[3];
		    tmp = (((double) mom[j][0]*x)/geo.L[1] + ((double) mom[j][1]*y)/geo.L[2] + ((double) mom[j][2]*z)/geo.L[3])*2*M_PI;
		    C2=(qcd_complex_16) {cos(tmp), -sin(tmp)}; //TABULATE FOR LARGE SPEEDUP!!!
		    bcorrsum32[0]=qcd_CADD(bcorrsum32[0], qcd_CMUL(blocknp[mu][nu][v3],C2));
		    bcorrsum32[1]=qcd_CADD(bcorrsum32[1], qcd_CMUL(block12[mu][nu][v3],C2));               
		    bcorrsum32[2]=qcd_CADD(bcorrsum32[2], qcd_CMUL(block32[mu][nu][v3],C2));
		  }
	      MPI_Reduce(&(bcorrsum32[0].re), &(bcorr32[iarr][0][pc][t][j].re), 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	      MPI_Reduce(&(bcorrsum32[1].re), &(bcorr32[iarr][1][pc][t][j].re), 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);            
	      MPI_Reduce(&(bcorrsum32[2].re), &(bcorr32[iarr][2][pc][t][j].re), 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	      iarr++;
	    }//-nu	
	  }//-mu
	}//-j
         
      }//end lt inside local block condition  
	
    }//end t-loop  
  }//-for pc   
    
  //#####################################################################    

  //-write up mesons
  if(myid==0) printf("\nWrite up...\n");

  if(mesyes){
    if(myid==0){
      for(i=0;i<19;i++){
	for(j=0;j<nmom;j++){
	  for(t=t_start;t<=t_stop;t++){
	    fprintf(fp_corr_mes,"%s %+i %+i %+i %i \t%+e %+e\t%+e %+e\n",meson_names[i+1],mom[j][0],mom[j][1],mom[j][2],t,
		    mcorrp[i][t][j].re,mcorrp[i][t][j].im, mcorra[i][t][j].re,mcorra[i][t][j].im);
	  }
	}
	fprintf(fp_corr_mes,"\n\n");
      }
      printf("Mesons write up finished\n");
    }
  }

  if(myid==0){
    //-write up spin 1/2 baryons		
    for(pc=0;pc<p12;pc++){
      p_id=p_arr12[pc];
      for(j=0;j<nmom;j++){
	for(t=t_start;t<=t_stop;t++){
	  for(mu=0;mu<4;mu++){
	    fprintf(fp_corr,"%s 0 %+i %+i %+i %i %i     \t%+e %+e\t%+e %+e \t%+e %+e\t%+e %+e\n",particle_names[p_id],mom[j][0],mom[j][1],mom[j][2],t,mu,	
		    bcorr12[4*mu][pc][t][j].re,bcorr12[4*mu][pc][t][j].im,
		    bcorr12[4*mu+1][pc][t][j].re,bcorr12[4*mu+1][pc][t][j].im,
		    bcorr12[4*mu+2][pc][t][j].re,bcorr12[4*mu+2][pc][t][j].im,
		    bcorr12[4*mu+3][pc][t][j].re,bcorr12[4*mu+3][pc][t][j].im);
	  }
	}
      }
      fprintf(fp_corr,"\n\n");
    }
		
    if(p12) printf("Spin-1/2 baryons write up finished\n");
				
    //-write up spin 3/2 baryons		
    for(pc=0;pc<p32;pc++){
      p_id=p_arr32[pc];
      for(k=0;k<3;k++){
	for(j=0;j<nmom;j++){
	  for(t=t_start;t<=t_stop;t++){
	    for(mu=0;mu<4;mu++){
	      fprintf(fp_corr,"%s %i %+i %+i %+i %i %i     \t%+e %+e\t%+e %+e \t%+e %+e\t%+e %+e\n",particle_names[p_id],k,mom[j][0],mom[j][1],mom[j][2],t,mu,	
		      bcorr32[4*mu][k][pc][t][j].re,bcorr32[4*mu][k][pc][t][j].im,
		      bcorr32[4*mu+1][k][pc][t][j].re,bcorr32[4*mu+1][k][pc][t][j].im,
		      bcorr32[4*mu+2][k][pc][t][j].re,bcorr32[4*mu+2][k][pc][t][j].im,
		      bcorr32[4*mu+3][k][pc][t][j].re,bcorr32[4*mu+3][k][pc][t][j].im);
	    }
	  }
	}
      }
      fprintf(fp_corr,"\n\n");
    }
		
    if(p32) printf("Spin-3/2 baryons write up finished\n");
		   
  }//-myid
      
  // #####################################################################   
  // clean up
  if(myid==0) printf("Cleaning up...\n");
   
  if(myid==0){
    fclose(fp_corr);
    fclose(fp_corr_mes);
  }
   
  for(mu=0;mu<4;mu++)
    for(nu=0;nu<4;nu++){
      free(block12[mu][nu]);
      free(block32[mu][nu]);
      free(blocknp[mu][nu]);  
    }
	
  free(udblock);
  free(dublock);
  free(sublock);
  free(sdblock);
  free(cublock);
  free(cdblock);

  free(a0ublock);
  free(a0dblock);
  free(rho1ublock);
  free(rho1dblock);
  free(rho2ublock); 
  free(rho2dblock);
  free(rho3ublock);
  free(rho3dblock);
  free(a11ublock);
  free(a11dblock);
  free(a12ublock);
  free(a12dblock);
  free(a13ublock);
  free(a13dblock);
 
  if(mesyes){
    for(i=0;i<19;i++){	
      for(j=0;j<tslices;j++){
	free(mcorrp[i][j]);
	free(mcorra[i][j]);		
      }				
    }	
    for(i=0;i<19;i++){
      free(mcorrp[i]);
      free(mcorra[i]);
    }
  }

  if(p12){
    for(i=0;i<16;i++){
      for(k=0;k<p12;k++){					
	for(t=0;t<tslices;t++){
	  free(bcorr12[i][k][t]);
	}}}
    for(i=0;i<16;i++){
      for(k=0;k<p12;k++){					
	free(bcorr12[i][k]);
      }}		
    for(i=0;i<16;i++){				
      free(bcorr12[i]);
    }		
  }

  if(p32){
    for(i=0;i<16;i++){
      for(j=0;j<3;j++){	
	for(k=0;k<p32;k++){					
	  for(t=0;t<tslices;t++){
	    free(bcorr32[i][j][k][t]);
	  }}}}
    for(i=0;i<16;i++){
      for(j=0;j<3;j++){	
	for(k=0;k<p32;k++){					
	  free(bcorr32[i][j][k]);
	}}}		
    for(i=0;i<16;i++){				
      for(j=0;j<3;j++){
	free(bcorr32[i][j]);
      }}		
  }
   
  free(mom);
  free(p_arr);
  if(p12) free(p_arr12);
  if(p32) free(p_arr32);
	
  qcd_destroyPropagator(&uprop_pb);
  qcd_destroyPropagator(&dprop_pb);
  qcd_destroyPropagator(&sprop_pb);
  qcd_destroyPropagator(&cprop_pb);  
  qcd_destroyGeometry(&geo);
	
	
  MPI_Finalize();
}//end main
