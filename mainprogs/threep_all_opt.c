/* threep_all.c
 *
 * reads forward propagators and sequential propagators
 * and creates three point functions for the 40 particles
 * Here the inversion through the current is used.
 *
 * Christos Kallidonis
 *
 * June 2012
 *
 ****************************************/
 
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
  qcd_uint_4 i,j,k, v,lx,ly,lz,ip1,im1,v3; 
  qcd_int_4 x,y,z;
  qcd_uint_2 ic1,ic2,ic3;                    //
  qcd_uint_4 x_src[4];                       // source and sink coordinates
  qcd_uint_4 t_start, t_stop, t,lt;
  qcd_real_8 tmp;                            // general purpuse
  qcd_uint_4 t_curr,ctype,t_slice,proj,t_int,t_try,t_first;
  qcd_int_4 pc,p_id,p_num,*p_arr,nmom,p12,p32,partno,relt,nthreads,np;
   
  qcd_uint_2 projlist[5] = {3,4,13,15,16};
   
  FILE *fp_momlist,*fp_parlist;
   
  int params_len;                            // needed to read inputfiles
  char *params;                              // needed to read inputfiles

  char gauge_name[qcd_MAX_STRING_LENGTH];      // name of gauge-configuration file
  char corr_p_name[qcd_MAX_STRING_LENGTH];     // name of output file proton 2pt function
   
  char param_name[qcd_MAX_STRING_LENGTH];      // name of parameter file  
  char momlist_name[qcd_MAX_STRING_LENGTH];    // name of momenta-list file
  char particle_list[qcd_MAX_STRING_LENGTH];
  char uprop_name[qcd_MAX_STRING_LENGTH];      
  char dprop_name[qcd_MAX_STRING_LENGTH];      
  char sprop_name[qcd_MAX_STRING_LENGTH];      // file names of forward quark propagators    
  char cprop_name[qcd_MAX_STRING_LENGTH];
   
  char seq_uprop_name[qcd_MAX_STRING_LENGTH];      
  char seq_dprop_name[qcd_MAX_STRING_LENGTH];      
  char seq_sprop_name[qcd_MAX_STRING_LENGTH];      // file names of sequential quark propagators    
  char seq_cprop_name[qcd_MAX_STRING_LENGTH];   
   
   
  qcd_geometry geo;                            // geometry structure
  qcd_propagator uprop,dprop,uprop_pb,dprop_pb;                        // forward propagators
  qcd_propagator sprop,cprop,sprop_pb,cprop_pb;                        

  qcd_propagator seq_uprop,seq_dprop,seq_uprop_pb,seq_dprop_pb;       // sequential propagators
  qcd_propagator seq_sprop,seq_cprop,seq_sprop_pb,seq_cprop_pb;       


  qcd_vector vec,vec_u,vec_d,vec_s,vec_c,vec_su,vec_sd,vec_ss,vec_sc;                              // needed when smearing
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
  qcd_real_8 plaq;
  qcd_int_4 ctr, ctr2;
  qcd_int_2 cg5cg5_ind[16*16][4];
  qcd_complex_16 cg5cg5_val[16*16];
  qcd_complex_16 one_plus_ig5[4],one_minus_ig5[4]; //-for transformation purposes
   
  qcd_complex_16 *block[5];                       // to store the block (3pt function before FT)

  qcd_complex_16 ***thrp[5][3], thrp_sum;


  qcd_int_4 (*mom)[3];                         // momenta-list

  int myid,numprocs, namelen;    
  char processor_name[MPI_MAX_PROCESSOR_NAME]; 
             
  qcd_complex_16 i_im ; // Imaginary i         
             
  //set up MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);         // num. of processes taking part in the calculation
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);             // each process gets its ID
  MPI_Get_processor_name(processor_name,&namelen); // 


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
      

                     	  
  sscanf(qcd_getParam("<source_pos_txyz>",params,params_len),"%d %d %d %d",&x_src[0],&x_src[1],&x_src[2],&x_src[3]);
  if(myid==0) printf("Got source coords: %d %d %d %d\n",x_src[0],x_src[1],x_src[2],x_src[3]);
 
  sscanf(qcd_getParam("<n_slices>",params,params_len),"%d",&t_slice);
  if(myid==0) printf("Got number of time slices for sink: %d\n", t_slice);	

  sscanf(qcd_getParam("<t_interval>",params,params_len),"%d",&t_int);
  if(myid==0) printf("Got time interval between source and current: %d\n", t_int);   
   
  sscanf(qcd_getParam("<t_first_sink>",params,params_len),"%d",&t_first);
  if(myid==0) printf("Got time slice for first sink point: %d\n", t_first);


  t_try = x_src[0] + t_int;
	
  t_curr = ( t_try >= L[0] ) ? (t_try-L[0]) : t_try ;
      
  if(myid==0) printf("Time slice for current is: %d\n", t_curr);

  t_start = t_int + t_first;
  t_stop = t_start + t_slice - 1;
     
  if(myid==0) printf("Sink times are: %d up to %d\n", t_start,t_stop);   
     
  strcpy(uprop_name,qcd_getParam("<propagator_u>",params,params_len));
  if(myid==0) printf("Got propagator file name: %s\n",uprop_name);
  strcpy(dprop_name,qcd_getParam("<propagator_d>",params,params_len));
  if(myid==0) printf("Got propagator file name: %s\n",dprop_name);
  strcpy(sprop_name,qcd_getParam("<propagator_s>",params,params_len));
  if(myid==0) printf("Got propagator file name: %s\n",sprop_name);   
  strcpy(cprop_name,qcd_getParam("<propagator_c>",params,params_len));
  if(myid==0) printf("Got propagator file name: %s\n",cprop_name);

  strcpy(seq_uprop_name,qcd_getParam("<seq_propagator_u>",params,params_len));
  if(myid==0) printf("Got sequential up propagator file name: %s\n",seq_uprop_name);
  strcpy(seq_dprop_name,qcd_getParam("<seq_propagator_d>",params,params_len));
  if(myid==0) printf("Got sequential down propagator file name: %s\n",seq_dprop_name);
  strcpy(seq_sprop_name,qcd_getParam("<seq_propagator_s>",params,params_len));
  if(myid==0) printf("Got sequential strange propagator file name: %s\n",seq_sprop_name);
  strcpy(seq_cprop_name,qcd_getParam("<seq_propagator_c>",params,params_len));
  if(myid==0) printf("Got sequential charm propagator file name: %s\n",seq_cprop_name);

   
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
   
  //  sscanf(qcd_getParam("<particles>",params,params_len),"%d %d",&p_ini,&p_fin);
  //  if(myid==0) printf(" Got particles: %d %d\n",p_ini,p_fin); 

  strcpy(particle_list,qcd_getParam("<particle_list>",params,params_len));
  if(myid==0) printf("Got particle-list file name: %s\n",particle_list);

   
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

  if(myid==0) printf("gauge-field APE smeared\n");
  plaq = qcd_calculatePlaquette(&uAPE);
  if(myid==0) printf("new plaquette = %e\n",plaq);  

  //##############################################################################
   
  //-allocate memory for propagators
  //-load propagators and transform them to physical basis

  i_im.re = 0;
  i_im.im = 1;
 
  for(mu=0;mu<4;mu++){   
    one_plus_ig5[mu]  = qcd_CADD( qcd_ONE[mu][mu],qcd_CMUL(i_im,qcd_GAMMA[5][mu][mu]) );	   
    one_minus_ig5[mu] = qcd_CSUB( qcd_ONE[mu][mu],qcd_CMUL(i_im,qcd_GAMMA[5][mu][mu]) );
  }

  //-forward props

  j = 0;
  j += qcd_initPropagator(&uprop, &geo);
  j += qcd_initPropagator(&dprop, &geo);
  j += qcd_initPropagator(&sprop, &geo);
  j += qcd_initPropagator(&cprop, &geo); 
      
  j += qcd_initPropagator(&uprop_pb, &geo);
  j += qcd_initPropagator(&dprop_pb, &geo);
  j += qcd_initPropagator(&sprop_pb, &geo);
  j += qcd_initPropagator(&cprop_pb, &geo);   
   

  MPI_Allreduce(&j, &k, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if(k>0)
    {
      if(myid==0) printf("not enough memory\n");
      exit(EXIT_FAILURE);
    }
  if(myid==0) printf("memory for forward propagators allocated\n");
    
   
  if(qcd_getPropagator(uprop_name,qcd_PROP_LIME, &uprop)) exit(EXIT_FAILURE);
  if(myid==0) printf("up propagator loaded\n");
  if(qcd_getPropagator(dprop_name,qcd_PROP_LIME, &dprop)) exit(EXIT_FAILURE);
  if(myid==0) printf("down propagator loaded\n");  
  if(qcd_getPropagator(sprop_name,qcd_PROP_LIME, &sprop)) exit(EXIT_FAILURE);
  if(myid==0) printf("strange propagator loaded\n");
  if(qcd_getPropagator(cprop_name,qcd_PROP_LIME, &cprop)) exit(EXIT_FAILURE);
  if(myid==0) printf("charm propagator loaded\n");   
   
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
	     
		  uprop_pb.D[v][mu][nu][c1][c2] = qcd_CSCALE(qcd_CMUL(uprop.D[v][mu][nu][c1][c2],
								      qcd_CMUL(one_plus_ig5[mu],one_plus_ig5[nu])
								      ), 0.5);

		  dprop_pb.D[v][mu][nu][c1][c2] = qcd_CSCALE(qcd_CMUL(dprop.D[v][mu][nu][c1][c2],
								      qcd_CMUL(one_minus_ig5[mu],one_minus_ig5[nu])
								      ), 0.5);
	     
		  sprop_pb.D[v][mu][nu][c1][c2] = qcd_CSCALE(qcd_CMUL(sprop.D[v][mu][nu][c1][c2],
								      qcd_CMUL(one_plus_ig5[mu],one_plus_ig5[nu])
								      ), 0.5);
	     
		  cprop_pb.D[v][mu][nu][c1][c2] = qcd_CSCALE(qcd_CMUL(cprop.D[v][mu][nu][c1][c2],
								      qcd_CMUL(one_plus_ig5[mu],one_plus_ig5[nu])
								      ), 0.5);
		}//mu,nu,c2,c1	  
	}//-space
  }//-time

  qcd_destroyPropagator(&uprop);
  qcd_destroyPropagator(&dprop);   
  qcd_destroyPropagator(&sprop);
  qcd_destroyPropagator(&cprop);

  if(myid==0) printf("forward propagators transformed to physical basis\n");

  //-sequential props

  j = 0;
  j += qcd_initPropagator(&seq_uprop, &geo);
  j += qcd_initPropagator(&seq_dprop, &geo);
  j += qcd_initPropagator(&seq_sprop, &geo);
  j += qcd_initPropagator(&seq_cprop, &geo);
   
  j += qcd_initPropagator(&seq_uprop_pb, &geo);
  j += qcd_initPropagator(&seq_dprop_pb, &geo);
  j += qcd_initPropagator(&seq_sprop_pb, &geo);
  j += qcd_initPropagator(&seq_cprop_pb, &geo);   

  MPI_Allreduce(&j, &k, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if(k>0)
    {
      if(myid==0) printf("not enough memory\n");
      exit(EXIT_FAILURE);
    }
  if(myid==0) printf("memory for sequential propagators allocated\n");



  if(qcd_getPropagator(seq_uprop_name,qcd_PROP_LIME, &seq_uprop)) exit(EXIT_FAILURE);
  if(myid==0) printf("sequential up propagator loaded\n");
  if(qcd_getPropagator(seq_dprop_name,qcd_PROP_LIME, &seq_dprop)) exit(EXIT_FAILURE);
  if(myid==0) printf("sequential down propagator loaded\n");
  if(qcd_getPropagator(seq_sprop_name,qcd_PROP_LIME, &seq_sprop)) exit(EXIT_FAILURE);
  if(myid==0) printf("sequential strange propagator loaded\n");
  if(qcd_getPropagator(seq_cprop_name,qcd_PROP_LIME, &seq_cprop)) exit(EXIT_FAILURE);
  if(myid==0) printf("sequential charm propagator loaded\n");
 
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

		  seq_uprop_pb.D[v][mu][nu][c1][c2] = qcd_CSCALE(qcd_CMUL(seq_uprop.D[v][mu][nu][c1][c2],
									  qcd_CMUL(one_plus_ig5[mu],one_plus_ig5[nu])
									  ), 0.5);

		  seq_dprop_pb.D[v][mu][nu][c1][c2] = qcd_CSCALE(qcd_CMUL(seq_dprop.D[v][mu][nu][c1][c2],
									  qcd_CMUL(one_minus_ig5[mu],one_minus_ig5[nu])
									  ), 0.5);
	     
		  seq_sprop_pb.D[v][mu][nu][c1][c2] = qcd_CSCALE(qcd_CMUL(seq_sprop.D[v][mu][nu][c1][c2],
									  qcd_CMUL(one_plus_ig5[mu],one_plus_ig5[nu])
									  ), 0.5);
	     
		  seq_cprop_pb.D[v][mu][nu][c1][c2] = qcd_CSCALE(qcd_CMUL(seq_cprop.D[v][mu][nu][c1][c2],
									  qcd_CMUL(one_plus_ig5[mu],one_plus_ig5[nu])
									  ), 0.5);
		}//mu,nu,c2,c1	  
	}//-space
  }//-time

  qcd_destroyPropagator(&seq_uprop);
  qcd_destroyPropagator(&seq_dprop);   
  qcd_destroyPropagator(&seq_sprop);
  qcd_destroyPropagator(&seq_cprop);

  if(myid==0) printf("sequential propagators transformed to physcial basis\n");

  //################################################################################ 
  // transform propagators to basis with theta-periodic boundaries in the temporal direction
  for(lt=0; lt<geo.lL[0]; lt++)
    {
      t = lt + geo.Pos[0] * geo.lL[0];
      phase_factor   = (qcd_complex_16) {cos(theta[0]*t/geo.L[0]),sin(theta[0]*t/geo.L[0])};
      qcd_mulPropagatorC3d(&uprop_pb, phase_factor, (t+x_src[0]) % geo.L[0]);
      qcd_mulPropagatorC3d(&dprop_pb, phase_factor, (t+x_src[0]) % geo.L[0]);
      qcd_mulPropagatorC3d(&sprop_pb, phase_factor, (t+x_src[0]) % geo.L[0]);
      qcd_mulPropagatorC3d(&cprop_pb, phase_factor, (t+x_src[0]) % geo.L[0]);
      
      qcd_mulPropagatorC3d(&seq_uprop_pb, phase_factor, (t+x_src[0]) % geo.L[0]);
      qcd_mulPropagatorC3d(&seq_dprop_pb, phase_factor, (t+x_src[0]) % geo.L[0]);
      qcd_mulPropagatorC3d(&seq_sprop_pb, phase_factor, (t+x_src[0]) % geo.L[0]);
      qcd_mulPropagatorC3d(&seq_cprop_pb, phase_factor, (t+x_src[0]) % geo.L[0]);                  
    }
  if(myid==0) printf("propagators transformed to basis with theta-periodic boundary conditions\n");   
  //################################################################################   
  //-Gaussian smearing new

  if(myid==0) printf("Smearing forward propagators...\n");
                
  if( qcd_gaussIteration3d_opt(&uprop_pb,&dprop_pb,&sprop_pb,&cprop_pb,&geo,&uAPE,nsmear,alpha,t_start,t_stop,x_src[0]) )
    {
      fprintf(stderr,"process %i: Error while smearing!\n",geo.myid);
      exit(EXIT_FAILURE);
    }

  if(myid==0) printf("forward propagators smeared\n");

  if(myid==0) printf("Smearing sequential propagators...\n");
                
  if( qcd_gaussIteration3d_opt(&seq_uprop_pb,&seq_dprop_pb,&seq_sprop_pb,&seq_cprop_pb,&geo,&uAPE,nsmear,alpha,t_start,t_stop,x_src[0]) )
    {
      fprintf(stderr,"process %i: Error while smearing!\n",geo.myid);
      exit(EXIT_FAILURE);
    }

  if(myid==0) printf("sequential propagators smeared\n");

  qcd_destroyGaugeField(&uAPE);

  //################################################################################

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

  //################################################################################

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


  //################################################################################

  //-allocate memory for the blocks
   
  for(i=0;i<5;i++){
    block[i] = (qcd_complex_16*) malloc(geo.lV3*sizeof(qcd_complex_16));
		
    if(block[i]==NULL){
      if(myid==0) printf("Block %d not properly initialized\n",i);
      exit(EXIT_FAILURE);
    }	
  }
  if(myid==0) printf("Blocks allocated properly\n");  
  
  //-allocate memory for the correlators
  for(i=0;i<5;i++){
    for(j=0;j<3;j++){
      thrp[i][j]  = malloc(p_num*sizeof(qcd_complex_16**));
      if(thrp[i][j]==NULL){
	if(myid==0) printf("thrp[%d][%d] not properly initialized\n",i,j);
	exit(EXIT_FAILURE);
      }
      
      for(k=0;k<p_num;k++){
	thrp[i][j][k] = malloc(t_slice*sizeof(qcd_complex_16*));
	if(thrp[i][j][k]==NULL){
	  if(myid==0) printf("thrp[%d][%d][%d] not properly initialized\n",i,j,k);
	  exit(EXIT_FAILURE);
	}
	for(t=0;t<t_slice;t++){
	  thrp[i][j][k][t] = malloc(nmom*sizeof(qcd_complex_16));
	  if(thrp[i][j][k][t]==NULL){
	    if(myid==0) printf("thrp[%d][%d][%d][%d] not properly initialized\n",i,j,k,t);
	    exit(EXIT_FAILURE);
	  }
	}
      }
    }
  }
  if(myid==0) printf("Correlator allocated properly\n");  

  //################################################################################
     
  //-open output files to write in
  char corr_f_name[qcd_MAX_STRING_LENGTH];
  FILE *fp_corr;

  k=0;
  if(myid==0){
    sprintf(corr_f_name,"%s.out",corr_p_name);
    
    fp_corr = fopen(corr_f_name,"w");   
    if(fp_corr==NULL){
      printf("failed to open %s for writing\n",corr_f_name);
      k=1;
    }				
  }//-myid
  MPI_Bcast(&k,1,MPI_INT, 0, MPI_COMM_WORLD);
  if(k==1) exit(EXIT_FAILURE);	

  if(myid==0) printf("File opened properly\n");

  //################################################################################
  //################################################################################
  //################################################################################                                                                                                                                                          

  //--------------------------- C O N T R A C T I O N S


  if(myid==0) printf("Performing the contractions...\n");


  for(pc=0;pc<p_num;pc++){
    p_id = p_arr[pc];
    
    partno=-1;
    for(np=0;np<4;np++){ // if the particle in interest has non zero part
      if(particles_pnum[p_id][np]){	
	partno++;
	relt=-1;	    
	for(t=t_start; t<=t_stop; t++){
	  lt = ((t+x_src[0])%geo.L[0]) - geo.Pos[0]*geo.lL[0];
	  if(lt>=0 && lt<geo.lL[0]){  //inside the local lattice, otherwise nothing to calculate
	    relt++;
	    for(v3=0; v3<geo.lV3; v3++)
	      for(proj=0;proj<5;proj++)  block[proj][v3] = (qcd_complex_16) {0,0};   //set blocks to zero

	    qcd_contractions3pt_new(p_id, np, block, &uprop_pb, &dprop_pb, &sprop_pb, &cprop_pb, 
				    &seq_uprop_pb, &seq_dprop_pb, &seq_sprop_pb, &seq_cprop_pb,&geo, lt);
			
	    //-Fourier transform time-slice for all projectors
	    for(proj=0;proj<5;proj++){ 
	      for(j=0; j<nmom; j++){					
		
		thrp_sum = (qcd_complex_16) {0,0};
      
		for(lx=0; lx<geo.lL[1]; lx++)
		  for(ly=0; ly<geo.lL[2]; ly++)
		    for(lz=0; lz<geo.lL[3]; lz++){
		      v3 = qcd_LEXIC0(lx,ly,lz,geo.lL);
		      x=lx+geo.Pos[1]*geo.lL[1] - x_src[1];
		      y=ly+geo.Pos[2]*geo.lL[2] - x_src[2];
		      z=lz+geo.Pos[3]*geo.lL[3] - x_src[3];
		      tmp = (((double) mom[j][0]*x)/geo.L[1] + ((double) mom[j][1]*y)/geo.L[2] + ((double) mom[j][2]*z)/geo.L[3])*2*M_PI;
		      C2=(qcd_complex_16) {cos(tmp), -sin(tmp)}; //TABULATE FOR LARGE SPEEDUP!!!

		      thrp_sum = qcd_CADD(thrp_sum, qcd_CMUL(block[proj][v3],C2));
		    }
		MPI_Reduce(&(thrp_sum.re), &(thrp[proj][partno][pc][relt][j].re), 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	      }//-mom	
	    }//-projector
		
	  }//end lt inside local block condition
	}//end t-loop
      }
    }//-part check	
  }//-particles
  
  if(myid==0) printf("Contractions finished successfully\n");
   
  //#####################################################################   
  // write up


  if(myid==0){
    printf("Write up\n");
 
    for(pc=0;pc<p_num;pc++){
      p_id=p_arr[pc];
      partno=-1;
      for(np=0;np<4;np++){ // if the particle in interest has non zero part
	if(particles_pnum[p_id][np]){	  
	  partno++;
	  for(j=0;j<nmom;j++){
	    for(proj=0;proj<5;proj++){
	      relt=-1;
	      for(t=t_start;t<=t_stop;t++){
		relt++;
		fprintf(fp_corr,"%s %s %+i %+i %+i %i %i %i     \t%+e  %+e\n",particle_names[p_id],particles_parts[p_id][np],mom[j][0],mom[j][1],mom[j][2],projlist[proj],t,relt+1,
			thrp[proj][partno][pc][relt][j].re,thrp[proj][partno][pc][relt][j].im);
	      }//-t
	    }//-proj
	  }//-mom
	  fprintf(fp_corr,"\n");
	}//-if
      }//-np
      fprintf(fp_corr,"\n\n");
    }//-pc
    printf("Write up finished\n");
  }//-myid

  //#####################################################################   
  // clean up
  if(myid==0) printf("cleaning up...\n");

  if(myid==0){
    fclose(fp_corr);
  }    
   
   
  for(proj=0;proj<5;proj++) free(block[proj]);
  //   free(mom);
  
  for(i=0;i<5;i++){
    for(j=0;j<3;j++){
      for(k=0;k<p_num;k++){
	for(t=0;t<t_slice;t++){
	  free(thrp[i][j][k][t]);
	}}}
  }
  for(i=0;i<5;i++){
    for(j=0;j<3;j++){
      for(k=0;k<p_num;k++){
	free(thrp[i][j][k]);
      }}
  }
  for(i=0;i<5;i++){
    for(j=0;j<3;j++){
      free(thrp[i][j]);
    }
  }

  free(p_arr);

 
  qcd_destroyPropagator(&uprop_pb);
  qcd_destroyPropagator(&dprop_pb);   
  qcd_destroyPropagator(&sprop_pb);
  qcd_destroyPropagator(&cprop_pb);  

  qcd_destroyPropagator(&seq_uprop_pb);
  qcd_destroyPropagator(&seq_dprop_pb);   
  qcd_destroyPropagator(&seq_sprop_pb);
  qcd_destroyPropagator(&seq_cprop_pb);

  if(myid==0) printf("Done.\n");
   
  qcd_destroyGeometry(&geo);
  MPI_Finalize();
}//end main
