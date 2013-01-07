/* twop_routine.c
 *
 * reads forward propagators
 * and creates two point functions for the 40 particles
 *
 * Christos Kallidonis
 *
 * April 2012
 *
 ****************************************/
 
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <qcd.h>

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
   qcd_uint_4 p_id,p_ini,p_fin;
   
   FILE *fp_momlist;
   
   int params_len;                            // needed to read inputfiles
   char *params;                              // needed to read inputfiles

   char gauge_name[qcd_MAX_STRING_LENGTH];      // name of gauge-configuration file
   char corr_p_name[qcd_MAX_STRING_LENGTH];     // name of output file proton 2pt function
   
   char param_name[qcd_MAX_STRING_LENGTH];      // name of parameter file  
   char momlist_name[qcd_MAX_STRING_LENGTH];    // name of momenta-list file
   char uprop_name[qcd_MAX_STRING_LENGTH];      // file names of up and down quark propagators
   char dprop_name[qcd_MAX_STRING_LENGTH];      
   char sprop_name[qcd_MAX_STRING_LENGTH];      // file names of up and down quark propagators
   char cprop_name[qcd_MAX_STRING_LENGTH];
   
   
   qcd_geometry geo;                            // geometry structure
   qcd_propagator uprop,sprop,uprop_pb,dprop_pb;                        // propagator
   qcd_propagator dprop,cprop,sprop_pb,cprop_pb;                        // propagator
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
   qcd_complex_16 corr, corr2,corr12,corr32,corr12_2,corr32_2;
   qcd_real_8 plaq;
   qcd_int_4 ctr, ctr2;
   qcd_int_2 cg5cg5_ind[16*16][4];
   qcd_complex_16 cg5cg5_val[16*16];
   qcd_complex_16 one_plus_ig5[4],one_minus_ig5[4]; //-for transformation purposes
   
   qcd_complex_16 *block12[4][4],*block32[4][4];                       // to store the block (2pt function before FT)

   qcd_int_4 (*mom)[3];                         // momenta-list

   int myid,numprocs, namelen;    
   char processor_name[MPI_MAX_PROCESSOR_NAME];
   				 
				 
             
    qcd_complex_16 i_im ; // Imaginary i         
             
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
   
   sscanf(qcd_getParam("<particles>",params,params_len),"%d %d",&p_ini,&p_fin);
   if(myid==0) printf(" Got particles: %d %d\n",p_ini,p_fin);
   
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

   j = 0;
   j += qcd_initPropagator(&uprop, &geo);
   j += qcd_initPropagator(&dprop, &geo);
   j += qcd_initPropagator(&sprop, &geo);
   j += qcd_initPropagator(&cprop, &geo);
   
   j += qcd_initPropagator(&uprop_pb, &geo);
   j += qcd_initPropagator(&dprop_pb, &geo);
   j += qcd_initPropagator(&sprop_pb, &geo);
   j += qcd_initPropagator(&cprop_pb, &geo);
   j += qcd_initVector(&vec, &geo);
   
   
   for(mu=0;mu<4;mu++){
	   for(nu=0;nu<4;nu++){
			block12[mu][nu] = (qcd_complex_16*) malloc(geo.lV3*sizeof(qcd_complex_16));
			block32[mu][nu] = (qcd_complex_16*) malloc(geo.lV3*sizeof(qcd_complex_16)); 
			
			if(block12[mu][nu]==NULL){
				if(myid==0) printf("Block12 %d %d not properly initialized\n",mu,nu);
				exit(EXIT_FAILURE);
			}
			if(block32[mu][nu]==NULL){
				if(myid==0) printf("Block32 %d %d not properly initialized\n",mu,nu);
				exit(EXIT_FAILURE);
			}
		}
	}
  
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
   if(qcd_getPropagator(sprop_name,qcd_PROP_LIME, &sprop)) exit(EXIT_FAILURE);
   if(myid==0) printf("strange propagator loaded\n");
   if(qcd_getPropagator(cprop_name,qcd_PROP_LIME, &cprop)) exit(EXIT_FAILURE);
   if(myid==0) printf("charm propagator loaded\n");   

   //################################################################################
   
   //transform propagators to the physical basis
   
   i_im.re = 0;
   i_im.im = 1;
 
    for(mu=0;mu<4;mu++){   
	   one_plus_ig5[mu]  = qcd_CADD( qcd_ONE[mu][mu],qcd_CMUL(i_im,qcd_GAMMA[5][mu][mu]) );	   
	   one_minus_ig5[mu] = qcd_CSUB( qcd_ONE[mu][mu],qcd_CMUL(i_im,qcd_GAMMA[5][mu][mu]) );
    }
   
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
   
    if(myid==0) printf("propagators transformed to physcial basis\n");  
   
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
   }
   if(myid==0) printf("propagators transformed to basis with theta-periodic boundary conditions\n");
   
   
   //################################################################################
   // gaussian smearing of propagators (only the time-slices that will be used)
   for(mu=0;mu<4;mu++)
   for(c1=0;c1<3;c1++)
   {
      qcd_copyVectorPropagator(&vec,&uprop_pb,mu,c1);
      for(i=0; i<nsmear; i++)
      {
         if(qcd_gaussIteration3dAll(&vec,&uAPE,alpha,i==0))
         {
            fprintf(stderr,"process %i: Error while smearing!\n",geo.myid);
            exit(EXIT_FAILURE);
         }
      }
      qcd_copyPropagatorVector(&uprop_pb,&vec,mu,c1);
      
      qcd_copyVectorPropagator(&vec,&dprop_pb,mu,c1);
      for(i=0; i<nsmear; i++)
      {
         if(qcd_gaussIteration3dAll(&vec,&uAPE,alpha,i==0))
         {
            fprintf(stderr,"process %i: Error while smearing!\n",geo.myid);
            exit(EXIT_FAILURE);
         }
      }
      qcd_copyPropagatorVector(&dprop_pb,&vec,mu,c1);
      
      qcd_copyVectorPropagator(&vec,&sprop_pb,mu,c1);
      for(i=0; i<nsmear; i++)
      {
         if(qcd_gaussIteration3dAll(&vec,&uAPE,alpha,i==0))
         {
            fprintf(stderr,"process %i: Error while smearing!\n",geo.myid);
            exit(EXIT_FAILURE);
         }
      }
      qcd_copyPropagatorVector(&sprop_pb,&vec,mu,c1);
      
      qcd_copyVectorPropagator(&vec,&cprop_pb,mu,c1);
      for(i=0; i<nsmear; i++)
      {
         if(qcd_gaussIteration3dAll(&vec,&uAPE,alpha,i==0))
         {
            fprintf(stderr,"process %i: Error while smearing!\n",geo.myid);
            exit(EXIT_FAILURE);
         }
      }
      qcd_copyPropagatorVector(&cprop_pb,&vec,mu,c1);      
   }
   qcd_destroyGaugeField(&uAPE);
   qcd_destroyVector(&vec);   
   if(myid==0) printf("propagators smeared\n");
   
   //################################################################################   
     
//-open output files to write in
	qcd_uint_4 partno = p_fin-p_ini+1,counter,counter32,partno32,partnorest;
	char corr_f_name[qcd_MAX_STRING_LENGTH],corr_f_name32[qcd_MAX_STRING_LENGTH],corr_f_name12[qcd_MAX_STRING_LENGTH];
	FILE *fp_corr[24],*fp_corr32[28],*fp_corr12[28];

	if(myid==0){
		counter=0;
		counter32=0; partno32 = 0;partnorest = 0;
		for(p_id=p_ini;p_id<=p_fin;p_id++){
			if(!particles32[p_id]){			
				sprintf(corr_f_name,"%s_%d_%s_nop.out",corr_p_name,p_id,particle_names[p_id]);
				fp_corr[counter] = fopen(corr_f_name,"w");   
				if(fp_corr[counter]==NULL){
					printf("failed to open %s for writing\n",corr_f_name);
					k=1;
				}
				counter++;
				partnorest++;
			}
			
			if(particles32[p_id]){
				sprintf(corr_f_name32,"%s_%d_%s_pr32_nop.out",corr_p_name,p_id,particle_names[p_id]);
				sprintf(corr_f_name12,"%s_%d_%s_pr12_nop.out",corr_p_name,p_id,particle_names[p_id]);				
				fp_corr32[counter32] = fopen(corr_f_name32,"w");
				fp_corr12[counter32] = fopen(corr_f_name12,"w");				   
				if(fp_corr32[counter32]==NULL){
					printf("failed to open %s for writing\n",corr_f_name32);
					k=1;
				}
				if(fp_corr12[counter32]==NULL){
					printf("failed to open %s for writing\n",corr_f_name12);
					k=1;
				}							
				counter32++;
				partno32++;
			}
	
		}
	}
	MPI_Bcast(&k,1,MPI_INT, 0, MPI_COMM_WORLD);
	if(k>=1) exit(EXIT_FAILURE);	


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


    
   counter=0; 
   for(p_id=p_ini;p_id<=p_fin;p_id++){ 
            
  	if(!particles32[p_id]){          
             
   for(t=t_start; t<=t_stop; t++){
      lt = ((t+x_src[0])%geo.L[0]) - geo.Pos[0]*geo.lL[0];
      if(lt>=0 && lt<geo.lL[0])  //inside the local lattice, otherwise nothing to calculate
      {
 
         //if(myid==0) printf("t=%i\n",t);
         for(v3=0; v3<geo.lV3; v3++)   //set blocks to zero
         for(mu=0;mu<4;mu++)
         for(nu=0;nu<4;nu++){
	   block12[mu][nu][v3]= (qcd_complex_16) {0,0};
	 }

	 qcd_contractions2pt(p_id, block12, &uprop_pb, &dprop_pb, &sprop_pb, &cprop_pb, &geo, lt);
         
//Fourier transform time-slice

         for(j=0; j<i; j++)
         {
	   
	   for(mu=0;mu<4;mu++){
	     if(myid==0)
	       {
		 fprintf(fp_corr[counter],"%i %+i %+i %+i %d ",t,mom[j][0],mom[j][1],mom[j][2],mu);
	       }
	     
	     for(nu=0;nu<4;nu++){		
	       
	       
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
               corr=qcd_CADD(corr, qcd_CMUL(block12[mu][nu][v3],C2));
            }
            //printf("process %i: corr = %f %+fi\n",myid,0.5*corr.re,0.5*corr.im);
            MPI_Reduce(&(corr.re), &(corr2.re), 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            if(myid==0) 
            {
               fprintf(fp_corr[counter],"\t%+e %+e",corr2.re,corr2.im);
            }
		}//-nu
            if(myid==0) 
            {
               fprintf(fp_corr[counter],"\n");
            }		
		}//-mu
       }//-j
             
      
      }//end lt inside local block condition
	  }//end t-loop      
		counter++;
	}//-particles if
	}//-particles
   
   
   counter32 = 0;  
   for(p_id=p_ini;p_id<=p_fin;p_id++){ 
             
 	if(particles32[p_id]){            
             
   for(t=t_start; t<=t_stop; t++){
      lt = ((t+x_src[0])%geo.L[0]) - geo.Pos[0]*geo.lL[0];
      if(lt>=0 && lt<geo.lL[0])  //inside the local lattice, otherwise nothing to calculate
      {   
   
        for(v3=0; v3<geo.lV3; v3++){   //set blocks to zero (again...)
			for(mu=0;mu<4;mu++)
			for(nu=0;nu<4;nu++){    
				block12[mu][nu][v3]= (qcd_complex_16) {0,0};
				block32[mu][nu][v3]= (qcd_complex_16) {0,0};           
			}
		}
			qcd_contractions2pt_pr(p_id, block12, block32, &uprop_pb, &dprop_pb, &sprop_pb, &cprop_pb, &geo, lt);
     
			//Fourier transform time-slice
         
         for(j=0; j<i; j++)
         {
			 
			for(mu=0;mu<4;mu++){
			  if(myid==0)
			    {
			      fprintf(fp_corr32[counter32],"%i %+i %+i %+i %d ",t,mom[j][0],mom[j][1],mom[j][2],mu);
			      fprintf(fp_corr12[counter32],"%i %+i %+i %+i %d ",t,mom[j][0],mom[j][1],mom[j][2],mu);
			    }
			for(nu=0;nu<4;nu++){			 
			 
            corr12 = (qcd_complex_16) {0,0};
            corr32 = (qcd_complex_16) {0,0};
                        
            for(lx=0; lx<geo.lL[1]; lx++)
            for(ly=0; ly<geo.lL[2]; ly++)
            for(lz=0; lz<geo.lL[3]; lz++)
            {
               v3 = qcd_LEXIC0(lx,ly,lz,geo.lL);
               x=lx+geo.Pos[1]*geo.lL[1] - x_src[1];
               y=ly+geo.Pos[2]*geo.lL[2] - x_src[2];
               z=lz+geo.Pos[3]*geo.lL[3] - x_src[3];
 //              if(myid==2) printf("%d %d %e\n",lt,v3,block[v3].re);
               tmp = (((double) mom[j][0]*x)/geo.L[1] + ((double) mom[j][1]*y)/geo.L[2] + ((double) mom[j][2]*z)/geo.L[3])*2*M_PI;
               C2=(qcd_complex_16) {cos(tmp), -sin(tmp)}; //TABULATE FOR LARGE SPEEDUP!!!
               corr12=qcd_CADD(corr12, qcd_CMUL(block12[mu][nu][v3],C2));
               corr32=qcd_CADD(corr32, qcd_CMUL(block32[mu][nu][v3],C2));               
            }
            //printf("process %i: corr = %f %+fi\n",myid,0.5*corr.re,0.5*corr.im);
            MPI_Reduce(&(corr12.re), &(corr12_2.re), 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce(&(corr32.re), &(corr32_2.re), 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);            
            if(myid==0) 
            {
               fprintf(fp_corr12[counter32],"\t%+e %+e",corr12_2.re,corr12_2.im);
               fprintf(fp_corr32[counter32],"\t%+e %+e",corr32_2.re,corr32_2.im);               
            }
		}//-nu
            if(myid==0) 
            {
               fprintf(fp_corr12[counter32],"\n");
               fprintf(fp_corr32[counter32],"\n");               
            }		
		}//-mu
         }//-j
      }//end lt inside local block condition  
	 }//end t-loop  
			counter32++;
	 }//-particles 32 if            
		
	}//-particles   
    
   if(myid==0)
   {
      for(i=0;i<partnorest;i++){
		  fclose(fp_corr[i]);
	  }
      for(i=0;i<partno32;i++){
		  fclose(fp_corr32[i]);
		  fclose(fp_corr12[i]);		  
	  }	  
   }   
   
                  
   
   //#####################################################################   
   // clean up
   if(myid==0) printf("cleaning up...\n");
   
   for(mu=0;mu<4;mu++)
   for(nu=0;nu<4;nu++){
		free(block12[mu][nu]);
		free(block32[mu][nu]);  
   }
   
   free(mom);
   qcd_destroyPropagator(&uprop_pb);
   qcd_destroyPropagator(&dprop_pb);
   qcd_destroyPropagator(&sprop_pb);
   qcd_destroyPropagator(&cprop_pb);  
   qcd_destroyGeometry(&geo);
   MPI_Finalize();
}//end main
