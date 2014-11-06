/* Christos Kallidonis
 * February 2012
 * 
 * This program reads forward propagators and creates Sequential Sources using fixed current method.
 *
 * 
 *** Sample Input file ***
  
 <processors_txyz>1 2 4 4</processors_txyz>
 <lattice_txyz>24 48 48 48</lattice_txyz>
 <source_pos_txyz>0 0 0 0</source_pos_txyz>
 <uprop_list>prop_u.txt</uprop_list>
 <dprop_list>prop_d.txt</dprop_list>
 <sprop_list>prop_s.txt</sprop_list>
 <usource_list>source_u</usource_list>
 <dsource_list>source_d</dsource_list>
 <ssource_list>source_s</ssource_list>
 <nsources>12</nsources>
 <t_current>5</t_current>
 <current_type>1</current_type>
 <momentum_xyz>0 0 0</momentum_xyz>
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <math.h>
#include <qcd.h>
 
int main(int argc,char* argv[])
{  
  qcd_uint_4 i,mu,nu,ku,lu,c1,c2,v,x,y,z;   // loop variables
  qcd_uint_4 c3,c4,ctr,ctr2,v3,lx,ly,lz,t;
  qcd_uint_4 cc1,cc2,a,b,j,gu,ju,k;
  qcd_uint_4 isource;                 // ..

  qcd_uint_4 nsources=12,ctype;                // number of different sources

  qcd_uint_4 nsmearAPE;       // APE smearing: n
  qcd_real_8 alphaAPE;         // APE smearing: alpha
  int params_len;                     // needed to read inputfiles
  char *params;                       // needed to read inputfiles
  char tmp_string[qcd_MAX_STRING_LENGTH]; // general purpuse
  char param_name[qcd_MAX_STRING_LENGTH];
  double tmp;                         // general purpuse

  char gauge_name[qcd_MAX_STRING_LENGTH]; // name of gauge-config file
  char momlist_name[qcd_MAX_STRING_LENGTH];    // name of momenta-list file
      
  qcd_uint_2 L[4], P[4],t_current,lt,t_int,t_try;     // lattice size and subdivision of lattice, time position of current insertion
  qcd_uint_4 x_src[4];					   // original source position
  qcd_int_4  mom[3];                // momentum
  qcd_geometry geo;                   // geometry structure   
  qcd_propagator uprop;              // u-propagator
  qcd_propagator dprop;              // d-propagator
  qcd_propagator sprop;			   // s-propagator	 
  qcd_propagator cprop;                           // c-propagator
  qcd_propagator prop_tmp;            // needed when rotating etc.
  qcd_vector usource;          // the seq_source for up
  qcd_vector dsource;          // the seq_source for down
  qcd_vector ssource;          // the seq_source for strange
  qcd_vector csource;          // the seq_source for charm   
  qcd_vector tmpvec; 			// general purposes
  qcd_gaugeField u;                   // gauge field
  qcd_gaugeField uAPE;                // APE smeared gaugeField
  qcd_gaugeField *u_ptr, *uAPE_ptr, *utmp_ptr;
   
  char uprop_name[qcd_MAX_STRING_LENGTH]; // file names of up and down and strange quark propagators
  char dprop_name[qcd_MAX_STRING_LENGTH];
  char sprop_name[qcd_MAX_STRING_LENGTH];
  char cprop_name[qcd_MAX_STRING_LENGTH];   

  char uout_name[qcd_MAX_STRING_LENGTH];     
  char dout_name[qcd_MAX_STRING_LENGTH];
  char sout_name[qcd_MAX_STRING_LENGTH]; 
  char cout_name[qcd_MAX_STRING_LENGTH];   

  qcd_real_8 theta[4] = {M_PI, 0.0, 0.0, 0.0}; // antiperiodic b.c. in time
   
  qcd_complex_16 z1, z2;               // temp variables
  qcd_complex_16 phase_factor;
  qcd_complex_16 C, factor;          
  qcd_real_8 plaq;
  qcd_complex_16 C2;
   
  int myid,numprocs, namelen;    
  char processor_name[MPI_MAX_PROCESSOR_NAME];
   
  FILE *fp_momlist, *fpusource,*fpdsource,*fpssource,*fpcsource;

  /*-----------------------------------------------------------------------------------------------------*/

  //set up MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);         // num. of processes taking part in the calculation
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);             // each process gets its ID
  MPI_Get_processor_name(processor_name,&namelen); // 
      
  //------------------------------ READ INPUT FILE ------------------------------//  
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
  
  sscanf(qcd_getParam("<source_pos_txyz>",params,params_len),"%d %d %d %d",&x_src[0],&x_src[1],&x_src[2],&x_src[3]);
  if(myid==0) printf("Got source coords: %d %d %d %d\n",x_src[0],x_src[1],x_src[2],x_src[3]);
  
  //   strcpy(gauge_name,qcd_getParam("<cfg_name>",params,params_len));
  //  if(myid==0) printf(" Got conf name: %s\n",gauge_name); 

  //-Propagator and source names  
  strcpy(uprop_name,qcd_getParam("<uprop_list>",params,params_len));
  strcpy(dprop_name,qcd_getParam("<dprop_list>",params,params_len));    
  strcpy(sprop_name,qcd_getParam("<sprop_list>",params,params_len));
  strcpy(cprop_name,qcd_getParam("<cprop_list>",params,params_len));
     
  if(myid==0) printf(" Got u-propagator list file: %s\n",uprop_name);
  if(myid==0) printf(" Got d-propagator list file: %s\n",dprop_name);
  if(myid==0) printf(" Got s-propagator list file: %s\n",sprop_name);
  if(myid==0) printf(" Got c-propagator list file: %s\n",cprop_name);   

  strcpy(uout_name,qcd_getParam("<usource_list>",params,params_len));
  strcpy(dout_name,qcd_getParam("<dsource_list>",params,params_len));  
  strcpy(sout_name,qcd_getParam("<ssource_list>",params,params_len));
  strcpy(cout_name,qcd_getParam("<csource_list>",params,params_len));   
   
  if(myid==0) printf(" Got sequential source list file for up %s\n",uout_name);    
  if(myid==0) printf(" Got sequential source list file for down %s\n",dout_name);
  if(myid==0) printf(" Got sequential source list file for strange %s\n",sout_name);
  if(myid==0) printf(" Got sequential source list file for charm %s\n",cout_name);
  //---  
    
  sscanf(qcd_getParam("<t_interval>",params,params_len),"%hd",&t_int);
   
  if(myid==0) printf("Got time interval between source-current: %hd\n", t_int);

  t_try = x_src[0] + t_int;
   

  //  t_current = ( t_try >= L[0] ) ? (t_try-L[0]) : t_try ;



  t_current = ((t_int+x_src[0])%geo.L[0]) - geo.Pos[0]*geo.lL[0];


  if(myid==0) printf("The time slice for current is: %hd\n", t_current);


  sscanf(qcd_getParam("<current_type>",params,params_len),"%d",&ctype);
  if(myid==0) printf(" Got current type: %d\n",ctype);

  sscanf(qcd_getParam("<momentum_xyz>",params,params_len),"%d %d %d",&mom[0],&mom[1],&mom[2]);
  if(myid==0) printf("Got momentum (%d,%d,%d)\n",mom[0],mom[1],mom[2]);
   
  free(params);
  //------------------------------------------------------------------------------------------  
  //-Read source names
  char usource_names[nsources][qcd_MAX_STRING_LENGTH];
  char dsource_names[nsources][qcd_MAX_STRING_LENGTH];
  char ssource_names[nsources][qcd_MAX_STRING_LENGTH];	
  char csource_names[nsources][qcd_MAX_STRING_LENGTH];
	
  if((fpusource=fopen(uout_name,"r")) == NULL){
    if(myid==0) fprintf(stderr,"Error! Cannot open %s for reading.\n",uout_name);
    exit(EXIT_FAILURE);
  }
  if((fpdsource=fopen(dout_name,"r")) == NULL){
    if(myid==0) fprintf(stderr,"Error! Cannot open %s for reading.\n",dout_name);
    exit(EXIT_FAILURE);
  }
  if((fpssource=fopen(sout_name,"r")) == NULL){
    if(myid==0) fprintf(stderr,"Error! Cannot open %s for reading.\n",sout_name);
    exit(EXIT_FAILURE);
  }			   
  if((fpcsource=fopen(cout_name,"r")) == NULL){
    if(myid==0) fprintf(stderr,"Error! Cannot open %s for reading.\n",cout_name);
    exit(EXIT_FAILURE);
  }


  for(i=0;i<nsources;i++){
    if(feof(fpusource)){
      if(myid==0) printf(" %s: not enough entries\n", uout_name);
      exit(EXIT_FAILURE);
    }
    if(feof(fpdsource)){
      if(myid==0) printf(" %s: not enough entries\n", dout_name);
      exit(EXIT_FAILURE);
    }
    if(feof(fpssource)){
      if(myid==0) printf(" %s: not enough entries\n", sout_name);
      exit(EXIT_FAILURE);
    }
    if(feof(fpcsource)){
      if(myid==0) printf(" %s: not enough entries\n", cout_name);
      exit(EXIT_FAILURE);
    }

    fscanf(fpusource,"%s\n",usource_names[i]);
    fscanf(fpdsource,"%s\n",dsource_names[i]);
    fscanf(fpssource,"%s\n",ssource_names[i]);	
    fscanf(fpcsource,"%s\n",csource_names[i]);						
  }
  //---

  //################################################################################

  //-Allocate memory for Props
  j = 0;
  j += qcd_initPropagator(&uprop, &geo);
  j += qcd_initPropagator(&dprop, &geo);
  j += qcd_initPropagator(&sprop, &geo);
  j += qcd_initPropagator(&cprop, &geo);
  j += qcd_initVector(&usource,&geo);
  j += qcd_initVector(&dsource,&geo);
  j += qcd_initVector(&ssource,&geo);      
  j += qcd_initVector(&csource,&geo);            

  MPI_Allreduce(&j, &k, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if(k>0)
    {
      if(myid==0) printf("not enough memory\n");
      exit(EXIT_FAILURE);
    }

  //-Load Props
  if(qcd_getPropagator(uprop_name,qcd_PROP_LIME, &uprop)) exit(EXIT_FAILURE);
  if(myid==0) printf("up propagator loaded\n");
  if(qcd_getPropagator(dprop_name,qcd_PROP_LIME, &dprop)) exit(EXIT_FAILURE);
  if(myid==0) printf("down propagator loaded\n");
  if(qcd_getPropagator(sprop_name,qcd_PROP_LIME, &sprop)) exit(EXIT_FAILURE);
  if(myid==0) printf("strange propagator loaded\n");
  if(qcd_getPropagator(cprop_name,qcd_PROP_LIME, &cprop)) exit(EXIT_FAILURE);
  if(myid==0) printf("charm propagator loaded\n");

  //################################################################################

  //--Main kernel of the program

  lt=t_current;
        
  for(mu=0;mu<4;mu++){	
    for(c1=0;c1<3;c1++){
      qcd_zeroVector(&usource); 
      qcd_zeroVector(&dsource);
      qcd_zeroVector(&ssource);
      qcd_zeroVector(&csource);
      for(nu=0;nu<4;nu++)
	for(c2=0;c2<3;c2++){
	  for(lx=0; lx<geo.lL[1]; lx++)
	    for(ly=0; ly<geo.lL[2]; ly++)
	      for(lz=0; lz<geo.lL[3]; lz++){
		v =  qcd_LEXIC(lt,lx,ly,lz,geo.lL);
					
		x=lx+geo.Pos[1]*geo.lL[1] - x_src[1];
		y=ly+geo.Pos[2]*geo.lL[2] - x_src[2];
		z=lz+geo.Pos[3]*geo.lL[3] - x_src[3];
	
		for(ku=0;ku<4;ku++){
		  usource.D[v][nu][c2] = qcd_CADD(usource.D[v][nu][c2],
						  qcd_CMUL(qcd_CURRENT_U[ctype][nu][ku],uprop.D[v][ku][mu][c2][c1]));
				
		  dsource.D[v][nu][c2] = qcd_CADD(dsource.D[v][nu][c2],
						  qcd_CMUL(qcd_CURRENT_D[ctype][nu][ku],dprop.D[v][ku][mu][c2][c1]));
											   
		  ssource.D[v][nu][c2] = qcd_CADD(ssource.D[v][nu][c2],
						  qcd_CMUL(qcd_CURRENT_S[ctype][nu][ku],sprop.D[v][ku][mu][c2][c1]));

		  csource.D[v][nu][c2] = qcd_CADD(csource.D[v][nu][c2],
						  qcd_CMUL(qcd_CURRENT_C[ctype][nu][ku],cprop.D[v][ku][mu][c2][c1]));
					
		}

		if( mom[0]!=0 && mom[1]!=0 && mom[2]!=0 ){
		  tmp = ( (double)mom[0]*x/geo.L[1] + (double)mom[1]*y/geo.L[2] + (double)mom[2]*z/geo.L[3] )*2*M_PI;
		  C2=(qcd_complex_16) {cos(tmp), sin(tmp)};
					
		  usource.D[v][nu][c2] = qcd_CMUL( usource.D[v][nu][c2] ,C2);
		  dsource.D[v][nu][c2] = qcd_CMUL( dsource.D[v][nu][c2] ,C2);
		  ssource.D[v][nu][c2] = qcd_CMUL( ssource.D[v][nu][c2] ,C2);
		  csource.D[v][nu][c2] = qcd_CMUL( csource.D[v][nu][c2] ,C2);				
		}
	      }//-space,lx,ly,lz
	}//-nu,c2
			
      qcd_writeVectorLime(usource_names[mu*3+c1],qcd_SOURCE_LIME,&usource);
      qcd_writeVectorLime(dsource_names[mu*3+c1],qcd_SOURCE_LIME,&dsource);
      qcd_writeVectorLime(ssource_names[mu*3+c1],qcd_SOURCE_LIME,&ssource);
      qcd_writeVectorLime(csource_names[mu*3+c1],qcd_SOURCE_LIME,&csource);

      if(myid==0) printf("Vector %d / 12 done\n",mu*3+c1+1);

    }//-c1
  }//-mu

  qcd_destroyVector(&usource);
  qcd_destroyVector(&dsource);
  qcd_destroyVector(&ssource);
  qcd_destroyVector(&csource);
  qcd_destroyPropagator(&uprop);
  qcd_destroyPropagator(&dprop);
  qcd_destroyPropagator(&sprop);
  qcd_destroyPropagator(&cprop);
  //	qcd_destroyGaugeField(&u);
	
  qcd_destroyGeometry(&geo);
  MPI_Finalize();
}//-main



