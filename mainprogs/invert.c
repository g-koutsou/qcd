/* invert.c
 *
 * solves Dx=b for x
 * where D is the twisted mass Wilson Dirac 
 * operator and b a source read from disk
 *
 *
 * Tomasz Korzec 2010
 ****************************************/
 
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <qcd.h>



qcd_int_4 gcr(qcd_vector *sol, qcd_vector *src, qcd_gaugeField *u, qcd_real_8 kappa, qcd_real_8 muTM, qcd_uint_4 n, qcd_geometry *geo)
{
   /* perform n iterations of gcr algorithm */
   qcd_vector chi[n];  //orthonormal basis of Krylov space
   qcd_vector res[n];  //residues
   qcd_complex_16 a[n][n]; //coefficients 
   qcd_complex_16 c[n];    //coefficients (chi,b)
   qcd_complex_16 tmp;
   qcd_int_4 i,j=0,k,l;
   
   /* get memory for basis vectors and residue */
   for(i=0; i<n; i++)
   {
      j += qcd_initVector(&(chi[i]), geo);
      j += qcd_initVector(&(res[i]), geo);
   }
   MPI_Allreduce(&j, &k, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   if(k>0)
   {
      if(geo->myid==0) printf("not enough memory\n");
      return(EXIT_FAILURE);
   }
   
   
   /* first iteration */
   qcd_copyVector(&(res[0]), src);
   qcd_applyQTMOp(&(chi[0]),&(res[0]),u, 1.0/(2.0*kappa)-4.0 ,muTM);
   a[0][0] = (qcd_complex_16) {1.0/qcd_normVector(&(chi[0])), 0};
   qcd_scaleVector(&(chi[0]), a[0][0].re);
   c[0] = qcd_mulAdjointVectorVector(&(chi[0]), src);
   qcd_zeroVector(sol);
   
   //qcd_copyVector(sol, &(res[0]));
   //qcd_mulVectorC(sol, qcd_CMUL(c[0],a[0][0]));   
   
   /* iterations 2 to n */
   for(k=1; k<n; k++)
   {
      /* calculate new residuum */
      qcd_axpyVector(&(res[k]), (qcd_complex_16) {-c[k-1].re, -c[k-1].im}, &(chi[k-1]), &(res[k-1]));
      
      tmp.re = qcd_normVector(&(res[k]));
      if(geo->myid==0) printf("iteration %i: |b-Dx| = %e\n",k, tmp.re);
      
      /* increase Krylov space by */
      qcd_applyQTMOp(&(chi[k]), &(res[k]), u, 1.0/(2.0*kappa)-4.0 ,muTM);
      
      /* orthonormalization */
      for(j=0; j<k; j++)
      {
         tmp = qcd_mulAdjointVectorVector(&(chi[j]), &(chi[k]));
         a[k][j] = (qcd_complex_16) { -tmp.re, -tmp.im};
         qcd_axpyVector(&(chi[k]), a[k][j], &(chi[j]), &(chi[k]));
      }
      a[k][k] = (qcd_complex_16) {1.0/qcd_normVector(&(chi[k])), 0};
      if(isnan(a[k][k].re))
      {
         if(geo->myid==0) fprintf(stderr,"Error! Krylov space is sick\n");
         return(EXIT_FAILURE);
      }
      qcd_scaleVector(&(chi[k]), a[k][k].re);
      
      /* get the coefficients */
      for(j=0; j<k; j++)
      {
         a[k][j] = qcd_CSCALE(a[k][j],a[k][k].re);
      }
      for(j=0; j<k; j++)
      {
         a[k][j] = qcd_CSCALE(a[k][j],a[j][j].re);
         for(l=j+1; l<k; l++)
         {
            a[k][j] = qcd_CADD(a[k][j], qcd_CMUL(a[k][l],a[l][j]));
         }
      }
      
      /* calculate next c */
      c[k] = qcd_mulAdjointVectorVector(&(chi[k]),src);
      
      /* update solution */
      for(j=0; j<=k; j++)
      {
         qcd_axpyVector(sol, qcd_CMUL(c[k],a[k][j]), &(res[j]), sol);
      }   
   }
   
   /* free memory */
   for(i=0; i<n; i++)
   {
      qcd_destroyVector(&(chi[i]));
      qcd_destroyVector(&(res[i]));
   }
   return(0);
}//end gcr



int main(int argc,char* argv[])
{   
   qcd_int_4 i,j,k;              // various loop variables
   qcd_uint_2 mu,col;
   
   int params_len;               // needed to read inputfiles
   char *params = NULL;                 // needed to read inputfiles

   
   qcd_uint_4 Nrestart;                         // restart GCR every Nrestart iterations
   char gauge_name[qcd_MAX_STRING_LENGTH];      // name of gauge configuration
   char param_name[qcd_MAX_STRING_LENGTH];      // name of parameter file  
   char sol_name[qcd_MAX_STRING_LENGTH];        // name of solution file
   char src_name[qcd_MAX_STRING_LENGTH];        // name of source file
   char src_type[qcd_MAX_STRING_LENGTH];        // source type
   char sol_type[qcd_MAX_STRING_LENGTH];        // solution type
   qcd_real_8 kappa;                            // hopping parameter
   qcd_real_8 muTM;                             // twisted mass parameter
   qcd_real_8 normsrc,normres;                  // norm of source, norm of residue
   qcd_uint_4 maxIter = 10000;
   qcd_uint_4 iter;
   qcd_real_8 maxRes = 1e-8;
 
   qcd_geometry geo;                            // geometry structure
 
   qcd_real_8 theta[4] = {M_PI,0.0,0.0,0.0};    // antiperiodic b.c. in time
   qcd_uint_2 L[4];
   qcd_uint_2 P[4];

   qcd_vector src;
   qcd_vector sol;
   qcd_vector res;
   qcd_vector correction;
   qcd_gaugeField u;

   
   int myid,numprocs, namelen;    
   char processor_name[MPI_MAX_PROCESSOR_NAME];
   				 
				 
             
             
             
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
   if(qcd_initGeometry(&geo,L,P, theta, myid, numprocs)) exit(EXIT_FAILURE);
   
   if(myid==0) printf(" Local lattice: %i x %i x %i x %i\n",geo.lL[0],geo.lL[1],geo.lL[2],geo.lL[3]);
  
   strcpy(src_type,qcd_getParam("<source_type>",params,params_len));
   if(myid==0) printf("Got source type: %s\n",src_type);
   
   /* src_type == "HMC_PROPAGATOR", propagator with 12 vectors */
   strcpy(src_name,qcd_getParam("<source>",params,params_len));
   if(myid==0) printf("Got source file name: %s\n",src_name);

   strcpy(sol_type,qcd_getParam("<solution_type>",params,params_len));
   if(myid==0) printf("Got solution type: %s\n",sol_type);
  
   strcpy(sol_name,qcd_getParam("<solution>",params,params_len));
   if(myid==0) printf("Got solution file name: %s\n",sol_name);
   
   strcpy(gauge_name,qcd_getParam("<cfg_name>",params,params_len));
   if(myid==0) printf("Got conf name: %s\n",gauge_name);
   
   sscanf(qcd_getParam("<N_restart>",params,params_len),"%u",&Nrestart);
   if(myid==0) printf("Got N_restart: %u\n",Nrestart);

   sscanf(qcd_getParam("<kappa>",params,params_len),"%lf",&kappa);
   if(myid==0) printf("Got kappa: %e\n",kappa);
   
   sscanf(qcd_getParam("<mu>",params,params_len),"%lf",&muTM);
   if(myid==0) printf("Got mu: %e\n",muTM);
              
   free(params);



         
   //#####################################################################   
   // allocate memory
  
   /* src_type == HMC */
   j = 0;
   j += qcd_initVector(&src, &geo);
   j += qcd_initVector(&sol, &geo);
   j += qcd_initVector(&res, &geo);
   j += qcd_initVector(&correction, &geo);
   j += qcd_initGaugeField(&u, &geo);
   
   MPI_Allreduce(&j, &k, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   if(k>0)
   {
      if(myid==0) printf("not enough memory\n");
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
   }
   if(myid==0) printf("memory for propagators and gauge-field allocated\n");
         
   
   //##############################################################################
   // load gauge-field
   if(qcd_getGaugeField(gauge_name,qcd_GF_LIME,&u)) exit(EXIT_FAILURE);
   if(myid==0) printf("gauge-field loaded\n");   
   
   
   
   for(mu=0; mu<4; mu++)
   for(col=0; col<3; col++)
   {
      if(myid==0) printf("------------ vector: mu = %hi,  col = %hi ------------\n",mu,col);
      
      iter = 1;
      if(qcd_getVector(src_name,qcd_PROP_HMC, mu, col, &src)) exit(EXIT_FAILURE);
      if(myid==0) printf("vector from %s loaded\n",src_name);   
      
      normsrc = qcd_normVector(&src);
      if(myid==0) printf("Norm of source: %e\n",normsrc);
      
      gcr(&sol, &src, &u, kappa, muTM, Nrestart, &geo);      
      
      /* calculate true residue */
      qcd_applyQTMOp(&res, &sol, &u, 1.0/(2.0*kappa)-4.0 ,muTM);
      
      qcd_subVector(&res, &src, &res);
      normres = qcd_normVector(&res);
      if(myid==0) printf("True norm of residue: %e\n",normres);
      normres /= normsrc;
      if(myid==0) printf("Relative residue: %e\n",normres);
      
      /* iterative improvement until precision reached */
      while(normres>maxRes && iter < maxIter)
      {
         /* solve D correction = residue */
         /* and set solution <- solution - correction */
         gcr(&correction, &res, &u, kappa, muTM, Nrestart, &geo);
         qcd_addVector(&sol,&sol,&correction);
         
         /* calculate true residue */
         qcd_applyQTMOp(&res, &sol, &u, 1.0/(2.0*kappa)-4.0 ,muTM);
         qcd_subVector(&res, &src, &res);
         normres = qcd_normVector(&res)/normsrc;
         if(myid==0) printf("True relative residue: %e\n",normres);
         iter++;
      }
      if(myid==0) printf("Converged after %i x %i iterations.\n\n",iter,Nrestart);
   }
   
   
   //#####################################################################   
   // clean up
   if(myid==0) printf("cleaning up...\n");
   qcd_destroyVector(&src);
   qcd_destroyVector(&sol);
   qcd_destroyVector(&res);
   qcd_destroyVector(&correction);
   qcd_destroyGaugeField(&u);
   qcd_destroyGeometry(&geo);
   MPI_Finalize();
}//end main
