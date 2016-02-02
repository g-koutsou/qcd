/* Christos Kallidonis
 * 
 * Modifications to the original qcd_gaussIteration3dAll to include all four propagators
 * and doing the iterations in the routine to save time in allocating-deallocating-calling.
 * 
 * qcd_smearing_opt.c
 *
 * gaussian smearing for all propagators for time slices needed (modification of original)
 *
 **********************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <qcd.h>
 
 /* perform nsmear iterations of gaussian smearing 
 * on time slices t_start to t_stop with
 * parameter alpha.
 *
 * v-> (v+alpha*Hv) / (1+6*alpha)
 */
int qcd_gaussIteration3d_opt(qcd_propagator *uprop,qcd_propagator *dprop,qcd_propagator *sprop,qcd_propagator *cprop,
                                    qcd_geometry *geo, qcd_gaugeField *u, qcd_uint_4 nsmear, qcd_real_8 alpha,
                                    qcd_uint_4 t_start, qcd_uint_4 t_stop, qcd_uint_4 t_src)
{
   qcd_uint_8 i,j;
   qcd_uint_2 c1,mu,nu; 
   qcd_uint_4 x,y,z,t,tt,k;
   qcd_real_8 nrm = 1.0/(1.0 + alpha*6.0);
   qcd_vector v2_u,v2_d,v2_s,v2_c;
   qcd_real_8 *uu;
   qcd_real_8 *psi;
   qcd_real_8 upsi_u[24],upsi_d[24],upsi_s[24],upsi_c[24];
   qcd_real_8 udaggerpsi_u[24],udaggerpsi_d[24],udaggerpsi_s[24],udaggerpsi_c[24];
   qcd_real_8 *total_u,*total_d,*total_s,*total_c;
   qcd_vector vec_u,vec_d,vec_s,vec_c;

   j=0;
   j += qcd_initVector(&vec_u, geo);
   j += qcd_initVector(&vec_d, geo);  
   j += qcd_initVector(&vec_s, geo);
   j += qcd_initVector(&vec_c, geo);
   j += qcd_initVector(&v2_u, geo);
   j += qcd_initVector(&v2_d, geo);
   j += qcd_initVector(&v2_s, geo);
   j += qcd_initVector(&v2_c, geo);
            
   MPI_Allreduce(&j, &k, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

   if(k>0)
   {
      if(geo->myid==0) printf("not enough memory for smearing\n");
      exit(EXIT_FAILURE);
   }


	for(mu=0;mu<4;mu++)
	for(c1=0;c1<3;c1++)
	{	
	
	qcd_copyVectorPropagator(&vec_u,uprop,mu,c1);
	qcd_copyVectorPropagator(&vec_d,dprop,mu,c1);
	qcd_copyVectorPropagator(&vec_s,sprop,mu,c1);
	qcd_copyVectorPropagator(&vec_c,cprop,mu,c1);

	for(int ismear = 0 ; ismear < nsmear ; ismear++){

		if(!ismear){
			qcd_communicateGaugeP(u);
			qcd_waitall(u->geo);
		}
		
		qcd_communicateVectorPM(&vec_u);
		qcd_waitall(vec_u.geo);
		qcd_communicateVectorPM(&vec_d);
		qcd_waitall(vec_d.geo);
		qcd_communicateVectorPM(&vec_s);
		qcd_waitall(vec_s.geo);		
		qcd_communicateVectorPM(&vec_c);				
		qcd_waitall(vec_c.geo);
 
		qcd_zeroVector(&v2_u);
		qcd_zeroVector(&v2_d);
		qcd_zeroVector(&v2_s);
		qcd_zeroVector(&v2_c);				

		for(z=0;z<geo->lL[3];z++)
		for(y=0;y<geo->lL[2];y++)
		for(x=0;x<geo->lL[1];x++)
 								
		for(t=t_start;t<=t_stop;t++){
				
			tt = ((t+t_src)%geo->L[0]) - geo->Pos[0]*geo->lL[0];
			if(tt>=0 && tt<geo->lL[0]){  //inside the local lattice, otherwise nothing to calculate
					
				i=qcd_LEXIC(tt,x,y,z,geo->lL);
				for(nu=1;nu<4;nu++){
	
					uu  = (qcd_real_8*) &(u->D[i][nu][0][0].re);
					psi = (qcd_real_8*) &(vec_u.D[geo->plus[i][nu]][0][0].re);
					qcd_APPLY_U(uu,upsi_u,psi);
					psi = (qcd_real_8*) &(vec_d.D[geo->plus[i][nu]][0][0].re);
					qcd_APPLY_U(uu,upsi_d,psi);
					psi = (qcd_real_8*) &(vec_s.D[geo->plus[i][nu]][0][0].re);
					qcd_APPLY_U(uu,upsi_s,psi);
					psi = (qcd_real_8*) &(vec_c.D[geo->plus[i][nu]][0][0].re);
					qcd_APPLY_U(uu,upsi_c,psi);
					
					uu  = (qcd_real_8*) &(u->D[geo->minus[i][nu]][nu][0][0].re);
					psi = (qcd_real_8*) &(vec_u.D[geo->minus[i][nu]][0][0].re); 
					qcd_APPLY_U_DAGGER(uu,udaggerpsi_u,psi);
					psi = (qcd_real_8*) &(vec_d.D[geo->minus[i][nu]][0][0].re); 
					qcd_APPLY_U_DAGGER(uu,udaggerpsi_d,psi);
					psi = (qcd_real_8*) &(vec_s.D[geo->minus[i][nu]][0][0].re); 
					qcd_APPLY_U_DAGGER(uu,udaggerpsi_s,psi);
					psi = (qcd_real_8*) &(vec_c.D[geo->minus[i][nu]][0][0].re); 
					qcd_APPLY_U_DAGGER(uu,udaggerpsi_c,psi);

					total_u = (qcd_real_8*) &(v2_u.D[i][0][0].re);
					total_d = (qcd_real_8*) &(v2_d.D[i][0][0].re);
					total_s = (qcd_real_8*) &(v2_s.D[i][0][0].re);
					total_c = (qcd_real_8*) &(v2_c.D[i][0][0].re);
					
					qcd_SUM_UP_HOPP(total_u,upsi_u,udaggerpsi_u);
					qcd_SUM_UP_HOPP(total_d,upsi_d,udaggerpsi_d);
					qcd_SUM_UP_HOPP(total_s,upsi_s,udaggerpsi_s);
					qcd_SUM_UP_HOPP(total_c,upsi_c,udaggerpsi_c);
																									
				}//-nu
			}//-tt condition
		}//-outer loop (space-time)
		
		qcd_scaleVector(&v2_u,alpha);
		qcd_scaleVector(&v2_d,alpha);
		qcd_scaleVector(&v2_s,alpha);
		qcd_scaleVector(&v2_c,alpha);	
					
		qcd_addVector(&vec_u,&vec_u,&v2_u);
		qcd_addVector(&vec_d,&vec_d,&v2_d);
		qcd_addVector(&vec_s,&vec_s,&v2_s);
		qcd_addVector(&vec_c,&vec_c,&v2_c);
		
		qcd_scaleVector(&vec_u,nrm); 
		qcd_scaleVector(&vec_d,nrm);
		qcd_scaleVector(&vec_s,nrm);
		qcd_scaleVector(&vec_c,nrm);				

	
	}//-smear iters
	
	qcd_copyPropagatorVector(uprop,&vec_u,mu,c1);
	qcd_copyPropagatorVector(dprop,&vec_d,mu,c1);
	qcd_copyPropagatorVector(sprop,&vec_s,mu,c1);
	qcd_copyPropagatorVector(cprop,&vec_c,mu,c1);		

 	if(geo->myid==0) printf("Smearing: Indices mu=%d,c1=%d smeared\n",mu,c1);

	}//-mu,c1

   qcd_destroyVector(&v2_u);
   qcd_destroyVector(&v2_d);
   qcd_destroyVector(&v2_s);   
   qcd_destroyVector(&v2_c);      
   qcd_destroyVector(&vec_u);
   qcd_destroyVector(&vec_d);	
   qcd_destroyVector(&vec_s);
   qcd_destroyVector(&vec_c);
         
   return 0;
   
}//end qcd_gaussIteration3d_opt


