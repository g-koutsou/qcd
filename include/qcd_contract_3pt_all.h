/* 
 * Christos Kallidonis
 * June 2012
 * 
 * Header file for qcd_contract_3pt_all.c
 * 
 */

int qcd_projector_3pt_spin12(qcd_complex_16 *block[5], qcd_complex_16 *block_pr[9][4][4], qcd_geometry *geo);

int qcd_projector_pr32_3pt_spin32(qcd_complex_16 *block[5], qcd_complex_16 *block_pr[9][4][4],
                                  qcd_complex_16 gamma12[4][4], qcd_complex_16 gamma13[4][4], qcd_complex_16 gamma23[4][4], qcd_geometry *geo);					       
//-----------------------------------

int qcd_f1f1f1_3pt(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[9][4][4],
		   qcd_propagator *prop, qcd_propagator *seqprop, qcd_geometry *geo, qcd_uint_4 lt,qcd_uint_4 elem);

int qcd_f1f2f1_f1_3pt(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[9][4][4],
		      qcd_propagator *propf1, qcd_propagator *propf2, qcd_propagator *seqpropf1, qcd_geometry *geo, qcd_uint_4 lt, qcd_uint_4 elem);

int qcd_pnextra_f1_3pt(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[9][4][4],
		       qcd_propagator *propf1, qcd_propagator *propf2, qcd_propagator *seqpropf1, qcd_geometry *geo, qcd_uint_4 lt, qcd_uint_4 elem);
	    						
int qcd_f1f2f1_f2_3pt(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[9][4][4],
		      qcd_propagator *propf1, qcd_propagator *propf2, qcd_propagator *seqpropf2, qcd_geometry *geo, qcd_uint_4 lt, qcd_uint_4 elem);	    			

int qcd_pnextra_f2_3pt(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[9][4][4],
		       qcd_propagator *propf1, qcd_propagator *propf2, qcd_propagator *seqpropf2, qcd_geometry *geo, qcd_uint_4 lt, qcd_uint_4 elem);
	    			
int qcd_deltas_f1_3pt(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[9][4][4],
		      qcd_propagator *propf1, qcd_propagator *propf2, qcd_propagator *seqpropf1, qcd_geometry *geo, qcd_uint_4 lt, qcd_uint_4 elem);
						
int qcd_deltas_f2_3pt(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[9][4][4],
		      qcd_propagator *propf1, qcd_propagator *propf2, qcd_propagator *seqpropf2, qcd_geometry *geo, qcd_uint_4 lt, qcd_uint_4 elem);

int qcd_xistar_extra11_f1_3pt(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[9][4][4],
			      qcd_propagator *propf1, qcd_propagator *propf2, qcd_propagator *seqpropf1, qcd_geometry *geo, qcd_uint_4 lt, qcd_uint_4 elem);

int qcd_xistar_extra21_f1_3pt(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[9][4][4],
			      qcd_propagator *propf1, qcd_propagator *propf2, qcd_propagator *seqpropf1, qcd_geometry *geo, qcd_uint_4 lt, qcd_uint_4 elem);

int qcd_xistar_extra12_f1_3pt(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[9][4][4],
			      qcd_propagator *propf1, qcd_propagator *propf2, qcd_propagator *seqpropf1, qcd_geometry *geo, qcd_uint_4 lt, qcd_uint_4 elem);

int qcd_xistar_extra22_f1_3pt(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[9][4][4],
			      qcd_propagator *propf1, qcd_propagator *propf2, qcd_propagator *seqpropf1, qcd_geometry *geo, qcd_uint_4 lt, qcd_uint_4 elem);

int qcd_xistar_extra11_f2_3pt(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[9][4][4],
			      qcd_propagator *propf1, qcd_propagator *propf2, qcd_propagator *seqpropf2, qcd_geometry *geo, qcd_uint_4 lt, qcd_uint_4 elem);

int qcd_xistar_extra12_f2_3pt(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[9][4][4],
			      qcd_propagator *propf1, qcd_propagator *propf2, qcd_propagator *seqpropf2, qcd_geometry *geo, qcd_uint_4 lt, qcd_uint_4 elem);

int qcd_xistar_extra21_f2_3pt(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[9][4][4],
			      qcd_propagator *propf1, qcd_propagator *propf2, qcd_propagator *seqpropf2, qcd_geometry *geo, qcd_uint_4 lt, qcd_uint_4 elem);

int qcd_xistar_extra22_f2_3pt(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[9][4][4],
			      qcd_propagator *propf1, qcd_propagator *propf2, qcd_propagator *seqpropf2, qcd_geometry *geo, qcd_uint_4 lt, qcd_uint_4 elem);

int qcd_f123f321_f1_3pt(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[9][4][4],
			qcd_propagator *propf2, qcd_propagator *propf3, qcd_propagator *seqpropf1, qcd_geometry *geo, qcd_uint_4 lt, qcd_uint_4 elem);
					    
int qcd_f123f321_f2_3pt(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[9][4][4],
			qcd_propagator *propf1, qcd_propagator *propf3, qcd_propagator *seqpropf2, qcd_geometry *geo, qcd_uint_4 lt, qcd_uint_4 elem);
					    
int qcd_f123f321_f3_3pt(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[9][4][4],
			qcd_propagator *propf1, qcd_propagator *propf2, qcd_propagator *seqpropf3, qcd_geometry *geo, qcd_uint_4 lt, qcd_uint_4 elem);					   

int qcd_f1f2f3_3pt(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[9][4][4],
		   qcd_propagator *propf1, qcd_propagator *propf2, qcd_propagator *propf3, qcd_geometry *geo, qcd_uint_4 lt, qcd_uint_4 elem);

int qcd_sigmas4_f1_3pt(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[9][4][4],
		       qcd_propagator *propf1, qcd_propagator *propf2, qcd_propagator *seqpropf1, qcd_geometry *geo, qcd_uint_4 lt, qcd_uint_4 elem);

int qcd_sigmas4_f2_3pt(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[9][4][4],
		       qcd_propagator *propf1, qcd_propagator *propf2, qcd_propagator *seqpropf2, qcd_geometry *geo, qcd_uint_4 lt, qcd_uint_4 elem);

int qcd_lambdas_f2_3pt(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[9][4][4],
		       qcd_propagator *propf1, qcd_propagator *propf3, qcd_propagator *seqpropf2, qcd_geometry *geo, qcd_uint_4 lt, qcd_uint_4 elem);
		                
int qcd_lambdas_f3_3pt(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[9][4][4],
		       qcd_propagator *propf1, qcd_propagator *propf2, qcd_propagator *seqpropf3, qcd_geometry *geo, qcd_uint_4 lt, qcd_uint_4 elem);
		                
int qcd_sigmas2_f1_3pt(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[9][4][4],
		       qcd_propagator *propf2, qcd_propagator *propf3, qcd_propagator *seqpropf1, qcd_geometry *geo, qcd_uint_4 lt, qcd_uint_4 elem, qcd_uint_4 xis);	
		                
int qcd_sigmas2_f2_3pt(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[9][4][4],
		       qcd_propagator *propf1, qcd_propagator *propf3, qcd_propagator *seqpropf2, qcd_geometry *geo, qcd_uint_4 lt, qcd_uint_4 elem, qcd_uint_4 xis);		                	                		                
int qcd_sigmas2_f3_3pt(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[9][4][4],
		       qcd_propagator *propf1, qcd_propagator *propf2, qcd_propagator *seqpropf3, qcd_geometry *geo, qcd_uint_4 lt, qcd_uint_4 elem, qcd_uint_4 xis);		
		                                					   						
//=====================================================================

void qcd_contractions3pt_new(qcd_uint_4 p_id, qcd_uint_4 np, qcd_complex_16 *block[5],
			     qcd_propagator *uprop, qcd_propagator *dprop, qcd_propagator *sprop, qcd_propagator *cprop,
			     qcd_propagator *sequprop, qcd_propagator *seqdprop, qcd_propagator *seqsprop, qcd_propagator *seqcprop,			 
			     qcd_geometry *geo, qcd_uint_4 lt);
