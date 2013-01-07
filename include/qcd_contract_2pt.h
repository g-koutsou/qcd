/* 
 * Christos Kallidonis
 * April 2012
 * 
 * Header file for qcd_contract_2pt.c
 * 
 */
 

int qcd_f1f2f1_2pt(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[4][4],
						qcd_propagator *propf1,qcd_propagator *propf2, qcd_geometry *geo, qcd_uint_4 lt);

int qcd_pnextra_2pt(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[4][4],
						qcd_propagator *propf1,qcd_propagator *propf2, qcd_geometry *geo, qcd_uint_4 lt);

int qcd_f1f2f3_2pt(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[4][4],
						qcd_propagator *propf1, qcd_propagator *propf2, qcd_propagator *propf3, qcd_geometry *geo, qcd_uint_4 lt);

int qcd_f123f321_2pt(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[4][4],
						qcd_propagator *propf1, qcd_propagator *propf2, qcd_propagator *propf3, qcd_geometry *geo, qcd_uint_4 lt);

int qcd_lambdas_xis_2pt(qcd_int_4 ctr, qcd_int_2 cg5cg5_ind[16*16][4],qcd_complex_16 cg5cg5_val[16*16], qcd_complex_16 *block[4][4],
						qcd_propagator *propf1, qcd_propagator *propf2, qcd_propagator *propf3, qcd_geometry *geo, qcd_uint_4 lt);

 //=====================================================================
void qcd_contractions2pt(qcd_uint_4 particle_id, qcd_complex_16 *block[4][4], qcd_propagator *uprop,qcd_propagator *dprop,qcd_propagator *sprop,qcd_propagator *cprop,
			 qcd_geometry *geo, qcd_uint_4 lt);
