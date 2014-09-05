/* 
 * Christos Kallidonis
 * April 2012
 * 
 * Header file for qcd_contract_2pt_pr.c
 * 
 */

int qcd_projector12_2pt(qcd_complex_16 gamma12[4][4], qcd_complex_16 gamma13[4][4], qcd_complex_16 gamma23[4][4], 
						qcd_complex_16 *block[4][4], qcd_complex_16 *block_pr[9][4][4], qcd_geometry *geo);
 
 
int qcd_projector32_2pt(qcd_complex_16 gamma12[4][4], qcd_complex_16 gamma13[4][4], qcd_complex_16 gamma23[4][4], 
						qcd_complex_16 *block[4][4], qcd_complex_16 *block_pr[9][4][4], qcd_geometry *geo);

int qcd_noprojector_2pt(qcd_complex_16 *block[4][4], qcd_complex_16 *block_pr[9][4][4], qcd_geometry *geo);

int qcd_f1f1f1_2pt_pr(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[9][4][4],
						qcd_propagator *prop, qcd_geometry *geo, qcd_uint_4 lt,qcd_uint_4 elem);

int qcd_deltas_xistar_omegastar_2pt_pr(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[9][4][4],
						qcd_propagator *propf1, qcd_propagator *propf2, qcd_geometry *geo, qcd_uint_4 lt, qcd_uint_4 elem);
						
int qcd_f1f2f1_2pt_pr(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[9][4][4],
						qcd_propagator *propf1, qcd_propagator *propf2, qcd_geometry *geo, qcd_uint_4 lt,qcd_uint_4 elem);						

int qcd_f1f2f3_2pt_pr(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[9][4][4],
						qcd_propagator *propf1, qcd_propagator *propf2, qcd_propagator *propf3, qcd_geometry *geo, qcd_uint_4 lt, qcd_uint_4 elem);
						
int qcd_sigmas4_2pt_pr(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[9][4][4],
						qcd_propagator *propf1, qcd_propagator *propf2, qcd_geometry *geo, qcd_uint_4 lt, qcd_uint_4 elem);
		    
int qcd_sigmas2_xistar_2pt_pr(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[9][4][4],
						qcd_propagator *propf1, qcd_propagator *propf2, qcd_propagator *propf3, qcd_geometry *geo, qcd_uint_4 lt, qcd_uint_4 elem, qcd_uint_4 sflag);	    						
						
 //=====================================================================
void qcd_contractions2pt_pr(qcd_uint_4 particle_id, qcd_complex_16 *block12[4][4], qcd_complex_16 *block32[4][4], qcd_complex_16 *blocknp[4][4],
							qcd_propagator *uprop,qcd_propagator *dprop,qcd_propagator *sprop,qcd_propagator *cprop,
							qcd_geometry *geo, qcd_uint_4 lt);
