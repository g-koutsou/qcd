/*
 * Christos Kallidonis
 * June 2012
 * 
 * This program contains contractions for the 3/2 projected 2-point functions 
 * of the spin-3/2 particles 
 * 
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <qcd.h>

int qcd_projector12_2pt(qcd_complex_16 gamma12[4][4], qcd_complex_16 gamma13[4][4], qcd_complex_16 gamma23[4][4], 
			qcd_complex_16 *block[4][4], qcd_complex_16 *block_pr[9][4][4], qcd_geometry *geo){
	
  qcd_uint_2 al,be,ga,j,i;
  qcd_uint_4 lx,ly,lz,v,v3;
  qcd_complex_16 Pr12[9][4][4];
  qcd_real_8 fac = 1.0/3.0;
  qcd_complex_16 cfac;
  
  for(i=0;i<9;i++){
    for(al=0;al<4;al++){
      for(be=0;be<4;be++){
	Pr12[i][al][be] = (qcd_complex_16) {0.0,0.0};
      }
    }
  }
  
  // Define the digonal elements of the projector to 1/2
  // elements of [0][al][be], [4][al][be] [8][al][be] for al!=be are zero by definition
  cfac  = (qcd_complex_16) {1.0/3.0,0.0};
  for(al=0;al<4;al++){
    Pr12[0][al][al] = cfac; // Pr_11
    Pr12[4][al][al] = cfac; // Pr_22
    Pr12[8][al][al] = cfac; // Pr_33
  }
  
  // Define the rest elements
  for(al=0;al<4;al++){
    for(be=0;be<4;be++){
      Pr12[1][al][be] = qcd_CSCALE(gamma12[al][be],fac); // Pr_12
      Pr12[2][al][be] = qcd_CSCALE(gamma13[al][be],fac); // Pr_13
      Pr12[5][al][be] = qcd_CSCALE(gamma23[al][be],fac); // Pr_23
      
      Pr12[3][al][be] = qcd_CSCALE(qcd_CSCALE(gamma12[al][be],-1.0),fac); // Pr_21
      Pr12[6][al][be] = qcd_CSCALE(qcd_CSCALE(gamma13[al][be],-1.0),fac); // Pr_31
      Pr12[7][al][be] = qcd_CSCALE(qcd_CSCALE(gamma23[al][be],-1.0),fac); // Pr_32
    }
  }
  
  for(lx=0; lx<geo->lL[1]; lx++) 
    for(ly=0; ly<geo->lL[2]; ly++) 
      for(lz=0; lz<geo->lL[3]; lz++){ 	
	v3 = qcd_LEXIC0(lx,ly,lz,geo->lL);

	for(j=0;j<9;j++){
	  i = ( (j < 8) ? (j*3)%8 : j );
	  
	  for(al=0;al<4;al++){
	    for(ga=0;ga<4;ga++){
	      for(be=0;be<4;be++){
		
		block[al][ga][v3] = qcd_CADD(block[al][ga][v3],qcd_CMUL(Pr12[j][al][be],block_pr[i][be][ga][v3]));
	      }
	    }
	  }
	}

      }//-space

  return (1);
} //-routine
//======================================================================
int qcd_projector32_2pt(qcd_complex_16 gamma12[4][4], qcd_complex_16 gamma13[4][4], qcd_complex_16 gamma23[4][4], 
			qcd_complex_16 *block[4][4], qcd_complex_16 *block_pr[9][4][4], qcd_geometry *geo){
	
  qcd_uint_2 al,be,ga,j,i;
  qcd_uint_4 lx,ly,lz,v,v3;
  qcd_complex_16 Pr32[9][4][4],delta;
  qcd_real_8 fac = 1.0/3.0;
  qcd_complex_16 cfac;

  for(i=0;i<9;i++){
    for(al=0;al<4;al++){
      for(be=0;be<4;be++){
	Pr32[i][al][be] = (qcd_complex_16) {0.0,0.0};
      }
    }
  }

  // Define the digonal elements of the projector to 3/2
  // elements of [0][al][be], [4][al][be] [8][al][be] for al!=be are zero by definition 
  delta = (qcd_complex_16) {1.0,0.0};
  cfac  = (qcd_complex_16) {1.0/3.0,0.0};
  for(al=0;al<4;al++){
    Pr32[0][al][al] = qcd_CSUB(delta,cfac); // Pr_11
    Pr32[4][al][al] = qcd_CSUB(delta,cfac); // Pr_22
    Pr32[8][al][al] = qcd_CSUB(delta,cfac); // Pr_33
  }

  // Define the rest elements
  delta = (qcd_complex_16) {0.0,0.0};
  for(al=0;al<4;al++){
    for(be=0;be<4;be++){
      Pr32[1][al][be] = qcd_CSUB(delta,qcd_CSCALE(gamma12[al][be],fac)); // Pr_12
      Pr32[2][al][be] = qcd_CSUB(delta,qcd_CSCALE(gamma13[al][be],fac)); // Pr_13
      Pr32[5][al][be] = qcd_CSUB(delta,qcd_CSCALE(gamma23[al][be],fac)); // Pr_23
      
      Pr32[3][al][be] = qcd_CSUB(delta,qcd_CSCALE(qcd_CSCALE(gamma12[al][be],-1.0),fac)); // Pr_21
      Pr32[6][al][be] = qcd_CSUB(delta,qcd_CSCALE(qcd_CSCALE(gamma13[al][be],-1.0),fac)); // Pr_31
      Pr32[7][al][be] = qcd_CSUB(delta,qcd_CSCALE(qcd_CSCALE(gamma23[al][be],-1.0),fac)); // Pr_32
    }
  }
   
  for(lx=0; lx<geo->lL[1]; lx++) 
    for(ly=0; ly<geo->lL[2]; ly++) 
      for(lz=0; lz<geo->lL[3]; lz++){ 	
	v3 = qcd_LEXIC0(lx,ly,lz,geo->lL);

	for(j=0;j<9;j++){
	  i = ( (j < 8) ? (j*3)%8 : j );
	  
	  for(al=0;al<4;al++){
	    for(ga=0;ga<4;ga++){
	      for(be=0;be<4;be++){
		
		block[al][ga][v3] = qcd_CADD(block[al][ga][v3],qcd_CMUL(Pr32[j][al][be],block_pr[i][be][ga][v3]));
	      }
	    }
	  }
	}

      }//-space

  return (1);
} //-routine
//======================================================================
int qcd_noprojector_2pt(qcd_complex_16 *block[4][4], qcd_complex_16 *block_pr[9][4][4], qcd_geometry *geo){
  
  qcd_uint_2 ga,gap,j;
  qcd_uint_4 lx,ly,lz,v,v3;
  int lv3 = geo->lL[1]*geo->lL[2]*geo->lL[3];
    
  for(lx=0; lx<geo->lL[1]; lx++) 
    for(ly=0; ly<geo->lL[2]; ly++) 
      for(lz=0; lz<geo->lL[3]; lz++){ 
	
	v3 = qcd_LEXIC0(lx,ly,lz,geo->lL);
  
	for(ga =0; ga <4; ga ++) 
	  for(gap=0; gap<4; gap++){
	    block[ga][gap][v3] = qcd_CSCALE(qcd_CADD(block_pr[0][ga][gap][v3],qcd_CADD(block_pr[4][ga][gap][v3],block_pr[8][ga][gap][v3])),1.0/3.0); //- only diagonal terms 
	  }//-ga,gap
      }//-space

  return (1);
} //-routine
//======================================================================

int qcd_f1f1f1_2pt_pr(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[9][4][4],
		      qcd_propagator *prop, qcd_geometry *geo, qcd_uint_4 lt,qcd_uint_4 elem){
							
  qcd_int_4 ctr2;
  qcd_uint_2 a,b,c,ap,bp,cp,al,be,ga,alp,bep,gap,cc1,cc2;
  qcd_uint_4 lx,ly,lz,v,v3;
  int lv3 = geo->lL[1]*geo->lL[2]*geo->lL[3];

#pragma omp parallel for private(lz,ly,lx,v,ga,gap,ctr2,al,be,bep,alp,cc1,a,b,c,cc2,ap,bp,cp)
  for(v3=0; v3<lv3; v3++) {
    /* for(lx=0; lx<geo->lL[1]; lx++) */
    /*   for(ly=0; ly<geo->lL[2]; ly++) */
    /*     for(lz=0; lz<geo->lL[3]; lz++){ */
    //v3 = qcd_LEXIC0(lx,ly,lz,geo->lL);
    //qcd_LEXIC(lt,lx,ly,lz,geo->lL);
    /*     lz = v3 % geo->lL[3]; */
    /*     ly = ((v3 - lz)/geo->lL[3]) % geo->lL[2]; */
    /*     lx = ((v3 - lz)/geo->lL[3] - ly) / geo->lL[2]; */
    /*     v =  qcd_LEXIC(lt,lx,ly,lz,geo->lL); */
    v = lt + v3*geo->lL[0];
	
    for(ga =0; ga <4; ga ++) 
      for(gap=0; gap<4; gap++){
       
	for(ctr2=0; ctr2<ctr; ctr2++){          
	  al = cgcg_ind[ctr2][0];
	  be = cgcg_ind[ctr2][1];
	  bep= cgcg_ind[ctr2][2];
	  alp= cgcg_ind[ctr2][3];                              
            
	  for(cc1=0;cc1<6;cc1++){
	    a=qcd_EPS[cc1][0];
	    b=qcd_EPS[cc1][1];
	    c=qcd_EPS[cc1][2];
			
	    for(cc2=0;cc2<6;cc2++){          
	      ap=qcd_EPS[cc2][0];
	      bp=qcd_EPS[cc2][1];
	      cp=qcd_EPS[cc2][2];
         

	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(prop->D[v][al][bep][a][bp],
													qcd_CMUL(prop->D[v][be][gap][b][cp],
														 prop->D[v][ga][alp][c][ap]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2])
						  );
	    
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(prop->D[v][al][gap][a][cp],
													qcd_CMUL(prop->D[v][be][alp][b][ap],
														 prop->D[v][ga][bep][c][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2])
						  );

	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(prop->D[v][al][alp][a][ap],
													qcd_CMUL(prop->D[v][be][bep][b][bp],
														 prop->D[v][ga][gap][c][cp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2])
						  );	    
										    
	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(prop->D[v][al][gap][a][cp],
													qcd_CMUL(prop->D[v][be][bep][b][bp],
														 prop->D[v][ga][alp][c][ap]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2])
						  );
										    
	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(prop->D[v][al][alp][a][ap],
													qcd_CMUL(prop->D[v][be][gap][b][cp],
														 prop->D[v][ga][bep][c][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2])
						  );
										    
	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(prop->D[v][al][bep][a][bp],
													qcd_CMUL(prop->D[v][be][alp][b][ap],
														 prop->D[v][ga][gap][c][cp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2])
						  );	    
	    
	    }//space loop
	  }//color2 loop    
	}//color1 loop
      }//nonvanishing cgcg loop						
  }//nonvanishing projector condition
	
  return (1);
}

//======================================================================

int qcd_deltas_xistar_omegastar_2pt_pr(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[9][4][4],
				       qcd_propagator *propf1, qcd_propagator *propf2, qcd_geometry *geo, qcd_uint_4 lt, qcd_uint_4 elem){
							
  qcd_int_4 ctr2;
  qcd_uint_2 a,b,c,ap,bp,cp,al,be,ga,alp,bep,gap,cc1,cc2;
  qcd_uint_4 lx,ly,lz,v,v3;
  int lv3 = geo->lL[1]*geo->lL[2]*geo->lL[3];
  
#pragma omp parallel for private(lz,ly,lx,v,ga,gap,ctr2,al,be,bep,alp,cc1,a,b,c,cc2,ap,bp,cp)
  for(v3=0; v3<lv3; v3++) {
    /* for(lx=0; lx<geo->lL[1]; lx++) */
    /*   for(ly=0; ly<geo->lL[2]; ly++) */
    /*     for(lz=0; lz<geo->lL[3]; lz++){ */
    //v3 = qcd_LEXIC0(lx,ly,lz,geo->lL);
    //qcd_LEXIC(lt,lx,ly,lz,geo->lL);
/*     lz = v3 % geo->lL[3]; */
/*     ly = ((v3 - lz)/geo->lL[3]) % geo->lL[2]; */
/*     lx = ((v3 - lz)/geo->lL[3] - ly) / geo->lL[2]; */
/*     v =  qcd_LEXIC(lt,lx,ly,lz,geo->lL); */
    v = lt + v3*geo->lL[0];
	
    for(ga =0; ga <4; ga ++) 
      for(gap=0; gap<4; gap++){
	
	for(ctr2=0; ctr2<ctr; ctr2++){          
	  al = cgcg_ind[ctr2][0];
	  be = cgcg_ind[ctr2][1];
	  bep= cgcg_ind[ctr2][2];
	  alp= cgcg_ind[ctr2][3];                              
            
	  for(cc1=0;cc1<6;cc1++){
	    a=qcd_EPS[cc1][0];
	    b=qcd_EPS[cc1][1];
	    c=qcd_EPS[cc1][2];
			
	    for(cc2=0;cc2<6;cc2++){          
	      ap=qcd_EPS[cc2][0];
	      bp=qcd_EPS[cc2][1];
	      cp=qcd_EPS[cc2][2];

	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][be][bep][b][bp],
													qcd_CMUL(propf1->D[v][al][alp][a][ap],
														 propf1->D[v][ga][gap][c][cp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*4.0/3.0)
						  ); //-ok

	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][be][bep][b][bp],
													qcd_CMUL(propf1->D[v][al][gap][a][cp],
														 propf1->D[v][ga][alp][c][ap]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*4.0/3.0)
						  );//-ok
	    
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][be][gap][b][cp],
													qcd_CMUL(propf1->D[v][al][bep][a][bp],
														 propf1->D[v][ga][alp][c][ap]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*2.0/3.0)
						  );//-ok	    
	    
	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][be][gap][b][cp],
													qcd_CMUL(propf1->D[v][al][alp][a][ap],
														 propf1->D[v][ga][bep][c][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*2.0/3.0)
						  );//-ok	    
	    
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][ga][bep][c][bp],
													qcd_CMUL(propf1->D[v][al][gap][a][cp],
														 propf1->D[v][be][alp][b][ap]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*2.0/3.0)
						  );//-ok

	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][ga][bep][c][bp],
													qcd_CMUL(propf1->D[v][al][alp][a][ap],
														 propf1->D[v][be][gap][b][cp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*2.0/3.0)
						  );//-ok	    
	    
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][ga][gap][c][cp],
													qcd_CMUL(propf1->D[v][al][alp][a][ap],
														 propf1->D[v][be][bep][b][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]/3.0)
						  );//-ok	    

	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][ga][gap][c][cp],
													qcd_CMUL(propf1->D[v][al][bep][a][bp],
														 propf1->D[v][be][alp][b][ap]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]/3.0)
						  );//-ok
				 
	    }//color2 loop     
	  }//color1 loop                                                                               
	}//spin loop (cg_cg)
      }//ga gap
  }//space loop                                                                                   

  return (1);
}

//======================================================================

int qcd_xistar_extra11_2pt_pr(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[9][4][4],
			      qcd_propagator *propf1, qcd_propagator *propf2, qcd_geometry *geo, qcd_uint_4 lt, qcd_uint_4 elem){
							
  qcd_int_4 ctr2;
  qcd_uint_2 a,b,c,ap,bp,cp,al,be,ga,alp,bep,gap,cc1,cc2;
  qcd_uint_4 lx,ly,lz,v,v3;
  int lv3 = geo->lL[1]*geo->lL[2]*geo->lL[3];
  
#pragma omp parallel for private(lz,ly,lx,v,ga,gap,ctr2,al,be,bep,alp,cc1,a,b,c,cc2,ap,bp,cp)
  for(v3=0; v3<lv3; v3++) {
    /* for(lx=0; lx<geo->lL[1]; lx++) */
    /*   for(ly=0; ly<geo->lL[2]; ly++) */
    /*     for(lz=0; lz<geo->lL[3]; lz++){ */
    //v3 = qcd_LEXIC0(lx,ly,lz,geo->lL);
    //qcd_LEXIC(lt,lx,ly,lz,geo->lL);
/*     lz = v3 % geo->lL[3]; */
/*     ly = ((v3 - lz)/geo->lL[3]) % geo->lL[2]; */
/*     lx = ((v3 - lz)/geo->lL[3] - ly) / geo->lL[2]; */
/*     v =  qcd_LEXIC(lt,lx,ly,lz,geo->lL); */
    v = lt + v3*geo->lL[0];
	
    for(ga =0; ga <4; ga ++) 
      for(gap=0; gap<4; gap++){
	
	for(ctr2=0; ctr2<ctr; ctr2++){          
	  al = cgcg_ind[ctr2][0];
	  be = cgcg_ind[ctr2][1];
	  bep= cgcg_ind[ctr2][2];
	  alp= cgcg_ind[ctr2][3];                              
            
	  for(cc1=0;cc1<6;cc1++){
	    a=qcd_EPS[cc1][0];
	    b=qcd_EPS[cc1][1];
	    c=qcd_EPS[cc1][2];
			
	    for(cc2=0;cc2<6;cc2++){          
	      ap=qcd_EPS[cc2][0];
	      bp=qcd_EPS[cc2][1];
	      cp=qcd_EPS[cc2][2];

	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][be][bep][b][bp],
													qcd_CMUL(propf1->D[v][al][alp][a][ap],
														 propf1->D[v][ga][gap][c][cp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2])
						  ); //-ok

	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][be][bep][b][bp],
													qcd_CMUL(propf1->D[v][al][gap][a][cp],
														 propf1->D[v][ga][alp][c][ap]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2])
						  );//-ok
				 
	    }//color2 loop     
	  }//color1 loop                                                                               
	}//spin loop (cg_cg)
      }//ga gap
  }//space loop                                                                                   

  return (1);
}

//======================================================================

int qcd_xistar_extra12_2pt_pr(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[9][4][4],
			      qcd_propagator *propf1, qcd_propagator *propf2, qcd_geometry *geo, qcd_uint_4 lt, qcd_uint_4 elem){
							
  qcd_int_4 ctr2;
  qcd_uint_2 a,b,c,ap,bp,cp,al,be,ga,alp,bep,gap,cc1,cc2;
  qcd_uint_4 lx,ly,lz,v,v3;
  int lv3 = geo->lL[1]*geo->lL[2]*geo->lL[3];
  
#pragma omp parallel for private(lz,ly,lx,v,ga,gap,ctr2,al,be,bep,alp,cc1,a,b,c,cc2,ap,bp,cp)
  for(v3=0; v3<lv3; v3++) {
    /* for(lx=0; lx<geo->lL[1]; lx++) */
    /*   for(ly=0; ly<geo->lL[2]; ly++) */
    /*     for(lz=0; lz<geo->lL[3]; lz++){ */
    //v3 = qcd_LEXIC0(lx,ly,lz,geo->lL);
    //qcd_LEXIC(lt,lx,ly,lz,geo->lL);
/*     lz = v3 % geo->lL[3]; */
/*     ly = ((v3 - lz)/geo->lL[3]) % geo->lL[2]; */
/*     lx = ((v3 - lz)/geo->lL[3] - ly) / geo->lL[2]; */
/*     v =  qcd_LEXIC(lt,lx,ly,lz,geo->lL); */
    v = lt + v3*geo->lL[0];
	
    for(ga =0; ga <4; ga ++) 
      for(gap=0; gap<4; gap++){
	
	for(ctr2=0; ctr2<ctr; ctr2++){          
	  al = cgcg_ind[ctr2][0];
	  be = cgcg_ind[ctr2][1];
	  bep= cgcg_ind[ctr2][2];
	  alp= cgcg_ind[ctr2][3];                              
            
	  for(cc1=0;cc1<6;cc1++){
	    a=qcd_EPS[cc1][0];
	    b=qcd_EPS[cc1][1];
	    c=qcd_EPS[cc1][2];
			
	    for(cc2=0;cc2<6;cc2++){          
	      ap=qcd_EPS[cc2][0];
	      bp=qcd_EPS[cc2][1];
	      cp=qcd_EPS[cc2][2];

	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][be][gap][b][cp],
													qcd_CMUL(propf1->D[v][al][bep][a][bp],
														 propf1->D[v][ga][alp][c][ap]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2])
						  );//-ok	    
	    
	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][be][gap][b][cp],
													qcd_CMUL(propf1->D[v][al][alp][a][ap],
														 propf1->D[v][ga][bep][c][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2])
						  );//-ok	    
				 
	    }//color2 loop     
	  }//color1 loop                                                                               
	}//spin loop (cg_cg)
      }//ga gap
  }//space loop                                                                                   

  return (1);
}

//======================================================================

int qcd_xistar_extra21_2pt_pr(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[9][4][4],
			      qcd_propagator *propf1, qcd_propagator *propf2, qcd_geometry *geo, qcd_uint_4 lt, qcd_uint_4 elem){
							
  qcd_int_4 ctr2;
  qcd_uint_2 a,b,c,ap,bp,cp,al,be,ga,alp,bep,gap,cc1,cc2;
  qcd_uint_4 lx,ly,lz,v,v3;
  int lv3 = geo->lL[1]*geo->lL[2]*geo->lL[3];
  
#pragma omp parallel for private(lz,ly,lx,v,ga,gap,ctr2,al,be,bep,alp,cc1,a,b,c,cc2,ap,bp,cp)
  for(v3=0; v3<lv3; v3++) {
    /* for(lx=0; lx<geo->lL[1]; lx++) */
    /*   for(ly=0; ly<geo->lL[2]; ly++) */
    /*     for(lz=0; lz<geo->lL[3]; lz++){ */
    //v3 = qcd_LEXIC0(lx,ly,lz,geo->lL);
    //qcd_LEXIC(lt,lx,ly,lz,geo->lL);
/*     lz = v3 % geo->lL[3]; */
/*     ly = ((v3 - lz)/geo->lL[3]) % geo->lL[2]; */
/*     lx = ((v3 - lz)/geo->lL[3] - ly) / geo->lL[2]; */
/*     v =  qcd_LEXIC(lt,lx,ly,lz,geo->lL); */
    v = lt + v3*geo->lL[0];
	
    for(ga =0; ga <4; ga ++) 
      for(gap=0; gap<4; gap++){
	
	for(ctr2=0; ctr2<ctr; ctr2++){          
	  al = cgcg_ind[ctr2][0];
	  be = cgcg_ind[ctr2][1];
	  bep= cgcg_ind[ctr2][2];
	  alp= cgcg_ind[ctr2][3];                              
            
	  for(cc1=0;cc1<6;cc1++){
	    a=qcd_EPS[cc1][0];
	    b=qcd_EPS[cc1][1];
	    c=qcd_EPS[cc1][2];
			
	    for(cc2=0;cc2<6;cc2++){          
	      ap=qcd_EPS[cc2][0];
	      bp=qcd_EPS[cc2][1];
	      cp=qcd_EPS[cc2][2];

	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][ga][bep][c][bp],
													qcd_CMUL(propf1->D[v][al][gap][a][cp],
														 propf1->D[v][be][alp][b][ap]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2])
						  );//-ok

	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][ga][bep][c][bp],
													qcd_CMUL(propf1->D[v][al][alp][a][ap],
														 propf1->D[v][be][gap][b][cp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2])
						  );//-ok	    
				 
	    }//color2 loop     
	  }//color1 loop                                                                               
	}//spin loop (cg_cg)
      }//ga gap
  }//space loop                                                                                   

  return (1);
}

//======================================================================

int qcd_xistar_extra22_2pt_pr(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[9][4][4],
			      qcd_propagator *propf1, qcd_propagator *propf2, qcd_geometry *geo, qcd_uint_4 lt, qcd_uint_4 elem){
							
  qcd_int_4 ctr2;
  qcd_uint_2 a,b,c,ap,bp,cp,al,be,ga,alp,bep,gap,cc1,cc2;
  qcd_uint_4 lx,ly,lz,v,v3;
  int lv3 = geo->lL[1]*geo->lL[2]*geo->lL[3];
  
#pragma omp parallel for private(lz,ly,lx,v,ga,gap,ctr2,al,be,bep,alp,cc1,a,b,c,cc2,ap,bp,cp)
  for(v3=0; v3<lv3; v3++) {
    /* for(lx=0; lx<geo->lL[1]; lx++) */
    /*   for(ly=0; ly<geo->lL[2]; ly++) */
    /*     for(lz=0; lz<geo->lL[3]; lz++){ */
    //v3 = qcd_LEXIC0(lx,ly,lz,geo->lL);
    //qcd_LEXIC(lt,lx,ly,lz,geo->lL);
/*     lz = v3 % geo->lL[3]; */
/*     ly = ((v3 - lz)/geo->lL[3]) % geo->lL[2]; */
/*     lx = ((v3 - lz)/geo->lL[3] - ly) / geo->lL[2]; */
/*     v =  qcd_LEXIC(lt,lx,ly,lz,geo->lL); */
    v = lt + v3*geo->lL[0];
	
    for(ga =0; ga <4; ga ++) 
      for(gap=0; gap<4; gap++){
	
	for(ctr2=0; ctr2<ctr; ctr2++){          
	  al = cgcg_ind[ctr2][0];
	  be = cgcg_ind[ctr2][1];
	  bep= cgcg_ind[ctr2][2];
	  alp= cgcg_ind[ctr2][3];                              
            
	  for(cc1=0;cc1<6;cc1++){
	    a=qcd_EPS[cc1][0];
	    b=qcd_EPS[cc1][1];
	    c=qcd_EPS[cc1][2];
			
	    for(cc2=0;cc2<6;cc2++){          
	      ap=qcd_EPS[cc2][0];
	      bp=qcd_EPS[cc2][1];
	      cp=qcd_EPS[cc2][2];

	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][ga][gap][c][cp],
													qcd_CMUL(propf1->D[v][al][alp][a][ap],
														 propf1->D[v][be][bep][b][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2])
						  );//-ok	    

	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][ga][gap][c][cp],
													qcd_CMUL(propf1->D[v][al][bep][a][bp],
														 propf1->D[v][be][alp][b][ap]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2])
						  );//-ok
				 
	    }//color2 loop     
	  }//color1 loop                                                                               
	}//spin loop (cg_cg)
      }//ga gap
  }//space loop                                                                                   

  return (1);
}

//======================================================================

int qcd_f1f2f1_2pt_pr(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[9][4][4],
		      qcd_propagator *propf1, qcd_propagator *propf2, qcd_geometry *geo, qcd_uint_4 lt,qcd_uint_4 elem){
							
  qcd_int_4 ctr2;
  qcd_uint_2 a,b,c,ap,bp,cp,al,be,ga,alp,bep,gap,cc1,cc2;
  qcd_uint_4 lx,ly,lz,v,v3;
  int lv3 = geo->lL[1]*geo->lL[2]*geo->lL[3];

#pragma omp parallel for private(lz,ly,lx,v,ga,gap,ctr2,al,be,bep,alp,cc1,a,b,c,cc2,ap,bp,cp)
  for(v3=0; v3<lv3; v3++) {
    /* for(lx=0; lx<geo->lL[1]; lx++) */
    /*   for(ly=0; ly<geo->lL[2]; ly++) */
    /*     for(lz=0; lz<geo->lL[3]; lz++){ */
    //v3 = qcd_LEXIC0(lx,ly,lz,geo->lL);
    //qcd_LEXIC(lt,lx,ly,lz,geo->lL);
/*     lz = v3 % geo->lL[3]; */
/*     ly = ((v3 - lz)/geo->lL[3]) % geo->lL[2]; */
/*     lx = ((v3 - lz)/geo->lL[3] - ly) / geo->lL[2]; */
/*     v =  qcd_LEXIC(lt,lx,ly,lz,geo->lL); */
    v = lt + v3*geo->lL[0];
    
    for(ga =0; ga <4; ga ++) 
      for(gap=0; gap<4; gap++){         
	
	for(ctr2=0; ctr2<ctr; ctr2++){          
	  al = cgcg_ind[ctr2][0];
	  be = cgcg_ind[ctr2][1];
	  bep= cgcg_ind[ctr2][2];
	  alp= cgcg_ind[ctr2][3];                              
            
	  for(cc1=0;cc1<6;cc1++){
	    a=qcd_EPS[cc1][0];
	    b=qcd_EPS[cc1][1];
	    c=qcd_EPS[cc1][2];
			
	    for(cc2=0;cc2<6;cc2++){          
	      ap=qcd_EPS[cc2][0];
	      bp=qcd_EPS[cc2][1];
	      cp=qcd_EPS[cc2][2];
                
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][be][bep][b][bp],
													qcd_CMUL(propf1->D[v][al][alp][a][ap],
														 propf1->D[v][ga][gap][c][cp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2])
						  );
										    
	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][be][bep][b][bp],
													qcd_CMUL(propf1->D[v][al][gap][a][cp],
														 propf1->D[v][ga][alp][c][ap]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2])
						  );
	    }//space loop
	  }//color2 loop    
	}//color1 loop
      }//nonvanishing cgcg loop						
  }// ga gap indices
	
  return (1);
}

//======================================================================

int qcd_f1f2f3_2pt_pr(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[9][4][4],
		      qcd_propagator *propf1, qcd_propagator *propf2, qcd_propagator *propf3, qcd_geometry *geo, qcd_uint_4 lt, qcd_uint_4 elem){
							
  qcd_int_4 ctr2;
  qcd_uint_2 a,b,c,ap,bp,cp,al,be,ga,alp,bep,gap,cc1,cc2;
  qcd_uint_4 lx,ly,lz,v,v3;
  int lv3 = geo->lL[1]*geo->lL[2]*geo->lL[3];

#pragma omp parallel for private(lz,ly,lx,v,ga,gap,ctr2,al,be,bep,alp,cc1,a,b,c,cc2,ap,bp,cp)
  for(v3=0; v3<lv3; v3++) {
    /* for(lx=0; lx<geo->lL[1]; lx++) */
    /*   for(ly=0; ly<geo->lL[2]; ly++) */
    /*     for(lz=0; lz<geo->lL[3]; lz++){ */
    //v3 = qcd_LEXIC0(lx,ly,lz,geo->lL);
    //qcd_LEXIC(lt,lx,ly,lz,geo->lL);
/*     lz = v3 % geo->lL[3]; */
/*     ly = ((v3 - lz)/geo->lL[3]) % geo->lL[2]; */
/*     lx = ((v3 - lz)/geo->lL[3] - ly) / geo->lL[2]; */
/*     v =  qcd_LEXIC(lt,lx,ly,lz,geo->lL); */
    v = lt + v3*geo->lL[0];  	
  	
  	
    for(ga =0; ga <4; ga ++) 
      for(gap=0; gap<4; gap++){        
	
	for(ctr2=0; ctr2<ctr; ctr2++){          
	  al = cgcg_ind[ctr2][0];
	  be = cgcg_ind[ctr2][1];
	  bep= cgcg_ind[ctr2][2];
	  alp= cgcg_ind[ctr2][3];                              
            
	  for(cc1=0;cc1<6;cc1++){
	    a=qcd_EPS[cc1][0];
	    b=qcd_EPS[cc1][1];
	    c=qcd_EPS[cc1][2];
			
	    for(cc2=0;cc2<6;cc2++){          
	      ap=qcd_EPS[cc2][0];
	      bp=qcd_EPS[cc2][1];
	      cp=qcd_EPS[cc2][2];
                
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][al][alp][a][ap],
													qcd_CMUL(propf2->D[v][be][bep][b][bp],
														 propf3->D[v][ga][gap][c][cp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2])
						  );	    
	    
	    }//space loop
	  }//color2 loop    
	}//color1 loop
      }//nonvanishing cgcg loop						
  }//ga gap indices
	
  return (1);
}

//======================================================================

int qcd_sigmas4_2pt_pr(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[9][4][4],
		       qcd_propagator *propf1, qcd_propagator *propf2, qcd_geometry *geo, qcd_uint_4 lt, qcd_uint_4 elem){
							
  qcd_int_4 ctr2;
  qcd_uint_2 a,b,c,ap,bp,cp,al,be,ga,alp,bep,gap,cc1,cc2;
  qcd_uint_4 lx,ly,lz,v,v3;
  int lv3 = geo->lL[1]*geo->lL[2]*geo->lL[3];
  
#pragma omp parallel for private(lz,ly,lx,v,ga,gap,ctr2,al,be,bep,alp,cc1,a,b,c,cc2,ap,bp,cp)
  for(v3=0; v3<lv3; v3++) {
    /* for(lx=0; lx<geo->lL[1]; lx++) */
    /*   for(ly=0; ly<geo->lL[2]; ly++) */
    /*     for(lz=0; lz<geo->lL[3]; lz++){ */
    //v3 = qcd_LEXIC0(lx,ly,lz,geo->lL);
    //qcd_LEXIC(lt,lx,ly,lz,geo->lL);
/*     lz = v3 % geo->lL[3]; */
/*     ly = ((v3 - lz)/geo->lL[3]) % geo->lL[2]; */
/*     lx = ((v3 - lz)/geo->lL[3] - ly) / geo->lL[2]; */
/*     v =  qcd_LEXIC(lt,lx,ly,lz,geo->lL);	 */
    v = lt + v3*geo->lL[0];

    for(ga =0; ga <4; ga ++) 
      for(gap=0; gap<4; gap++){
	
	for(ctr2=0; ctr2<ctr; ctr2++){          
	  al = cgcg_ind[ctr2][0];
	  be = cgcg_ind[ctr2][1];
	  bep= cgcg_ind[ctr2][2];
	  alp= cgcg_ind[ctr2][3];                              
            
	  for(cc1=0;cc1<6;cc1++){
	    a=qcd_EPS[cc1][0];
	    b=qcd_EPS[cc1][1];
	    c=qcd_EPS[cc1][2];
			
	    for(cc2=0;cc2<6;cc2++){          
	      ap=qcd_EPS[cc2][0];
	      bp=qcd_EPS[cc2][1];
	      cp=qcd_EPS[cc2][2];
  
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][al][alp][a][ap],
													qcd_CMUL(propf1->D[v][be][bep][b][bp],
														 propf1->D[v][ga][gap][c][cp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*4.0/3.0)
						  );//-ok

	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][al][alp][a][ap],
													qcd_CMUL(propf1->D[v][be][gap][b][cp],
														 propf1->D[v][ga][bep][c][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*4.0/3.0)
						  );//-ok
	    
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][al][gap][a][cp],
													qcd_CMUL(propf1->D[v][be][alp][b][ap],
														 propf1->D[v][ga][bep][c][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*2.0/3.0)
						  );//-ok	    
	    
	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][al][gap][a][cp],
													qcd_CMUL(propf1->D[v][be][bep][b][bp],
														 propf1->D[v][ga][alp][c][ap]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*2.0/3.0)
						  );//-ok	    
	    
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][ga][alp][c][ap],
													qcd_CMUL(propf1->D[v][al][bep][a][bp],
														 propf1->D[v][be][gap][b][cp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*2.0/3.0)
						  );//-ok

	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][ga][alp][c][ap],
													qcd_CMUL(propf1->D[v][al][gap][a][cp],
														 propf1->D[v][be][bep][b][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*2.0/3.0)
						  );//-ok	    
	    
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][ga][gap][c][cp],
													qcd_CMUL(propf1->D[v][al][alp][a][ap],
														 propf1->D[v][be][bep][b][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]/3.0)
						  );//-ok	    

	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][ga][gap][c][cp],
													qcd_CMUL(propf1->D[v][al][bep][a][bp],
														 propf1->D[v][be][alp][b][ap]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]/3.0)
						  );//-ok
	    }//color2 loop                                            
	  }//color1 loop                                                                 
	}//spin loop (cg_cg)                                                          
      }//ga gap indices
  }//space loop                                                                                                                                                       
	
  return (1);
}

//======================================================================

int qcd_sigmas2_xistar_2pt_pr(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[9][4][4],
			      qcd_propagator *propf1, qcd_propagator *propf2, qcd_propagator *propf3, qcd_geometry *geo, qcd_uint_4 lt, qcd_uint_4 elem, qcd_uint_4 sflag){
							
  qcd_int_4 ctr2;
  qcd_uint_2 a,b,c,ap,bp,cp,al,be,ga,alp,bep,gap,cc1,cc2;
  qcd_uint_4 lx,ly,lz,v,v3;
  qcd_real_8 fact;
  int lv3 = geo->lL[1]*geo->lL[2]*geo->lL[3];
      
    
  if (sflag) fact = 2.0/3.0;
  else fact = 1.0/3.0;
	
#pragma omp parallel for private(lz,ly,lx,v,ga,gap,ctr2,al,be,bep,alp,cc1,a,b,c,cc2,ap,bp,cp)
  for(v3=0; v3<lv3; v3++) {
    /* for(lx=0; lx<geo->lL[1]; lx++) */
    /*   for(ly=0; ly<geo->lL[2]; ly++) */
    /*     for(lz=0; lz<geo->lL[3]; lz++){ */
    //v3 = qcd_LEXIC0(lx,ly,lz,geo->lL);
    //qcd_LEXIC(lt,lx,ly,lz,geo->lL);
/*     lz = v3 % geo->lL[3]; */
/*     ly = ((v3 - lz)/geo->lL[3]) % geo->lL[2]; */
/*     lx = ((v3 - lz)/geo->lL[3] - ly) / geo->lL[2]; */
/*     v =  qcd_LEXIC(lt,lx,ly,lz,geo->lL); */
    v = lt + v3*geo->lL[0];

    for(ga =0; ga <4; ga ++) 
      for(gap=0; gap<4; gap++){         
	
	for(ctr2=0; ctr2<ctr; ctr2++){          
	  al = cgcg_ind[ctr2][0];
	  be = cgcg_ind[ctr2][1];
	  bep= cgcg_ind[ctr2][2];
	  alp= cgcg_ind[ctr2][3];                              
            
	  for(cc1=0;cc1<6;cc1++){
	    a=qcd_EPS[cc1][0];
	    b=qcd_EPS[cc1][1];
	    c=qcd_EPS[cc1][2];
			
	    for(cc2=0;cc2<6;cc2++){          
	      ap=qcd_EPS[cc2][0];
	      bp=qcd_EPS[cc2][1];
	      cp=qcd_EPS[cc2][2];

	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][al][alp][a][ap],
													qcd_CMUL(propf2->D[v][be][bep][b][bp],
														 propf3->D[v][ga][gap][c][cp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*fact)
						  );//-ok

	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][al][gap][a][cp],
													qcd_CMUL(propf2->D[v][be][alp][b][ap],
														 propf3->D[v][ga][bep][c][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*fact)
						  );//-ok
	    
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][al][bep][a][bp],
													qcd_CMUL(propf2->D[v][be][gap][b][cp],
														 propf3->D[v][ga][alp][c][ap]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*fact)
						  );//-ok	    
	    
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][ga][alp][c][ap],
													qcd_CMUL(propf2->D[v][al][bep][a][bp],
														 propf3->D[v][be][gap][b][cp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*fact)
						  );//-ok	    
	    
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][ga][gap][c][cp],
													qcd_CMUL(propf2->D[v][al][alp][a][ap],
														 propf3->D[v][be][bep][b][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*fact)
						  );//-ok

	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][ga][bep][c][bp],
													qcd_CMUL(propf2->D[v][al][gap][a][cp],
														 propf3->D[v][be][alp][b][ap]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*fact)
						  );//-ok	    
	    
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][be][alp][b][ap],
													qcd_CMUL(propf2->D[v][ga][bep][c][bp],
														 propf3->D[v][al][gap][a][cp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*fact)
						  );//-ok	    

	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][be][gap][b][cp],
													qcd_CMUL(propf2->D[v][ga][alp][c][ap],
														 propf3->D[v][al][bep][a][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*fact)
						  );//-ok
	    
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][be][bep][b][bp],
													qcd_CMUL(propf2->D[v][ga][gap][c][cp],
														 propf3->D[v][al][alp][a][ap]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*fact)
						  );//-ok
				 	    
	    }//color2 loop        
	  }//color1 loop                                                                                                   
	}//spin loop (cg_cg)                                           
      }//ga gap
  }//space loop                                                                     
	
  return (1);
}

//======================================================================
//======================================================================
//======================================================================

void qcd_contractions2pt_pr(qcd_uint_4 particle_id, qcd_complex_16 *block12[4][4],qcd_complex_16 *block32[4][4],qcd_complex_16 *blocknp[4][4],
			    qcd_propagator *uprop,qcd_propagator *dprop,qcd_propagator *sprop,qcd_propagator *cprop,
			    qcd_geometry *geo, qcd_uint_4 lt){

  qcd_int_4 ctr1,ctr2,ctr3,ctr12,ctr21,ctr13,ctr31,ctr23,ctr32,v3;
  qcd_uint_2 al,be,alp,bep,check;
  qcd_complex_16 C,C1,C2,C3,C21,C12,C31,C13,C23,C32;
  qcd_complex_16 gamma12[4][4],gamma13[4][4],gamma23[4][4];

  qcd_int_2 cg1cg1_ind[16*16][4],cg2cg2_ind[16*16][4],cg3cg3_ind[16*16][4];
  qcd_int_2 cg1cg2_ind[16*16][4],cg2cg1_ind[16*16][4],cg3cg2_ind[16*16][4],cg2cg3_ind[16*16][4];
  qcd_int_2 cg1cg3_ind[16*16][4],cg3cg1_ind[16*16][4]; 
  
  qcd_complex_16 cg1cg1_val[16*16],cg2cg2_val[16*16],cg3cg3_val[16*16];
  qcd_complex_16 cg1cg2_val[16*16],cg2cg1_val[16*16],cg3cg2_val[16*16],cg2cg3_val[16*16];
  qcd_complex_16 cg1cg3_val[16*16],cg3cg1_val[16*16];
  
  qcd_complex_16 *block_pr[9][4][4];
  
  qcd_uint_4 i,j,k,det,sflag;
  
  for(i=0;i<9;i++)
    for(j=0;j<4;j++)
      for(k=0;k<4;k++){	
	block_pr[i][j][k] = (qcd_complex_16*) malloc(geo->lV3*sizeof(qcd_complex_16));
			
	if(block_pr[i][j][k]==NULL){
	  printf("Block %d %d %d not properly initialized\n",i,j,k);
	  exit(EXIT_FAILURE);
	}
			
	for(v3=0;v3<(geo->lV3);v3++){
	  block_pr[i][j][k][v3] = (qcd_complex_16) {0,0};			
	}
      }
    
  //-- Non zero elements of all combinations of gamma_{1,2,3} matrices
  ctr1 = 0; ctr2 = 0; ctr3 = 0;
  ctr13 = 0; ctr31 = 0; ctr12 = 0; ctr21 = 0; ctr23 = 0; ctr32 = 0;
    
  for(al=0;al<4;al++)
    for(be=0;be<4;be++)
      for(alp=0;alp<4;alp++)
	for(bep=0;bep<4;bep++)
	  {  
	    C1  = qcd_CMUL(qcd_CGAMMA[1][al][be],qcd_BAR_CGAMMA[1][bep][alp]);
	    C12 = qcd_CMUL(qcd_CGAMMA[1][al][be],qcd_BAR_CGAMMA[2][bep][alp]);    
	    C13 = qcd_CMUL(qcd_CGAMMA[1][al][be],qcd_BAR_CGAMMA[3][bep][alp]);
        
	    C21 = qcd_CMUL(qcd_CGAMMA[2][al][be],qcd_BAR_CGAMMA[1][bep][alp]);    
	    C2  = qcd_CMUL(qcd_CGAMMA[2][al][be],qcd_BAR_CGAMMA[2][bep][alp]);
	    C23 = qcd_CMUL(qcd_CGAMMA[2][al][be],qcd_BAR_CGAMMA[3][bep][alp]);    
    
	    C31 = qcd_CMUL(qcd_CGAMMA[3][al][be],qcd_BAR_CGAMMA[1][bep][alp]);
	    C32 = qcd_CMUL(qcd_CGAMMA[3][al][be],qcd_BAR_CGAMMA[2][bep][alp]);   
	    C3  = qcd_CMUL(qcd_CGAMMA[3][al][be],qcd_BAR_CGAMMA[3][bep][alp]);    
        
	    if(qcd_NORM(C1)>1e-3)
	      {
		cg1cg1_val[ctr1].re = C1.re;
		cg1cg1_val[ctr1].im = C1.im;
		cg1cg1_ind[ctr1][0] = al;
		cg1cg1_ind[ctr1][1] = be;
		cg1cg1_ind[ctr1][2] = bep;
		cg1cg1_ind[ctr1][3] = alp;                                                            
		ctr1++;
	      }
	    if(qcd_NORM(C2)>1e-3)
	      {
		cg2cg2_val[ctr2].re = C2.re;
		cg2cg2_val[ctr2].im = C2.im;
		cg2cg2_ind[ctr2][0] = al;
		cg2cg2_ind[ctr2][1] = be;
		cg2cg2_ind[ctr2][2] = bep;
		cg2cg2_ind[ctr2][3] = alp;                                                            
		ctr2++;
	      }
	    if(qcd_NORM(C3)>1e-3)
	      {
		cg3cg3_val[ctr3].re = C3.re;
		cg3cg3_val[ctr3].im = C3.im;
		cg3cg3_ind[ctr3][0] = al;
		cg3cg3_ind[ctr3][1] = be;
		cg3cg3_ind[ctr3][2] = bep;
		cg3cg3_ind[ctr3][3] = alp;                                                            
		ctr3++;
	      }
	    if(qcd_NORM(C23)>2e-3)
	      {
		cg2cg3_val[ctr23].re = C23.re;
		cg2cg3_val[ctr23].im = C23.im;
		cg2cg3_ind[ctr23][0] = al;
		cg2cg3_ind[ctr23][1] = be;
		cg2cg3_ind[ctr23][2] = bep;
		cg2cg3_ind[ctr23][3] = alp;                                                            
		ctr23++;
	      }
	    if(qcd_NORM(C32)>2e-3)
	      {
		cg3cg2_val[ctr32].re = C32.re;
		cg3cg2_val[ctr32].im = C32.im;
		cg3cg2_ind[ctr32][0] = al;
		cg3cg2_ind[ctr32][1] = be;
		cg3cg2_ind[ctr32][2] = bep;
		cg3cg2_ind[ctr32][3] = alp;                                                            
		ctr32++;
	      }    
	    if(qcd_NORM(C13)>1e-3)
	      {
		cg1cg3_val[ctr13].re = C13.re;
		cg1cg3_val[ctr13].im = C13.im;
		cg1cg3_ind[ctr13][0] = al;
		cg1cg3_ind[ctr13][1] = be;
		cg1cg3_ind[ctr13][2] = bep;
		cg1cg3_ind[ctr13][3] = alp;                                                            
		ctr13++;
	      }
	    if(qcd_NORM(C31)>1e-3)
	      {
		cg3cg1_val[ctr31].re = C31.re;
		cg3cg1_val[ctr31].im = C31.im;
		cg3cg1_ind[ctr31][0] = al;
		cg3cg1_ind[ctr31][1] = be;
		cg3cg1_ind[ctr31][2] = bep;
		cg3cg1_ind[ctr31][3] = alp;                                                            
		ctr31++;
	      }
	    if(qcd_NORM(C12)>1e-2)
	      {
		cg1cg2_val[ctr12].re = C12.re;
		cg1cg2_val[ctr12].im = C12.im;
		cg1cg2_ind[ctr12][0] = al;
		cg1cg2_ind[ctr12][1] = be;
		cg1cg2_ind[ctr12][2] = bep;
		cg1cg2_ind[ctr12][3] = alp;                                                            
		ctr12++;
	      }
	    if(qcd_NORM(C21)>1e-2)
	      {
		cg2cg1_val[ctr21].re = C21.re;
		cg2cg1_val[ctr21].im = C21.im;
		cg2cg1_ind[ctr21][0] = al;
		cg2cg1_ind[ctr21][1] = be;
		cg2cg1_ind[ctr21][2] = bep;
		cg2cg1_ind[ctr21][3] = alp;                                                            
		ctr21++;
	      }                  
	  }  

  //-- Multiplication of gamma_1,gamma_2,gamma_3 between each other

  for(i=0;i<4;i++){
    for(j=0;j<4;j++){
      gamma12[i][j] = (qcd_complex_16) {0,0};
      gamma13[i][j] = (qcd_complex_16) {0,0};
      gamma23[i][j] = (qcd_complex_16) {0,0};
		
      for(k=0;k<4;k++){
	gamma12[i][j] = qcd_CADD(gamma12[i][j],qcd_CMUL(qcd_GAMMA[1][i][k],qcd_GAMMA[2][k][j]));
	gamma13[i][j] = qcd_CADD(gamma13[i][j],qcd_CMUL(qcd_GAMMA[1][i][k],qcd_GAMMA[3][k][j]));
	gamma23[i][j] = qcd_CADD(gamma23[i][j],qcd_CMUL(qcd_GAMMA[2][i][k],qcd_GAMMA[3][k][j]));								
      }
    }
  }	
	
  
  switch(particle_id){
  case 9:
    check =  qcd_f1f1f1_2pt_pr(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, uprop, geo, lt,0);
    check += qcd_f1f1f1_2pt_pr(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, uprop, geo, lt,1);  
    check += qcd_f1f1f1_2pt_pr(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, uprop, geo, lt,2);
    check += qcd_f1f1f1_2pt_pr(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, uprop, geo, lt,3);
    check += qcd_f1f1f1_2pt_pr(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, uprop, geo, lt,4);  
    check += qcd_f1f1f1_2pt_pr(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, uprop, geo, lt,5);
    check += qcd_f1f1f1_2pt_pr(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, uprop, geo, lt,6);
    check += qcd_f1f1f1_2pt_pr(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, uprop, geo, lt,7);  
    check += qcd_f1f1f1_2pt_pr(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, uprop, geo, lt,8); 	  	  
    break;

  case 10:
    check =  qcd_deltas_xistar_omegastar_2pt_pr(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, uprop, dprop, geo, lt,0);
    check += qcd_deltas_xistar_omegastar_2pt_pr(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, uprop, dprop, geo, lt,1);  
    check += qcd_deltas_xistar_omegastar_2pt_pr(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, uprop, dprop, geo, lt,2);
    check += qcd_deltas_xistar_omegastar_2pt_pr(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, uprop, dprop, geo, lt,3);
    check += qcd_deltas_xistar_omegastar_2pt_pr(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, uprop, dprop, geo, lt,4);  
    check += qcd_deltas_xistar_omegastar_2pt_pr(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, uprop, dprop, geo, lt,5);
    check += qcd_deltas_xistar_omegastar_2pt_pr(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, uprop, dprop, geo, lt,6);
    check += qcd_deltas_xistar_omegastar_2pt_pr(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, uprop, dprop, geo, lt,7);  
    check += qcd_deltas_xistar_omegastar_2pt_pr(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, uprop, dprop, geo, lt,8); 
    break;
	   
  case 11:
    check =  qcd_deltas_xistar_omegastar_2pt_pr(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, dprop, uprop, geo, lt,0);
    check += qcd_deltas_xistar_omegastar_2pt_pr(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, dprop, uprop, geo, lt,1);  
    check += qcd_deltas_xistar_omegastar_2pt_pr(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, dprop, uprop, geo, lt,2);
    check += qcd_deltas_xistar_omegastar_2pt_pr(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, dprop, uprop, geo, lt,3);
    check += qcd_deltas_xistar_omegastar_2pt_pr(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, dprop, uprop, geo, lt,4);  
    check += qcd_deltas_xistar_omegastar_2pt_pr(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, dprop, uprop, geo, lt,5);
    check += qcd_deltas_xistar_omegastar_2pt_pr(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, dprop, uprop, geo, lt,6);
    check += qcd_deltas_xistar_omegastar_2pt_pr(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, dprop, uprop, geo, lt,7);  
    check += qcd_deltas_xistar_omegastar_2pt_pr(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, dprop, uprop, geo, lt,8);
    break;

  case 12:
    check =  qcd_f1f1f1_2pt_pr(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, dprop, geo, lt,0);
    check += qcd_f1f1f1_2pt_pr(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, dprop, geo, lt,1);  
    check += qcd_f1f1f1_2pt_pr(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, dprop, geo, lt,2);
 
    check += qcd_f1f1f1_2pt_pr(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, dprop, geo, lt,3);
    check += qcd_f1f1f1_2pt_pr(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, dprop, geo, lt,4);  
    check += qcd_f1f1f1_2pt_pr(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, dprop, geo, lt,5);

    check += qcd_f1f1f1_2pt_pr(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, dprop, geo, lt,6);
    check += qcd_f1f1f1_2pt_pr(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, dprop, geo, lt,7);  
    check += qcd_f1f1f1_2pt_pr(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, dprop, geo, lt,8); 
    break;
	  
  case 13:
    check =  qcd_sigmas4_2pt_pr(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, uprop, sprop, geo, lt,0);
    check += qcd_sigmas4_2pt_pr(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, uprop, sprop, geo, lt,1);  
    check += qcd_sigmas4_2pt_pr(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, uprop, sprop, geo, lt,2);
    check += qcd_sigmas4_2pt_pr(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, uprop, sprop, geo, lt,3);
    check += qcd_sigmas4_2pt_pr(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, uprop, sprop, geo, lt,4);  
    check += qcd_sigmas4_2pt_pr(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, uprop, sprop, geo, lt,5);
    check += qcd_sigmas4_2pt_pr(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, uprop, sprop, geo, lt,6);
    check += qcd_sigmas4_2pt_pr(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, uprop, sprop, geo, lt,7);  
    check += qcd_sigmas4_2pt_pr(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, uprop, sprop, geo, lt,8); 
    break;
	  
  case 14:
    check =  qcd_sigmas2_xistar_2pt_pr(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, uprop, dprop, sprop, geo, lt,0,1);
    check += qcd_sigmas2_xistar_2pt_pr(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, uprop, dprop, sprop, geo, lt,1,1);  
    check += qcd_sigmas2_xistar_2pt_pr(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, uprop, dprop, sprop, geo, lt,2,1);
    check += qcd_sigmas2_xistar_2pt_pr(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, uprop, dprop, sprop, geo, lt,3,1);
    check += qcd_sigmas2_xistar_2pt_pr(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, uprop, dprop, sprop, geo, lt,4,1);  
    check += qcd_sigmas2_xistar_2pt_pr(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, uprop, dprop, sprop, geo, lt,5,1);
    check += qcd_sigmas2_xistar_2pt_pr(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, uprop, dprop, sprop, geo, lt,6,1);
    check += qcd_sigmas2_xistar_2pt_pr(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, uprop, dprop, sprop, geo, lt,7,1);  
    check += qcd_sigmas2_xistar_2pt_pr(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, uprop, dprop, sprop, geo, lt,8,1); 
    break;
 	  
  case 15:
    check =  qcd_sigmas4_2pt_pr(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, dprop, sprop, geo, lt,0);
    check += qcd_sigmas4_2pt_pr(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, dprop, sprop, geo, lt,1);  
    check += qcd_sigmas4_2pt_pr(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, dprop, sprop, geo, lt,2);
    check += qcd_sigmas4_2pt_pr(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, dprop, sprop, geo, lt,3);
    check += qcd_sigmas4_2pt_pr(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, dprop, sprop, geo, lt,4);  
    check += qcd_sigmas4_2pt_pr(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, dprop, sprop, geo, lt,5);
    check += qcd_sigmas4_2pt_pr(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, dprop, sprop, geo, lt,6);
    check += qcd_sigmas4_2pt_pr(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, dprop, sprop, geo, lt,7);  
    check += qcd_sigmas4_2pt_pr(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, dprop, sprop, geo, lt,8); 
    break;
 
  case 16:
    check =  qcd_f1f2f1_2pt_pr(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, sprop, uprop, geo, lt,0);
    check += qcd_f1f2f1_2pt_pr(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, sprop, uprop, geo, lt,1);  
    check += qcd_f1f2f1_2pt_pr(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, sprop, uprop, geo, lt,2);
    check += qcd_f1f2f1_2pt_pr(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, sprop, uprop, geo, lt,3);
    check += qcd_f1f2f1_2pt_pr(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, sprop, uprop, geo, lt,4);  
    check += qcd_f1f2f1_2pt_pr(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, sprop, uprop, geo, lt,5);
    check += qcd_f1f2f1_2pt_pr(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, sprop, uprop, geo, lt,6);
    check += qcd_f1f2f1_2pt_pr(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, sprop, uprop, geo, lt,7);  
    check += qcd_f1f2f1_2pt_pr(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, sprop, uprop, geo, lt,8); 
    break;
	
  case 17:
    check =  qcd_f1f2f1_2pt_pr(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, sprop, dprop, geo, lt,0);
    check += qcd_f1f2f1_2pt_pr(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, sprop, dprop, geo, lt,1);  
    check += qcd_f1f2f1_2pt_pr(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, sprop, dprop, geo, lt,2);
    check += qcd_f1f2f1_2pt_pr(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, sprop, dprop, geo, lt,3);
    check += qcd_f1f2f1_2pt_pr(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, sprop, dprop, geo, lt,4);  
    check += qcd_f1f2f1_2pt_pr(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, sprop, dprop, geo, lt,5);
    check += qcd_f1f2f1_2pt_pr(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, sprop, dprop, geo, lt,6);
    check += qcd_f1f2f1_2pt_pr(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, sprop, dprop, geo, lt,7);  
    check += qcd_f1f2f1_2pt_pr(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, sprop, dprop, geo, lt,8); 
    break;
 	  
  case 18:
    check =  qcd_f1f1f1_2pt_pr(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, sprop, geo, lt,0);
    check += qcd_f1f1f1_2pt_pr(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, sprop, geo, lt,1);  
    check += qcd_f1f1f1_2pt_pr(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, sprop, geo, lt,2);
    check += qcd_f1f1f1_2pt_pr(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, sprop, geo, lt,3);
    check += qcd_f1f1f1_2pt_pr(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, sprop, geo, lt,4);  
    check += qcd_f1f1f1_2pt_pr(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, sprop, geo, lt,5);
    check += qcd_f1f1f1_2pt_pr(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, sprop, geo, lt,6);
    check += qcd_f1f1f1_2pt_pr(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, sprop, geo, lt,7);  
    check += qcd_f1f1f1_2pt_pr(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, sprop, geo, lt,8); 
    break;
    //---------------------------------------------------------------------	  
  case 28:
    check =  qcd_sigmas4_2pt_pr(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, uprop, cprop, geo, lt,0);
    check += qcd_sigmas4_2pt_pr(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, uprop, cprop, geo, lt,1);  
    check += qcd_sigmas4_2pt_pr(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, uprop, cprop, geo, lt,2);
    check += qcd_sigmas4_2pt_pr(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, uprop, cprop, geo, lt,3);
    check += qcd_sigmas4_2pt_pr(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, uprop, cprop, geo, lt,4);  
    check += qcd_sigmas4_2pt_pr(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, uprop, cprop, geo, lt,5);
    check += qcd_sigmas4_2pt_pr(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, uprop, cprop, geo, lt,6);
    check += qcd_sigmas4_2pt_pr(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, uprop, cprop, geo, lt,7);  
    check += qcd_sigmas4_2pt_pr(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, uprop, cprop, geo, lt,8); 
    break;
	  
  case 29:
    check =  qcd_sigmas2_xistar_2pt_pr(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, uprop, dprop, cprop, geo, lt,0,1);
    check += qcd_sigmas2_xistar_2pt_pr(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, uprop, dprop, cprop, geo, lt,1,1);  
    check += qcd_sigmas2_xistar_2pt_pr(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, uprop, dprop, cprop, geo, lt,2,1);
    check += qcd_sigmas2_xistar_2pt_pr(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, uprop, dprop, cprop, geo, lt,3,1);
    check += qcd_sigmas2_xistar_2pt_pr(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, uprop, dprop, cprop, geo, lt,4,1);  
    check += qcd_sigmas2_xistar_2pt_pr(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, uprop, dprop, cprop, geo, lt,5,1);
    check += qcd_sigmas2_xistar_2pt_pr(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, uprop, dprop, cprop, geo, lt,6,1);
    check += qcd_sigmas2_xistar_2pt_pr(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, uprop, dprop, cprop, geo, lt,7,1);  
    check += qcd_sigmas2_xistar_2pt_pr(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, uprop, dprop, cprop, geo, lt,8,1);
    break;
	  
  case 30:
    check =  qcd_sigmas4_2pt_pr(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, dprop, cprop, geo, lt,0);
    check += qcd_sigmas4_2pt_pr(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, dprop, cprop, geo, lt,1);  
    check += qcd_sigmas4_2pt_pr(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, dprop, cprop, geo, lt,2);
    check += qcd_sigmas4_2pt_pr(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, dprop, cprop, geo, lt,3);
    check += qcd_sigmas4_2pt_pr(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, dprop, cprop, geo, lt,4);  
    check += qcd_sigmas4_2pt_pr(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, dprop, cprop, geo, lt,5);
    check += qcd_sigmas4_2pt_pr(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, dprop, cprop, geo, lt,6);
    check += qcd_sigmas4_2pt_pr(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, dprop, cprop, geo, lt,7);  
    check += qcd_sigmas4_2pt_pr(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, dprop, cprop, geo, lt,8); 
    break;
	  
  case 31:
    check =  qcd_f1f2f3_2pt_pr(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, sprop, uprop, cprop, geo, lt,0);
    check += qcd_f1f2f3_2pt_pr(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, sprop, uprop, cprop, geo, lt,1);  
    check += qcd_f1f2f3_2pt_pr(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, sprop, uprop, cprop, geo, lt,2);
    check += qcd_f1f2f3_2pt_pr(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, sprop, uprop, cprop, geo, lt,3);
    check += qcd_f1f2f3_2pt_pr(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, sprop, uprop, cprop, geo, lt,4);  
    check += qcd_f1f2f3_2pt_pr(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, sprop, uprop, cprop, geo, lt,5);
    check += qcd_f1f2f3_2pt_pr(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, sprop, uprop, cprop, geo, lt,6);
    check += qcd_f1f2f3_2pt_pr(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, sprop, uprop, cprop, geo, lt,7);  
    check += qcd_f1f2f3_2pt_pr(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, sprop, uprop, cprop, geo, lt,8); 
    break;
	
  case 32:
    check =  qcd_f1f2f3_2pt_pr(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, sprop, dprop, cprop, geo, lt,0);
    check += qcd_f1f2f3_2pt_pr(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, sprop, dprop, cprop, geo, lt,1);  
    check += qcd_f1f2f3_2pt_pr(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, sprop, dprop, cprop, geo, lt,2);
    check += qcd_f1f2f3_2pt_pr(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, sprop, dprop, cprop, geo, lt,3);
    check += qcd_f1f2f3_2pt_pr(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, sprop, dprop, cprop, geo, lt,4);  
    check += qcd_f1f2f3_2pt_pr(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, sprop, dprop, cprop, geo, lt,5);
    check += qcd_f1f2f3_2pt_pr(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, sprop, dprop, cprop, geo, lt,6);
    check += qcd_f1f2f3_2pt_pr(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, sprop, dprop, cprop, geo, lt,7);  
    check += qcd_f1f2f3_2pt_pr(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, sprop, dprop, cprop, geo, lt,8); 
    break;
	  
  case 33:
    check =  qcd_f1f2f1_2pt_pr(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, sprop, cprop, geo, lt,0);
    check += qcd_f1f2f1_2pt_pr(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, sprop, cprop, geo, lt,1);  
    check += qcd_f1f2f1_2pt_pr(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, sprop, cprop, geo, lt,2);
    check += qcd_f1f2f1_2pt_pr(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, sprop, cprop, geo, lt,3);
    check += qcd_f1f2f1_2pt_pr(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, sprop, cprop, geo, lt,4);  
    check += qcd_f1f2f1_2pt_pr(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, sprop, cprop, geo, lt,5);
    check += qcd_f1f2f1_2pt_pr(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, sprop, cprop, geo, lt,6);
    check += qcd_f1f2f1_2pt_pr(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, sprop, cprop, geo, lt,7);  
    check += qcd_f1f2f1_2pt_pr(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, sprop, cprop, geo, lt,8); 
    break;
    //-----------------------------------------------------------------------------------       
  case 37:
    check =  qcd_f1f2f1_2pt_pr(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, cprop, uprop, geo, lt,0);
    check += qcd_f1f2f1_2pt_pr(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, cprop, uprop, geo, lt,1);  
    check += qcd_f1f2f1_2pt_pr(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, cprop, uprop, geo, lt,2);
    check += qcd_f1f2f1_2pt_pr(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, cprop, uprop, geo, lt,3);
    check += qcd_f1f2f1_2pt_pr(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, cprop, uprop, geo, lt,4);  
    check += qcd_f1f2f1_2pt_pr(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, cprop, uprop, geo, lt,5);
    check += qcd_f1f2f1_2pt_pr(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, cprop, uprop, geo, lt,6);
    check += qcd_f1f2f1_2pt_pr(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, cprop, uprop, geo, lt,7);  
    check += qcd_f1f2f1_2pt_pr(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, cprop, uprop, geo, lt,8); 
    break;
    
  case 38:
    check =  qcd_f1f2f1_2pt_pr(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, cprop, dprop, geo, lt,0);
    check += qcd_f1f2f1_2pt_pr(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, cprop, dprop, geo, lt,1);  
    check += qcd_f1f2f1_2pt_pr(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, cprop, dprop, geo, lt,2);
    check += qcd_f1f2f1_2pt_pr(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, cprop, dprop, geo, lt,3);
    check += qcd_f1f2f1_2pt_pr(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, cprop, dprop, geo, lt,4);  
    check += qcd_f1f2f1_2pt_pr(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, cprop, dprop, geo, lt,5);
    check += qcd_f1f2f1_2pt_pr(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, cprop, dprop, geo, lt,6);
    check += qcd_f1f2f1_2pt_pr(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, cprop, dprop, geo, lt,7);  
    check += qcd_f1f2f1_2pt_pr(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, cprop, dprop, geo, lt,8); 
    break;	  
    
  case 39:
    check =  qcd_f1f2f1_2pt_pr(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, cprop, sprop, geo, lt,0);
    check += qcd_f1f2f1_2pt_pr(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, cprop, sprop, geo, lt,1);  
    check += qcd_f1f2f1_2pt_pr(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, cprop, sprop, geo, lt,2);
    check += qcd_f1f2f1_2pt_pr(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, cprop, sprop, geo, lt,3);
    check += qcd_f1f2f1_2pt_pr(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, cprop, sprop, geo, lt,4);  
    check += qcd_f1f2f1_2pt_pr(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, cprop, sprop, geo, lt,5);
    check += qcd_f1f2f1_2pt_pr(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, cprop, sprop, geo, lt,6);
    check += qcd_f1f2f1_2pt_pr(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, cprop, sprop, geo, lt,7);  
    check += qcd_f1f2f1_2pt_pr(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, cprop, sprop, geo, lt,8); 
    break;
    //-----------------------------------------------------------------------------------	  
  case 40:
    check =  qcd_f1f1f1_2pt_pr(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, cprop, geo, lt,0);
    check += qcd_f1f1f1_2pt_pr(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, cprop, geo, lt,1);  
    check += qcd_f1f1f1_2pt_pr(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, cprop, geo, lt,2);
    check += qcd_f1f1f1_2pt_pr(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, cprop, geo, lt,3);
    check += qcd_f1f1f1_2pt_pr(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, cprop, geo, lt,4);  
    check += qcd_f1f1f1_2pt_pr(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, cprop, geo, lt,5);
    check += qcd_f1f1f1_2pt_pr(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, cprop, geo, lt,6);
    check += qcd_f1f1f1_2pt_pr(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, cprop, geo, lt,7);  
    check += qcd_f1f1f1_2pt_pr(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, cprop, geo, lt,8); 
    break;
	
    //-----------------------------------------------------------------------------------	

  case 43:
    check =  qcd_deltas_xistar_omegastar_2pt_pr(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, sprop, uprop, geo, lt,0);
    check += qcd_deltas_xistar_omegastar_2pt_pr(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, sprop, uprop, geo, lt,1);  
    check += qcd_deltas_xistar_omegastar_2pt_pr(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, sprop, uprop, geo, lt,2);
    check += qcd_deltas_xistar_omegastar_2pt_pr(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, sprop, uprop, geo, lt,3);
    check += qcd_deltas_xistar_omegastar_2pt_pr(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, sprop, uprop, geo, lt,4);  
    check += qcd_deltas_xistar_omegastar_2pt_pr(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, sprop, uprop, geo, lt,5);
    check += qcd_deltas_xistar_omegastar_2pt_pr(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, sprop, uprop, geo, lt,6);
    check += qcd_deltas_xistar_omegastar_2pt_pr(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, sprop, uprop, geo, lt,7);  
    check += qcd_deltas_xistar_omegastar_2pt_pr(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, sprop, uprop, geo, lt,8); 
    break;
	   
  case 44:
    check =  qcd_deltas_xistar_omegastar_2pt_pr(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, sprop, dprop, geo, lt,0);
    check += qcd_deltas_xistar_omegastar_2pt_pr(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, sprop, dprop, geo, lt,1);  
    check += qcd_deltas_xistar_omegastar_2pt_pr(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, sprop, dprop, geo, lt,2);
    check += qcd_deltas_xistar_omegastar_2pt_pr(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, sprop, dprop, geo, lt,3);
    check += qcd_deltas_xistar_omegastar_2pt_pr(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, sprop, dprop, geo, lt,4);  
    check += qcd_deltas_xistar_omegastar_2pt_pr(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, sprop, dprop, geo, lt,5);
    check += qcd_deltas_xistar_omegastar_2pt_pr(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, sprop, dprop, geo, lt,6);
    check += qcd_deltas_xistar_omegastar_2pt_pr(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, sprop, dprop, geo, lt,7);  
    check += qcd_deltas_xistar_omegastar_2pt_pr(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, sprop, dprop, geo, lt,8);
    break;

  case 47:
    check =  qcd_sigmas2_xistar_2pt_pr(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, uprop, sprop, cprop, geo, lt,0,0);
    check += qcd_sigmas2_xistar_2pt_pr(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, uprop, sprop, cprop, geo, lt,1,0);  
    check += qcd_sigmas2_xistar_2pt_pr(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, uprop, sprop, cprop, geo, lt,2,0);
    check += qcd_sigmas2_xistar_2pt_pr(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, uprop, sprop, cprop, geo, lt,3,0);
    check += qcd_sigmas2_xistar_2pt_pr(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, uprop, sprop, cprop, geo, lt,4,0);  
    check += qcd_sigmas2_xistar_2pt_pr(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, uprop, sprop, cprop, geo, lt,5,0);
    check += qcd_sigmas2_xistar_2pt_pr(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, uprop, sprop, cprop, geo, lt,6,0);
    check += qcd_sigmas2_xistar_2pt_pr(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, uprop, sprop, cprop, geo, lt,7,0);  
    check += qcd_sigmas2_xistar_2pt_pr(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, uprop, sprop, cprop, geo, lt,8,0); 
    break;

  case 48:
    check =  qcd_sigmas2_xistar_2pt_pr(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, dprop, sprop, cprop, geo, lt,0,0);
    check += qcd_sigmas2_xistar_2pt_pr(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, dprop, sprop, cprop, geo, lt,1,0);  
    check += qcd_sigmas2_xistar_2pt_pr(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, dprop, sprop, cprop, geo, lt,2,0);
    check += qcd_sigmas2_xistar_2pt_pr(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, dprop, sprop, cprop, geo, lt,3,0);
    check += qcd_sigmas2_xistar_2pt_pr(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, dprop, sprop, cprop, geo, lt,4,0);  
    check += qcd_sigmas2_xistar_2pt_pr(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, dprop, sprop, cprop, geo, lt,5,0);
    check += qcd_sigmas2_xistar_2pt_pr(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, dprop, sprop, cprop, geo, lt,6,0);
    check += qcd_sigmas2_xistar_2pt_pr(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, dprop, sprop, cprop, geo, lt,7,0);  
    check += qcd_sigmas2_xistar_2pt_pr(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, dprop, sprop, cprop, geo, lt,8,0); 
    break;

  case 49:
    check =  qcd_deltas_xistar_omegastar_2pt_pr(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, sprop, cprop, geo, lt,0);
    check += qcd_deltas_xistar_omegastar_2pt_pr(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, sprop, cprop, geo, lt,1);  
    check += qcd_deltas_xistar_omegastar_2pt_pr(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, sprop, cprop, geo, lt,2);
    check += qcd_deltas_xistar_omegastar_2pt_pr(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, sprop, cprop, geo, lt,3);
    check += qcd_deltas_xistar_omegastar_2pt_pr(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, sprop, cprop, geo, lt,4);  
    check += qcd_deltas_xistar_omegastar_2pt_pr(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, sprop, cprop, geo, lt,5);
    check += qcd_deltas_xistar_omegastar_2pt_pr(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, sprop, cprop, geo, lt,6);
    check += qcd_deltas_xistar_omegastar_2pt_pr(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, sprop, cprop, geo, lt,7);  
    check += qcd_deltas_xistar_omegastar_2pt_pr(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, sprop, cprop, geo, lt,8);
    break;

  case 50:
    check =  qcd_deltas_xistar_omegastar_2pt_pr(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, cprop, uprop, geo, lt,0);
    check += qcd_deltas_xistar_omegastar_2pt_pr(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, cprop, uprop, geo, lt,1);  
    check += qcd_deltas_xistar_omegastar_2pt_pr(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, cprop, uprop, geo, lt,2);
    check += qcd_deltas_xistar_omegastar_2pt_pr(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, cprop, uprop, geo, lt,3);
    check += qcd_deltas_xistar_omegastar_2pt_pr(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, cprop, uprop, geo, lt,4);  
    check += qcd_deltas_xistar_omegastar_2pt_pr(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, cprop, uprop, geo, lt,5);
    check += qcd_deltas_xistar_omegastar_2pt_pr(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, cprop, uprop, geo, lt,6);
    check += qcd_deltas_xistar_omegastar_2pt_pr(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, cprop, uprop, geo, lt,7);  
    check += qcd_deltas_xistar_omegastar_2pt_pr(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, cprop, uprop, geo, lt,8);
    break;

  case 51:
    check =  qcd_deltas_xistar_omegastar_2pt_pr(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, cprop, dprop, geo, lt,0);
    check += qcd_deltas_xistar_omegastar_2pt_pr(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, cprop, dprop, geo, lt,1);  
    check += qcd_deltas_xistar_omegastar_2pt_pr(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, cprop, dprop, geo, lt,2);
    check += qcd_deltas_xistar_omegastar_2pt_pr(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, cprop, dprop, geo, lt,3);
    check += qcd_deltas_xistar_omegastar_2pt_pr(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, cprop, dprop, geo, lt,4);  
    check += qcd_deltas_xistar_omegastar_2pt_pr(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, cprop, dprop, geo, lt,5);
    check += qcd_deltas_xistar_omegastar_2pt_pr(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, cprop, dprop, geo, lt,6);
    check += qcd_deltas_xistar_omegastar_2pt_pr(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, cprop, dprop, geo, lt,7);  
    check += qcd_deltas_xistar_omegastar_2pt_pr(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, cprop, dprop, geo, lt,8);
    break;	
	
  case 52:
    check =  qcd_deltas_xistar_omegastar_2pt_pr(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, cprop, sprop, geo, lt,0);
    check += qcd_deltas_xistar_omegastar_2pt_pr(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, cprop, sprop, geo, lt,1);  
    check += qcd_deltas_xistar_omegastar_2pt_pr(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, cprop, sprop, geo, lt,2);
    check += qcd_deltas_xistar_omegastar_2pt_pr(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, cprop, sprop, geo, lt,3);
    check += qcd_deltas_xistar_omegastar_2pt_pr(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, cprop, sprop, geo, lt,4);  
    check += qcd_deltas_xistar_omegastar_2pt_pr(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, cprop, sprop, geo, lt,5);
    check += qcd_deltas_xistar_omegastar_2pt_pr(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, cprop, sprop, geo, lt,6);
    check += qcd_deltas_xistar_omegastar_2pt_pr(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, cprop, sprop, geo, lt,7);  
    check += qcd_deltas_xistar_omegastar_2pt_pr(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, cprop, sprop, geo, lt,8);
    break;	

    //-----------------------------------------------------------------------------------	

  case 53:
    check =  qcd_xistar_extra11_2pt_pr(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, sprop, uprop, geo, lt,0);
    check += qcd_xistar_extra11_2pt_pr(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, sprop, uprop, geo, lt,1);  
    check += qcd_xistar_extra11_2pt_pr(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, sprop, uprop, geo, lt,2);
    check += qcd_xistar_extra11_2pt_pr(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, sprop, uprop, geo, lt,3);
    check += qcd_xistar_extra11_2pt_pr(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, sprop, uprop, geo, lt,4);  
    check += qcd_xistar_extra11_2pt_pr(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, sprop, uprop, geo, lt,5);
    check += qcd_xistar_extra11_2pt_pr(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, sprop, uprop, geo, lt,6);
    check += qcd_xistar_extra11_2pt_pr(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, sprop, uprop, geo, lt,7);  
    check += qcd_xistar_extra11_2pt_pr(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, sprop, uprop, geo, lt,8); 
    break;
    
  case 54:
    check =  qcd_xistar_extra12_2pt_pr(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, sprop, uprop, geo, lt,0);
    check += qcd_xistar_extra12_2pt_pr(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, sprop, uprop, geo, lt,1);  
    check += qcd_xistar_extra12_2pt_pr(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, sprop, uprop, geo, lt,2);
    check += qcd_xistar_extra12_2pt_pr(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, sprop, uprop, geo, lt,3);
    check += qcd_xistar_extra12_2pt_pr(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, sprop, uprop, geo, lt,4);  
    check += qcd_xistar_extra12_2pt_pr(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, sprop, uprop, geo, lt,5);
    check += qcd_xistar_extra12_2pt_pr(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, sprop, uprop, geo, lt,6);
    check += qcd_xistar_extra12_2pt_pr(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, sprop, uprop, geo, lt,7);  
    check += qcd_xistar_extra12_2pt_pr(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, sprop, uprop, geo, lt,8); 
    break;
    
  case 55:
    check =  qcd_xistar_extra21_2pt_pr(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, sprop, uprop, geo, lt,0);
    check += qcd_xistar_extra21_2pt_pr(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, sprop, uprop, geo, lt,1);  
    check += qcd_xistar_extra21_2pt_pr(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, sprop, uprop, geo, lt,2);
    check += qcd_xistar_extra21_2pt_pr(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, sprop, uprop, geo, lt,3);
    check += qcd_xistar_extra21_2pt_pr(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, sprop, uprop, geo, lt,4);  
    check += qcd_xistar_extra21_2pt_pr(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, sprop, uprop, geo, lt,5);
    check += qcd_xistar_extra21_2pt_pr(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, sprop, uprop, geo, lt,6);
    check += qcd_xistar_extra21_2pt_pr(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, sprop, uprop, geo, lt,7);  
    check += qcd_xistar_extra21_2pt_pr(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, sprop, uprop, geo, lt,8); 
    break;
    
  case 56:
    check =  qcd_xistar_extra22_2pt_pr(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, sprop, uprop, geo, lt,0);
    check += qcd_xistar_extra22_2pt_pr(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, sprop, uprop, geo, lt,1);  
    check += qcd_xistar_extra22_2pt_pr(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, sprop, uprop, geo, lt,2);
    check += qcd_xistar_extra22_2pt_pr(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, sprop, uprop, geo, lt,3);
    check += qcd_xistar_extra22_2pt_pr(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, sprop, uprop, geo, lt,4);  
    check += qcd_xistar_extra22_2pt_pr(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, sprop, uprop, geo, lt,5);
    check += qcd_xistar_extra22_2pt_pr(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, sprop, uprop, geo, lt,6);
    check += qcd_xistar_extra22_2pt_pr(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, sprop, uprop, geo, lt,7);  
    check += qcd_xistar_extra22_2pt_pr(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, sprop, uprop, geo, lt,8); 
    break;
    
    //-----------------------------------------------------------------------------------	

  case 57:
    check =  qcd_xistar_extra11_2pt_pr(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, sprop, dprop, geo, lt,0);
    check += qcd_xistar_extra11_2pt_pr(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, sprop, dprop, geo, lt,1);  
    check += qcd_xistar_extra11_2pt_pr(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, sprop, dprop, geo, lt,2);
    check += qcd_xistar_extra11_2pt_pr(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, sprop, dprop, geo, lt,3);
    check += qcd_xistar_extra11_2pt_pr(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, sprop, dprop, geo, lt,4);  
    check += qcd_xistar_extra11_2pt_pr(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, sprop, dprop, geo, lt,5);
    check += qcd_xistar_extra11_2pt_pr(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, sprop, dprop, geo, lt,6);
    check += qcd_xistar_extra11_2pt_pr(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, sprop, dprop, geo, lt,7);  
    check += qcd_xistar_extra11_2pt_pr(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, sprop, dprop, geo, lt,8); 
    break;
    
  case 58:
    check =  qcd_xistar_extra12_2pt_pr(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, sprop, dprop, geo, lt,0);
    check += qcd_xistar_extra12_2pt_pr(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, sprop, dprop, geo, lt,1);  
    check += qcd_xistar_extra12_2pt_pr(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, sprop, dprop, geo, lt,2);
    check += qcd_xistar_extra12_2pt_pr(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, sprop, dprop, geo, lt,3);
    check += qcd_xistar_extra12_2pt_pr(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, sprop, dprop, geo, lt,4);  
    check += qcd_xistar_extra12_2pt_pr(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, sprop, dprop, geo, lt,5);
    check += qcd_xistar_extra12_2pt_pr(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, sprop, dprop, geo, lt,6);
    check += qcd_xistar_extra12_2pt_pr(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, sprop, dprop, geo, lt,7);  
    check += qcd_xistar_extra12_2pt_pr(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, sprop, dprop, geo, lt,8); 
    break;
    
  case 59:
    check =  qcd_xistar_extra21_2pt_pr(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, sprop, dprop, geo, lt,0);
    check += qcd_xistar_extra21_2pt_pr(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, sprop, dprop, geo, lt,1);  
    check += qcd_xistar_extra21_2pt_pr(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, sprop, dprop, geo, lt,2);
    check += qcd_xistar_extra21_2pt_pr(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, sprop, dprop, geo, lt,3);
    check += qcd_xistar_extra21_2pt_pr(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, sprop, dprop, geo, lt,4);  
    check += qcd_xistar_extra21_2pt_pr(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, sprop, dprop, geo, lt,5);
    check += qcd_xistar_extra21_2pt_pr(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, sprop, dprop, geo, lt,6);
    check += qcd_xistar_extra21_2pt_pr(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, sprop, dprop, geo, lt,7);  
    check += qcd_xistar_extra21_2pt_pr(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, sprop, dprop, geo, lt,8); 
    break;
    
  case 60:
    check =  qcd_xistar_extra22_2pt_pr(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, sprop, dprop, geo, lt,0);
    check += qcd_xistar_extra22_2pt_pr(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, sprop, dprop, geo, lt,1);  
    check += qcd_xistar_extra22_2pt_pr(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, sprop, dprop, geo, lt,2);
    check += qcd_xistar_extra22_2pt_pr(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, sprop, dprop, geo, lt,3);
    check += qcd_xistar_extra22_2pt_pr(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, sprop, dprop, geo, lt,4);  
    check += qcd_xistar_extra22_2pt_pr(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, sprop, dprop, geo, lt,5);
    check += qcd_xistar_extra22_2pt_pr(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, sprop, dprop, geo, lt,6);
    check += qcd_xistar_extra22_2pt_pr(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, sprop, dprop, geo, lt,7);  
    check += qcd_xistar_extra22_2pt_pr(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, sprop, dprop, geo, lt,8); 
    break;
		   
  }//-switch

  check += qcd_projector12_2pt(gamma12,gamma13,gamma23,block12,block_pr,geo);  
  check += qcd_projector32_2pt(gamma12,gamma13,gamma23,block32,block_pr,geo);
  check += qcd_noprojector_2pt(blocknp,block_pr,geo);

  if ( check==12 ){	
    if(geo->myid==0) printf("Particle: %s projector 1/2, 3/2 and diagonal t=%d\n",particle_names[particle_id],lt);
  }
  
  for(i=0;i<9;i++)
    for(j=0;j<4;j++)
      for(k=0;k<4;k++){
	free(block_pr[i][j][k]);
      }
						  
}//--main function closes
