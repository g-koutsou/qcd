/*
 * Christos Kallidonis
 * June 2012
 * 
 * This program contains contractions for the 40 particles
 * to obtain the 3-point functions using the fixed current method.
 * 
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <qcd.h>
#include <projectors_pb.h>

int qcd_projector_3pt_spin12(qcd_complex_16 *block[5], qcd_complex_16 *block_pr[9][4][4], qcd_geometry *geo){
	
  qcd_uint_2 ga,gap,prid;
  qcd_uint_4 lx,ly,lz,v,v3;
  
  qcd_uint_4 projlist[5] = {3,4,13,15,16},proj;
  
  for(prid=0;prid<5;prid++){
    proj = projlist[prid];

    for(lx=0; lx<geo->lL[1]; lx++)
      for(ly=0; ly<geo->lL[2]; ly++)
	for(lz=0; lz<geo->lL[3]; lz++){
	  v3 = qcd_LEXIC0(lx,ly,lz,geo->lL);
	  
	  for(ga =0; ga <4; ga ++) 
	    for(gap=0; gap<4; gap++)	
		block[prid][v3] = qcd_CADD( block[prid][v3],qcd_CMUL(PROJECTOR[proj][ga][gap],block_pr[8][gap][ga][v3]) );

	}//-space
  }//-prid

  return (1);
} //-routine
//======================================================================

int qcd_projector_pr32_3pt_spin32(qcd_complex_16 *block[5], qcd_complex_16 *block_pr[9][4][4], 
				  qcd_complex_16 gamma12[4][4], qcd_complex_16 gamma13[4][4], qcd_complex_16 gamma23[4][4], qcd_geometry *geo){
	
  qcd_uint_2 al,be,ga,de,j,i,k,i1,i2,prid,proj;
  qcd_uint_4 lx,ly,lz,v,v3;
  qcd_complex_16 Pr32[9][4][4],delta;
  qcd_complex_16 Fin_proj[5][4][4][3][3];
  qcd_real_8 fac = 1.0/3.0;
  qcd_complex_16 cfac;

  qcd_uint_4 projlist[5] = {3,4,13,15,16};

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

  /* //- Print the projector to 3/2 */
  /* for(i=0;i<3;i++){ */
  /*   for(j=0;j<3;j++){ */
  /*     for(al=0;al<4;al++){ */
  /* 	for(be=0;be<4;be++){ */
  /* 	  i1 = j+3*i; */
  /* 	  printf("Pr %d %d %d %d = %lf %lf\n",al+1,be+1,i+1,j+1,Pr32[i1][al][be].re,Pr32[i1][al][be].im); */
  /* 	} */
  /*     } */
  /*   } */
  /* } */


  //- Print the unprojected 3pt-function
  /* v3=0; */
  /* for(i=0;i<3;i++){ */
  /*   for(j=0;j<3;j++){ */
  /*     for(al=0;al<4;al++){ */
  /* 	for(be=0;be<4;be++){ */
  /* 	  i1 = j+3*i; */
  /* 	  printf("Rank %d: 3pt_unproj[v=%d] %d %d %d %d = %+e %+e\n",geo->myid,v3,al+1,be+1,i+1,j+1,block_pr[i1][al][be][v3].re,block_pr[i1][al][be][v3].im); */
  /* 	} */
  /*     } */
  /*   } */
  /* } */

  //-Define the projector 
  for(prid=0;prid<5;prid++){
    proj = projlist[prid];

    for(al=0;al<4;al++){
      for(de=0;de<4;de++){
  	for(k=0;k<3;k++){
  	  for(j=0;j<3;j++){
  	    Fin_proj[prid][al][de][k][j] = (qcd_complex_16) {0.0,0.0};

  	    for(i=0;i<3;i++){
  	      i1 = i+3*k;
  	      i2 = j+3*i;
  	      for(be=0;be<4;be++){
  		for(ga=0;ga<4;ga++){

  		  Fin_proj[prid][al][de][k][j] = qcd_CADD(Fin_proj[prid][al][de][k][j],qcd_CMUL(Pr32[i1][al][be],qcd_CMUL(PROJECTOR[proj][be][ga],Pr32[i2][ga][de])));

  		}//-ga
  	      }//-be
  	    }//-i
  	  }//-j
  	}//-k
      }//-de
    }//-al
  }//-prid

  //-Project
  for(prid=0;prid<5;prid++){
    proj = projlist[prid];
    for(lx=0; lx<geo->lL[1]; lx++)
      for(ly=0; ly<geo->lL[2]; ly++)
	for(lz=0; lz<geo->lL[3]; lz++){
	    v3 = qcd_LEXIC0(lx,ly,lz,geo->lL);

	    for(k=0;k<3;k++){
	      for(j=0;j<3;j++){
		i2 = k+j*3;
		for(al=0;al<4;al++){
		  for(de=0;de<4;de++){
		    block[prid][v3] = qcd_CADD(block[prid][v3],qcd_CMUL(Fin_proj[prid][al][de][k][j],block_pr[i2][de][al][v3]));
		  }//-de
		}//-al
	      }//-k
	    }//-j
	}//-volume	
  }//-prid	    

  //- Print the projected 3pt-function
  /* v3=0; */
  /* printf("Rank %d: 3pt_g5g1[v=%d] = %+e %+e\n",geo->myid,v3,block[4][v3].re,block[4][v3].im); */
  /* printf("Rank %d: 3pt_g5g2[v=%d] = %+e %+e\n",geo->myid,v3,block[3][v3].re,block[3][v3].im); */
  /* printf("Rank %d: 3pt_g5g3[v=%d] = %+e %+e\n",geo->myid,v3,block[1][v3].re,block[1][v3].im); */

  return (1);
} //-routine
//======================================================================

int qcd_f1f1f1_3pt(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[9][4][4],
		   qcd_propagator *prop, qcd_propagator *seqprop, qcd_geometry *geo, qcd_uint_4 lt,qcd_uint_4 elem){
							
  qcd_int_4 ctr2;
  qcd_uint_2 a,b,c,ap,bp,cp,al,be,ga,alp,bep,gap,cc1,cc2;
  qcd_uint_4 lx,ly,lz,v,v3;
  int lv3 = geo->lL[1]*geo->lL[2]*geo->lL[3];
  	
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
                
#pragma omp parallel for private(lz,ly,lx,v)
	    for(v3=0; v3<lv3; v3++) {
	      /* for(lx=0; lx<geo->lL[1]; lx++) */
	      /*   for(ly=0; ly<geo->lL[2]; ly++) */
	      /*     for(lz=0; lz<geo->lL[3]; lz++){ */
	      //v3 = qcd_LEXIC0(lx,ly,lz,geo->lL);
	      //qcd_LEXIC(lt,lx,ly,lz,geo->lL);
	      lz = v3 % geo->lL[3];
	      ly = ((v3 - lz)/geo->lL[3]) % geo->lL[2];
	      lx = ((v3 - lz)/geo->lL[3] - ly) / geo->lL[2];
	      v =  qcd_LEXIC(lt,lx,ly,lz,geo->lL);

	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(prop->D[v][be][bep][b][bp],
													qcd_CMUL(prop->D[v][ga][alp][c][ap],
														 seqprop->D[v][al][gap][a][cp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2])
						  ); //-1
	    
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(prop->D[v][be][gap][b][cp],
													qcd_CMUL(prop->D[v][ga][alp][c][ap],
														 seqprop->D[v][al][bep][a][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2])
						  );//-2

	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(prop->D[v][be][alp][b][ap],
													qcd_CMUL(prop->D[v][ga][bep][c][bp],
														 seqprop->D[v][al][gap][a][cp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2])
						  );//-3	    
										    
	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(prop->D[v][be][gap][b][cp],
													qcd_CMUL(prop->D[v][ga][bep][c][bp],
														 seqprop->D[v][al][alp][a][ap]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2])
						  );//-4
										    
	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(prop->D[v][be][alp][b][ap],
													qcd_CMUL(prop->D[v][ga][gap][c][cp],
														 seqprop->D[v][al][bep][a][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2])
						  );//-5
										    
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(prop->D[v][be][bep][b][bp],
													qcd_CMUL(prop->D[v][ga][gap][c][cp],
														 seqprop->D[v][al][alp][a][ap]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2])
						  );//-6
				
	      //-------------------------				
				
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(prop->D[v][al][bep][a][bp],
													qcd_CMUL(prop->D[v][ga][alp][c][ap],
														 seqprop->D[v][be][gap][b][cp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2])
						  ); //-7
	    
	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(prop->D[v][al][gap][a][cp],
													qcd_CMUL(prop->D[v][ga][alp][c][ap],
														 seqprop->D[v][be][bep][b][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2])
						  );//-8

	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(prop->D[v][al][alp][a][ap],
													qcd_CMUL(prop->D[v][ga][bep][c][bp],
														 seqprop->D[v][be][gap][b][cp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2])
						  );//-9	    
										    
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(prop->D[v][al][gap][a][cp],
													qcd_CMUL(prop->D[v][ga][bep][c][bp],
														 seqprop->D[v][be][alp][b][ap]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2])
						  );//-10
										    
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(prop->D[v][al][alp][a][ap],
													qcd_CMUL(prop->D[v][ga][gap][c][cp],
														 seqprop->D[v][be][bep][b][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2])
						  );//-11
										    
	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(prop->D[v][al][bep][a][bp],
													qcd_CMUL(prop->D[v][ga][gap][c][cp],
														 seqprop->D[v][be][alp][b][ap]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2])
						  );//-12				
				
	      //-------------------------				
				
	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(prop->D[v][al][bep][a][bp],
													qcd_CMUL(prop->D[v][be][alp][b][ap],
														 seqprop->D[v][ga][gap][c][cp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2])
						  ); //-13
	    
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(prop->D[v][al][gap][a][cp],
													qcd_CMUL(prop->D[v][be][alp][b][ap],
														 seqprop->D[v][ga][bep][c][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2])
						  );//-14

	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(prop->D[v][al][alp][a][ap],
													qcd_CMUL(prop->D[v][be][bep][b][bp],
														 seqprop->D[v][ga][gap][c][cp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2])
						  );//-15	    
										    
	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(prop->D[v][al][gap][a][cp],
													qcd_CMUL(prop->D[v][be][bep][b][bp],
														 seqprop->D[v][ga][alp][c][ap]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2])
						  );//-16
										    
	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(prop->D[v][al][alp][a][ap],
													qcd_CMUL(prop->D[v][be][gap][b][cp],
														 seqprop->D[v][ga][bep][c][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2])
						  );//-17
										    
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(prop->D[v][al][bep][a][bp],
													qcd_CMUL(prop->D[v][be][gap][b][cp],
														 seqprop->D[v][ga][alp][c][ap]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2])
						  );//-18				
					    
	    
	    }//space loop
	  }//color2 loop    
	}//color1 loop
      }//nonvanishing cgcg loop						
    }//nonvanishing projector condition
	
  return (1);
}
//======================================================================

int qcd_f1f2f1_f1_3pt(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[9][4][4],
		      qcd_propagator *propf1, qcd_propagator *propf2, qcd_propagator *seqpropf1, qcd_geometry *geo, qcd_uint_4 lt, qcd_uint_4 elem){
							
  qcd_int_4 ctr2;
  qcd_uint_2 a,b,c,ap,bp,cp,al,be,ga,alp,bep,gap,cc1,cc2;
  qcd_uint_4 lx,ly,lz,v,v3;
  int lv3 = geo->lL[1]*geo->lL[2]*geo->lL[3];
    
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
                
#pragma omp parallel for private(lz,ly,lx,v)
	    for(v3=0; v3<lv3; v3++) {
	      /* for(lx=0; lx<geo->lL[1]; lx++) */
	      /*   for(ly=0; ly<geo->lL[2]; ly++) */
	      /*     for(lz=0; lz<geo->lL[3]; lz++){ */
	      //v3 = qcd_LEXIC0(lx,ly,lz,geo->lL);
	      //qcd_LEXIC(lt,lx,ly,lz,geo->lL);
	      lz = v3 % geo->lL[3];
	      ly = ((v3 - lz)/geo->lL[3]) % geo->lL[2];
	      lx = ((v3 - lz)/geo->lL[3] - ly) / geo->lL[2];
	      v =  qcd_LEXIC(lt,lx,ly,lz,geo->lL);

	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][be][bep][b][bp],
													qcd_CMUL(propf1->D[v][ga][gap][c][cp],
														 seqpropf1->D[v][al][alp][a][ap]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2])
						  );
				
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][be][bep][b][bp],
													qcd_CMUL(propf1->D[v][al][alp][a][ap],
														 seqpropf1->D[v][ga][gap][c][cp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2])
						  );				
										    
	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][be][bep][b][bp],
													qcd_CMUL(propf1->D[v][ga][alp][c][ap],
														 seqpropf1->D[v][al][gap][a][cp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2])
						  );

	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][be][bep][b][bp],
													qcd_CMUL(propf1->D[v][al][gap][a][cp],
														 seqpropf1->D[v][ga][alp][c][ap]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2])
						  );				
				
	    }//space loop
	  }//color2 loop    
	}//color1 loop
      }//nonvanishing cgcg loop						
    }// ga gap indices
	
  return (1);
}

//======================================================================

int qcd_pnextra_f1_3pt(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[9][4][4],
		       qcd_propagator *propf1, qcd_propagator *propf2, qcd_propagator *seqpropf1, qcd_geometry *geo, qcd_uint_4 lt, qcd_uint_4 elem){
							
  qcd_int_4 ctr2;
  qcd_uint_2 a,b,c,ap,bp,cp,al,be,ga,alp,bep,gap,cc1,cc2,de,dep;
  qcd_uint_4 lx,ly,lz,v,v3;
  qcd_complex_16 term;
  int lv3 = geo->lL[1]*geo->lL[2]*geo->lL[3];
    
  for(de =0; de <4; de ++) 
    for(dep=0; dep<4; dep++){ 
    
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

		term = qcd_CSCALE( qcd_CMUL(qcd_GAMMA[5][de][ga],qcd_GAMMA[5][gap][dep]),-qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2] );
                
#pragma omp parallel for private(lz,ly,lx,v)
		for(v3=0; v3<lv3; v3++) {
		  /* for(lx=0; lx<geo->lL[1]; lx++) */
		  /*   for(ly=0; ly<geo->lL[2]; ly++) */
		  /*     for(lz=0; lz<geo->lL[3]; lz++){ */
		  //v3 = qcd_LEXIC0(lx,ly,lz,geo->lL);
		  //qcd_LEXIC(lt,lx,ly,lz,geo->lL);
		  lz = v3 % geo->lL[3];
		  ly = ((v3 - lz)/geo->lL[3]) % geo->lL[2];
		  lx = ((v3 - lz)/geo->lL[3] - ly) / geo->lL[2];
		  v =  qcd_LEXIC(lt,lx,ly,lz,geo->lL);

		  block[elem][de][dep][v3] = qcd_CADD(block[elem][de][dep][v3],qcd_CMUL(
											qcd_CMUL(cgcg_val[ctr2],
												 qcd_CMUL(propf2->D[v][be][bep][b][bp],
													  qcd_CMUL(propf1->D[v][ga][gap][c][cp],
														   seqpropf1->D[v][al][alp][a][ap]))),term)
						      );
				
		  block[elem][de][dep][v3] = qcd_CADD(block[elem][de][dep][v3],qcd_CMUL(
											qcd_CMUL(cgcg_val[ctr2],
												 qcd_CMUL(propf2->D[v][be][bep][b][bp],
													  qcd_CMUL(propf1->D[v][al][alp][a][ap],
														   seqpropf1->D[v][ga][gap][c][cp]))),term)
						      );				
										    
		  block[elem][de][dep][v3] = qcd_CSUB(block[elem][de][dep][v3],qcd_CMUL(
											qcd_CMUL(cgcg_val[ctr2],
												 qcd_CMUL(propf2->D[v][be][bep][b][bp],
													  qcd_CMUL(propf1->D[v][ga][alp][c][ap],
														   seqpropf1->D[v][al][gap][a][cp]))),term)
						      );

		  block[elem][de][dep][v3] = qcd_CSUB(block[elem][de][dep][v3],qcd_CMUL(
											qcd_CMUL(cgcg_val[ctr2],
												 qcd_CMUL(propf2->D[v][be][bep][b][bp],
													  qcd_CMUL(propf1->D[v][al][gap][a][cp],
														   seqpropf1->D[v][ga][alp][c][ap]))),term)
						      );				
				
		}//space loop
	      }//color2 loop    
	    }//color1 loop
	  }//nonvanishing cgcg loop						
	}// ga gap indices
    }// de dep indices
	
  return (1);
}
//======================================================================

int qcd_f1f2f1_f2_3pt(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[9][4][4],
		      qcd_propagator *propf1, qcd_propagator *propf2, qcd_propagator *seqpropf2, qcd_geometry *geo, qcd_uint_4 lt, qcd_uint_4 elem){
							
  qcd_int_4 ctr2;
  qcd_uint_2 a,b,c,ap,bp,cp,al,be,ga,alp,bep,gap,cc1,cc2;
  qcd_uint_4 lx,ly,lz,v,v3;
  int lv3 = geo->lL[1]*geo->lL[2]*geo->lL[3];
    
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
                
#pragma omp parallel for private(lz,ly,lx,v)
	    for(v3=0; v3<lv3; v3++) {
	      /* for(lx=0; lx<geo->lL[1]; lx++) */
	      /*   for(ly=0; ly<geo->lL[2]; ly++) */
	      /*     for(lz=0; lz<geo->lL[3]; lz++){ */
	      //v3 = qcd_LEXIC0(lx,ly,lz,geo->lL);
	      //qcd_LEXIC(lt,lx,ly,lz,geo->lL);
	      lz = v3 % geo->lL[3];
	      ly = ((v3 - lz)/geo->lL[3]) % geo->lL[2];
	      lx = ((v3 - lz)/geo->lL[3] - ly) / geo->lL[2];
	      v =  qcd_LEXIC(lt,lx,ly,lz,geo->lL);

	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][al][alp][a][ap],
													qcd_CMUL(propf1->D[v][ga][gap][c][cp],
														 seqpropf2->D[v][be][bep][b][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2])
						  );
				

	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][al][gap][a][cp],
													qcd_CMUL(propf1->D[v][ga][alp][c][ap],
														 seqpropf2->D[v][be][bep][b][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2])
						  );				
				
	    }//space loop
	  }//color2 loop    
	}//color1 loop
      }//nonvanishing cgcg loop						
    }// ga gap indices
	
  return (1);
}
//======================================================================

int qcd_pnextra_f2_3pt(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[9][4][4],
		       qcd_propagator *propf1, qcd_propagator *propf2, qcd_propagator *seqpropf2, qcd_geometry *geo, qcd_uint_4 lt, qcd_uint_4 elem){
							
  qcd_int_4 ctr2;
  qcd_uint_2 a,b,c,ap,bp,cp,al,be,ga,alp,bep,gap,cc1,cc2,de,dep;
  qcd_uint_4 lx,ly,lz,v,v3;
  
  qcd_complex_16 term;
  int lv3 = geo->lL[1]*geo->lL[2]*geo->lL[3];
    
  for(de =0; de <4; de ++) 
    for(dep=0; dep<4; dep++){  
  
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

		term = qcd_CSCALE( qcd_CMUL(qcd_GAMMA[5][de][ga],qcd_GAMMA[5][gap][dep]),-qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2] );
                
#pragma omp parallel for private(lz,ly,lx,v)
		for(v3=0; v3<lv3; v3++) {
		  /* for(lx=0; lx<geo->lL[1]; lx++) */
		  /*   for(ly=0; ly<geo->lL[2]; ly++) */
		  /*     for(lz=0; lz<geo->lL[3]; lz++){ */
		  //v3 = qcd_LEXIC0(lx,ly,lz,geo->lL);
		  //qcd_LEXIC(lt,lx,ly,lz,geo->lL);
		  lz = v3 % geo->lL[3];
		  ly = ((v3 - lz)/geo->lL[3]) % geo->lL[2];
		  lx = ((v3 - lz)/geo->lL[3] - ly) / geo->lL[2];
		  v =  qcd_LEXIC(lt,lx,ly,lz,geo->lL);

		  block[elem][de][dep][v3] = qcd_CADD(block[elem][de][dep][v3],qcd_CMUL(
											qcd_CMUL(cgcg_val[ctr2],
												 qcd_CMUL(propf1->D[v][al][alp][a][ap],
													  qcd_CMUL(propf1->D[v][ga][gap][c][cp],
														   seqpropf2->D[v][be][bep][b][bp]))),term)
						      );
				

		  block[elem][de][dep][v3] = qcd_CSUB(block[elem][de][dep][v3],qcd_CMUL(
											qcd_CMUL(cgcg_val[ctr2],
												 qcd_CMUL(propf1->D[v][al][gap][a][cp],
													  qcd_CMUL(propf1->D[v][ga][alp][c][ap],
														   seqpropf2->D[v][be][bep][b][bp]))),term)
						      );				
				
		}//space loop
	      }//color2 loop    
	    }//color1 loop
	  }//nonvanishing cgcg loop						
	}// ga gap indices
    }// de dep indices
	
  return (1);
}

//======================================================================

int qcd_deltas_f1_3pt(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[9][4][4],
		      qcd_propagator *propf1, qcd_propagator *propf2, qcd_propagator *seqpropf1, qcd_geometry *geo, qcd_uint_4 lt, qcd_uint_4 elem){
							
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
    lz = v3 % geo->lL[3];
    ly = ((v3 - lz)/geo->lL[3]) % geo->lL[2];
    lx = ((v3 - lz)/geo->lL[3] - ly) / geo->lL[2];
    v =  qcd_LEXIC(lt,lx,ly,lz,geo->lL);

	
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

	      //-- 1-1
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][be][bep][b][bp],
													qcd_CMUL(propf1->D[v][ga][gap][c][cp],
														 seqpropf1->D[v][al][alp][a][ap]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*4.0/3.0)
						  );
				
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][be][bep][b][bp],
													qcd_CMUL(propf1->D[v][al][alp][a][ap],
														 seqpropf1->D[v][ga][gap][c][cp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*4.0/3.0)
						  );				
										    
	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][be][bep][b][bp],
													qcd_CMUL(propf1->D[v][ga][alp][c][ap],
														 seqpropf1->D[v][al][gap][a][cp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*4.0/3.0)
						  );

	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][be][bep][b][bp],
													qcd_CMUL(propf1->D[v][al][gap][a][cp],
														 seqpropf1->D[v][ga][alp][c][ap]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*4.0/3.0)
						  );			
	      //-- 1-2								
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][be][gap][b][cp],
													qcd_CMUL(propf1->D[v][ga][alp][c][ap],
														 seqpropf1->D[v][al][bep][a][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*2.0/3.0)
						  );
								
	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][be][gap][b][cp],
													qcd_CMUL(propf1->D[v][ga][bep][c][bp],
														 seqpropf1->D[v][al][alp][a][ap]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*2.0/3.0)
						  );				
										    
	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][be][gap][b][cp],
													qcd_CMUL(propf1->D[v][al][alp][a][ap],
														 seqpropf1->D[v][ga][bep][c][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*2.0/3.0)
						  );

	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][be][gap][b][cp],
													qcd_CMUL(propf1->D[v][al][bep][a][bp],
														 seqpropf1->D[v][ga][alp][c][ap]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*2.0/3.0)
						  );			
	      //-- 2-1								
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][ga][bep][c][bp],
													qcd_CMUL(propf1->D[v][be][alp][b][ap],
														 seqpropf1->D[v][al][gap][a][cp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*2.0/3.0)
						  );
								
	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][ga][bep][c][bp],
													qcd_CMUL(propf1->D[v][be][gap][b][cp],
														 seqpropf1->D[v][al][alp][a][ap]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*2.0/3.0)
						  );				
										    
	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][ga][bep][c][bp],
													qcd_CMUL(propf1->D[v][al][alp][a][ap],
														 seqpropf1->D[v][be][gap][b][cp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*2.0/3.0)
						  );

	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][ga][bep][c][bp],
													qcd_CMUL(propf1->D[v][al][gap][a][cp],
														 seqpropf1->D[v][be][alp][b][ap]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*2.0/3.0)
						  );				
	      //-- 2-2
	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][ga][gap][c][cp],
													qcd_CMUL(propf1->D[v][be][alp][b][ap],
														 seqpropf1->D[v][al][bep][a][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]/3.0)
						  );
								
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][ga][gap][c][cp],
													qcd_CMUL(propf1->D[v][be][bep][b][bp],
														 seqpropf1->D[v][al][alp][a][ap]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]/3.0)
						  );				
										    
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][ga][gap][c][cp],
													qcd_CMUL(propf1->D[v][al][alp][a][ap],
														 seqpropf1->D[v][be][bep][b][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]/3.0)
						  );

	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][ga][gap][c][cp],
													qcd_CMUL(propf1->D[v][al][bep][a][bp],
														 seqpropf1->D[v][be][alp][b][ap]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]/3.0)
						  );				
				 
	    }//color2 loop     
	  }//color1 loop                                                                               
	}//spin loop (cg_cg)
      }//ga gap
  }//space loop                                                                                   

  return (1);
}
//======================================================================

int qcd_xistar_extra11_f1_3pt(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[9][4][4],
		      qcd_propagator *propf1, qcd_propagator *propf2, qcd_propagator *seqpropf1, qcd_geometry *geo, qcd_uint_4 lt, qcd_uint_4 elem){
							
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
    lz = v3 % geo->lL[3];
    ly = ((v3 - lz)/geo->lL[3]) % geo->lL[2];
    lx = ((v3 - lz)/geo->lL[3] - ly) / geo->lL[2];
    v =  qcd_LEXIC(lt,lx,ly,lz,geo->lL);

	
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

	      //-- 1-1
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][be][bep][b][bp],
													qcd_CMUL(propf1->D[v][ga][gap][c][cp],
														 seqpropf1->D[v][al][alp][a][ap]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2])
						  );
				
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][be][bep][b][bp],
													qcd_CMUL(propf1->D[v][al][alp][a][ap],
														 seqpropf1->D[v][ga][gap][c][cp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2])
						  );				
										    
	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][be][bep][b][bp],
													qcd_CMUL(propf1->D[v][ga][alp][c][ap],
														 seqpropf1->D[v][al][gap][a][cp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2])
						  );

	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][be][bep][b][bp],
													qcd_CMUL(propf1->D[v][al][gap][a][cp],
														 seqpropf1->D[v][ga][alp][c][ap]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2])
						  );			
				 
	    }//color2 loop     
	  }//color1 loop                                                                               
	}//spin loop (cg_cg)
      }//ga gap
  }//space loop                                                                                   

  return (1);
}
//======================================================================


int qcd_xistar_extra12_f1_3pt(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[9][4][4],
		      qcd_propagator *propf1, qcd_propagator *propf2, qcd_propagator *seqpropf1, qcd_geometry *geo, qcd_uint_4 lt, qcd_uint_4 elem){
							
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
    lz = v3 % geo->lL[3];
    ly = ((v3 - lz)/geo->lL[3]) % geo->lL[2];
    lx = ((v3 - lz)/geo->lL[3] - ly) / geo->lL[2];
    v =  qcd_LEXIC(lt,lx,ly,lz,geo->lL);

	
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

	      //-- 1-2								
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][be][gap][b][cp],
													qcd_CMUL(propf1->D[v][ga][alp][c][ap],
														 seqpropf1->D[v][al][bep][a][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2])
						  );
								
	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][be][gap][b][cp],
													qcd_CMUL(propf1->D[v][ga][bep][c][bp],
														 seqpropf1->D[v][al][alp][a][ap]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2])
						  );				
										    
	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][be][gap][b][cp],
													qcd_CMUL(propf1->D[v][al][alp][a][ap],
														 seqpropf1->D[v][ga][bep][c][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2])
						  );

	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][be][gap][b][cp],
													qcd_CMUL(propf1->D[v][al][bep][a][bp],
														 seqpropf1->D[v][ga][alp][c][ap]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2])
						  );			
				 
	    }//color2 loop     
	  }//color1 loop                                                                               
	}//spin loop (cg_cg)
      }//ga gap
  }//space loop                                                                                   

  return (1);
}
//======================================================================


int qcd_xistar_extra21_f1_3pt(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[9][4][4],
		      qcd_propagator *propf1, qcd_propagator *propf2, qcd_propagator *seqpropf1, qcd_geometry *geo, qcd_uint_4 lt, qcd_uint_4 elem){
							
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
    lz = v3 % geo->lL[3];
    ly = ((v3 - lz)/geo->lL[3]) % geo->lL[2];
    lx = ((v3 - lz)/geo->lL[3] - ly) / geo->lL[2];
    v =  qcd_LEXIC(lt,lx,ly,lz,geo->lL);

	
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

	      //-- 2-1								
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][ga][bep][c][bp],
													qcd_CMUL(propf1->D[v][be][alp][b][ap],
														 seqpropf1->D[v][al][gap][a][cp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2])
						  );
								
	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][ga][bep][c][bp],
													qcd_CMUL(propf1->D[v][be][gap][b][cp],
														 seqpropf1->D[v][al][alp][a][ap]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2])
						  );				
										    
	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][ga][bep][c][bp],
													qcd_CMUL(propf1->D[v][al][alp][a][ap],
														 seqpropf1->D[v][be][gap][b][cp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2])
						  );

	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][ga][bep][c][bp],
													qcd_CMUL(propf1->D[v][al][gap][a][cp],
														 seqpropf1->D[v][be][alp][b][ap]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2])
						  );				
				 
	    }//color2 loop     
	  }//color1 loop                                                                               
	}//spin loop (cg_cg)
      }//ga gap
  }//space loop                                                                                   

  return (1);
}
//======================================================================

int qcd_xistar_extra22_f1_3pt(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[9][4][4],
		      qcd_propagator *propf1, qcd_propagator *propf2, qcd_propagator *seqpropf1, qcd_geometry *geo, qcd_uint_4 lt, qcd_uint_4 elem){
							
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
    lz = v3 % geo->lL[3];
    ly = ((v3 - lz)/geo->lL[3]) % geo->lL[2];
    lx = ((v3 - lz)/geo->lL[3] - ly) / geo->lL[2];
    v =  qcd_LEXIC(lt,lx,ly,lz,geo->lL);

	
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

	      //-- 2-2
	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][ga][gap][c][cp],
													qcd_CMUL(propf1->D[v][be][alp][b][ap],
														 seqpropf1->D[v][al][bep][a][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2])
						  );
								
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][ga][gap][c][cp],
													qcd_CMUL(propf1->D[v][be][bep][b][bp],
														 seqpropf1->D[v][al][alp][a][ap]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2])
						  );				
										    
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][ga][gap][c][cp],
													qcd_CMUL(propf1->D[v][al][alp][a][ap],
														 seqpropf1->D[v][be][bep][b][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2])
						  );

	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][ga][gap][c][cp],
													qcd_CMUL(propf1->D[v][al][bep][a][bp],
														 seqpropf1->D[v][be][alp][b][ap]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2])
						  );				
				 
	    }//color2 loop     
	  }//color1 loop                                                                               
	}//spin loop (cg_cg)
      }//ga gap
  }//space loop                                                                                   

  return (1);
}
//======================================================================

int qcd_deltas_f2_3pt(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[9][4][4],
		      qcd_propagator *propf1, qcd_propagator *propf2, qcd_propagator *seqpropf2, qcd_geometry *geo, qcd_uint_4 lt, qcd_uint_4 elem){
							
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
    lz = v3 % geo->lL[3];
    ly = ((v3 - lz)/geo->lL[3]) % geo->lL[2];
    lx = ((v3 - lz)/geo->lL[3] - ly) / geo->lL[2];
    v =  qcd_LEXIC(lt,lx,ly,lz,geo->lL);

	
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

	      //-- 1-1
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][al][alp][a][ap],
													qcd_CMUL(propf1->D[v][ga][gap][c][cp],
														 seqpropf2->D[v][be][bep][b][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*4.0/3.0)
						  );
				

	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][al][gap][a][cp],
													qcd_CMUL(propf1->D[v][ga][alp][c][ap],
														 seqpropf2->D[v][be][bep][b][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*4.0/3.0)
						  );			
	      //-- 1-2								
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][al][bep][a][bp],
													qcd_CMUL(propf1->D[v][ga][alp][c][ap],
														 seqpropf2->D[v][be][gap][b][cp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*2.0/3.0)
						  );
				

	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][al][alp][a][ap],
													qcd_CMUL(propf1->D[v][ga][bep][c][bp],
														 seqpropf2->D[v][be][gap][b][cp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*2.0/3.0)
						  );			
	      //-- 2-1								
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][al][gap][a][cp],
													qcd_CMUL(propf1->D[v][be][alp][b][ap],
														 seqpropf2->D[v][ga][bep][c][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*2.0/3.0)
						  );
				

	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][al][alp][a][ap],
													qcd_CMUL(propf1->D[v][be][gap][b][cp],
														 seqpropf2->D[v][ga][bep][c][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*2.0/3.0)
						  );				
	      //-- 2-2
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][al][alp][a][ap],
													qcd_CMUL(propf1->D[v][be][bep][b][bp],
														 seqpropf2->D[v][ga][gap][c][cp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]/3.0)
						  );
				

	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][al][bep][a][bp],
													qcd_CMUL(propf1->D[v][be][alp][b][ap],
														 seqpropf2->D[v][ga][gap][c][cp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]/3.0)
						  );				
				 
	    }//color2 loop     
	  }//color1 loop                                                                               
	}//spin loop (cg_cg)
      }//ga gap
  }//space loop                                                                                   

  return (1);
}
//======================================================================

int qcd_xistar_extra11_f2_3pt(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[9][4][4],
		      qcd_propagator *propf1, qcd_propagator *propf2, qcd_propagator *seqpropf2, qcd_geometry *geo, qcd_uint_4 lt, qcd_uint_4 elem){
							
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
    lz = v3 % geo->lL[3];
    ly = ((v3 - lz)/geo->lL[3]) % geo->lL[2];
    lx = ((v3 - lz)/geo->lL[3] - ly) / geo->lL[2];
    v =  qcd_LEXIC(lt,lx,ly,lz,geo->lL);

	
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

	      //-- 1-1
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][al][alp][a][ap],
													qcd_CMUL(propf1->D[v][ga][gap][c][cp],
														 seqpropf2->D[v][be][bep][b][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2])
						  );
				

	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][al][gap][a][cp],
													qcd_CMUL(propf1->D[v][ga][alp][c][ap],
														 seqpropf2->D[v][be][bep][b][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2])
						  );			
				 
	    }//color2 loop     
	  }//color1 loop                                                                               
	}//spin loop (cg_cg)
      }//ga gap
  }//space loop                                                                                   

  return (1);
}
//======================================================================

int qcd_xistar_extra12_f2_3pt(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[9][4][4],
		      qcd_propagator *propf1, qcd_propagator *propf2, qcd_propagator *seqpropf2, qcd_geometry *geo, qcd_uint_4 lt, qcd_uint_4 elem){
							
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
    lz = v3 % geo->lL[3];
    ly = ((v3 - lz)/geo->lL[3]) % geo->lL[2];
    lx = ((v3 - lz)/geo->lL[3] - ly) / geo->lL[2];
    v =  qcd_LEXIC(lt,lx,ly,lz,geo->lL);

	
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

	      //-- 1-2								
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][al][bep][a][bp],
													qcd_CMUL(propf1->D[v][ga][alp][c][ap],
														 seqpropf2->D[v][be][gap][b][cp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2])
						  );
				

	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][al][alp][a][ap],
													qcd_CMUL(propf1->D[v][ga][bep][c][bp],
														 seqpropf2->D[v][be][gap][b][cp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2])
						  );			
				 
	    }//color2 loop     
	  }//color1 loop                                                                               
	}//spin loop (cg_cg)
      }//ga gap
  }//space loop                                                                                   

  return (1);
}
//======================================================================

int qcd_xistar_extra21_f2_3pt(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[9][4][4],
		      qcd_propagator *propf1, qcd_propagator *propf2, qcd_propagator *seqpropf2, qcd_geometry *geo, qcd_uint_4 lt, qcd_uint_4 elem){
							
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
    lz = v3 % geo->lL[3];
    ly = ((v3 - lz)/geo->lL[3]) % geo->lL[2];
    lx = ((v3 - lz)/geo->lL[3] - ly) / geo->lL[2];
    v =  qcd_LEXIC(lt,lx,ly,lz,geo->lL);

	
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

	      //-- 2-1								
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][al][gap][a][cp],
													qcd_CMUL(propf1->D[v][be][alp][b][ap],
														 seqpropf2->D[v][ga][bep][c][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2])
						  );
				

	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][al][alp][a][ap],
													qcd_CMUL(propf1->D[v][be][gap][b][cp],
														 seqpropf2->D[v][ga][bep][c][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2])
						  );				
				 
	    }//color2 loop     
	  }//color1 loop                                                                               
	}//spin loop (cg_cg)
      }//ga gap
  }//space loop                                                                                   

  return (1);
}
//======================================================================

int qcd_xistar_extra22_f2_3pt(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[9][4][4],
		      qcd_propagator *propf1, qcd_propagator *propf2, qcd_propagator *seqpropf2, qcd_geometry *geo, qcd_uint_4 lt, qcd_uint_4 elem){
							
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
    lz = v3 % geo->lL[3];
    ly = ((v3 - lz)/geo->lL[3]) % geo->lL[2];
    lx = ((v3 - lz)/geo->lL[3] - ly) / geo->lL[2];
    v =  qcd_LEXIC(lt,lx,ly,lz,geo->lL);

	
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

	      //-- 2-2
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][al][alp][a][ap],
													qcd_CMUL(propf1->D[v][be][bep][b][bp],
														 seqpropf2->D[v][ga][gap][c][cp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2])
						  );
				

	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][al][bep][a][bp],
													qcd_CMUL(propf1->D[v][be][alp][b][ap],
														 seqpropf2->D[v][ga][gap][c][cp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2])
						  );				
				 
	    }//color2 loop     
	  }//color1 loop                                                                               
	}//spin loop (cg_cg)
      }//ga gap
  }//space loop                                                                                   

  return (1);
}
//======================================================================

int qcd_f123f321_f1_3pt(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[9][4][4],
			qcd_propagator *propf2, qcd_propagator *propf3, qcd_propagator *seqpropf1, qcd_geometry *geo, qcd_uint_4 lt, qcd_uint_4 elem){
							
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
    lz = v3 % geo->lL[3];
    ly = ((v3 - lz)/geo->lL[3]) % geo->lL[2];
    lx = ((v3 - lz)/geo->lL[3] - ly) / geo->lL[2];
    v =  qcd_LEXIC(lt,lx,ly,lz,geo->lL);

	
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

	      //-- 1-1
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][be][bep][b][bp],
													qcd_CMUL(propf3->D[v][ga][gap][c][cp],
														 seqpropf1->D[v][al][alp][a][ap]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*0.5)
						  );		
	      //-- 1-2								
	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][be][bep][b][bp],
													qcd_CMUL(propf3->D[v][ga][alp][c][ap],
														 seqpropf1->D[v][al][gap][a][cp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*0.5)
						  );			
	      //-- 2-1								
	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][be][bep][b][bp],
													qcd_CMUL(propf3->D[v][al][gap][a][cp],
														 seqpropf1->D[v][ga][alp][c][ap]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*0.5)
						  );				
	      //-- 2-2
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][be][bep][b][bp],
													qcd_CMUL(propf3->D[v][al][alp][a][ap],
														 seqpropf1->D[v][ga][gap][c][cp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*0.5)
						  );				
				 
	    }//color2 loop     
	  }//color1 loop                                                                               
	}//spin loop (cg_cg)
      }//ga gap
  }//space loop
                                                                              

  return (1);
}
//======================================================================

int qcd_f123f321_f2_3pt(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[9][4][4],
			qcd_propagator *propf1, qcd_propagator *propf3, qcd_propagator *seqpropf2, qcd_geometry *geo, qcd_uint_4 lt, qcd_uint_4 elem){
							
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
    lz = v3 % geo->lL[3];
    ly = ((v3 - lz)/geo->lL[3]) % geo->lL[2];
    lx = ((v3 - lz)/geo->lL[3] - ly) / geo->lL[2];
    v =  qcd_LEXIC(lt,lx,ly,lz,geo->lL);

	
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

	      //-- 1-1
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][al][alp][a][ap],
													qcd_CMUL(propf3->D[v][ga][gap][c][cp],
														 seqpropf2->D[v][be][bep][b][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*0.5)
						  );		
	      //-- 1-2								
	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][al][gap][a][cp],
													qcd_CMUL(propf3->D[v][ga][alp][c][ap],
														 seqpropf2->D[v][be][bep][b][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*0.5)
						  );			
	      //-- 2-1								
	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][ga][alp][c][ap],
													qcd_CMUL(propf3->D[v][al][gap][a][cp],
														 seqpropf2->D[v][be][bep][b][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*0.5)
						  );				
	      //-- 2-2
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][ga][gap][c][cp],
													qcd_CMUL(propf3->D[v][al][alp][a][ap],
														 seqpropf2->D[v][be][bep][b][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*0.5)
						  );				
				 
	    }//color2 loop     
	  }//color1 loop                                                                               
	}//spin loop (cg_cg)
      }//ga gap
  }//space loop 
                                                                                  

  return (1);
}
//======================================================================

int qcd_f123f321_f3_3pt(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[9][4][4],
			qcd_propagator *propf1, qcd_propagator *propf2, qcd_propagator *seqpropf3, qcd_geometry *geo, qcd_uint_4 lt, qcd_uint_4 elem){
							
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
    lz = v3 % geo->lL[3];
    ly = ((v3 - lz)/geo->lL[3]) % geo->lL[2];
    lx = ((v3 - lz)/geo->lL[3] - ly) / geo->lL[2];
    v =  qcd_LEXIC(lt,lx,ly,lz,geo->lL);

	
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

	      //-- 1-1
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][al][alp][a][ap],
													qcd_CMUL(propf2->D[v][be][bep][b][bp],
														 seqpropf3->D[v][ga][gap][c][cp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*0.5)
						  );		
	      //-- 1-2								
	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][al][gap][a][cp],
													qcd_CMUL(propf2->D[v][be][bep][b][bp],
														 seqpropf3->D[v][ga][alp][c][ap]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*0.5)
						  );			
	      //-- 2-1								
	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][ga][alp][c][ap],
													qcd_CMUL(propf2->D[v][be][bep][b][bp],
														 seqpropf3->D[v][al][gap][a][cp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*0.5)
						  );				
	      //-- 2-2
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][ga][gap][c][cp],
													qcd_CMUL(propf2->D[v][be][bep][b][bp],
														 seqpropf3->D[v][al][alp][a][ap]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*0.5)
						  );				
				 
	    }//color2 loop     
	  }//color1 loop                                                                               
	}//spin loop (cg_cg)
      }//ga gap
  }//space loop                                                                                  

  return (1);
}
//======================================================================

int qcd_f1f2f3_3pt(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[9][4][4],
		   qcd_propagator *propf1, qcd_propagator *propf2, qcd_propagator *propf3, qcd_geometry *geo, qcd_uint_4 lt, qcd_uint_4 elem){
							
  qcd_int_4 ctr2;
  qcd_uint_2 a,b,c,ap,bp,cp,al,be,ga,alp,bep,gap,cc1,cc2;
  qcd_uint_4 lx,ly,lz,v,v3;
  int lv3 = geo->lL[1]*geo->lL[2]*geo->lL[3];
  	
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
                
#pragma omp parallel for private(lz,ly,lx,v)
	    for(v3=0; v3<lv3; v3++) {
	      /* for(lx=0; lx<geo->lL[1]; lx++) */
	      /*   for(ly=0; ly<geo->lL[2]; ly++) */
	      /*     for(lz=0; lz<geo->lL[3]; lz++){ */
	      //v3 = qcd_LEXIC0(lx,ly,lz,geo->lL);
	      //qcd_LEXIC(lt,lx,ly,lz,geo->lL);
	      lz = v3 % geo->lL[3];
	      ly = ((v3 - lz)/geo->lL[3]) % geo->lL[2];
	      lx = ((v3 - lz)/geo->lL[3] - ly) / geo->lL[2];
	      v =  qcd_LEXIC(lt,lx,ly,lz,geo->lL);

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

int qcd_sigmas4_f1_3pt(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[9][4][4],
		       qcd_propagator *propf1, qcd_propagator *propf2, qcd_propagator *seqpropf1, qcd_geometry *geo, qcd_uint_4 lt, qcd_uint_4 elem){
							
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
    lz = v3 % geo->lL[3];
    ly = ((v3 - lz)/geo->lL[3]) % geo->lL[2];
    lx = ((v3 - lz)/geo->lL[3] - ly) / geo->lL[2];
    v =  qcd_LEXIC(lt,lx,ly,lz,geo->lL);


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

	      //-- 1-1  
	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][be][alp][b][ap],
													qcd_CMUL(propf2->D[v][ga][gap][c][cp],
														 seqpropf1->D[v][al][bep][a][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]/3.0)
						  );

	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][be][bep][b][bp],
													qcd_CMUL(propf2->D[v][ga][gap][c][cp],
														 seqpropf1->D[v][al][alp][a][ap]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]/3.0)
						  );
				 
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][al][alp][a][ap],
													qcd_CMUL(propf2->D[v][ga][gap][c][cp],
														 seqpropf1->D[v][be][bep][b][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]/3.0)
						  );				 

	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][al][bep][a][bp],
													qcd_CMUL(propf2->D[v][ga][gap][c][cp],
														 seqpropf1->D[v][be][alp][b][ap]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]/3.0)
						  );
	      //-- 1-2
	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][be][bep][b][bp],
													qcd_CMUL(propf2->D[v][ga][alp][c][ap],
														 seqpropf1->D[v][al][gap][a][cp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*2.0/3.0)
						  );				 
				
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][be][gap][b][cp],
													qcd_CMUL(propf2->D[v][ga][alp][c][ap],
														 seqpropf1->D[v][al][bep][a][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*2.0/3.0)
						  );				 

	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][al][bep][a][bp],
													qcd_CMUL(propf2->D[v][ga][alp][c][ap],
														 seqpropf1->D[v][be][gap][b][cp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*2.0/3.0)
						  );

	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][al][gap][a][cp],
													qcd_CMUL(propf2->D[v][ga][alp][c][ap],
														 seqpropf1->D[v][be][bep][b][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*2.0/3.0)
						  );
	      //-- 2-1
	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][ga][alp][c][ap],
													qcd_CMUL(propf2->D[v][al][gap][a][cp],
														 seqpropf1->D[v][be][bep][b][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*2.0/3.0)
						  );
				 
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][ga][bep][c][bp],
													qcd_CMUL(propf2->D[v][al][gap][a][cp],
														 seqpropf1->D[v][be][alp][b][ap]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*2.0/3.0)
						  );
				 
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][be][alp][b][ap],
													qcd_CMUL(propf2->D[v][al][gap][a][cp],
														 seqpropf1->D[v][ga][bep][c][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*2.0/3.0)
						  );				 			 				 

	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][be][bep][b][bp],
													qcd_CMUL(propf2->D[v][al][gap][a][cp],
														 seqpropf1->D[v][ga][alp][c][ap]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*2.0/3.0)
						  );
				 
	      //-- 2-2				 
	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][ga][bep][c][bp],
													qcd_CMUL(propf2->D[v][al][alp][a][ap],
														 seqpropf1->D[v][be][gap][b][cp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*4.0/3.0)
						  );
				 
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][ga][gap][c][cp],
													qcd_CMUL(propf2->D[v][al][alp][a][ap],
														 seqpropf1->D[v][be][bep][b][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*4.0/3.0)
						  );				 

	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][be][bep][b][bp],
													qcd_CMUL(propf2->D[v][al][alp][a][ap],
														 seqpropf1->D[v][ga][gap][c][cp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*4.0/3.0)
						  );	

	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][be][gap][b][cp],
													qcd_CMUL(propf2->D[v][al][alp][a][ap],
														 seqpropf1->D[v][ga][bep][c][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*4.0/3.0)
						  );
				 				 
	    }//color2 loop                                            
	  }//color1 loop                                                                 
	}//spin loop (cg_cg)                                                          
      }//ga gap indices
  }//space loop                                                                                                                                                       
	
  return (1);
}
//======================================================================

int qcd_sigmas4_f2_3pt(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[9][4][4],
		       qcd_propagator *propf1, qcd_propagator *propf2, qcd_propagator *seqpropf2, qcd_geometry *geo, qcd_uint_4 lt, qcd_uint_4 elem){
							
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
    lz = v3 % geo->lL[3];
    ly = ((v3 - lz)/geo->lL[3]) % geo->lL[2];
    lx = ((v3 - lz)/geo->lL[3] - ly) / geo->lL[2];
    v =  qcd_LEXIC(lt,lx,ly,lz,geo->lL);	


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

	      //-- 1-1  
	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][al][bep][a][bp],
													qcd_CMUL(propf1->D[v][be][alp][b][ap],
														 seqpropf2->D[v][ga][gap][c][cp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]/3.0)
						  );

	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][al][alp][a][ap],
													qcd_CMUL(propf1->D[v][be][bep][b][bp],
														 seqpropf2->D[v][ga][gap][c][cp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]/3.0)
						  );
	      //-- 1-2
	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][al][gap][a][cp],
													qcd_CMUL(propf1->D[v][be][bep][b][bp],
														 seqpropf2->D[v][ga][alp][c][ap]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*2.0/3.0)
						  );

	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][al][bep][a][bp],
													qcd_CMUL(propf1->D[v][be][gap][b][cp],
														 seqpropf2->D[v][ga][alp][c][ap]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*2.0/3.0)
						  );
	      //-- 2-1
	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][be][bep][b][bp],
													qcd_CMUL(propf1->D[v][ga][alp][c][ap],
														 seqpropf2->D[v][al][gap][a][cp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*2.0/3.0)
						  );

	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][be][alp][b][ap],
													qcd_CMUL(propf1->D[v][ga][bep][c][bp],
														 seqpropf2->D[v][al][gap][a][cp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*2.0/3.0)
						  );
				 
	      //-- 2-2				 
	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][be][gap][b][cp],
													qcd_CMUL(propf1->D[v][ga][bep][c][bp],
														 seqpropf2->D[v][al][alp][a][ap]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*4.0/3.0)
						  );

	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][be][bep][b][bp],
													qcd_CMUL(propf1->D[v][ga][gap][c][cp],
														 seqpropf2->D[v][al][alp][a][ap]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*4.0/3.0)
						  );
				 				 
	    }//color2 loop                                            
	  }//color1 loop                                                                 
	}//spin loop (cg_cg)                                                          
      }//ga gap indices
  }//space loop                                                                                                                                                       
	
  return (1);
}
//======================================================================

int qcd_lambdas_f1_3pt(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[9][4][4],
		       qcd_propagator *propf2, qcd_propagator *propf3, qcd_propagator *seqpropf1, qcd_geometry *geo, qcd_uint_4 lt, qcd_uint_4 elem){
							
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
    lz = v3 % geo->lL[3];
    ly = ((v3 - lz)/geo->lL[3]) % geo->lL[2];
    lx = ((v3 - lz)/geo->lL[3] - ly) / geo->lL[2];
    v =  qcd_LEXIC(lt,lx,ly,lz,geo->lL);


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

	      //-- 1-1  
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][be][bep][b][bp],
													qcd_CMUL(propf3->D[v][ga][gap][c][cp],
														 seqpropf1->D[v][al][alp][a][ap]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*4.0/6.0)
						  );
	      //-- 1-2
	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][be][gap][b][cp],
													qcd_CMUL(propf3->D[v][ga][bep][c][bp],
														 seqpropf1->D[v][al][alp][a][ap]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*2.0/6.0)
						  );
	      //-- 1-3
	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][be][alp][b][ap],
													qcd_CMUL(propf3->D[v][ga][bep][c][bp],
														 seqpropf1->D[v][al][gap][a][cp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*2.0/6.0)
						  );				 
	      //-- 2-1
	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][ga][bep][c][bp],
													qcd_CMUL(propf3->D[v][be][gap][b][cp],
														 seqpropf1->D[v][al][alp][a][ap]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*2.0/6.0)
						  );				 
	      //-- 2-2				 
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][ga][gap][c][cp],
													qcd_CMUL(propf3->D[v][be][bep][b][bp],
														 seqpropf1->D[v][al][alp][a][ap]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]/6.0)
						  );
	      //-- 2-3				 
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][ga][alp][c][ap],
													qcd_CMUL(propf3->D[v][be][bep][b][bp],
														 seqpropf1->D[v][al][gap][a][cp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]/6.0)
						  );
	      //-- 3-1				 
	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][al][bep][a][bp],
													qcd_CMUL(propf3->D[v][be][gap][b][cp],
														 seqpropf1->D[v][ga][alp][c][ap]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*2.0/6.0)
						  );				
	      //-- 3-2				 
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][al][gap][a][cp],
													qcd_CMUL(propf3->D[v][be][bep][b][bp],
														 seqpropf1->D[v][ga][alp][c][ap]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]/6.0)
						  );
	      //-- 3-3				 
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][al][alp][a][ap],
													qcd_CMUL(propf3->D[v][be][bep][b][bp],
														 seqpropf1->D[v][ga][gap][c][cp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]/6.0)
						  );								
				 				 
	    }//color2 loop                                            
	  }//color1 loop                                                                 
	}//spin loop (cg_cg)                                                          
      }//ga gap indices
  }//space loop                                                                                                                                                       
	
  return (1);
}
//======================================================================

int qcd_lambdas_f2_3pt(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[9][4][4],
		       qcd_propagator *propf1, qcd_propagator *propf3, qcd_propagator *seqpropf2, qcd_geometry *geo, qcd_uint_4 lt, qcd_uint_4 elem){
							
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
    lz = v3 % geo->lL[3];
    ly = ((v3 - lz)/geo->lL[3]) % geo->lL[2];
    lx = ((v3 - lz)/geo->lL[3] - ly) / geo->lL[2];
    v =  qcd_LEXIC(lt,lx,ly,lz,geo->lL);	


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

	      //-- 1-1  
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][al][alp][a][ap],
													qcd_CMUL(propf3->D[v][ga][gap][c][cp],
														 seqpropf2->D[v][be][bep][b][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*4.0/6.0)
						  );
	      //-- 1-2
	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][al][alp][a][ap],
													qcd_CMUL(propf3->D[v][ga][bep][c][bp],
														 seqpropf2->D[v][be][gap][b][cp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*2.0/6.0)
						  );
	      //-- 1-3
	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][al][gap][a][cp],
													qcd_CMUL(propf3->D[v][ga][bep][c][bp],
														 seqpropf2->D[v][be][alp][b][ap]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*2.0/6.0)
						  );				 
	      //-- 2-1
	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][al][alp][a][ap],
													qcd_CMUL(propf3->D[v][be][gap][b][cp],
														 seqpropf2->D[v][ga][bep][c][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*2.0/6.0)
						  );				 
	      //-- 2-2				 
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][al][alp][a][ap],
													qcd_CMUL(propf3->D[v][be][bep][b][bp],
														 seqpropf2->D[v][ga][gap][c][cp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]/6.0)
						  );
	      //-- 2-3				 
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][al][gap][a][cp],
													qcd_CMUL(propf3->D[v][be][bep][b][bp],
														 seqpropf2->D[v][ga][alp][c][ap]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]/6.0)
						  );
	      //-- 3-1				 
	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][ga][alp][c][ap],
													qcd_CMUL(propf3->D[v][be][gap][b][cp],
														 seqpropf2->D[v][al][bep][a][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*2.0/6.0)
						  );				
	      //-- 3-2				 
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][ga][alp][c][ap],
													qcd_CMUL(propf3->D[v][be][bep][b][bp],
														 seqpropf2->D[v][al][gap][a][cp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]/6.0)
						  );
	      //-- 3-3				 
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][ga][gap][c][cp],
													qcd_CMUL(propf3->D[v][be][bep][b][bp],
														 seqpropf2->D[v][al][alp][a][ap]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]/6.0)
						  );								
				 				 
	    }//color2 loop                                            
	  }//color1 loop                                                                 
	}//spin loop (cg_cg)                                                          
      }//ga gap indices
  }//space loop                                                                                                                                                       
	
  return (1);
}
//======================================================================

int qcd_lambdas_f3_3pt(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[9][4][4],
		       qcd_propagator *propf1, qcd_propagator *propf2, qcd_propagator *seqpropf3, qcd_geometry *geo, qcd_uint_4 lt, qcd_uint_4 elem){
							
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
    lz = v3 % geo->lL[3];
    ly = ((v3 - lz)/geo->lL[3]) % geo->lL[2];
    lx = ((v3 - lz)/geo->lL[3] - ly) / geo->lL[2];
    v =  qcd_LEXIC(lt,lx,ly,lz,geo->lL);	


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

	      //-- 1-1  
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][al][alp][a][ap],
													qcd_CMUL(propf2->D[v][be][bep][b][bp],
														 seqpropf3->D[v][ga][gap][c][cp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*4.0/6.0)
						  );
	      //-- 1-2
	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][al][alp][a][ap],
													qcd_CMUL(propf2->D[v][be][gap][b][cp],
														 seqpropf3->D[v][ga][bep][c][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*2.0/6.0)
						  );
	      //-- 1-3
	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][al][gap][a][cp],
													qcd_CMUL(propf2->D[v][be][alp][b][ap],
														 seqpropf3->D[v][ga][bep][c][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*2.0/6.0)
						  );				 
	      //-- 2-1
	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][al][alp][a][ap],
													qcd_CMUL(propf2->D[v][ga][bep][c][bp],
														 seqpropf3->D[v][be][gap][b][cp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*2.0/6.0)
						  );				 
	      //-- 2-2				 
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][al][alp][a][ap],
													qcd_CMUL(propf2->D[v][ga][gap][c][cp],
														 seqpropf3->D[v][be][bep][b][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]/6.0)
						  );
	      //-- 2-3				 
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][al][gap][a][cp],
													qcd_CMUL(propf2->D[v][ga][alp][c][ap],
														 seqpropf3->D[v][be][bep][b][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]/6.0)
						  );
	      //-- 3-1				 
	      block[elem][ga][gap][v3] = qcd_CSUB(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][ga][alp][c][ap],
													qcd_CMUL(propf2->D[v][al][bep][a][bp],
														 seqpropf3->D[v][be][gap][b][cp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*2.0/6.0)
						  );				
	      //-- 3-2				 
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][ga][alp][c][ap],
													qcd_CMUL(propf2->D[v][al][gap][a][cp],
														 seqpropf3->D[v][be][bep][b][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]/6.0)
						  );
	      //-- 3-3				 
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][ga][gap][c][cp],
													qcd_CMUL(propf2->D[v][al][alp][a][ap],
														 seqpropf3->D[v][be][bep][b][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]/6.0)
						  );								
				 				 
	    }//color2 loop                                            
	  }//color1 loop                                                                 
	}//spin loop (cg_cg)                                                          
      }//ga gap indices
  }//space loop                                                                                                                                                       
	
  return (1);
}
//======================================================================

int qcd_sigmas2_f1_3pt(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[9][4][4],
		       qcd_propagator *propf2, qcd_propagator *propf3, qcd_propagator *seqpropf1, qcd_geometry *geo, qcd_uint_4 lt, qcd_uint_4 elem, qcd_uint_4 xis){
							
  qcd_int_4 ctr2;
  qcd_uint_2 a,b,c,ap,bp,cp,al,be,ga,alp,bep,gap,cc1,cc2;
  qcd_uint_4 lx,ly,lz,v,v3;
  qcd_real_8 fact;
  int lv3 = geo->lL[1]*geo->lL[2]*geo->lL[3];
    	
  if(xis) fact = 1.0/3.0;
  else fact = 2.0/3.0;	
	
#pragma omp parallel for private(lz,ly,lx,v,ga,gap,ctr2,al,be,bep,alp,cc1,a,b,c,cc2,ap,bp,cp)
  for(v3=0; v3<lv3; v3++) {
    /* for(lx=0; lx<geo->lL[1]; lx++) */
    /*   for(ly=0; ly<geo->lL[2]; ly++) */
    /*     for(lz=0; lz<geo->lL[3]; lz++){ */
    //v3 = qcd_LEXIC0(lx,ly,lz,geo->lL);
    //qcd_LEXIC(lt,lx,ly,lz,geo->lL);
    lz = v3 % geo->lL[3];
    ly = ((v3 - lz)/geo->lL[3]) % geo->lL[2];
    lx = ((v3 - lz)/geo->lL[3] - ly) / geo->lL[2];
    v =  qcd_LEXIC(lt,lx,ly,lz,geo->lL);

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

	      //-- 1-1  
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][be][bep][b][bp],
													qcd_CMUL(propf3->D[v][ga][gap][c][cp],
														 seqpropf1->D[v][al][alp][a][ap]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*fact)
						  );
	      //-- 1-2
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][be][alp][b][ap],
													qcd_CMUL(propf3->D[v][ga][bep][c][bp],
														 seqpropf1->D[v][al][gap][a][cp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*fact)
						  );
	      //-- 1-3
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][be][gap][b][cp],
													qcd_CMUL(propf3->D[v][ga][alp][c][ap],
														 seqpropf1->D[v][al][bep][a][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*fact)
						  );				 
	      //-- 2-1
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][al][bep][a][bp],
													qcd_CMUL(propf3->D[v][be][gap][b][cp],
														 seqpropf1->D[v][ga][alp][c][ap]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*fact)
						  );				 
	      //-- 2-2				 
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][al][alp][a][ap],
													qcd_CMUL(propf3->D[v][be][bep][b][bp],
														 seqpropf1->D[v][ga][gap][c][cp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*fact)
						  );
	      //-- 2-3				 
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][al][gap][a][cp],
													qcd_CMUL(propf3->D[v][be][alp][b][ap],
														 seqpropf1->D[v][ga][bep][c][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*fact)
						  );
	      //-- 3-1				 
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][ga][bep][c][bp],
													qcd_CMUL(propf3->D[v][al][gap][a][cp],
														 seqpropf1->D[v][be][alp][b][ap]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*fact)
						  );				
	      //-- 3-2				 
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][ga][alp][c][ap],
													qcd_CMUL(propf3->D[v][al][bep][a][bp],
														 seqpropf1->D[v][be][gap][b][cp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*fact)
						  );
	      //-- 3-3				 
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf2->D[v][ga][gap][c][cp],
													qcd_CMUL(propf3->D[v][al][alp][a][ap],
														 seqpropf1->D[v][be][bep][b][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*fact)
						  );

				 	    
	    }//color2 loop        
	  }//color1 loop                                                                                                   
	}//spin loop (cg_cg)                                           
      }//ga gap
  }//space loop                                                                     
	
  return (1);
}
//======================================================================

int qcd_sigmas2_f2_3pt(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[9][4][4],
		       qcd_propagator *propf1, qcd_propagator *propf3, qcd_propagator *seqpropf2, qcd_geometry *geo, qcd_uint_4 lt, qcd_uint_4 elem, qcd_uint_4 xis){
							
  qcd_int_4 ctr2;
  qcd_uint_2 a,b,c,ap,bp,cp,al,be,ga,alp,bep,gap,cc1,cc2;
  qcd_uint_4 lx,ly,lz,v,v3;
  qcd_real_8 fact;
  int lv3 = geo->lL[1]*geo->lL[2]*geo->lL[3];
    
  if(xis) fact = 1.0/3.0;
  else fact = 2.0/3.0;
	
#pragma omp parallel for private(lz,ly,lx,v,ga,gap,ctr2,al,be,bep,alp,cc1,a,b,c,cc2,ap,bp,cp)
  for(v3=0; v3<lv3; v3++) {
    /* for(lx=0; lx<geo->lL[1]; lx++) */
    /*   for(ly=0; ly<geo->lL[2]; ly++) */
    /*     for(lz=0; lz<geo->lL[3]; lz++){ */
    //v3 = qcd_LEXIC0(lx,ly,lz,geo->lL);
    //qcd_LEXIC(lt,lx,ly,lz,geo->lL);
    lz = v3 % geo->lL[3];
    ly = ((v3 - lz)/geo->lL[3]) % geo->lL[2];
    lx = ((v3 - lz)/geo->lL[3] - ly) / geo->lL[2];
    v =  qcd_LEXIC(lt,lx,ly,lz,geo->lL);

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

	      //-- 1-1  
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][al][alp][a][ap],
													qcd_CMUL(propf3->D[v][ga][gap][c][cp],
														 seqpropf2->D[v][be][bep][b][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*fact)
						  );
	      //-- 1-2
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][al][gap][a][cp],
													qcd_CMUL(propf3->D[v][ga][bep][c][bp],
														 seqpropf2->D[v][be][alp][b][ap]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*fact)
						  );
	      //-- 1-3
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][al][bep][a][bp],
													qcd_CMUL(propf3->D[v][ga][alp][c][ap],
														 seqpropf2->D[v][be][gap][b][cp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*fact)
						  );				 
	      //-- 2-1
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][ga][alp][c][ap],
													qcd_CMUL(propf3->D[v][be][gap][b][cp],
														 seqpropf2->D[v][al][bep][a][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*fact)
						  );				 
	      //-- 2-2				 
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][ga][gap][c][cp],
													qcd_CMUL(propf3->D[v][be][bep][b][bp],
														 seqpropf2->D[v][al][alp][a][ap]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*fact)
						  );
	      //-- 2-3				 
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][ga][bep][c][bp],
													qcd_CMUL(propf3->D[v][be][alp][b][ap],
														 seqpropf2->D[v][al][gap][a][cp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*fact)
						  );
	      //-- 3-1				 
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][be][alp][b][ap],
													qcd_CMUL(propf3->D[v][al][gap][a][cp],
														 seqpropf2->D[v][ga][bep][c][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*fact)
						  );				
	      //-- 3-2				 
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][be][gap][b][cp],
													qcd_CMUL(propf3->D[v][al][bep][a][bp],
														 seqpropf2->D[v][ga][alp][c][ap]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*fact)
						  );
	      //-- 3-3				 
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][be][bep][b][bp],
													qcd_CMUL(propf3->D[v][al][alp][a][ap],
														 seqpropf2->D[v][ga][gap][c][cp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*fact)
						  );

				 	    
	    }//color2 loop        
	  }//color1 loop                                                                                                   
	}//spin loop (cg_cg)                                           
      }//ga gap
  }//space loop                                                                     
	
  return (1);
}
//======================================================================

int qcd_sigmas2_f3_3pt(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[9][4][4],
		       qcd_propagator *propf1, qcd_propagator *propf2, qcd_propagator *seqpropf3, qcd_geometry *geo, qcd_uint_4 lt, qcd_uint_4 elem, qcd_uint_4 xis){
							
  qcd_int_4 ctr2;
  qcd_uint_2 a,b,c,ap,bp,cp,al,be,ga,alp,bep,gap,cc1,cc2;
  qcd_uint_4 lx,ly,lz,v,v3;
  qcd_real_8 fact;
  int lv3 = geo->lL[1]*geo->lL[2]*geo->lL[3];
  
  if (xis) fact = 1.0/3.0;
  else fact = 2.0/3.0;
	
#pragma omp parallel for private(lz,ly,lx,v,ga,gap,ctr2,al,be,bep,alp,cc1,a,b,c,cc2,ap,bp,cp)
  for(v3=0; v3<lv3; v3++) {
    /* for(lx=0; lx<geo->lL[1]; lx++) */
    /*   for(ly=0; ly<geo->lL[2]; ly++) */
    /*     for(lz=0; lz<geo->lL[3]; lz++){ */
    //v3 = qcd_LEXIC0(lx,ly,lz,geo->lL);
    //qcd_LEXIC(lt,lx,ly,lz,geo->lL);
    lz = v3 % geo->lL[3];
    ly = ((v3 - lz)/geo->lL[3]) % geo->lL[2];
    lx = ((v3 - lz)/geo->lL[3] - ly) / geo->lL[2];
    v =  qcd_LEXIC(lt,lx,ly,lz,geo->lL);

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

	      //-- 1-1  
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][al][alp][a][ap],
													qcd_CMUL(propf2->D[v][be][bep][b][bp],
														 seqpropf3->D[v][ga][gap][c][cp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*fact)
						  );
	      //-- 1-2
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][al][gap][a][cp],
													qcd_CMUL(propf2->D[v][be][alp][b][ap],
														 seqpropf3->D[v][ga][bep][c][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*fact)
						  );
	      //-- 1-3
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][al][bep][a][bp],
													qcd_CMUL(propf2->D[v][be][gap][b][cp],
														 seqpropf3->D[v][ga][alp][c][ap]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*fact)
						  );				 
	      //-- 2-1
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][ga][alp][c][ap],
													qcd_CMUL(propf2->D[v][al][bep][a][bp],
														 seqpropf3->D[v][be][gap][b][cp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*fact)
						  );				 
	      //-- 2-2				 
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][ga][gap][c][cp],
													qcd_CMUL(propf2->D[v][al][alp][a][ap],
														 seqpropf3->D[v][be][bep][b][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*fact)
						  );
	      //-- 2-3				 
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][ga][bep][c][bp],
													qcd_CMUL(propf2->D[v][al][gap][a][cp],
														 seqpropf3->D[v][be][alp][b][ap]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*fact)
						  );
	      //-- 3-1				 
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][be][alp][b][ap],
													qcd_CMUL(propf2->D[v][ga][bep][c][bp],
														 seqpropf3->D[v][al][gap][a][cp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*fact)
						  );				
	      //-- 3-2				 
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][be][gap][b][cp],
													qcd_CMUL(propf2->D[v][ga][alp][c][ap],
														 seqpropf3->D[v][al][bep][a][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*fact)
						  );
	      //-- 3-3				 
	      block[elem][ga][gap][v3] = qcd_CADD(block[elem][ga][gap][v3],qcd_CSCALE(
										      qcd_CMUL(cgcg_val[ctr2],
											       qcd_CMUL(propf1->D[v][be][bep][b][bp],
													qcd_CMUL(propf2->D[v][ga][gap][c][cp],
														 seqpropf3->D[v][al][alp][a][ap]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*fact)
						  );

				 	    
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

void qcd_contractions3pt_new(qcd_uint_4 p_id, qcd_uint_4 np, qcd_complex_16 *block[5],
			     qcd_propagator *uprop, qcd_propagator *dprop, qcd_propagator *sprop, qcd_propagator *cprop,
			     qcd_propagator *sequprop, qcd_propagator *seqdprop, qcd_propagator *seqsprop, qcd_propagator *seqcprop,			 
			     qcd_geometry *geo, qcd_uint_4 lt){
				 	 

  qcd_uint_4 ctr,ctr5,ctr1,ctr2,ctr3,ctr12,ctr21,ctr13,ctr31,ctr23,ctr32,v3;
  qcd_uint_2 al,be,alp,bep,ga,check = 0;
  qcd_complex_16 C,C1,C2,C3,C21,C12,C31,C13,C23,C32,C5;
  qcd_complex_16 gamma12[4][4],gamma13[4][4],gamma23[4][4];

  qcd_int_2 cg1cg1_ind[16*16][4],cg2cg2_ind[16*16][4],cg3cg3_ind[16*16][4];
  qcd_int_2 cg1cg2_ind[16*16][4],cg2cg1_ind[16*16][4],cg3cg2_ind[16*16][4],cg2cg3_ind[16*16][4];
  qcd_int_2 cg1cg3_ind[16*16][4],cg3cg1_ind[16*16][4],cg5cg5_ind[16*16][4];

  qcd_int_2 cgcg_ind[16*16][4];
  
  qcd_complex_16 cg1cg1_val[16*16],cg2cg2_val[16*16],cg3cg3_val[16*16];
  qcd_complex_16 cg1cg2_val[16*16],cg2cg1_val[16*16],cg3cg2_val[16*16],cg2cg3_val[16*16];
  qcd_complex_16 cg1cg3_val[16*16],cg3cg1_val[16*16],cg5cg5_val[16*16];
  
  qcd_complex_16 cgcg_val[16*16];  
  
  
  qcd_complex_16 CG[4][4],CG_bar[4][4];
  
  qcd_complex_16 *block_pr[9][4][4];
  
  qcd_uint_4 i,j,k,p,det;
    
  //-- Allocation of local blocks  
  
  if( particles32[p_id] ){
    for(i=0;i<9;i++)
      for(j=0;j<4;j++)
	for(k=0;k<4;k++){	
	  block_pr[i][j][k] = (qcd_complex_16*) malloc(geo->lV3*sizeof(qcd_complex_16));
			
	  if(block_pr[i][j][k]==NULL){
	    printf("Process %d: Block_pr[%d][%d][%d] not properly initialized\n",i,j,k,geo->myid);
	    exit(EXIT_FAILURE);
	  }

	  for(v3=0;v3<(geo->lV3);v3++){
	    block_pr[i][j][k][v3] = (qcd_complex_16) {0,0};
	  }
	}							
  }
  else{
    for(j=0;j<4;j++)
      for(k=0;k<4;k++){	
	block_pr[8][j][k] = (qcd_complex_16*) malloc(geo->lV3*sizeof(qcd_complex_16));
			
	if(block_pr[8][j][k]==NULL){
	  printf("Process %d: Block_pr[8][%d][%d] not properly initialized\n",j,k,geo->myid);
	  exit(EXIT_FAILURE);
	}

	for(v3=0;v3<(geo->lV3);v3++){
	  block_pr[8][j][k][v3] = (qcd_complex_16) {0,0};
	}
      }				
  }

  //if(geo->myid==0) printf("Local blocks allocated properly for particle t=%d\n",lt);
    
  //-- Non zero elements of all combinations of gamma_{1,2,3} matrices and gamma_5

  if ( particles32[p_id] ){

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
	      if(qcd_NORM(C23)>1e-3)
		{
		  cg2cg3_val[ctr23].re = C23.re;
		  cg2cg3_val[ctr23].im = C23.im;
		  cg2cg3_ind[ctr23][0] = al;
		  cg2cg3_ind[ctr23][1] = be;
		  cg2cg3_ind[ctr23][2] = bep;
		  cg2cg3_ind[ctr23][3] = alp;                                                            
		  ctr23++;
		}
	      if(qcd_NORM(C32)>1e-3)
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
	      if(qcd_NORM(C12)>1e-3)
		{
		  cg1cg2_val[ctr12].re = C12.re;
		  cg1cg2_val[ctr12].im = C12.im;
		  cg1cg2_ind[ctr12][0] = al;
		  cg1cg2_ind[ctr12][1] = be;
		  cg1cg2_ind[ctr12][2] = bep;
		  cg1cg2_ind[ctr12][3] = alp;                                                            
		  ctr12++;
		}
	      if(qcd_NORM(C21)>1e-3)
		{
		  cg2cg1_val[ctr21].re = C21.re;
		  cg2cg1_val[ctr21].im = C21.im;
		  cg2cg1_ind[ctr21][0] = al;
		  cg2cg1_ind[ctr21][1] = be;
		  cg2cg1_ind[ctr21][2] = bep;
		  cg2cg1_ind[ctr21][3] = alp;                                                            
		  ctr21++;
		}                  
	    }//-for's  
  }//-particles condition
  else{
    ctr5 = 0;
    
    for(al=0;al<4;al++)
      for(be=0;be<4;be++)
	for(alp=0;alp<4;alp++)
	  for(bep=0;bep<4;bep++)
	    {  
	      C5  = qcd_CMUL(qcd_CGAMMA[5][al][be],qcd_BAR_CGAMMA[5][bep][alp]);
       
	      if(qcd_NORM(C5)>1e-3)
		{
		  cg5cg5_val[ctr5].re = C5.re;
		  cg5cg5_val[ctr5].im = C5.im;
		  cg5cg5_ind[ctr5][0] = al;
		  cg5cg5_ind[ctr5][1] = be;
		  cg5cg5_ind[ctr5][2] = bep;
		  cg5cg5_ind[ctr5][3] = alp;                                                            
		  ctr5++;
		}	  	  
	    }
  }


  //-- Calculate the Charge Conjugation and non zero elements 
  if( (p_id==41) || (p_id==42) ){  
    for(al=0;al<4;al++)
      for(be=0;be<4;be++){
	CG[al][be]     = (qcd_complex_16) {0.0,0.0};
	CG_bar[al][be] = (qcd_complex_16) {0.0,0.0};
	for(ga=0;ga<4;ga++){
	  CG[al][be] = qcd_CADD( CG[al][be],qcd_CMUL(qcd_GAMMA[4][al][ga],qcd_GAMMA[2][ga][be]) );		
	  CG_bar[al][be] = CG[al][be];
	}
      }  
  
    ctr = 0;
    for(al=0;al<4;al++)
      for(be=0;be<4;be++)
	for(alp=0;alp<4;alp++)
	  for(bep=0;bep<4;bep++)
	    {  
	      C  = qcd_CMUL(CG[al][be],CG_bar[bep][alp]);
       
	      if(qcd_NORM(C)>1e-3)
		{
		  cgcg_val[ctr].re = C.re;
		  cgcg_val[ctr].im = C.im;
		  cgcg_ind[ctr][0] = al;
		  cgcg_ind[ctr][1] = be;
		  cgcg_ind[ctr][2] = bep;
		  cgcg_ind[ctr][3] = alp;                                                            
		  ctr++;
		}	  	  
	    }
	
  }//-if

  //-- Multiplication of gamma_1,gamma_2,gamma_3 between each other
  if (particles32[p_id]){
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
  }
  //-----------------------------------------------  

  
  switch(p_id){
  case 1://-proton
    switch(np){
    case 0:
      check += qcd_f1f2f1_f1_3pt(ctr5, cg5cg5_ind, cg5cg5_val, block_pr, uprop, dprop, sequprop, geo, lt,8);				
      break;
    case 1:
      check += qcd_f1f2f1_f2_3pt(ctr5, cg5cg5_ind, cg5cg5_val, block_pr, uprop, dprop, seqdprop, geo, lt,8);
      break;
    }	 	  
    break;
	
  case 2://-neutron
    switch(np){
    case 0:
      check += qcd_f1f2f1_f2_3pt(ctr5, cg5cg5_ind, cg5cg5_val, block_pr, dprop, uprop, sequprop, geo, lt,8);
      break;
    case 1:
      check += qcd_f1f2f1_f1_3pt(ctr5, cg5cg5_ind, cg5cg5_val, block_pr, dprop, uprop, seqdprop, geo, lt,8);
      break;
    }	 	  
    break;
	
  case 3://-lambda
    switch(np){
    case 0:
      check += qcd_lambdas_f1_3pt(ctr5, cg5cg5_ind, cg5cg5_val, block_pr, dprop, sprop, sequprop, geo, lt,8);
      break;
    case 1:
      check += qcd_lambdas_f2_3pt(ctr5, cg5cg5_ind, cg5cg5_val, block_pr, uprop, sprop, seqdprop, geo, lt,8);
      break;
    case 2:
      check += qcd_lambdas_f3_3pt(ctr5, cg5cg5_ind, cg5cg5_val, block_pr, uprop, dprop, seqsprop, geo, lt,8);
      break;			
    }	 	  
    break;		
    
  case 4://-sigma plus
    switch(np){
    case 0:
      check += qcd_f1f2f1_f1_3pt(ctr5, cg5cg5_ind, cg5cg5_val, block_pr, uprop, sprop, sequprop, geo, lt,8);	
      break;
    case 2:
      check += qcd_f1f2f1_f2_3pt(ctr5, cg5cg5_ind, cg5cg5_val, block_pr, uprop, sprop, seqsprop, geo, lt,8);	
      break;
    }	 	  
    break;

  case 5://-sigma zero
    switch(np){
    case 0:
      check += qcd_f123f321_f1_3pt(ctr5, cg5cg5_ind, cg5cg5_val, block_pr, sprop, dprop, sequprop, geo, lt,8);	
      break;
    case 1:
      check += qcd_f123f321_f3_3pt(ctr5, cg5cg5_ind, cg5cg5_val, block_pr, uprop, sprop, seqdprop, geo, lt,8);	
      break;
    case 2:
      check += qcd_f123f321_f2_3pt(ctr5, cg5cg5_ind, cg5cg5_val, block_pr, uprop, dprop, seqsprop, geo, lt,8);	
      break;			
    }	 	  
    break;	

  case 6://-sigma minus
    switch(np){
    case 1:
      check += qcd_f1f2f1_f1_3pt(ctr5, cg5cg5_ind, cg5cg5_val, block_pr, dprop, sprop, seqdprop, geo, lt,8);
      break;
    case 2:
      check += qcd_f1f2f1_f2_3pt(ctr5, cg5cg5_ind, cg5cg5_val, block_pr, dprop, sprop, seqsprop, geo, lt,8);	
      break;
    }	 	  
    break;
	
  case 7://-xi zero
    switch(np){
    case 0:
      check += qcd_f1f2f1_f2_3pt(ctr5, cg5cg5_ind, cg5cg5_val, block_pr, sprop, uprop, sequprop, geo, lt,8);	
      break;
    case 2:
      check += qcd_f1f2f1_f1_3pt(ctr5, cg5cg5_ind, cg5cg5_val, block_pr, sprop, uprop, seqsprop, geo, lt,8);	
      break;
    }	 	  
    break;
	
  case 8://-xi minus
    switch(np){
    case 1:
      check += qcd_f1f2f1_f2_3pt(ctr5, cg5cg5_ind, cg5cg5_val, block_pr, sprop, dprop, seqdprop, geo, lt,8);	
      break;
    case 2:
      check += qcd_f1f2f1_f1_3pt(ctr5, cg5cg5_ind, cg5cg5_val, block_pr, sprop, dprop, seqsprop, geo, lt,8);	
      break;
    } 	  
    break;			
    //-----------------------------------------------------------------------------------	
  case 9://-delta plus plus
    check =  qcd_f1f1f1_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, uprop, sequprop, geo, lt,0);
    check += qcd_f1f1f1_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, uprop, sequprop, geo, lt,1);  
    check += qcd_f1f1f1_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, uprop, sequprop, geo, lt,2);
    check += qcd_f1f1f1_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, uprop, sequprop, geo, lt,3);
    check += qcd_f1f1f1_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, uprop, sequprop, geo, lt,4);  
    check += qcd_f1f1f1_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, uprop, sequprop, geo, lt,5);
    check += qcd_f1f1f1_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, uprop, sequprop, geo, lt,6);
    check += qcd_f1f1f1_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, uprop, sequprop, geo, lt,7);  
    check += qcd_f1f1f1_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, uprop, sequprop, geo, lt,8); 
    break;

  case 10://-delta plus
    switch(np){
    case 0:
      check =  qcd_deltas_f1_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, uprop, dprop, sequprop, geo, lt,0);
      check += qcd_deltas_f1_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, uprop, dprop, sequprop, geo, lt,1);  
      check += qcd_deltas_f1_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, uprop, dprop, sequprop, geo, lt,2);
      check += qcd_deltas_f1_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, uprop, dprop, sequprop, geo, lt,3);
      check += qcd_deltas_f1_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, uprop, dprop, sequprop, geo, lt,4);  
      check += qcd_deltas_f1_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, uprop, dprop, sequprop, geo, lt,5);
      check += qcd_deltas_f1_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, uprop, dprop, sequprop, geo, lt,6);
      check += qcd_deltas_f1_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, uprop, dprop, sequprop, geo, lt,7);  
      check += qcd_deltas_f1_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, uprop, dprop, sequprop, geo, lt,8); 
      break;
    case 1:
      check =  qcd_deltas_f2_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, uprop, dprop, seqdprop, geo, lt,0);
      check += qcd_deltas_f2_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, uprop, dprop, seqdprop, geo, lt,1);  
      check += qcd_deltas_f2_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, uprop, dprop, seqdprop, geo, lt,2);
      check += qcd_deltas_f2_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, uprop, dprop, seqdprop, geo, lt,3);
      check += qcd_deltas_f2_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, uprop, dprop, seqdprop, geo, lt,4);  
      check += qcd_deltas_f2_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, uprop, dprop, seqdprop, geo, lt,5);
      check += qcd_deltas_f2_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, uprop, dprop, seqdprop, geo, lt,6);
      check += qcd_deltas_f2_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, uprop, dprop, seqdprop, geo, lt,7);  
      check += qcd_deltas_f2_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, uprop, dprop, seqdprop, geo, lt,8);  
      break;
    }	
    break;
	   
  case 11://-delta zero
    switch(np){
    case 0:
      check =  qcd_deltas_f2_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, dprop, uprop, sequprop, geo, lt,0);
      check += qcd_deltas_f2_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, dprop, uprop, sequprop, geo, lt,1);  
      check += qcd_deltas_f2_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, dprop, uprop, sequprop, geo, lt,2);
      check += qcd_deltas_f2_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, dprop, uprop, sequprop, geo, lt,3);
      check += qcd_deltas_f2_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, dprop, uprop, sequprop, geo, lt,4);  
      check += qcd_deltas_f2_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, dprop, uprop, sequprop, geo, lt,5);
      check += qcd_deltas_f2_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, dprop, uprop, sequprop, geo, lt,6);
      check += qcd_deltas_f2_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, dprop, uprop, sequprop, geo, lt,7);  
      check += qcd_deltas_f2_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, dprop, uprop, sequprop, geo, lt,8); 
      break;
    case 1:
      check =  qcd_deltas_f1_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, dprop, uprop, seqdprop, geo, lt,0);
      check += qcd_deltas_f1_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, dprop, uprop, seqdprop, geo, lt,1);  
      check += qcd_deltas_f1_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, dprop, uprop, seqdprop, geo, lt,2);
      check += qcd_deltas_f1_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, dprop, uprop, seqdprop, geo, lt,3);
      check += qcd_deltas_f1_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, dprop, uprop, seqdprop, geo, lt,4);  
      check += qcd_deltas_f1_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, dprop, uprop, seqdprop, geo, lt,5);
      check += qcd_deltas_f1_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, dprop, uprop, seqdprop, geo, lt,6);
      check += qcd_deltas_f1_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, dprop, uprop, seqdprop, geo, lt,7);  
      check += qcd_deltas_f1_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, dprop, uprop, seqdprop, geo, lt,8);  
      break;
    }
    break;

  case 12://-delta minus
    check =  qcd_f1f1f1_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, dprop, seqdprop, geo, lt,0);
    check += qcd_f1f1f1_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, dprop, seqdprop, geo, lt,1);  
    check += qcd_f1f1f1_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, dprop, seqdprop, geo, lt,2);
    check += qcd_f1f1f1_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, dprop, seqdprop, geo, lt,3);
    check += qcd_f1f1f1_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, dprop, seqdprop, geo, lt,4);  
    check += qcd_f1f1f1_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, dprop, seqdprop, geo, lt,5);
    check += qcd_f1f1f1_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, dprop, seqdprop, geo, lt,6);
    check += qcd_f1f1f1_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, dprop, seqdprop, geo, lt,7);  
    check += qcd_f1f1f1_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, dprop, seqdprop, geo, lt,8); 
    break;
	  
  case 13://-sigma star plus
    switch(np){
    case 0:
      check =  qcd_sigmas4_f1_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, uprop, sprop, sequprop, geo, lt,0);
      check += qcd_sigmas4_f1_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, uprop, sprop, sequprop, geo, lt,1);  
      check += qcd_sigmas4_f1_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, uprop, sprop, sequprop, geo, lt,2);
      check += qcd_sigmas4_f1_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, uprop, sprop, sequprop, geo, lt,3);
      check += qcd_sigmas4_f1_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, uprop, sprop, sequprop, geo, lt,4);  
      check += qcd_sigmas4_f1_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, uprop, sprop, sequprop, geo, lt,5);
      check += qcd_sigmas4_f1_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, uprop, sprop, sequprop, geo, lt,6);
      check += qcd_sigmas4_f1_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, uprop, sprop, sequprop, geo, lt,7);  
      check += qcd_sigmas4_f1_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, uprop, sprop, sequprop, geo, lt,8); 
      break;
    case 2:
      check =  qcd_sigmas4_f2_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, uprop, sprop, seqsprop, geo, lt,0);
      check += qcd_sigmas4_f2_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, uprop, sprop, seqsprop, geo, lt,1);  
      check += qcd_sigmas4_f2_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, uprop, sprop, seqsprop, geo, lt,2);
      check += qcd_sigmas4_f2_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, uprop, sprop, seqsprop, geo, lt,3);
      check += qcd_sigmas4_f2_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, uprop, sprop, seqsprop, geo, lt,4);  
      check += qcd_sigmas4_f2_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, uprop, sprop, seqsprop, geo, lt,5);
      check += qcd_sigmas4_f2_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, uprop, sprop, seqsprop, geo, lt,6);
      check += qcd_sigmas4_f2_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, uprop, sprop, seqsprop, geo, lt,7);  
      check += qcd_sigmas4_f2_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, uprop, sprop, seqsprop, geo, lt,8);  
      break;
    }
    break;
	  
  case 14://-sigma star zero
    switch(np){
    case 0:
      check =  qcd_sigmas2_f1_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, dprop, sprop, sequprop, geo, lt,0,0);
      check += qcd_sigmas2_f1_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, dprop, sprop, sequprop, geo, lt,1,0);  
      check += qcd_sigmas2_f1_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, dprop, sprop, sequprop, geo, lt,2,0);
      check += qcd_sigmas2_f1_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, dprop, sprop, sequprop, geo, lt,3,0);
      check += qcd_sigmas2_f1_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, dprop, sprop, sequprop, geo, lt,4,0);  
      check += qcd_sigmas2_f1_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, dprop, sprop, sequprop, geo, lt,5,0);
      check += qcd_sigmas2_f1_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, dprop, sprop, sequprop, geo, lt,6,0);
      check += qcd_sigmas2_f1_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, dprop, sprop, sequprop, geo, lt,7,0);  
      check += qcd_sigmas2_f1_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, dprop, sprop, sequprop, geo, lt,8,0); 
      break;
    case 1:
      check =  qcd_sigmas2_f2_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, uprop, sprop, seqdprop, geo, lt,0,0);
      check += qcd_sigmas2_f2_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, uprop, sprop, seqdprop, geo, lt,1,0);  
      check += qcd_sigmas2_f2_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, uprop, sprop, seqdprop, geo, lt,2,0);
      check += qcd_sigmas2_f2_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, uprop, sprop, seqdprop, geo, lt,3,0);
      check += qcd_sigmas2_f2_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, uprop, sprop, seqdprop, geo, lt,4,0);  
      check += qcd_sigmas2_f2_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, uprop, sprop, seqdprop, geo, lt,5,0);
      check += qcd_sigmas2_f2_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, uprop, sprop, seqdprop, geo, lt,6,0);
      check += qcd_sigmas2_f2_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, uprop, sprop, seqdprop, geo, lt,7,0);  
      check += qcd_sigmas2_f2_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, uprop, sprop, seqdprop, geo, lt,8,0);  
      break;
    case 2:
      check =  qcd_sigmas2_f3_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, uprop, dprop, seqsprop, geo, lt,0,0);
      check += qcd_sigmas2_f3_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, uprop, dprop, seqsprop, geo, lt,1,0);  
      check += qcd_sigmas2_f3_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, uprop, dprop, seqsprop, geo, lt,2,0);
      check += qcd_sigmas2_f3_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, uprop, dprop, seqsprop, geo, lt,3,0);
      check += qcd_sigmas2_f3_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, uprop, dprop, seqsprop, geo, lt,4,0);  
      check += qcd_sigmas2_f3_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, uprop, dprop, seqsprop, geo, lt,5,0);
      check += qcd_sigmas2_f3_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, uprop, dprop, seqsprop, geo, lt,6,0);
      check += qcd_sigmas2_f3_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, uprop, dprop, seqsprop, geo, lt,7,0);  
      check += qcd_sigmas2_f3_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, uprop, dprop, seqsprop, geo, lt,8,0);  
      break;			
    }
    break;
 	  
  case 15://-sigma star minus
    switch(np){
    case 1:
      check =  qcd_sigmas4_f1_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, dprop, sprop, seqdprop, geo, lt,0);
      check += qcd_sigmas4_f1_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, dprop, sprop, seqdprop, geo, lt,1);  
      check += qcd_sigmas4_f1_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, dprop, sprop, seqdprop, geo, lt,2);
      check += qcd_sigmas4_f1_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, dprop, sprop, seqdprop, geo, lt,3);
      check += qcd_sigmas4_f1_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, dprop, sprop, seqdprop, geo, lt,4);  
      check += qcd_sigmas4_f1_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, dprop, sprop, seqdprop, geo, lt,5);
      check += qcd_sigmas4_f1_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, dprop, sprop, seqdprop, geo, lt,6);
      check += qcd_sigmas4_f1_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, dprop, sprop, seqdprop, geo, lt,7);  
      check += qcd_sigmas4_f1_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, dprop, sprop, seqdprop, geo, lt,8); 
      break;
    case 2:
      check =  qcd_sigmas4_f2_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, dprop, sprop, seqsprop, geo, lt,0);
      check += qcd_sigmas4_f2_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, dprop, sprop, seqsprop, geo, lt,1);  
      check += qcd_sigmas4_f2_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, dprop, sprop, seqsprop, geo, lt,2);
      check += qcd_sigmas4_f2_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, dprop, sprop, seqsprop, geo, lt,3);
      check += qcd_sigmas4_f2_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, dprop, sprop, seqsprop, geo, lt,4);  
      check += qcd_sigmas4_f2_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, dprop, sprop, seqsprop, geo, lt,5);
      check += qcd_sigmas4_f2_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, dprop, sprop, seqsprop, geo, lt,6);
      check += qcd_sigmas4_f2_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, dprop, sprop, seqsprop, geo, lt,7);  
      check += qcd_sigmas4_f2_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, dprop, sprop, seqsprop, geo, lt,8);  
      break;
    }
    break;
 
  case 16://-xi star zero
    switch(np){
    case 0:
      check =  qcd_f1f2f1_f2_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, sprop, uprop, sequprop, geo, lt,0);
      check += qcd_f1f2f1_f2_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, sprop, uprop, sequprop, geo, lt,1);  
      check += qcd_f1f2f1_f2_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, sprop, uprop, sequprop, geo, lt,2);
      check += qcd_f1f2f1_f2_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, sprop, uprop, sequprop, geo, lt,3);
      check += qcd_f1f2f1_f2_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, sprop, uprop, sequprop, geo, lt,4);  
      check += qcd_f1f2f1_f2_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, sprop, uprop, sequprop, geo, lt,5);
      check += qcd_f1f2f1_f2_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, sprop, uprop, sequprop, geo, lt,6);
      check += qcd_f1f2f1_f2_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, sprop, uprop, sequprop, geo, lt,7);  
      check += qcd_f1f2f1_f2_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, sprop, uprop, sequprop, geo, lt,8); 
      break;
    case 2:
      check =  qcd_f1f2f1_f1_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, sprop, uprop, seqsprop, geo, lt,0);
      check += qcd_f1f2f1_f1_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, sprop, uprop, seqsprop, geo, lt,1);  
      check += qcd_f1f2f1_f1_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, sprop, uprop, seqsprop, geo, lt,2);
      check += qcd_f1f2f1_f1_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, sprop, uprop, seqsprop, geo, lt,3);
      check += qcd_f1f2f1_f1_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, sprop, uprop, seqsprop, geo, lt,4);  
      check += qcd_f1f2f1_f1_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, sprop, uprop, seqsprop, geo, lt,5);
      check += qcd_f1f2f1_f1_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, sprop, uprop, seqsprop, geo, lt,6);
      check += qcd_f1f2f1_f1_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, sprop, uprop, seqsprop, geo, lt,7);  
      check += qcd_f1f2f1_f1_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, sprop, uprop, seqsprop, geo, lt,8); 
      break;
    }
    break;
	
  case 17://-xi star minus
    switch(np){
    case 1:
      check =  qcd_f1f2f1_f2_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, sprop, dprop, seqdprop, geo, lt,0);
      check += qcd_f1f2f1_f2_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, sprop, dprop, seqdprop, geo, lt,1);  
      check += qcd_f1f2f1_f2_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, sprop, dprop, seqdprop, geo, lt,2);
      check += qcd_f1f2f1_f2_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, sprop, dprop, seqdprop, geo, lt,3);
      check += qcd_f1f2f1_f2_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, sprop, dprop, seqdprop, geo, lt,4);  
      check += qcd_f1f2f1_f2_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, sprop, dprop, seqdprop, geo, lt,5);
      check += qcd_f1f2f1_f2_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, sprop, dprop, seqdprop, geo, lt,6);
      check += qcd_f1f2f1_f2_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, sprop, dprop, seqdprop, geo, lt,7);  
      check += qcd_f1f2f1_f2_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, sprop, dprop, seqdprop, geo, lt,8); 
      break;
    case 2:
      check =  qcd_f1f2f1_f1_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, sprop, dprop, seqsprop, geo, lt,0);
      check += qcd_f1f2f1_f1_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, sprop, dprop, seqsprop, geo, lt,1);  
      check += qcd_f1f2f1_f1_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, sprop, dprop, seqsprop, geo, lt,2);
      check += qcd_f1f2f1_f1_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, sprop, dprop, seqsprop, geo, lt,3);
      check += qcd_f1f2f1_f1_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, sprop, dprop, seqsprop, geo, lt,4);  
      check += qcd_f1f2f1_f1_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, sprop, dprop, seqsprop, geo, lt,5);
      check += qcd_f1f2f1_f1_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, sprop, dprop, seqsprop, geo, lt,6);
      check += qcd_f1f2f1_f1_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, sprop, dprop, seqsprop, geo, lt,7);  
      check += qcd_f1f2f1_f1_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, sprop, dprop, seqsprop, geo, lt,8); 
      break;
    }
    break;
 	  
  case 18://-omega
    check =  qcd_f1f1f1_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, sprop, seqsprop, geo, lt,0);
    check += qcd_f1f1f1_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, sprop, seqsprop, geo, lt,1);  
    check += qcd_f1f1f1_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, sprop, seqsprop, geo, lt,2);
    check += qcd_f1f1f1_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, sprop, seqsprop, geo, lt,3);
    check += qcd_f1f1f1_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, sprop, seqsprop, geo, lt,4);  
    check += qcd_f1f1f1_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, sprop, seqsprop, geo, lt,5);
    check += qcd_f1f1f1_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, sprop, seqsprop, geo, lt,6);
    check += qcd_f1f1f1_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, sprop, seqsprop, geo, lt,7);  
    check += qcd_f1f1f1_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, sprop, seqsprop, geo, lt,8); 
    break;
    //--------------------------------------------------------------------------------	
  case 19://-lambda c plus
    switch(np){
    case 0:
      check += qcd_lambdas_f1_3pt(ctr5, cg5cg5_ind, cg5cg5_val, block_pr, dprop, cprop, sequprop, geo, lt,8);	
      break;
    case 1:
      check += qcd_lambdas_f2_3pt(ctr5, cg5cg5_ind, cg5cg5_val, block_pr, uprop, cprop, seqdprop, geo, lt,8);	
      break;
    case 3:
      check += qcd_lambdas_f3_3pt(ctr5, cg5cg5_ind, cg5cg5_val, block_pr, uprop, dprop, seqcprop, geo, lt,8);	
      break;			
    }	 	  
    break;

  case 20://-sigma c plus plus
    switch(np){
    case 0:
      check += qcd_f1f2f1_f1_3pt(ctr5, cg5cg5_ind, cg5cg5_val, block_pr, uprop, cprop, sequprop, geo, lt,8);	
      break;
    case 3:
      check += qcd_f1f2f1_f2_3pt(ctr5, cg5cg5_ind, cg5cg5_val, block_pr, uprop, cprop, seqcprop, geo, lt,8);	
      break;
    }	 	  
    break;
	
  case 21://-sigma c plus
    switch(np){
    case 0:
      check += qcd_f123f321_f1_3pt(ctr5, cg5cg5_ind, cg5cg5_val, block_pr, cprop, dprop, sequprop, geo, lt,8);	
      break;
    case 1:
      check += qcd_f123f321_f3_3pt(ctr5, cg5cg5_ind, cg5cg5_val, block_pr, uprop, cprop, seqdprop, geo, lt,8);
      break;	
    case 3:
      check += qcd_f123f321_f2_3pt(ctr5, cg5cg5_ind, cg5cg5_val, block_pr, uprop, dprop, seqcprop, geo, lt,8);					
      break;
    }	 	  
    break;	

  case 22://-sigma c zero
    switch(np){
    case 1:
      check += qcd_f1f2f1_f1_3pt(ctr5, cg5cg5_ind, cg5cg5_val, block_pr, dprop, cprop, seqdprop, geo, lt,8);
      break;
    case 3:
      check += qcd_f1f2f1_f2_3pt(ctr5, cg5cg5_ind, cg5cg5_val, block_pr, dprop, cprop, seqcprop, geo, lt,8);	
      break;
    }	 	  
    break;

  case 23://-xi c plus
    switch(np){
    case 0:
      check += qcd_f1f2f3_3pt(ctr5, cg5cg5_ind, cg5cg5_val, block_pr, sequprop, sprop, cprop, geo, lt,8);	
      break;			
    case 2:
      check += qcd_f1f2f3_3pt(ctr5, cg5cg5_ind, cg5cg5_val, block_pr, uprop, seqsprop, cprop, geo, lt,8);	
      break;
    case 3:
      check += qcd_f1f2f3_3pt(ctr5, cg5cg5_ind, cg5cg5_val, block_pr, uprop, sprop, seqcprop, geo, lt,8);	
      break;
    }	 	  
    break;
  	
  case 24://-xi prime c plus
    switch(np){
    case 0:
      check += qcd_f123f321_f1_3pt(ctr5, cg5cg5_ind, cg5cg5_val, block_pr, cprop, sprop, sequprop, geo, lt,8);	
      break;
    case 2:
      check += qcd_f123f321_f3_3pt(ctr5, cg5cg5_ind, cg5cg5_val, block_pr, uprop, cprop, seqsprop, geo, lt,8);
      break;
    case 3:
      check += qcd_f123f321_f2_3pt(ctr5, cg5cg5_ind, cg5cg5_val, block_pr, uprop, sprop, seqcprop, geo, lt,8);					
      break;
    }	 	  
    break;	
  
  case 25://-xi c zero
    switch(np){
    case 1:
      check += qcd_f1f2f3_3pt(ctr5, cg5cg5_ind, cg5cg5_val, block_pr, seqdprop, sprop, cprop, geo, lt,8);	
      break;			
    case 2:
      check += qcd_f1f2f3_3pt(ctr5, cg5cg5_ind, cg5cg5_val, block_pr, dprop, seqsprop, cprop, geo, lt,8);	
      break;
    case 3:
      check += qcd_f1f2f3_3pt(ctr5, cg5cg5_ind, cg5cg5_val, block_pr, dprop, sprop, seqcprop, geo, lt,8);	
      break;
    }	 	  
    break;

  case 26://-xi prime c zero
    switch(np){
    case 1:
      check += qcd_f123f321_f1_3pt(ctr5, cg5cg5_ind, cg5cg5_val, block_pr, cprop, sprop, seqdprop, geo, lt,8);	
      break;
    case 2:
      check += qcd_f123f321_f3_3pt(ctr5, cg5cg5_ind, cg5cg5_val, block_pr, dprop, cprop, seqsprop, geo, lt,8);
      break;	
    case 3:
      check += qcd_f123f321_f2_3pt(ctr5, cg5cg5_ind, cg5cg5_val, block_pr, dprop, sprop, seqcprop, geo, lt,8);					
      break;
    }	 	  
    break;
	
  case 27://-omega c zero
    switch(np){
    case 2:
      check += qcd_f1f2f1_f1_3pt(ctr5, cg5cg5_ind, cg5cg5_val, block_pr, sprop, cprop, seqsprop, geo, lt,8);	
      break;
    case 3:
      check += qcd_f1f2f1_f2_3pt(ctr5, cg5cg5_ind, cg5cg5_val, block_pr, sprop, cprop, seqcprop, geo, lt,8);	
      break;
    }	 	  
    break;	
    //--------------------------------------------------------------------------------	
  case 28://-sigma star c plus plus
    switch(np){
    case 0:
      check =  qcd_sigmas4_f1_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, uprop, cprop, sequprop, geo, lt,0);
      check += qcd_sigmas4_f1_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, uprop, cprop, sequprop, geo, lt,1);  
      check += qcd_sigmas4_f1_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, uprop, cprop, sequprop, geo, lt,2);
      check += qcd_sigmas4_f1_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, uprop, cprop, sequprop, geo, lt,3);
      check += qcd_sigmas4_f1_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, uprop, cprop, sequprop, geo, lt,4);  
      check += qcd_sigmas4_f1_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, uprop, cprop, sequprop, geo, lt,5);
      check += qcd_sigmas4_f1_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, uprop, cprop, sequprop, geo, lt,6);
      check += qcd_sigmas4_f1_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, uprop, cprop, sequprop, geo, lt,7);  
      check += qcd_sigmas4_f1_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, uprop, cprop, sequprop, geo, lt,8); 
      break;
    case 3:
      check =  qcd_sigmas4_f2_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, uprop, cprop, seqcprop, geo, lt,0);
      check += qcd_sigmas4_f2_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, uprop, cprop, seqcprop, geo, lt,1);  
      check += qcd_sigmas4_f2_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, uprop, cprop, seqcprop, geo, lt,2);
      check += qcd_sigmas4_f2_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, uprop, cprop, seqcprop, geo, lt,3);
      check += qcd_sigmas4_f2_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, uprop, cprop, seqcprop, geo, lt,4);  
      check += qcd_sigmas4_f2_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, uprop, cprop, seqcprop, geo, lt,5);
      check += qcd_sigmas4_f2_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, uprop, cprop, seqcprop, geo, lt,6);
      check += qcd_sigmas4_f2_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, uprop, cprop, seqcprop, geo, lt,7);  
      check += qcd_sigmas4_f2_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, uprop, cprop, seqcprop, geo, lt,8);  
      break;
    }
    break;
 
  case 29://-sigma star c plus
    switch(np){
    case 0:
      check =  qcd_sigmas2_f1_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, dprop, cprop, sequprop, geo, lt,0,0);
      check += qcd_sigmas2_f1_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, dprop, cprop, sequprop, geo, lt,1,0);  
      check += qcd_sigmas2_f1_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, dprop, cprop, sequprop, geo, lt,2,0);
      check += qcd_sigmas2_f1_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, dprop, cprop, sequprop, geo, lt,3,0);
      check += qcd_sigmas2_f1_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, dprop, cprop, sequprop, geo, lt,4,0);  
      check += qcd_sigmas2_f1_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, dprop, cprop, sequprop, geo, lt,5,0);
      check += qcd_sigmas2_f1_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, dprop, cprop, sequprop, geo, lt,6,0);
      check += qcd_sigmas2_f1_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, dprop, cprop, sequprop, geo, lt,7,0);  
      check += qcd_sigmas2_f1_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, dprop, cprop, sequprop, geo, lt,8,0); 
      break;
    case 1:
      check =  qcd_sigmas2_f2_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, uprop, cprop, seqdprop, geo, lt,0,0);
      check += qcd_sigmas2_f2_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, uprop, cprop, seqdprop, geo, lt,1,0);  
      check += qcd_sigmas2_f2_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, uprop, cprop, seqdprop, geo, lt,2,0);
      check += qcd_sigmas2_f2_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, uprop, cprop, seqdprop, geo, lt,3,0);
      check += qcd_sigmas2_f2_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, uprop, cprop, seqdprop, geo, lt,4,0);  
      check += qcd_sigmas2_f2_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, uprop, cprop, seqdprop, geo, lt,5,0);
      check += qcd_sigmas2_f2_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, uprop, cprop, seqdprop, geo, lt,6,0);
      check += qcd_sigmas2_f2_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, uprop, cprop, seqdprop, geo, lt,7,0);  
      check += qcd_sigmas2_f2_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, uprop, cprop, seqdprop, geo, lt,8,0);  
      break;
    case 3:
      check =  qcd_sigmas2_f3_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, uprop, dprop, seqcprop, geo, lt,0,0);
      check += qcd_sigmas2_f3_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, uprop, dprop, seqcprop, geo, lt,1,0);  
      check += qcd_sigmas2_f3_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, uprop, dprop, seqcprop, geo, lt,2,0);
      check += qcd_sigmas2_f3_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, uprop, dprop, seqcprop, geo, lt,3,0);
      check += qcd_sigmas2_f3_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, uprop, dprop, seqcprop, geo, lt,4,0);  
      check += qcd_sigmas2_f3_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, uprop, dprop, seqcprop, geo, lt,5,0);
      check += qcd_sigmas2_f3_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, uprop, dprop, seqcprop, geo, lt,6,0);
      check += qcd_sigmas2_f3_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, uprop, dprop, seqcprop, geo, lt,7,0);  
      check += qcd_sigmas2_f3_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, uprop, dprop, seqcprop, geo, lt,8,0);  
      break;			
    }
    break;
	  
  case 30://-sigma star c zero
    switch(np){
    case 1:
      check =  qcd_sigmas4_f1_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, dprop, cprop, seqdprop, geo, lt,0);
      check += qcd_sigmas4_f1_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, dprop, cprop, seqdprop, geo, lt,1);  
      check += qcd_sigmas4_f1_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, dprop, cprop, seqdprop, geo, lt,2);
      check += qcd_sigmas4_f1_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, dprop, cprop, seqdprop, geo, lt,3);
      check += qcd_sigmas4_f1_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, dprop, cprop, seqdprop, geo, lt,4);  
      check += qcd_sigmas4_f1_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, dprop, cprop, seqdprop, geo, lt,5);
      check += qcd_sigmas4_f1_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, dprop, cprop, seqdprop, geo, lt,6);
      check += qcd_sigmas4_f1_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, dprop, cprop, seqdprop, geo, lt,7);  
      check += qcd_sigmas4_f1_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, dprop, cprop, seqdprop, geo, lt,8); 
      break;
    case 3:
      check =  qcd_sigmas4_f2_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, dprop, cprop, seqcprop, geo, lt,0);
      check += qcd_sigmas4_f2_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, dprop, cprop, seqcprop, geo, lt,1);  
      check += qcd_sigmas4_f2_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, dprop, cprop, seqcprop, geo, lt,2);
      check += qcd_sigmas4_f2_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, dprop, cprop, seqcprop, geo, lt,3);
      check += qcd_sigmas4_f2_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, dprop, cprop, seqcprop, geo, lt,4);  
      check += qcd_sigmas4_f2_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, dprop, cprop, seqcprop, geo, lt,5);
      check += qcd_sigmas4_f2_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, dprop, cprop, seqcprop, geo, lt,6);
      check += qcd_sigmas4_f2_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, dprop, cprop, seqcprop, geo, lt,7);  
      check += qcd_sigmas4_f2_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, dprop, cprop, seqcprop, geo, lt,8);  
      break;
    }
    break;
	  
  case 31://-xi star c plus
    switch(np){
    case 0:
      check =  qcd_f1f2f3_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, sprop, sequprop, cprop, geo, lt,0);
      check += qcd_f1f2f3_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, sprop, sequprop, cprop, geo, lt,1);  
      check += qcd_f1f2f3_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, sprop, sequprop, cprop, geo, lt,2);
      check += qcd_f1f2f3_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, sprop, sequprop, cprop, geo, lt,3);
      check += qcd_f1f2f3_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, sprop, sequprop, cprop, geo, lt,4);  
      check += qcd_f1f2f3_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, sprop, sequprop, cprop, geo, lt,5);
      check += qcd_f1f2f3_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, sprop, sequprop, cprop, geo, lt,6);
      check += qcd_f1f2f3_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, sprop, sequprop, cprop, geo, lt,7);  
      check += qcd_f1f2f3_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, sprop, sequprop, cprop, geo, lt,8); 
      break;	
    case 2:
      check =  qcd_f1f2f3_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, seqsprop, uprop, cprop, geo, lt,0);
      check += qcd_f1f2f3_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, seqsprop, uprop, cprop, geo, lt,1);  
      check += qcd_f1f2f3_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, seqsprop, uprop, cprop, geo, lt,2);
      check += qcd_f1f2f3_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, seqsprop, uprop, cprop, geo, lt,3);
      check += qcd_f1f2f3_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, seqsprop, uprop, cprop, geo, lt,4);  
      check += qcd_f1f2f3_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, seqsprop, uprop, cprop, geo, lt,5);
      check += qcd_f1f2f3_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, seqsprop, uprop, cprop, geo, lt,6);
      check += qcd_f1f2f3_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, seqsprop, uprop, cprop, geo, lt,7);  
      check += qcd_f1f2f3_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, seqsprop, uprop, cprop, geo, lt,8); 
      break;
    case 3:
      check =  qcd_f1f2f3_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, sprop, uprop, seqcprop, geo, lt,0);
      check += qcd_f1f2f3_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, sprop, uprop, seqcprop, geo, lt,1);  
      check += qcd_f1f2f3_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, sprop, uprop, seqcprop, geo, lt,2);
      check += qcd_f1f2f3_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, sprop, uprop, seqcprop, geo, lt,3);
      check += qcd_f1f2f3_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, sprop, uprop, seqcprop, geo, lt,4);  
      check += qcd_f1f2f3_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, sprop, uprop, seqcprop, geo, lt,5);
      check += qcd_f1f2f3_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, sprop, uprop, seqcprop, geo, lt,6);
      check += qcd_f1f2f3_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, sprop, uprop, seqcprop, geo, lt,7);  
      check += qcd_f1f2f3_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, sprop, uprop, seqcprop, geo, lt,8); 
      break;							
    }
    break;
	
  case 32://-xi star c zero
    switch(np){
    case 1:
      check =  qcd_f1f2f3_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, sprop, seqdprop, cprop, geo, lt,0);
      check += qcd_f1f2f3_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, sprop, seqdprop, cprop, geo, lt,1);  
      check += qcd_f1f2f3_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, sprop, seqdprop, cprop, geo, lt,2);
      check += qcd_f1f2f3_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, sprop, seqdprop, cprop, geo, lt,3);
      check += qcd_f1f2f3_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, sprop, seqdprop, cprop, geo, lt,4);  
      check += qcd_f1f2f3_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, sprop, seqdprop, cprop, geo, lt,5);
      check += qcd_f1f2f3_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, sprop, seqdprop, cprop, geo, lt,6);
      check += qcd_f1f2f3_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, sprop, seqdprop, cprop, geo, lt,7);  
      check += qcd_f1f2f3_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, sprop, seqdprop, cprop, geo, lt,8); 
      break;	
    case 2:
      check =  qcd_f1f2f3_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, seqsprop, dprop, cprop, geo, lt,0);
      check += qcd_f1f2f3_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, seqsprop, dprop, cprop, geo, lt,1);  
      check += qcd_f1f2f3_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, seqsprop, dprop, cprop, geo, lt,2);
      check += qcd_f1f2f3_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, seqsprop, dprop, cprop, geo, lt,3);
      check += qcd_f1f2f3_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, seqsprop, dprop, cprop, geo, lt,4);  
      check += qcd_f1f2f3_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, seqsprop, dprop, cprop, geo, lt,5);
      check += qcd_f1f2f3_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, seqsprop, dprop, cprop, geo, lt,6);
      check += qcd_f1f2f3_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, seqsprop, dprop, cprop, geo, lt,7);  
      check += qcd_f1f2f3_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, seqsprop, dprop, cprop, geo, lt,8); 
      break;
    case 3:
      check =  qcd_f1f2f3_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, sprop, dprop, seqcprop, geo, lt,0);
      check += qcd_f1f2f3_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, sprop, dprop, seqcprop, geo, lt,1);  
      check += qcd_f1f2f3_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, sprop, dprop, seqcprop, geo, lt,2);
      check += qcd_f1f2f3_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, sprop, dprop, seqcprop, geo, lt,3);
      check += qcd_f1f2f3_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, sprop, dprop, seqcprop, geo, lt,4);  
      check += qcd_f1f2f3_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, sprop, dprop, seqcprop, geo, lt,5);
      check += qcd_f1f2f3_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, sprop, dprop, seqcprop, geo, lt,6);
      check += qcd_f1f2f3_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, sprop, dprop, seqcprop, geo, lt,7);  
      check += qcd_f1f2f3_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, sprop, dprop, seqcprop, geo, lt,8); 
      break;							
    }
    break;
	  
  case 33://-omega star c zero
    switch(np){
    case 2:
      check =  qcd_f1f2f1_f1_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, sprop, cprop, seqsprop, geo, lt,0);
      check += qcd_f1f2f1_f1_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, sprop, cprop, seqsprop, geo, lt,1);  
      check += qcd_f1f2f1_f1_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, sprop, cprop, seqsprop, geo, lt,2);
      check += qcd_f1f2f1_f1_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, sprop, cprop, seqsprop, geo, lt,3);
      check += qcd_f1f2f1_f1_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, sprop, cprop, seqsprop, geo, lt,4);  
      check += qcd_f1f2f1_f1_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, sprop, cprop, seqsprop, geo, lt,5);
      check += qcd_f1f2f1_f1_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, sprop, cprop, seqsprop, geo, lt,6);
      check += qcd_f1f2f1_f1_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, sprop, cprop, seqsprop, geo, lt,7);  
      check += qcd_f1f2f1_f1_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, sprop, cprop, seqsprop, geo, lt,8); 			
      break;
    case 3:
      check =  qcd_f1f2f1_f2_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, sprop, cprop, seqcprop, geo, lt,0);
      check += qcd_f1f2f1_f2_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, sprop, cprop, seqcprop, geo, lt,1);  
      check += qcd_f1f2f1_f2_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, sprop, cprop, seqcprop, geo, lt,2);
      check += qcd_f1f2f1_f2_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, sprop, cprop, seqcprop, geo, lt,3);
      check += qcd_f1f2f1_f2_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, sprop, cprop, seqcprop, geo, lt,4);  
      check += qcd_f1f2f1_f2_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, sprop, cprop, seqcprop, geo, lt,5);
      check += qcd_f1f2f1_f2_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, sprop, cprop, seqcprop, geo, lt,6);
      check += qcd_f1f2f1_f2_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, sprop, cprop, seqcprop, geo, lt,7);  
      check += qcd_f1f2f1_f2_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, sprop, cprop, seqcprop, geo, lt,8); 
      break;
    }
    break;
    //----------------------------------------------------------------------------------- 	
  case 34://-xi c c plus plus
    switch(np){
    case 0:
      check += qcd_f1f2f1_f2_3pt(ctr5, cg5cg5_ind, cg5cg5_val, block_pr, cprop, uprop, sequprop, geo, lt,8);	
      break;
    case 3:
      check += qcd_f1f2f1_f1_3pt(ctr5, cg5cg5_ind, cg5cg5_val, block_pr, cprop, uprop, seqcprop, geo, lt,8);	
      break;
    }	 	  
    break;
	
  case 35://-xi c c plus
    switch(np){
    case 1:
      check += qcd_f1f2f1_f2_3pt(ctr5, cg5cg5_ind, cg5cg5_val, block_pr, cprop, dprop, seqdprop, geo, lt,8);	
      break;
    case 3:
      check += qcd_f1f2f1_f1_3pt(ctr5, cg5cg5_ind, cg5cg5_val, block_pr, cprop, dprop, seqcprop, geo, lt,8);	
      break;
    } 	  
    break;

  case 36://-omega c c plus
    switch(np){
    case 2:
      check += qcd_f1f2f1_f2_3pt(ctr5, cg5cg5_ind, cg5cg5_val, block_pr, cprop, sprop, seqsprop, geo, lt,8);	
      break;
    case 3:
      check += qcd_f1f2f1_f1_3pt(ctr5, cg5cg5_ind, cg5cg5_val, block_pr, cprop, sprop, seqcprop, geo, lt,8);	
      break;
    } 	  
    break;				
    //-----------------------------------------------------------------------------------       
  case 37://-xi star c c plus plus 
    switch(np){
    case 0:
      check =  qcd_f1f2f1_f2_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, cprop, uprop, sequprop, geo, lt,0);
      check += qcd_f1f2f1_f2_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, cprop, uprop, sequprop, geo, lt,1);  
      check += qcd_f1f2f1_f2_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, cprop, uprop, sequprop, geo, lt,2);
      check += qcd_f1f2f1_f2_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, cprop, uprop, sequprop, geo, lt,3);
      check += qcd_f1f2f1_f2_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, cprop, uprop, sequprop, geo, lt,4);  
      check += qcd_f1f2f1_f2_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, cprop, uprop, sequprop, geo, lt,5);
      check += qcd_f1f2f1_f2_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, cprop, uprop, sequprop, geo, lt,6);
      check += qcd_f1f2f1_f2_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, cprop, uprop, sequprop, geo, lt,7);  
      check += qcd_f1f2f1_f2_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, cprop, uprop, sequprop, geo, lt,8); 
      break;
    case 3:
      check =  qcd_f1f2f1_f1_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, cprop, uprop, seqcprop, geo, lt,0);
      check += qcd_f1f2f1_f1_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, cprop, uprop, seqcprop, geo, lt,1);  
      check += qcd_f1f2f1_f1_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, cprop, uprop, seqcprop, geo, lt,2);
      check += qcd_f1f2f1_f1_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, cprop, uprop, seqcprop, geo, lt,3);
      check += qcd_f1f2f1_f1_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, cprop, uprop, seqcprop, geo, lt,4);  
      check += qcd_f1f2f1_f1_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, cprop, uprop, seqcprop, geo, lt,5);
      check += qcd_f1f2f1_f1_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, cprop, uprop, seqcprop, geo, lt,6);
      check += qcd_f1f2f1_f1_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, cprop, uprop, seqcprop, geo, lt,7);  
      check += qcd_f1f2f1_f1_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, cprop, uprop, seqcprop, geo, lt,8); 
      break;
    }
    break;
	
  case 38://-xi star c c plus
    switch(np){
    case 1:
      check =  qcd_f1f2f1_f2_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, cprop, dprop, seqdprop, geo, lt,0);
      check += qcd_f1f2f1_f2_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, cprop, dprop, seqdprop, geo, lt,1);  
      check += qcd_f1f2f1_f2_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, cprop, dprop, seqdprop, geo, lt,2);
      check += qcd_f1f2f1_f2_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, cprop, dprop, seqdprop, geo, lt,3);
      check += qcd_f1f2f1_f2_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, cprop, dprop, seqdprop, geo, lt,4);  
      check += qcd_f1f2f1_f2_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, cprop, dprop, seqdprop, geo, lt,5);
      check += qcd_f1f2f1_f2_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, cprop, dprop, seqdprop, geo, lt,6);
      check += qcd_f1f2f1_f2_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, cprop, dprop, seqdprop, geo, lt,7);  
      check += qcd_f1f2f1_f2_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, cprop, dprop, seqdprop, geo, lt,8); 
      break;
    case 3:
      check =  qcd_f1f2f1_f1_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, cprop, dprop, seqcprop, geo, lt,0);
      check += qcd_f1f2f1_f1_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, cprop, dprop, seqcprop, geo, lt,1);  
      check += qcd_f1f2f1_f1_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, cprop, dprop, seqcprop, geo, lt,2);
      check += qcd_f1f2f1_f1_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, cprop, dprop, seqcprop, geo, lt,3);
      check += qcd_f1f2f1_f1_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, cprop, dprop, seqcprop, geo, lt,4);  
      check += qcd_f1f2f1_f1_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, cprop, dprop, seqcprop, geo, lt,5);
      check += qcd_f1f2f1_f1_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, cprop, dprop, seqcprop, geo, lt,6);
      check += qcd_f1f2f1_f1_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, cprop, dprop, seqcprop, geo, lt,7);  
      check += qcd_f1f2f1_f1_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, cprop, dprop, seqcprop, geo, lt,8); 
      break;
    }
    break;
	
  case 39://-omega star c c plus
    switch(np){
    case 2:
      check =  qcd_f1f2f1_f2_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, cprop, sprop, seqsprop, geo, lt,0);
      check += qcd_f1f2f1_f2_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, cprop, sprop, seqsprop, geo, lt,1);  
      check += qcd_f1f2f1_f2_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, cprop, sprop, seqsprop, geo, lt,2);
      check += qcd_f1f2f1_f2_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, cprop, sprop, seqsprop, geo, lt,3);
      check += qcd_f1f2f1_f2_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, cprop, sprop, seqsprop, geo, lt,4);  
      check += qcd_f1f2f1_f2_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, cprop, sprop, seqsprop, geo, lt,5);
      check += qcd_f1f2f1_f2_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, cprop, sprop, seqsprop, geo, lt,6);
      check += qcd_f1f2f1_f2_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, cprop, sprop, seqsprop, geo, lt,7);  
      check += qcd_f1f2f1_f2_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, cprop, sprop, seqsprop, geo, lt,8); 
      break;
    case 3:
      check =  qcd_f1f2f1_f1_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, cprop, sprop, seqcprop, geo, lt,0);
      check += qcd_f1f2f1_f1_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, cprop, sprop, seqcprop, geo, lt,1);  
      check += qcd_f1f2f1_f1_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, cprop, sprop, seqcprop, geo, lt,2);
      check += qcd_f1f2f1_f1_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, cprop, sprop, seqcprop, geo, lt,3);
      check += qcd_f1f2f1_f1_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, cprop, sprop, seqcprop, geo, lt,4);  
      check += qcd_f1f2f1_f1_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, cprop, sprop, seqcprop, geo, lt,5);
      check += qcd_f1f2f1_f1_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, cprop, sprop, seqcprop, geo, lt,6);
      check += qcd_f1f2f1_f1_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, cprop, sprop, seqcprop, geo, lt,7);  
      check += qcd_f1f2f1_f1_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, cprop, sprop, seqcprop, geo, lt,8); 

      			
      break;
    }
    break;	
    //-----------------------------------------------------------------------------------	  
  case 40://-omega c c c plus plus
    check =  qcd_f1f1f1_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, cprop, seqcprop, geo, lt,0);
    check += qcd_f1f1f1_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, cprop, seqcprop, geo, lt,1);  
    check += qcd_f1f1f1_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, cprop, seqcprop, geo, lt,2);
    check += qcd_f1f1f1_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, cprop, seqcprop, geo, lt,3);
    check += qcd_f1f1f1_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, cprop, seqcprop, geo, lt,4);  
    check += qcd_f1f1f1_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, cprop, seqcprop, geo, lt,5);
    check += qcd_f1f1f1_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, cprop, seqcprop, geo, lt,6);
    check += qcd_f1f1f1_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, cprop, seqcprop, geo, lt,7);  
    check += qcd_f1f1f1_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, cprop, seqcprop, geo, lt,8); 
    break;
    //-----------------------------------------------------------------------------------	
  case 41://-proton extra
    switch(np){
    case 0:
      check += qcd_pnextra_f1_3pt(ctr, cgcg_ind, cgcg_val, block_pr, uprop, dprop, sequprop, geo, lt,8);				
      break;
    case 1:
      check += qcd_pnextra_f2_3pt(ctr, cgcg_ind, cgcg_val, block_pr, uprop, dprop, seqdprop, geo, lt,8);
      break;
    }	 	  
    break;
	
  case 42://-neutron extra
    switch(np){
    case 0:
      check += qcd_pnextra_f2_3pt(ctr, cgcg_ind, cgcg_val, block_pr, dprop, uprop, sequprop, geo, lt,8);
      break;
    case 1:
      check += qcd_pnextra_f1_3pt(ctr, cgcg_ind, cgcg_val, block_pr, dprop, uprop, seqdprop, geo, lt,8);
      break;
    }	 	  
    break;
    //-----------------------------------------------------------------------------------	
  case 43://-xistar_zero extra
    switch(np){
    case 0:
      check =  qcd_deltas_f2_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, sprop, uprop, sequprop, geo, lt,0);
      check += qcd_deltas_f2_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, sprop, uprop, sequprop, geo, lt,1);  
      check += qcd_deltas_f2_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, sprop, uprop, sequprop, geo, lt,2);
      check += qcd_deltas_f2_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, sprop, uprop, sequprop, geo, lt,3);
      check += qcd_deltas_f2_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, sprop, uprop, sequprop, geo, lt,4);  
      check += qcd_deltas_f2_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, sprop, uprop, sequprop, geo, lt,5);
      check += qcd_deltas_f2_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, sprop, uprop, sequprop, geo, lt,6);
      check += qcd_deltas_f2_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, sprop, uprop, sequprop, geo, lt,7);  
      check += qcd_deltas_f2_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, sprop, uprop, sequprop, geo, lt,8); 
      break;
    case 2:
      check =  qcd_deltas_f1_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, sprop, uprop, seqsprop, geo, lt,0);
      check += qcd_deltas_f1_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, sprop, uprop, seqsprop, geo, lt,1);  
      check += qcd_deltas_f1_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, sprop, uprop, seqsprop, geo, lt,2);
      check += qcd_deltas_f1_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, sprop, uprop, seqsprop, geo, lt,3);
      check += qcd_deltas_f1_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, sprop, uprop, seqsprop, geo, lt,4);  
      check += qcd_deltas_f1_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, sprop, uprop, seqsprop, geo, lt,5);
      check += qcd_deltas_f1_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, sprop, uprop, seqsprop, geo, lt,6);
      check += qcd_deltas_f1_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, sprop, uprop, seqsprop, geo, lt,7);  
      check += qcd_deltas_f1_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, sprop, uprop, seqsprop, geo, lt,8);  
      break;
    }
    break;

  case 44://-xistar_minus extra
    switch(np){
    case 1:
      check =  qcd_deltas_f2_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, sprop, dprop, seqdprop, geo, lt,0);
      check += qcd_deltas_f2_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, sprop, dprop, seqdprop, geo, lt,1);  
      check += qcd_deltas_f2_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, sprop, dprop, seqdprop, geo, lt,2);
      check += qcd_deltas_f2_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, sprop, dprop, seqdprop, geo, lt,3);
      check += qcd_deltas_f2_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, sprop, dprop, seqdprop, geo, lt,4);  
      check += qcd_deltas_f2_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, sprop, dprop, seqdprop, geo, lt,5);
      check += qcd_deltas_f2_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, sprop, dprop, seqdprop, geo, lt,6);
      check += qcd_deltas_f2_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, sprop, dprop, seqdprop, geo, lt,7);  
      check += qcd_deltas_f2_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, sprop, dprop, seqdprop, geo, lt,8); 
      break;
    case 2:
      check =  qcd_deltas_f1_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, sprop, dprop, seqsprop, geo, lt,0);
      check += qcd_deltas_f1_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, sprop, dprop, seqsprop, geo, lt,1);  
      check += qcd_deltas_f1_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, sprop, dprop, seqsprop, geo, lt,2);
      check += qcd_deltas_f1_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, sprop, dprop, seqsprop, geo, lt,3);
      check += qcd_deltas_f1_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, sprop, dprop, seqsprop, geo, lt,4);  
      check += qcd_deltas_f1_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, sprop, dprop, seqsprop, geo, lt,5);
      check += qcd_deltas_f1_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, sprop, dprop, seqsprop, geo, lt,6);
      check += qcd_deltas_f1_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, sprop, dprop, seqsprop, geo, lt,7);  
      check += qcd_deltas_f1_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, sprop, dprop, seqsprop, geo, lt,8);  
      break;
    }
    break;			
    //-----------------------------------------------------------------------------------
  case 45://-xi_cplus extra
    switch(np){
    case 0:
      check += qcd_lambdas_f2_3pt(ctr5, cg5cg5_ind, cg5cg5_val, block_pr, sprop, cprop, sequprop, geo, lt,8);
      break;
    case 2:
      check += qcd_lambdas_f1_3pt(ctr5, cg5cg5_ind, cg5cg5_val, block_pr, uprop, cprop, seqsprop, geo, lt,8);
      break;
    case 3:
      check += qcd_lambdas_f3_3pt(ctr5, cg5cg5_ind, cg5cg5_val, block_pr, sprop, uprop, seqcprop, geo, lt,8);
      break;			
    }	 	  
    break;
	
  case 46://-xi_czero extra
    switch(np){
    case 1:
      check += qcd_lambdas_f2_3pt(ctr5, cg5cg5_ind, cg5cg5_val, block_pr, sprop, cprop, seqdprop, geo, lt,8);
      break;
    case 2:
      check += qcd_lambdas_f1_3pt(ctr5, cg5cg5_ind, cg5cg5_val, block_pr, dprop, cprop, seqsprop, geo, lt,8);
      break;
    case 3:
      check += qcd_lambdas_f3_3pt(ctr5, cg5cg5_ind, cg5cg5_val, block_pr, sprop, dprop, seqcprop, geo, lt,8);
      break;			
    }	 	  
    break;
    //-----------------------------------------------------------------------------------
  case 47://-xistar_cplus extra
    switch(np){
    case 0:
      check =  qcd_sigmas2_f1_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, sprop, cprop, sequprop, geo, lt,0,1);
      check += qcd_sigmas2_f1_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, sprop, cprop, sequprop, geo, lt,1,1);  
      check += qcd_sigmas2_f1_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, sprop, cprop, sequprop, geo, lt,2,1);
      check += qcd_sigmas2_f1_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, sprop, cprop, sequprop, geo, lt,3,1);
      check += qcd_sigmas2_f1_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, sprop, cprop, sequprop, geo, lt,4,1);  
      check += qcd_sigmas2_f1_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, sprop, cprop, sequprop, geo, lt,5,1);
      check += qcd_sigmas2_f1_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, sprop, cprop, sequprop, geo, lt,6,1);
      check += qcd_sigmas2_f1_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, sprop, cprop, sequprop, geo, lt,7,1);  
      check += qcd_sigmas2_f1_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, sprop, cprop, sequprop, geo, lt,8,1); 
      break;
    case 2:
      check =  qcd_sigmas2_f2_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, uprop, cprop, seqsprop, geo, lt,0,1);
      check += qcd_sigmas2_f2_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, uprop, cprop, seqsprop, geo, lt,1,1);  
      check += qcd_sigmas2_f2_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, uprop, cprop, seqsprop, geo, lt,2,1);
      check += qcd_sigmas2_f2_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, uprop, cprop, seqsprop, geo, lt,3,1);
      check += qcd_sigmas2_f2_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, uprop, cprop, seqsprop, geo, lt,4,1);  
      check += qcd_sigmas2_f2_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, uprop, cprop, seqsprop, geo, lt,5,1);
      check += qcd_sigmas2_f2_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, uprop, cprop, seqsprop, geo, lt,6,1);
      check += qcd_sigmas2_f2_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, uprop, cprop, seqsprop, geo, lt,7,1);  
      check += qcd_sigmas2_f2_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, uprop, cprop, seqsprop, geo, lt,8,1);  
      break;
    case 3:
      check =  qcd_sigmas2_f3_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, uprop, sprop, seqcprop, geo, lt,0,1);
      check += qcd_sigmas2_f3_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, uprop, sprop, seqcprop, geo, lt,1,1);  
      check += qcd_sigmas2_f3_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, uprop, sprop, seqcprop, geo, lt,2,1);
      check += qcd_sigmas2_f3_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, uprop, sprop, seqcprop, geo, lt,3,1);
      check += qcd_sigmas2_f3_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, uprop, sprop, seqcprop, geo, lt,4,1);  
      check += qcd_sigmas2_f3_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, uprop, sprop, seqcprop, geo, lt,5,1);
      check += qcd_sigmas2_f3_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, uprop, sprop, seqcprop, geo, lt,6,1);
      check += qcd_sigmas2_f3_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, uprop, sprop, seqcprop, geo, lt,7,1);  
      check += qcd_sigmas2_f3_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, uprop, sprop, seqcprop, geo, lt,8,1);  
      break;			
    }
    break;
	
  case 48://-xistar_czero extra
    switch(np){
    case 1:
      check =  qcd_sigmas2_f1_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, sprop, cprop, seqdprop, geo, lt,0,1);
      check += qcd_sigmas2_f1_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, sprop, cprop, seqdprop, geo, lt,1,1);  
      check += qcd_sigmas2_f1_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, sprop, cprop, seqdprop, geo, lt,2,1);
      check += qcd_sigmas2_f1_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, sprop, cprop, seqdprop, geo, lt,3,1);
      check += qcd_sigmas2_f1_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, sprop, cprop, seqdprop, geo, lt,4,1);  
      check += qcd_sigmas2_f1_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, sprop, cprop, seqdprop, geo, lt,5,1);
      check += qcd_sigmas2_f1_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, sprop, cprop, seqdprop, geo, lt,6,1);
      check += qcd_sigmas2_f1_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, sprop, cprop, seqdprop, geo, lt,7,1);  
      check += qcd_sigmas2_f1_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, sprop, cprop, seqdprop, geo, lt,8,1); 

      
      break;
    case 2:
      check =  qcd_sigmas2_f2_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, dprop, cprop, seqsprop, geo, lt,0,1);
      check += qcd_sigmas2_f2_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, dprop, cprop, seqsprop, geo, lt,1,1);  
      check += qcd_sigmas2_f2_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, dprop, cprop, seqsprop, geo, lt,2,1);
      check += qcd_sigmas2_f2_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, dprop, cprop, seqsprop, geo, lt,3,1);
      check += qcd_sigmas2_f2_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, dprop, cprop, seqsprop, geo, lt,4,1);  
      check += qcd_sigmas2_f2_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, dprop, cprop, seqsprop, geo, lt,5,1);
      check += qcd_sigmas2_f2_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, dprop, cprop, seqsprop, geo, lt,6,1);
      check += qcd_sigmas2_f2_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, dprop, cprop, seqsprop, geo, lt,7,1);  
      check += qcd_sigmas2_f2_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, dprop, cprop, seqsprop, geo, lt,8,1);  
      break;
    case 3:
      check =  qcd_sigmas2_f3_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, dprop, sprop, seqcprop, geo, lt,0,1);
      check += qcd_sigmas2_f3_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, dprop, sprop, seqcprop, geo, lt,1,1);  
      check += qcd_sigmas2_f3_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, dprop, sprop, seqcprop, geo, lt,2,1);
      check += qcd_sigmas2_f3_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, dprop, sprop, seqcprop, geo, lt,3,1);
      check += qcd_sigmas2_f3_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, dprop, sprop, seqcprop, geo, lt,4,1);  
      check += qcd_sigmas2_f3_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, dprop, sprop, seqcprop, geo, lt,5,1);
      check += qcd_sigmas2_f3_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, dprop, sprop, seqcprop, geo, lt,6,1);
      check += qcd_sigmas2_f3_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, dprop, sprop, seqcprop, geo, lt,7,1);  
      check += qcd_sigmas2_f3_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, dprop, sprop, seqcprop, geo, lt,8,1);  
      break;			
    }
    break;
    //-----------------------------------------------------------------------------------
  case 49://-omegastar_czero extra
    switch(np){
    case 2:
      check =  qcd_deltas_f1_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, sprop, cprop, seqsprop, geo, lt,0);
      check += qcd_deltas_f1_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, sprop, cprop, seqsprop, geo, lt,1);  
      check += qcd_deltas_f1_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, sprop, cprop, seqsprop, geo, lt,2);
      check += qcd_deltas_f1_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, sprop, cprop, seqsprop, geo, lt,3);
      check += qcd_deltas_f1_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, sprop, cprop, seqsprop, geo, lt,4);  
      check += qcd_deltas_f1_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, sprop, cprop, seqsprop, geo, lt,5);
      check += qcd_deltas_f1_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, sprop, cprop, seqsprop, geo, lt,6);
      check += qcd_deltas_f1_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, sprop, cprop, seqsprop, geo, lt,7);  
      check += qcd_deltas_f1_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, sprop, cprop, seqsprop, geo, lt,8);			
      break;
    case 3:
      check =  qcd_deltas_f2_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, sprop, cprop, seqcprop, geo, lt,0);
      check += qcd_deltas_f2_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, sprop, cprop, seqcprop, geo, lt,1);  
      check += qcd_deltas_f2_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, sprop, cprop, seqcprop, geo, lt,2);
      check += qcd_deltas_f2_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, sprop, cprop, seqcprop, geo, lt,3);
      check += qcd_deltas_f2_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, sprop, cprop, seqcprop, geo, lt,4);  
      check += qcd_deltas_f2_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, sprop, cprop, seqcprop, geo, lt,5);
      check += qcd_deltas_f2_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, sprop, cprop, seqcprop, geo, lt,6);
      check += qcd_deltas_f2_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, sprop, cprop, seqcprop, geo, lt,7);  
      check += qcd_deltas_f2_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, sprop, cprop, seqcprop, geo, lt,8);   
      break;
    }
    break;	

  case 50://-xistar_ccplusplus extra
    switch(np){
    case 0:
      check =  qcd_deltas_f2_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, cprop, uprop, sequprop, geo, lt,0);
      check += qcd_deltas_f2_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, cprop, uprop, sequprop, geo, lt,1);  
      check += qcd_deltas_f2_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, cprop, uprop, sequprop, geo, lt,2);
      check += qcd_deltas_f2_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, cprop, uprop, sequprop, geo, lt,3);
      check += qcd_deltas_f2_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, cprop, uprop, sequprop, geo, lt,4);  
      check += qcd_deltas_f2_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, cprop, uprop, sequprop, geo, lt,5);
      check += qcd_deltas_f2_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, cprop, uprop, sequprop, geo, lt,6);
      check += qcd_deltas_f2_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, cprop, uprop, sequprop, geo, lt,7);  
      check += qcd_deltas_f2_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, cprop, uprop, sequprop, geo, lt,8); 
      break;
    case 3:
      check =  qcd_deltas_f1_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, cprop, uprop, seqcprop, geo, lt,0);
      check += qcd_deltas_f1_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, cprop, uprop, seqcprop, geo, lt,1);  
      check += qcd_deltas_f1_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, cprop, uprop, seqcprop, geo, lt,2);
      check += qcd_deltas_f1_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, cprop, uprop, seqcprop, geo, lt,3);
      check += qcd_deltas_f1_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, cprop, uprop, seqcprop, geo, lt,4);  
      check += qcd_deltas_f1_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, cprop, uprop, seqcprop, geo, lt,5);
      check += qcd_deltas_f1_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, cprop, uprop, seqcprop, geo, lt,6);
      check += qcd_deltas_f1_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, cprop, uprop, seqcprop, geo, lt,7);  
      check += qcd_deltas_f1_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, cprop, uprop, seqcprop, geo, lt,8);  
      break;
    }
    break;

  case 51://-xistar_ccplus extra
    switch(np){
    case 1:
      check =  qcd_deltas_f2_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, cprop, dprop, seqdprop, geo, lt,0);
      check += qcd_deltas_f2_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, cprop, dprop, seqdprop, geo, lt,1);  
      check += qcd_deltas_f2_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, cprop, dprop, seqdprop, geo, lt,2);
      check += qcd_deltas_f2_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, cprop, dprop, seqdprop, geo, lt,3);
      check += qcd_deltas_f2_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, cprop, dprop, seqdprop, geo, lt,4);  
      check += qcd_deltas_f2_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, cprop, dprop, seqdprop, geo, lt,5);
      check += qcd_deltas_f2_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, cprop, dprop, seqdprop, geo, lt,6);
      check += qcd_deltas_f2_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, cprop, dprop, seqdprop, geo, lt,7);  
      check += qcd_deltas_f2_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, cprop, dprop, seqdprop, geo, lt,8); 
      break;
    case 3:
      check =  qcd_deltas_f1_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, cprop, dprop, seqcprop, geo, lt,0);
      check += qcd_deltas_f1_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, cprop, dprop, seqcprop, geo, lt,1);  
      check += qcd_deltas_f1_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, cprop, dprop, seqcprop, geo, lt,2);
      check += qcd_deltas_f1_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, cprop, dprop, seqcprop, geo, lt,3);
      check += qcd_deltas_f1_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, cprop, dprop, seqcprop, geo, lt,4);  
      check += qcd_deltas_f1_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, cprop, dprop, seqcprop, geo, lt,5);
      check += qcd_deltas_f1_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, cprop, dprop, seqcprop, geo, lt,6);
      check += qcd_deltas_f1_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, cprop, dprop, seqcprop, geo, lt,7);  
      check += qcd_deltas_f1_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, cprop, dprop, seqcprop, geo, lt,8);  
      break;
    }
    break;

  case 52://-omegastar_ccplus extra
    switch(np){
    case 2:
      check =  qcd_deltas_f2_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, cprop, sprop, seqsprop, geo, lt,0);
      check += qcd_deltas_f2_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, cprop, sprop, seqsprop, geo, lt,1);  
      check += qcd_deltas_f2_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, cprop, sprop, seqsprop, geo, lt,2);
      check += qcd_deltas_f2_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, cprop, sprop, seqsprop, geo, lt,3);
      check += qcd_deltas_f2_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, cprop, sprop, seqsprop, geo, lt,4);  
      check += qcd_deltas_f2_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, cprop, sprop, seqsprop, geo, lt,5);
      check += qcd_deltas_f2_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, cprop, sprop, seqsprop, geo, lt,6);
      check += qcd_deltas_f2_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, cprop, sprop, seqsprop, geo, lt,7);  
      check += qcd_deltas_f2_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, cprop, sprop, seqsprop, geo, lt,8); 
      break;
    case 3:
      check =  qcd_deltas_f1_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, cprop, sprop, seqcprop, geo, lt,0);
      check += qcd_deltas_f1_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, cprop, sprop, seqcprop, geo, lt,1);  
      check += qcd_deltas_f1_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, cprop, sprop, seqcprop, geo, lt,2);
      check += qcd_deltas_f1_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, cprop, sprop, seqcprop, geo, lt,3);
      check += qcd_deltas_f1_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, cprop, sprop, seqcprop, geo, lt,4);  
      check += qcd_deltas_f1_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, cprop, sprop, seqcprop, geo, lt,5);
      check += qcd_deltas_f1_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, cprop, sprop, seqcprop, geo, lt,6);
      check += qcd_deltas_f1_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, cprop, sprop, seqcprop, geo, lt,7);  
      check += qcd_deltas_f1_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, cprop, sprop, seqcprop, geo, lt,8);  
      break;
    }
    break;

    //-----------------------------------------------------------------------------------	   

  case 53://-xistar_zero extra 11
    switch(np){
    case 0:
      check =  qcd_xistar_extra11_f2_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, sprop, uprop, sequprop, geo, lt,0);
      check += qcd_xistar_extra11_f2_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, sprop, uprop, sequprop, geo, lt,1);  
      check += qcd_xistar_extra11_f2_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, sprop, uprop, sequprop, geo, lt,2);
      check += qcd_xistar_extra11_f2_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, sprop, uprop, sequprop, geo, lt,3);
      check += qcd_xistar_extra11_f2_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, sprop, uprop, sequprop, geo, lt,4);  
      check += qcd_xistar_extra11_f2_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, sprop, uprop, sequprop, geo, lt,5);
      check += qcd_xistar_extra11_f2_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, sprop, uprop, sequprop, geo, lt,6);
      check += qcd_xistar_extra11_f2_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, sprop, uprop, sequprop, geo, lt,7);  
      check += qcd_xistar_extra11_f2_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, sprop, uprop, sequprop, geo, lt,8); 
      break;
    case 2:
      check =  qcd_xistar_extra11_f1_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, sprop, uprop, seqsprop, geo, lt,0);
      check += qcd_xistar_extra11_f1_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, sprop, uprop, seqsprop, geo, lt,1);  
      check += qcd_xistar_extra11_f1_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, sprop, uprop, seqsprop, geo, lt,2);
      check += qcd_xistar_extra11_f1_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, sprop, uprop, seqsprop, geo, lt,3);
      check += qcd_xistar_extra11_f1_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, sprop, uprop, seqsprop, geo, lt,4);  
      check += qcd_xistar_extra11_f1_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, sprop, uprop, seqsprop, geo, lt,5);
      check += qcd_xistar_extra11_f1_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, sprop, uprop, seqsprop, geo, lt,6);
      check += qcd_xistar_extra11_f1_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, sprop, uprop, seqsprop, geo, lt,7);  
      check += qcd_xistar_extra11_f1_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, sprop, uprop, seqsprop, geo, lt,8);  
      break;
    }
    break;

  case 54://-xistar_zero extra 12
    switch(np){
    case 0:
      check =  qcd_xistar_extra12_f2_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, sprop, uprop, sequprop, geo, lt,0);
      check += qcd_xistar_extra12_f2_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, sprop, uprop, sequprop, geo, lt,1);  
      check += qcd_xistar_extra12_f2_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, sprop, uprop, sequprop, geo, lt,2);
      check += qcd_xistar_extra12_f2_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, sprop, uprop, sequprop, geo, lt,3);
      check += qcd_xistar_extra12_f2_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, sprop, uprop, sequprop, geo, lt,4);  
      check += qcd_xistar_extra12_f2_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, sprop, uprop, sequprop, geo, lt,5);
      check += qcd_xistar_extra12_f2_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, sprop, uprop, sequprop, geo, lt,6);
      check += qcd_xistar_extra12_f2_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, sprop, uprop, sequprop, geo, lt,7);  
      check += qcd_xistar_extra12_f2_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, sprop, uprop, sequprop, geo, lt,8); 
      break;
    case 2:
      check =  qcd_xistar_extra12_f1_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, sprop, uprop, seqsprop, geo, lt,0);
      check += qcd_xistar_extra12_f1_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, sprop, uprop, seqsprop, geo, lt,1);  
      check += qcd_xistar_extra12_f1_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, sprop, uprop, seqsprop, geo, lt,2);
      check += qcd_xistar_extra12_f1_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, sprop, uprop, seqsprop, geo, lt,3);
      check += qcd_xistar_extra12_f1_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, sprop, uprop, seqsprop, geo, lt,4);  
      check += qcd_xistar_extra12_f1_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, sprop, uprop, seqsprop, geo, lt,5);
      check += qcd_xistar_extra12_f1_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, sprop, uprop, seqsprop, geo, lt,6);
      check += qcd_xistar_extra12_f1_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, sprop, uprop, seqsprop, geo, lt,7);  
      check += qcd_xistar_extra12_f1_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, sprop, uprop, seqsprop, geo, lt,8);  
      break;
    }
    break;
  case 55://-xistar_zero extra 21
    switch(np){
    case 0:
      check =  qcd_xistar_extra21_f2_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, sprop, uprop, sequprop, geo, lt,0);
      check += qcd_xistar_extra21_f2_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, sprop, uprop, sequprop, geo, lt,1);  
      check += qcd_xistar_extra21_f2_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, sprop, uprop, sequprop, geo, lt,2);
      check += qcd_xistar_extra21_f2_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, sprop, uprop, sequprop, geo, lt,3);
      check += qcd_xistar_extra21_f2_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, sprop, uprop, sequprop, geo, lt,4);  
      check += qcd_xistar_extra21_f2_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, sprop, uprop, sequprop, geo, lt,5);
      check += qcd_xistar_extra21_f2_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, sprop, uprop, sequprop, geo, lt,6);
      check += qcd_xistar_extra21_f2_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, sprop, uprop, sequprop, geo, lt,7);  
      check += qcd_xistar_extra21_f2_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, sprop, uprop, sequprop, geo, lt,8); 
      break;
    case 2:
      check =  qcd_xistar_extra21_f1_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, sprop, uprop, seqsprop, geo, lt,0);
      check += qcd_xistar_extra21_f1_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, sprop, uprop, seqsprop, geo, lt,1);  
      check += qcd_xistar_extra21_f1_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, sprop, uprop, seqsprop, geo, lt,2);
      check += qcd_xistar_extra21_f1_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, sprop, uprop, seqsprop, geo, lt,3);
      check += qcd_xistar_extra21_f1_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, sprop, uprop, seqsprop, geo, lt,4);  
      check += qcd_xistar_extra21_f1_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, sprop, uprop, seqsprop, geo, lt,5);
      check += qcd_xistar_extra21_f1_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, sprop, uprop, seqsprop, geo, lt,6);
      check += qcd_xistar_extra21_f1_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, sprop, uprop, seqsprop, geo, lt,7);  
      check += qcd_xistar_extra21_f1_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, sprop, uprop, seqsprop, geo, lt,8);  
      break;
    }
    break;
  case 56://-xistar_zero extra 22
    switch(np){
    case 0:
      check =  qcd_xistar_extra22_f2_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, sprop, uprop, sequprop, geo, lt,0);
      check += qcd_xistar_extra22_f2_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, sprop, uprop, sequprop, geo, lt,1);  
      check += qcd_xistar_extra22_f2_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, sprop, uprop, sequprop, geo, lt,2);
      check += qcd_xistar_extra22_f2_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, sprop, uprop, sequprop, geo, lt,3);
      check += qcd_xistar_extra22_f2_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, sprop, uprop, sequprop, geo, lt,4);  
      check += qcd_xistar_extra22_f2_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, sprop, uprop, sequprop, geo, lt,5);
      check += qcd_xistar_extra22_f2_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, sprop, uprop, sequprop, geo, lt,6);
      check += qcd_xistar_extra22_f2_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, sprop, uprop, sequprop, geo, lt,7);  
      check += qcd_xistar_extra22_f2_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, sprop, uprop, sequprop, geo, lt,8); 
      break;
    case 2:
      check =  qcd_xistar_extra22_f1_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, sprop, uprop, seqsprop, geo, lt,0);
      check += qcd_xistar_extra22_f1_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, sprop, uprop, seqsprop, geo, lt,1);  
      check += qcd_xistar_extra22_f1_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, sprop, uprop, seqsprop, geo, lt,2);
      check += qcd_xistar_extra22_f1_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, sprop, uprop, seqsprop, geo, lt,3);
      check += qcd_xistar_extra22_f1_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, sprop, uprop, seqsprop, geo, lt,4);  
      check += qcd_xistar_extra22_f1_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, sprop, uprop, seqsprop, geo, lt,5);
      check += qcd_xistar_extra22_f1_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, sprop, uprop, seqsprop, geo, lt,6);
      check += qcd_xistar_extra22_f1_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, sprop, uprop, seqsprop, geo, lt,7);  
      check += qcd_xistar_extra22_f1_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, sprop, uprop, seqsprop, geo, lt,8);  
      break;
    }
    break;

    //-----------------------------------------------------------------------------------	   

  case 57://-xistar_minus extra 11
    switch(np){
    case 1:
      check =  qcd_xistar_extra11_f2_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, sprop, dprop, seqdprop, geo, lt,0);
      check += qcd_xistar_extra11_f2_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, sprop, dprop, seqdprop, geo, lt,1);  
      check += qcd_xistar_extra11_f2_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, sprop, dprop, seqdprop, geo, lt,2);
      check += qcd_xistar_extra11_f2_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, sprop, dprop, seqdprop, geo, lt,3);
      check += qcd_xistar_extra11_f2_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, sprop, dprop, seqdprop, geo, lt,4);  
      check += qcd_xistar_extra11_f2_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, sprop, dprop, seqdprop, geo, lt,5);
      check += qcd_xistar_extra11_f2_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, sprop, dprop, seqdprop, geo, lt,6);
      check += qcd_xistar_extra11_f2_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, sprop, dprop, seqdprop, geo, lt,7);  
      check += qcd_xistar_extra11_f2_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, sprop, dprop, seqdprop, geo, lt,8); 
      break;
    case 2:
      check =  qcd_xistar_extra11_f1_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, sprop, dprop, seqsprop, geo, lt,0);
      check += qcd_xistar_extra11_f1_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, sprop, dprop, seqsprop, geo, lt,1);  
      check += qcd_xistar_extra11_f1_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, sprop, dprop, seqsprop, geo, lt,2);
      check += qcd_xistar_extra11_f1_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, sprop, dprop, seqsprop, geo, lt,3);
      check += qcd_xistar_extra11_f1_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, sprop, dprop, seqsprop, geo, lt,4);  
      check += qcd_xistar_extra11_f1_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, sprop, dprop, seqsprop, geo, lt,5);
      check += qcd_xistar_extra11_f1_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, sprop, dprop, seqsprop, geo, lt,6);
      check += qcd_xistar_extra11_f1_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, sprop, dprop, seqsprop, geo, lt,7);  
      check += qcd_xistar_extra11_f1_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, sprop, dprop, seqsprop, geo, lt,8);  
      break;
    }
    break;

  case 58://-xistar_minus extra 12
    switch(np){
    case 1:
      check =  qcd_xistar_extra12_f2_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, sprop, dprop, seqdprop, geo, lt,0);
      check += qcd_xistar_extra12_f2_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, sprop, dprop, seqdprop, geo, lt,1);  
      check += qcd_xistar_extra12_f2_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, sprop, dprop, seqdprop, geo, lt,2);
      check += qcd_xistar_extra12_f2_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, sprop, dprop, seqdprop, geo, lt,3);
      check += qcd_xistar_extra12_f2_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, sprop, dprop, seqdprop, geo, lt,4);  
      check += qcd_xistar_extra12_f2_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, sprop, dprop, seqdprop, geo, lt,5);
      check += qcd_xistar_extra12_f2_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, sprop, dprop, seqdprop, geo, lt,6);
      check += qcd_xistar_extra12_f2_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, sprop, dprop, seqdprop, geo, lt,7);  
      check += qcd_xistar_extra12_f2_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, sprop, dprop, seqdprop, geo, lt,8); 
      break;
    case 2:
      check =  qcd_xistar_extra12_f1_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, sprop, dprop, seqsprop, geo, lt,0);
      check += qcd_xistar_extra12_f1_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, sprop, dprop, seqsprop, geo, lt,1);  
      check += qcd_xistar_extra12_f1_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, sprop, dprop, seqsprop, geo, lt,2);
      check += qcd_xistar_extra12_f1_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, sprop, dprop, seqsprop, geo, lt,3);
      check += qcd_xistar_extra12_f1_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, sprop, dprop, seqsprop, geo, lt,4);  
      check += qcd_xistar_extra12_f1_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, sprop, dprop, seqsprop, geo, lt,5);
      check += qcd_xistar_extra12_f1_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, sprop, dprop, seqsprop, geo, lt,6);
      check += qcd_xistar_extra12_f1_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, sprop, dprop, seqsprop, geo, lt,7);  
      check += qcd_xistar_extra12_f1_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, sprop, dprop, seqsprop, geo, lt,8);  
      break;
    }
    break;
  case 59://-xistar_minus extra 21
    switch(np){
    case 1:
      check =  qcd_xistar_extra21_f2_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, sprop, dprop, seqdprop, geo, lt,0);
      check += qcd_xistar_extra21_f2_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, sprop, dprop, seqdprop, geo, lt,1);  
      check += qcd_xistar_extra21_f2_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, sprop, dprop, seqdprop, geo, lt,2);
      check += qcd_xistar_extra21_f2_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, sprop, dprop, seqdprop, geo, lt,3);
      check += qcd_xistar_extra21_f2_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, sprop, dprop, seqdprop, geo, lt,4);  
      check += qcd_xistar_extra21_f2_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, sprop, dprop, seqdprop, geo, lt,5);
      check += qcd_xistar_extra21_f2_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, sprop, dprop, seqdprop, geo, lt,6);
      check += qcd_xistar_extra21_f2_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, sprop, dprop, seqdprop, geo, lt,7);  
      check += qcd_xistar_extra21_f2_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, sprop, dprop, seqdprop, geo, lt,8); 
      break;
    case 2:
      check =  qcd_xistar_extra21_f1_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, sprop, dprop, seqsprop, geo, lt,0);
      check += qcd_xistar_extra21_f1_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, sprop, dprop, seqsprop, geo, lt,1);  
      check += qcd_xistar_extra21_f1_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, sprop, dprop, seqsprop, geo, lt,2);
      check += qcd_xistar_extra21_f1_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, sprop, dprop, seqsprop, geo, lt,3);
      check += qcd_xistar_extra21_f1_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, sprop, dprop, seqsprop, geo, lt,4);  
      check += qcd_xistar_extra21_f1_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, sprop, dprop, seqsprop, geo, lt,5);
      check += qcd_xistar_extra21_f1_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, sprop, dprop, seqsprop, geo, lt,6);
      check += qcd_xistar_extra21_f1_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, sprop, dprop, seqsprop, geo, lt,7);  
      check += qcd_xistar_extra21_f1_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, sprop, dprop, seqsprop, geo, lt,8);  
      break;
    }
    break;
  case 60://-xistar_minus extra 22
    switch(np){
    case 1:
      check =  qcd_xistar_extra22_f2_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, sprop, dprop, seqdprop, geo, lt,0);
      check += qcd_xistar_extra22_f2_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, sprop, dprop, seqdprop, geo, lt,1);  
      check += qcd_xistar_extra22_f2_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, sprop, dprop, seqdprop, geo, lt,2);
      check += qcd_xistar_extra22_f2_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, sprop, dprop, seqdprop, geo, lt,3);
      check += qcd_xistar_extra22_f2_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, sprop, dprop, seqdprop, geo, lt,4);  
      check += qcd_xistar_extra22_f2_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, sprop, dprop, seqdprop, geo, lt,5);
      check += qcd_xistar_extra22_f2_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, sprop, dprop, seqdprop, geo, lt,6);
      check += qcd_xistar_extra22_f2_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, sprop, dprop, seqdprop, geo, lt,7);  
      check += qcd_xistar_extra22_f2_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, sprop, dprop, seqdprop, geo, lt,8); 
      break;
    case 2:
      check =  qcd_xistar_extra22_f1_3pt(ctr1 , cg1cg1_ind, cg1cg1_val, block_pr, sprop, dprop, seqsprop, geo, lt,0);
      check += qcd_xistar_extra22_f1_3pt(ctr12, cg1cg2_ind, cg1cg2_val, block_pr, sprop, dprop, seqsprop, geo, lt,1);  
      check += qcd_xistar_extra22_f1_3pt(ctr13, cg1cg3_ind, cg1cg3_val, block_pr, sprop, dprop, seqsprop, geo, lt,2);
      check += qcd_xistar_extra22_f1_3pt(ctr21, cg2cg1_ind, cg2cg1_val, block_pr, sprop, dprop, seqsprop, geo, lt,3);
      check += qcd_xistar_extra22_f1_3pt(ctr2,  cg2cg2_ind, cg2cg2_val, block_pr, sprop, dprop, seqsprop, geo, lt,4);  
      check += qcd_xistar_extra22_f1_3pt(ctr23, cg2cg3_ind, cg2cg3_val, block_pr, sprop, dprop, seqsprop, geo, lt,5);
      check += qcd_xistar_extra22_f1_3pt(ctr31, cg3cg1_ind, cg3cg1_val, block_pr, sprop, dprop, seqsprop, geo, lt,6);
      check += qcd_xistar_extra22_f1_3pt(ctr32, cg3cg2_ind, cg3cg2_val, block_pr, sprop, dprop, seqsprop, geo, lt,7);  
      check += qcd_xistar_extra22_f1_3pt(ctr3 , cg3cg3_ind, cg3cg3_val, block_pr, sprop, dprop, seqsprop, geo, lt,8);  
      break;
    }
    break;

  }//-main switch	

  if(particles32[p_id]){
    check += qcd_projector_pr32_3pt_spin32(block,block_pr,gamma12,gamma13,gamma23,geo);
    if(check==10){
      if(geo->myid==0) printf("Particle: %s part %s t=%d\n",particle_names[p_id],particles_parts[p_id][np], lt);
    }
  }
  else{
    check += qcd_projector_3pt_spin12(block,block_pr,geo);
    if(check==2){
      if(geo->myid==0) printf("Particle: %s part %s t=%d\n",particle_names[p_id],particles_parts[p_id][np], lt);
    }
  }
  
  if( particles32[p_id] ){
    for(i=0;i<9;i++)
      for(j=0;j<4;j++)
	for(k=0;k<4;k++){
	  free(block_pr[i][j][k]);
	}							
  }
  else{
    for(j=0;j<4;j++)
      for(k=0;k<4;k++){
	free(block_pr[8][j][k]);
      }				
  }        
	
    
}//--main function closes

