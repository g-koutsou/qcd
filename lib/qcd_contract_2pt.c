/*
 * Christos Kallidonis
 * April 2012
 * 
 * This program contains contractions for the 2-point functions
 * of the 40 particles.
 * 
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <qcd.h>

int qcd_f1f2f1_2pt(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[4][4],
						qcd_propagator *propf1,qcd_propagator *propf2, qcd_geometry *geo, qcd_uint_4 lt){
							
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
          

		
		block[ga][gap][v3] = qcd_CADD(block[ga][gap][v3],qcd_CSCALE(
									    qcd_CMUL(cgcg_val[ctr2],
										     qcd_CMUL(propf2->D[v][be][bep][b][bp],
											      qcd_CMUL(propf1->D[v][al][alp][a][ap],
												       propf1->D[v][ga][gap][c][cp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2])
					      );
		
		block[ga][gap][v3] = qcd_CSUB(block[ga][gap][v3],qcd_CSCALE(
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

/*====================================================================*/

int qcd_pnextra_2pt(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[4][4],
						qcd_propagator *propf1,qcd_propagator *propf2, qcd_geometry *geo, qcd_uint_4 lt){
							
  qcd_int_4 ctr2;
  qcd_uint_2 a,b,c,ap,bp,cp,al,be,ga,alp,bep,gap,cc1,cc2,de,dep;
  qcd_uint_4 lx,ly,lz,v,v3;
  qcd_complex_16 term;
  int lv3 = geo->lL[1]*geo->lL[2]*geo->lL[3];

#pragma omp parallel for private(lz,ly,lx,v,de,dep,ga,gap,ctr2,al,be,bep,alp,cc1,a,b,c,cc2,ap,bp,cp)
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



	    block[de][dep][v3] = qcd_CADD(block[de][dep][v3],qcd_CMUL(
					     qcd_CMUL(cgcg_val[ctr2],
					     qcd_CMUL(propf2->D[v][be][bep][b][bp],
					     qcd_CMUL(propf1->D[v][al][alp][a][ap],
					              propf1->D[v][ga][gap][c][cp]))),term)
				);
										    
	    block[de][dep][v3] = qcd_CSUB(block[de][dep][v3],qcd_CMUL(
					     qcd_CMUL(cgcg_val[ctr2],
					     qcd_CMUL(propf2->D[v][be][bep][b][bp],
					     qcd_CMUL(propf1->D[v][al][gap][a][cp],
					              propf1->D[v][ga][alp][c][ap]))),term)
				);
	  }//space loop
	}//color2 loop    
	}//color1 loop
      }//nonvanishing cgcg loop
	  }//ga, gap						
    }//de,dep

	
  return (1);
}

/*====================================================================*/

int qcd_f1f2f3_2pt(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[4][4],
						qcd_propagator *propf1, qcd_propagator *propf2, qcd_propagator *propf3, qcd_geometry *geo, qcd_uint_4 lt){
							
  qcd_int_4 ctr2;
  qcd_uint_2 a,b,c,ap,bp,cp,al,be,ga,alp,bep,gap,cc1,cc2;
  qcd_uint_4 lx,ly,lz,v,v3;
  int lv3 = geo->lL[1]*geo->lL[2]*geo->lL[3];

#pragma omp parallel for private(lz,ly,lx,v,ga,gap,ctr2,al,be,bep,alp,cc1,a,b,c,cc2,ap,bp,cp)
	  for(v3=0; v3<lv3; v3++) {
	  // for(lx=0; lx<geo->lL[1]; lx++) 
	  //   for(ly=0; ly<geo->lL[2]; ly++) 
	  //     for(lz=0; lz<geo->lL[3]; lz++){ 
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

	    block[ga][gap][v3] = qcd_CADD(block[ga][gap][v3],qcd_CSCALE(
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

/*====================================================================*/

int qcd_f123f321_2pt(qcd_int_4 ctr, qcd_int_2 cgcg_ind[16*16][4],qcd_complex_16 cgcg_val[16*16], qcd_complex_16 *block[4][4],
						qcd_propagator *propf1, qcd_propagator *propf2, qcd_propagator *propf3, qcd_geometry *geo, qcd_uint_4 lt){
							
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
  
	    block[ga][gap][v3] = qcd_CADD(block[ga][gap][v3],qcd_CSCALE(
					     qcd_CMUL(cgcg_val[ctr2],
					     qcd_CMUL(propf1->D[v][al][alp][a][ap],
					     qcd_CMUL(propf3->D[v][ga][gap][c][cp],
					              propf2->D[v][be][bep][b][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*0.5)
				);
	    
	    block[ga][gap][v3] = qcd_CADD(block[ga][gap][v3],qcd_CSCALE(
					     qcd_CMUL(cgcg_val[ctr2],
					     qcd_CMUL(propf1->D[v][ga][gap][c][cp],
					     qcd_CMUL(propf3->D[v][al][alp][a][ap],
					              propf2->D[v][be][bep][b][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*0.5)
				);	    
										    
	    block[ga][gap][v3] = qcd_CSUB(block[ga][gap][v3],qcd_CSCALE(
					     qcd_CMUL(cgcg_val[ctr2],
					     qcd_CMUL(propf1->D[v][al][gap][a][cp],
					     qcd_CMUL(propf3->D[v][ga][alp][c][ap],
					              propf2->D[v][be][bep][b][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*0.5)
				);
	    
	    block[ga][gap][v3] = qcd_CSUB(block[ga][gap][v3],qcd_CSCALE(
					     qcd_CMUL(cgcg_val[ctr2],
					     qcd_CMUL(propf1->D[v][ga][alp][c][ap],
					     qcd_CMUL(propf3->D[v][al][gap][a][cp],
					              propf2->D[v][be][bep][b][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*0.5)
				);
				
	}//color2 loop                                                                                                                                                                               
      }//color1 loop                                                                                                                                                                               
      }//spin loop (cg5cg5)                                                                                                                                                                         
      }//-ga,gap

      }//space loop                                                                                                                                                                                 
    
  return (1);
}

/*====================================================================*/

int qcd_lambdas_xis_2pt(qcd_int_4 ctr, qcd_int_2 cg5cg5_ind[16*16][4],qcd_complex_16 cg5cg5_val[16*16], qcd_complex_16 *block[4][4],
						qcd_propagator *propf1, qcd_propagator *propf2, qcd_propagator *propf3, qcd_geometry *geo, qcd_uint_4 lt){
							
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
	  al = cg5cg5_ind[ctr2][0];
	  be = cg5cg5_ind[ctr2][1];
	  bep= cg5cg5_ind[ctr2][2];
	  alp= cg5cg5_ind[ctr2][3];                              
            
          for(cc1=0;cc1<6;cc1++){
	    a=qcd_EPS[cc1][0];
	    b=qcd_EPS[cc1][1];
	    c=qcd_EPS[cc1][2];
			
	    for(cc2=0;cc2<6;cc2++){          
	      ap=qcd_EPS[cc2][0];
	      bp=qcd_EPS[cc2][1];
	      cp=qcd_EPS[cc2][2];
    
	    //--Diagonal terms

	    block[ga][gap][v3] = qcd_CADD(block[ga][gap][v3],qcd_CSCALE(
					     qcd_CMUL(cg5cg5_val[ctr2],
					     qcd_CMUL(propf1->D[v][al][alp][a][ap],
					     qcd_CMUL(propf2->D[v][be][bep][b][bp],
								 propf3->D[v][ga][gap][c][cp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*4.0/6.0)
				 ); //-ok

	    block[ga][gap][v3] = qcd_CADD(block[ga][gap][v3],qcd_CSCALE(
					     qcd_CMUL(cg5cg5_val[ctr2],
					     qcd_CMUL(propf1->D[v][al][alp][a][ap],
					     qcd_CMUL(propf2->D[v][ga][gap][c][cp],
					             propf3->D[v][be][bep][b][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]/6.0)
				 ); //-ok	    
	    
	    block[ga][gap][v3] = qcd_CADD(block[ga][gap][v3],qcd_CSCALE(
					     qcd_CMUL(cg5cg5_val[ctr2],
					     qcd_CMUL(propf1->D[v][ga][gap][c][cp],
					     qcd_CMUL(propf2->D[v][al][alp][a][ap],
					             propf3->D[v][be][bep][b][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]/6.0)
				 ); //-ok
	    //----------------------

	    block[ga][gap][v3] = qcd_CADD(block[ga][gap][v3],qcd_CSCALE(
					     qcd_CMUL(cg5cg5_val[ctr2],
					     qcd_CMUL(propf1->D[v][ga][alp][c][ap],
					     qcd_CMUL(propf2->D[v][al][gap][a][cp],
					             propf3->D[v][be][bep][b][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]/6.0)
				 ); //-ok	    
	    
	    block[ga][gap][v3] = qcd_CADD(block[ga][gap][v3],qcd_CSCALE(
					     qcd_CMUL(cg5cg5_val[ctr2],
					     qcd_CMUL(propf1->D[v][al][gap][a][cp],
					     qcd_CMUL(propf2->D[v][ga][alp][c][ap],
					             propf3->D[v][be][bep][b][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]/6.0)
				 ); //-ok
	    
	    
	    block[ga][gap][v3] = qcd_CSUB(block[ga][gap][v3],qcd_CSCALE(
					     qcd_CMUL(cg5cg5_val[ctr2],
					     qcd_CMUL(propf1->D[v][al][alp][a][ap],
					     qcd_CMUL(propf2->D[v][ga][bep][c][bp],
					             propf3->D[v][be][gap][b][cp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*2.0/6.0)
				 ); //-ok
	    
	    block[ga][gap][v3] = qcd_CSUB(block[ga][gap][v3],qcd_CSCALE(
					     qcd_CMUL(cg5cg5_val[ctr2],
					     qcd_CMUL(propf1->D[v][al][alp][a][ap],
					     qcd_CMUL(propf2->D[v][be][gap][b][cp],
					             propf3->D[v][ga][bep][c][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*2.0/6.0)
				 ); //-ok

	    block[ga][gap][v3] = qcd_CSUB(block[ga][gap][v3],qcd_CSCALE(
					     qcd_CMUL(cg5cg5_val[ctr2],
					     qcd_CMUL(propf1->D[v][al][gap][a][cp],
					     qcd_CMUL(propf2->D[v][be][alp][b][ap],
					             propf3->D[v][ga][bep][c][bp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*2.0/6.0)
				 ); //-ok	    

	    block[ga][gap][v3] = qcd_CSUB(block[ga][gap][v3],qcd_CSCALE(
					     qcd_CMUL(cg5cg5_val[ctr2],
					     qcd_CMUL(propf1->D[v][ga][alp][c][ap],
					     qcd_CMUL(propf2->D[v][al][bep][a][bp],
					             propf3->D[v][be][gap][b][cp]))),qcd_SGN_EPS[cc1]*qcd_SGN_EPS[cc2]*2.0/6.0)
				 ); //-ok   
	    
	}//color2 loop
	}//color1 loop    
	}//spin loop (cg5cg5)
    }//-ga,gap
    }//space loop

	
  return (1);
}


/*====================================================================*/
/*====================================================================*/
/*====================================================================*/

void qcd_contractions2pt(qcd_uint_4 particle_id, qcd_complex_16 *block[4][4],
			 qcd_propagator *uprop,qcd_propagator *dprop,qcd_propagator *sprop,qcd_propagator *cprop,
			 qcd_geometry *geo, qcd_uint_4 lt){

  qcd_int_4 ctr,ctr2 = 0;
  qcd_uint_2 al,be,alp,bep,ga;
  qcd_complex_16 C;

  qcd_int_2 cg5cg5_ind[16*16][4],cgcg_ind[16*16][4];
  qcd_complex_16 cg5cg5_val[16*16],cgcg_val[16*16];
  qcd_complex_16 CG[4][4],CG_bar[4][4];
  
//-- Calculate the Charge Conjugation and non zero elements for proton and neutron extra
  if( (particle_id==41) || (particle_id==42) ){  
    for(al=0;al<4;al++)
      for(be=0;be<4;be++){
	CG[al][be]     = (qcd_complex_16) {0.0,0.0};
	CG_bar[al][be] = (qcd_complex_16) {0.0,0.0};
	for(ga=0;ga<4;ga++){
	  CG[al][be] = qcd_CADD( CG[al][be],qcd_CMUL(qcd_GAMMA[4][al][ga],qcd_GAMMA[2][ga][be]) );		
	  CG_bar[al][be] = CG[al][be];
	}
      }  
    
    ctr2 = 0;
    for(al=0;al<4;al++)
      for(be=0;be<4;be++)
	for(alp=0;alp<4;alp++)
	  for(bep=0;bep<4;bep++)
	    {  
	      C  = qcd_CMUL(CG[al][be],CG_bar[bep][alp]);
	      
	      if(qcd_NORM(C)>1e-3)
		{
		  cgcg_val[ctr2].re = C.re;
		  cgcg_val[ctr2].im = C.im;
		  cgcg_ind[ctr2][0] = al;
		  cgcg_ind[ctr2][1] = be;
		  cgcg_ind[ctr2][2] = bep;
		  cgcg_ind[ctr2][3] = alp;                                                            
		  ctr2++;
		}	  	  
	    }	
  }//-if
  
//-- Non zero elements of (Cg5)*(Cg5)
  ctr = 0;
  for(al=0;al<4;al++)
  for(be=0;be<4;be++)
  for(alp=0;alp<4;alp++)
  for(bep=0;bep<4;bep++)
  {  
    C  = qcd_CMUL(qcd_CGAMMA[5][al][be],qcd_BAR_CGAMMA[5][bep][alp]);    
    if(qcd_NORM(C)>1e-3)
    {
        cg5cg5_val[ctr].re = C.re;
        cg5cg5_val[ctr].im = C.im;
        cg5cg5_ind[ctr][0] = al;
        cg5cg5_ind[ctr][1] = be;
        cg5cg5_ind[ctr][2] = bep;
        cg5cg5_ind[ctr][3] = alp;                                                            
        ctr++;
    }     
  }  
  
  
  switch(particle_id){
      case 1:
	  if ( qcd_f1f2f1_2pt(ctr, cg5cg5_ind, cg5cg5_val, block, uprop, dprop, geo, lt) ){						 
	    if(geo->myid==0) printf("proton t=%d\n",lt);
	  }
	  break;
      case 2:
	  if (qcd_f1f2f1_2pt(ctr, cg5cg5_ind, cg5cg5_val, block, dprop, uprop, geo, lt) ){						 
	    if(geo->myid==0) printf("neutron t=%d\n",lt);
	  }
	  break;
      case 3:
	  if ( qcd_lambdas_xis_2pt(ctr, cg5cg5_ind, cg5cg5_val, block, uprop, dprop, sprop, geo, lt) ){						 
	    if(geo->myid==0) printf("lambda t=%d\n",lt);
	  }
	  break;			
      case 4:
	  if ( qcd_f1f2f1_2pt(ctr, cg5cg5_ind, cg5cg5_val, block, uprop, sprop, geo, lt) ){						 
	    if(geo->myid==0) printf("sigmaplus t=%d\n",lt);
	  }
	  break;
      case 5:
	  if ( qcd_f123f321_2pt(ctr, cg5cg5_ind, cg5cg5_val, block, uprop, sprop, dprop, geo, lt) ){						 
	    if(geo->myid==0) printf("sigmazero t=%d\n",lt);
	  }
	  break;
      case 6:
	  if ( qcd_f1f2f1_2pt(ctr, cg5cg5_ind, cg5cg5_val, block, dprop, sprop, geo, lt) ){						 
	    if(geo->myid==0) printf("sigmaminus t=%d\n",lt);
	  }
	  break;			
      case 7:
	  if ( qcd_f1f2f1_2pt(ctr, cg5cg5_ind, cg5cg5_val, block, sprop, uprop, geo, lt) ){						 
	    if(geo->myid==0) printf("xizero t=%d\n",lt);
	  }
	  break;
      case 8:
	  if ( qcd_f1f2f1_2pt(ctr, cg5cg5_ind, cg5cg5_val, block, sprop, dprop, geo, lt) ){						 
	    if(geo->myid==0) printf("ximinus t=%d\n",lt);
	  }
	  break;
//---------------------------------------------------------------------	  
      case 19:
	  if ( qcd_lambdas_xis_2pt(ctr, cg5cg5_ind, cg5cg5_val, block, uprop, dprop, cprop, geo, lt) ){						 
	    if(geo->myid==0) printf("lambdac_plus t=%d\n",lt);
	  }
	  break;	  
      case 20:
	  if ( qcd_f1f2f1_2pt(ctr, cg5cg5_ind, cg5cg5_val, block, uprop, cprop, geo, lt) ){						 
	    if(geo->myid==0) printf("sigmac_plusplus t=%d\n",lt);
	  }
	  break;
      case 21:
	  if ( qcd_f123f321_2pt(ctr, cg5cg5_ind, cg5cg5_val, block, uprop, cprop, dprop, geo, lt) ){						 
	    if(geo->myid==0) printf("sigmac_plus t=%d\n",lt);
	  }
	  break;
      case 22:
	  if ( qcd_f1f2f1_2pt(ctr, cg5cg5_ind, cg5cg5_val, block, dprop, cprop, geo, lt) ){						 
	    if(geo->myid==0) printf("sigmac_zero t=%d\n",lt);
	  }
	  break;
      case 23:
	  if ( qcd_f1f2f3_2pt(ctr, cg5cg5_ind, cg5cg5_val, block, uprop, sprop, cprop, geo, lt) ){						 
	    if(geo->myid==0) printf("xicplus t=%d\n",lt);
	  }
	  break;
      case 24:
	  if ( qcd_f123f321_2pt(ctr, cg5cg5_ind, cg5cg5_val, block, uprop, cprop, sprop, geo, lt) ){						 
	    if(geo->myid==0) printf("xicplus_prime t=%d\n",lt);
	  }
	  break;
      case 25:
	  if ( qcd_f1f2f3_2pt(ctr, cg5cg5_ind, cg5cg5_val, block, dprop, sprop, cprop, geo, lt) ){						 
	    if(geo->myid==0) printf("xiczero t=%d\n",lt);
	  }
	  break;	  
      case 26:
	  if (qcd_f123f321_2pt(ctr, cg5cg5_ind, cg5cg5_val, block, dprop, cprop, sprop, geo, lt) ){						 
	    if(geo->myid==0) printf("xiczero_prime t=%d\n",lt);
	  }
	  break;	  
      case 27:
	  if (  qcd_f1f2f1_2pt(ctr, cg5cg5_ind, cg5cg5_val, block, sprop, cprop, geo, lt) ){						 
	    if(geo->myid==0) printf("omegac_zero t=%d\n",lt);
	  }
	  break;
//-----------------------------------------------------------------------------------       
      case 34:
	  if ( qcd_f1f2f1_2pt(ctr, cg5cg5_ind, cg5cg5_val, block, cprop, uprop, geo, lt) ){						 
	    if(geo->myid==0) printf("xiccplusplus t=%d\n",lt);
	  }
	  break;
      case 35:
	  if ( qcd_f1f2f1_2pt(ctr, cg5cg5_ind, cg5cg5_val, block, cprop, dprop, geo, lt) ){						 
	    if(geo->myid==0) printf("xiccplus t=%d\n",lt);
	  }
	  break;	  
      case 36:
	  if ( qcd_f1f2f1_2pt(ctr, cg5cg5_ind, cg5cg5_val, block, cprop, sprop, geo, lt) ){						 
	    if(geo->myid==0) printf("omegacc_plus t=%d\n",lt);
	  }
	  break;
//----------------------------------------------------------------------------------- 	  
      case 41:
	  if ( qcd_pnextra_2pt(ctr2, cgcg_ind, cgcg_val, block, uprop, dprop, geo, lt) ){						 
	    if(geo->myid==0) printf("proton_extra t=%d\n",lt);
	  }
	  break;
      case 42:
	  if ( qcd_pnextra_2pt(ctr2, cgcg_ind, cgcg_val, block, dprop, uprop, geo, lt) ){						 
	    if(geo->myid==0) printf("neutron_extra t=%d\n",lt);
	  }
	  break;
//----------------------------------------------------------------------------------- 
      case 45:
	  if ( qcd_lambdas_xis_2pt(ctr, cg5cg5_ind, cg5cg5_val, block, sprop, uprop, cprop, geo, lt) ){						 
	    if(geo->myid==0) printf("xicplus_extra t=%d\n",lt);
	  }
	  break;
      case 46:
	  if ( qcd_lambdas_xis_2pt(ctr, cg5cg5_ind, cg5cg5_val, block, sprop, dprop, cprop, geo, lt) ){						 
	    if(geo->myid==0) printf("xiczero_extra t=%d\n",lt);
	  }
	  break; 	  	  
//-----------------------------------------------------------------------------------	  	   
  }//-switch	

						  
}//--main function closes






